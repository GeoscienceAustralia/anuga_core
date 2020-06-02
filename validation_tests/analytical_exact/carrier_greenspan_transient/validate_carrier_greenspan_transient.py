"""Automatic verification of ANUGA flows.
See functions exercised by this wrapper for more details
"""

import unittest
import os
import numpy
import anuga

from math import sqrt

import warnings
warnings.simplefilter('ignore')

indent = anuga.indent

args = anuga.get_args()
verbose = args.verbose

class Test_results(unittest.TestCase):
    def setUp(self):
        for file in os.listdir('.'):    
            if file.endswith('.stdout') or\
                    file.endswith('.sww') or\
                    file.endswith('.msh') or\
                    file.endswith('.png'):
                os.remove(file)
                
        
    def tearDown(self):
        pass

    def test_carrier_greenspan_transient(self):
    

        if verbose:
            print()
            print(indent+'Running simulation script')

        s = 'numerical_cg_transient.py'
        res = anuga.run_anuga_script(s,args=args)
        
        # Test that script runs ok
        assert res == 0


        if verbose:
            print(indent+'Testing accuracy')


        import anuga.utilities.plot_utils as util
        from analytical_cg_transient import analytical_sol

        p_st = util.get_output('carrier_greenspan.sww')
        p2_st=util.get_centroids(p_st)

        v = p2_st.y[10]
        v2=(p2_st.y==v)

        x_n = p2_st.x[v2]


        ids = [1, 15, 30]
        use_absolute = [False, True, True]
        n = len(ids)

        W, UH, Z, H, U = [], [], [], [], []
        W_n, U_n, UH_n = [], [], []

        #Dimensional parameters
        L   = 5e4         # Length of channel (m)
        h0  = 5e2         # Height at origin when the water is still
        g   = 9.81        # Gravity

        for i, id in enumerate(ids):
            w, z, u = analytical_sol(x_n/L, p2_st.time[id]*sqrt(g*h0)/L)

            w = w*h0
            z = z*h0
            u = u*sqrt(g*h0)
            h = w-z
            uh = u*h

            W.append(w)
            Z.append(z)
            H.append(h)
            U.append(u)
            UH.append(uh)


            W_n.append(p2_st.stage[id,v2])
            UH_n.append(p2_st.xmom[id,v2])
            U_n.append(p2_st.xvel[id,v2])

            #print id, numpy.max(H[i]), numpy.min(H[i])
            #print id, numpy.max(U[i]), numpy.min(U[i])
            #print id, numpy.max(U_n[i]), numpy.min(U_n[i])

        #Test stages
        # Calculate L^1 error at times
        
        ew = numpy.zeros(n)
        for i, id in enumerate(ids):
            ew[i] = numpy.sum(numpy.abs(W_n[i]-W[i]))/numpy.sum(numpy.abs(W[i]))
 
        print() 
        print(indent+'L^1 Errors in stage: ', ew)



        #Test xmomenta
        # Calculate L^1 error at times
        euh = numpy.zeros(n)
        for i, id in enumerate(ids):
            euh[i] = numpy.sum(numpy.abs(UH_n[i]-U[i]*H[i]))/numpy.sum(numpy.abs(UH_n[i]))

        
        print(indent+'L^1 Errors in xmomentum: ', euh)    

        #Test xvelocity
        # Calculate L^1 error at times

        eu = numpy.zeros(n)
        for i, id in enumerate(ids):
            if use_absolute[i]:
                eu[i] = numpy.sum(numpy.abs(U_n[i]-U[i]))/len(U[i])
            else:
                eu[i] = numpy.sum(numpy.abs(U_n[i]-U[i]))/numpy.sum(numpy.abs(U[i]))

        print(indent+'L^1 Errors in xvelocity: ', eu)

        for i, id in enumerate(ids):
            assert ew[i] < 0.01,  'L^1 error %g greater than 1 percent'% ew[i]

        for i, id in enumerate(ids):
            assert euh[i] < 0.025,  'L^1 error %g greater than 2.5 percent'% euh[i]

        for i, id in enumerate(ids):
            assert eu[i] < 0.025,  'L^1 error %g greater than 2.5 percent'% eu[i]


#-------------------------------------------------------------
if __name__ == '__main__':
    suite = unittest.makeSuite(Test_results, 'test')
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
