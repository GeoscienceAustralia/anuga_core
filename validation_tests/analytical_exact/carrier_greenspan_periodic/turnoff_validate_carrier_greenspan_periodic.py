"""Automatic verification of ANUGA flows.
See functions exercised by this wrapper for more details
"""

import unittest
import os
import numpy
import anuga

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

    def test_carrier_greenspan_periodic(self):
    

        if verbose:
            print()
            print(indent+'Running simulation script')

        s = 'numerical_carrier_greenspan.py'
        res = anuga.run_anuga_script(s,args=args)

        # Test that script runs ok
        assert res == 0


        if verbose:
            print(indent+'Testing accuracy')

        import anuga.utilities.plot_utils as util
        from analytical_carrier_greenspan import analytic_cg

        p_st = util.get_output('carrier_greenspan.sww')
        p2_st=util.get_centroids(p_st)

        v = p2_st.y[10]
        v2=(p2_st.y==v)

        x_n = p2_st.x[v2]

        #ids = [288, 296, 304, 312, 320, 328]
        ids = [296, 304, 320, 328]

        # For the velocity and xmomentum calculations use
        # absolute error (as the exact solution is 0)
        use_absolute = [False, False, False, False]

        n = len(ids)

        W, P, Z, H, U = [], [], [], [], []

        W_n, U_n, UH_n = [], [], []

        for i, id in enumerate(ids):
            w0, p0, z0, h0, u0 = analytic_cg(x_n, p2_st.time[id], h0=5e2, L=5e4, a=1.0, Tp=900.0)

            W.append(w0)
            H.append(h0)
            U.append(u0)
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
            if use_absolute[i]:
                euh[i] = numpy.sum(numpy.abs(UH_n[i]-U[i]*H[i]))/len(U[i])
            else:
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
            assert eu[i] < 0.1,  'L^1 error %g greater than 10 percent'% eu[i]


#-------------------------------------------------------------
if __name__ == '__main__':
    suite = unittest.makeSuite(Test_results, 'test')
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
