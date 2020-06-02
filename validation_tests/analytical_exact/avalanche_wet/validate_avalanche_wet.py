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

    def test_avalanche_wet(self):
    

        if verbose:
            print()
            print(indent+'Running simulation script')

        
        s = 'numerical_avalanche_wet.py'
        res = anuga.run_anuga_script(s, args=args)


        # Test that script runs ok
        assert res == 0


        if verbose:
            print(indent+'Testing accuracy')
            

        import anuga.utilities.plot_utils as util
        from analytical_avalanche_wet import analytical_sol

        p_st = util.get_output('avalanche.sww')
        p2_st=util.get_centroids(p_st)

        v = p2_st.y[20]
        v2=(p2_st.y==v)

        x_n = p2_st.x[v2]

        u0,h0,w0,z0,p0 = analytical_sol(x_n, p2_st.time[0])
        u20,h20,w20,z20,p20 = analytical_sol(x_n, p2_st.time[20])
        u40,h40,w40,z40,p40 = analytical_sol(x_n, p2_st.time[40])



        w0_n  = p2_st.stage[0,v2]
        w20_n = p2_st.stage[20,v2]
        w40_n = p2_st.stage[40,v2]

        z_n = p2_st.elev[v2]

        uh0_n  = p2_st.xmom[0,v2]
        uh20_n = p2_st.xmom[20,v2]
        uh40_n = p2_st.xmom[40,v2]

        u0_n  = p2_st.xvel[0,v2]
        u20_n = p2_st.xvel[20,v2]
        u40_n = p2_st.xvel[40,v2]


        #Test stages
        # Calculate L^1 error at times corrsponding to slices 20, 40
        eh20 = numpy.sum(numpy.abs(w20_n-w20))/numpy.sum(numpy.abs(w20))
        eh40 = numpy.sum(numpy.abs(w40_n-w40))/numpy.sum(numpy.abs(w40))


        print() 
        print(indent+'Errors in stage: ', eh20, eh40)




        #Test xmomenta
        # Calculate L^1 error at times corrsponding to slices 20, 40
        euh20 = numpy.sum(numpy.abs(uh20_n-u20*h20))/numpy.sum(numpy.abs(u20*h20))
        euh40 = numpy.sum(numpy.abs(uh40_n-u40*h40))/numpy.sum(numpy.abs(u40*h40))


        print(indent+'Errors in xmomentum: ', euh20, euh40)

        #Test xvelocity
        # Calculate L^1 error at times corrsponding to slices 20, 40
        eu20 = numpy.sum(numpy.abs(u20_n-u20))/numpy.sum(numpy.abs(u20))
        eu40 = numpy.sum(numpy.abs(u40_n-u40))/numpy.sum(numpy.abs(u40))

        print(indent+'Errors in xvelocity: ', eu20, eu40)

        assert eh20 < 0.01,  'L^1 error %g greater than 1 percent'% eh20
        assert eh40 < 0.01,  'L^1 error %g greater than 1 percent'% eh40
        
        assert euh20 < 0.025,  'L^1 error %g greater than 2.5 percent'% euh20
        assert euh40 < 0.025,  'L^1 error %g greater than 2.5 percent'% euh40
       
        assert eu20 < 0.02,  'L^1 error %g greater than 2 percent'% eu20
        assert eu40 < 0.01,  'L^1 error %g greater than 2 percent'% eu40


#-------------------------------------------------------------
if __name__ == '__main__':
    suite = unittest.makeSuite(Test_results, 'test')
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
