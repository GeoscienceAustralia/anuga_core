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

    def test_avalanche_dry(self):
    

        if verbose:
            print()
            print(indent+'Running simulation script')

        # Run basic script (can be parallel if -np used in call
        # to this script
        s = 'numerical_avalanche_dry.py'
        res = anuga.run_anuga_script(s,args=args)


        # Test that script runs ok
        assert res == 0


        if verbose:
            print(indent+'Testing accuracy')



        import anuga.utilities.plot_utils as util
        from analytical_avalanche_dry import analytical_sol

        p_st = util.get_output('avalanche.sww')
        p2_st=util.get_centroids(p_st)

        v = p2_st.y[10]
        v2=(p2_st.y==v)

        x_n = p2_st.x[v2]



        u0,h0,w0,z0,p0 = analytical_sol(x_n, p2_st.time[0])
        u10,h10,w10,z10,p10 = analytical_sol(x_n, p2_st.time[10])
        u30,h30,w30,z30,p30 = analytical_sol(x_n, p2_st.time[30])


        w0_n  = p2_st.stage[0,v2]
        w10_n = p2_st.stage[10,v2]
        w30_n = p2_st.stage[30,v2]

        z_n = p2_st.elev[v2]

        uh0_n  = p2_st.xmom[0,v2]
        uh10_n = p2_st.xmom[10,v2]
        uh30_n = p2_st.xmom[30,v2]

        u0_n  = p2_st.xvel[0,v2]
        u10_n = p2_st.xvel[10,v2]
        u30_n = p2_st.xvel[30,v2]


        #Test stages
        # Calculate L^1 error at times corrsponding to slices 10, 50 and 100
        eh10 = numpy.sum(numpy.abs(w10_n-w10))/numpy.sum(numpy.abs(w10))
        eh30 = numpy.sum(numpy.abs(w30_n-w30))/numpy.sum(numpy.abs(w30))


        print() 
        print(indent+'Errors in stage: ', eh10, eh30)



        # Test xmomenta
        # Calculate L^1 error at times corrsponding to slices 10, 50 and 100
        euh10 = numpy.sum(numpy.abs(uh10_n-u10*h10))/numpy.sum(numpy.abs(u10*h10))
        euh30 = numpy.sum(numpy.abs(uh30_n-u30*h30))/numpy.sum(numpy.abs(u30*h30))
    

        print(indent+'Errors in xmomentum: ', euh10, euh30)

        #Test xvelocity
        # Calculate L^1 error at times corrsponding to slices 10, 50 and 100
        eu10 = numpy.sum(numpy.abs(u10_n-u10))/numpy.sum(numpy.abs(u10))
        eu30 = numpy.sum(numpy.abs(u30_n-u30))/numpy.sum(numpy.abs(u30))


        print(indent+'Errors in xvelocity: ', eu10, eu30)


        assert eh10 < 0.01,  'L^1 error %g greater than 1 percent'% eh10
        assert eh30 < 0.01,  'L^1 error %g greater than 1 percent'% eh30

        assert euh10 < 0.025,  'L^1 error %g greater than 2.5 percent'% euh10
        assert euh30 < 0.025,  'L^1 error %g greater than 2.5 percent'% euh30


        assert eu10 < 0.2,  'L^1 error %g greater than 20 percent'% eu10
        assert eu30 < 0.2,  'L^1 error %g greater than 20 percent'% eu30

 



#-------------------------------------------------------------
if __name__ == '__main__':
    suite = unittest.makeSuite(Test_results, 'test')
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
