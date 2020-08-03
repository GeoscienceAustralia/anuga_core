"""Automatic verification of ANUGA flows.
See functions exercised by this wrapper for more details
"""

import unittest
import os
import numpy
import anuga

args = anuga.get_args()


indent = anuga.indent

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

    def test_dam_break_dry(self):
    

        if verbose:
            print()
            print(indent+'Running simulation script')

        s = 'numerical_dam_break_dry.py'
        res = anuga.run_anuga_script(s,args=args)

        # Test that script runs ok
        assert res == 0


        if verbose:
            print(indent+'Testing accuracy')
            
        import anuga.utilities.plot_utils as util
        import analytical_dam_break_dry as analytic

        p_st = util.get_output('dam_break.sww')
        p2_st=util.get_centroids(p_st)

        v = p2_st.y[10]
        v2=(p2_st.y==v)


        h0 = 1.0
        h1 = 10.0

        # calculate analytic values at various time slices
        h10,u10 = analytic.vec_dam_break(p2_st.x[v2], p2_st.time[10], h0=h0, h1=h1)
        h50,u50 = analytic.vec_dam_break(p2_st.x[v2], p2_st.time[50], h0=h0, h1=h1)
        h100,u100 = analytic.vec_dam_break(p2_st.x[v2], p2_st.time[100], h0=h0, h1=h1)


        #Test stages
        # Calculate L^1 error at times corrsponding to slices 10, 50 and 100
        eh10 = numpy.sum(numpy.abs(p2_st.stage[10,v2]-h10))/numpy.sum(numpy.abs(h10))
        eh50 = numpy.sum(numpy.abs(p2_st.stage[50,v2]-h50))/numpy.sum(numpy.abs(h50))
        eh100 = numpy.sum(numpy.abs(p2_st.stage[100,v2]-h100))/numpy.sum(numpy.abs(h100))

        print() 
        print(indent+'Errors in stage: ', eh10, eh50, eh100)


        #Test xmomenta
        # Calculate L^1 error at times corrsponding to slices 10, 50 and 100
        euh10 = numpy.sum(numpy.abs(p2_st.xmom[10,v2]-u10*h10))/numpy.sum(numpy.abs(u10*h10))
        euh50 = numpy.sum(numpy.abs(p2_st.xmom[50,v2]-u50*h50))/numpy.sum(numpy.abs(u50*h50))
        euh100 = numpy.sum(numpy.abs(p2_st.xmom[100,v2]-u100*h100))/numpy.sum(numpy.abs(u100*h100))

        print(indent+'Errors in xmomentum: ', euh10, euh50, euh100)

        #Test xvelocity
        # Calculate L^1 error at times corrsponding to slices 10, 50 and 100
        eu10 = numpy.sum(numpy.abs(p2_st.xvel[10,v2]-u10))/numpy.sum(numpy.abs(u10))
        eu50 = numpy.sum(numpy.abs(p2_st.xvel[50,v2]-u50))/numpy.sum(numpy.abs(u50))
        eu100 = numpy.sum(numpy.abs(p2_st.xvel[100,v2]-u100))/numpy.sum(numpy.abs(u100))

        print(indent+'Errors in xvelocity: ', eu10, eu50, eu100)


        assert eh10 < 0.1,  'Relative L^1 error %g greater than 0.1'% eh10
        assert eh50 < 0.1,  'Relative L^1 error %g greater than 0.1'% eh50
        assert eh100 < 0.15, 'Relative L^1 error %g greater than 0.15'% eh100

        assert euh10 < 0.25,  'Relative L^1 error %g greater than 0.25'% euh10
        assert euh50 < 0.25,  'Relative L^1 error %g greater than 0.25'% euh50
        assert euh100 < 0.25, 'Relative L^1 error %g greater than 0.25'% euh100

        assert eu10 < 2.0,  'Relative L^1 error %g greater than 2.0'% eu10
        assert eu50 < 2.0,  'Relative L^1 error %g greater than 2.0'% eu50
        assert eu100 < 0.3, 'Relative L^1 error %g greater than 0.3'% eu100
 



#-------------------------------------------------------------
if __name__ == '__main__':
    suite = unittest.makeSuite(Test_results, 'test')
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
