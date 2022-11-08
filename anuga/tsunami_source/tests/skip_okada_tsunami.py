"""
Test the 'okada_tsunami' routine

We first check

"""

import unittest, os
import warnings
from os import sep

import anuga
import anuga.tsunami_source.okada_tsunami as okada_tsunami
import numpy
from anuga.utilities.system_tools import get_pathname_from_package

class Test_okada_tsunami(unittest.TestCase):
    def setUp(self):
        pass
   
    def tearDown(self):
        pass

    def test_reproduce_okadas_table(self):
        # Here we check that we can compute
        # the answers provided in Okada's Table 2
        # These refer to earthquakes with strike = 0 and origin 0,0
        for j in range(2,3):
            #print 'Testing Case ', j

            if(j==2):
                ##-------
                # Case 2:
                ##-------
                okada_x_origin=0.
                okada_y_origin=0
                # NOTE: In Okada's system, x = 2, y=3. However, 
                # this must be rotated backwards by 90 degrees
                # to convert to physical coordinates as used by
                # the routine. These physical coordinates correspond
                # to standard descriptions of earthquakes where
                # x,y = lon,lat
                x = -3. # =2 in okada's frame
                y = 2.  # =3 in okadas frame
                d = 4.
                strike=0.
                dip = 70.
                L = 3.
                W = 2.
                slip=1.0

                # uz for [strike, normal] faults
                okada_values = [-2.747E-3, -3.564E-2]
            elif(j==3):
                ##-------
                # Case 3:
                ##-------
                okada_x_origin=0.
                okada_y_origin=0
                x = 0.
                y = 0.
                d = 4.
                strike=0.
                dip = 90.
                L = 3.
                W = 2.
                slip=1.0

                # uz for [strike, normal] faults
                okada_values = [0., 0.]
           
            #import pdb
            #pdb.set_trace() 
            #Compute slip surface centroid
            x_cent, y_cent, d_cent = okada_tsunami.okada_origin_2_slip_centroid([okada_x_origin*1000., okada_y_origin*1000.],\
                                                    d, L, W, strike, dip)
            #y_cent = (okada_y_origin + (W/2.)*numpy.cos(dip/180.*numpy.pi))*1000.
            #x_cent = (okada_x_origin + (L/2.))*1000.
            #d_cent = d - (W/2.)*numpy.sin(dip/180.*numpy.pi)

            #print 'x_cent, y_cent, d_cent = ', x_cent, y_cent, d_cent
            
            x_wanted=x*1000. # Desired values of x, y in m
            y_wanted=y*1000.

            # Set up earthquakes, and check that the vertical deformation is the same as in okada
            for i, rake in enumerate([0., 90.]):
                dis1=slip*numpy.cos(rake/180.*numpy.pi)
                dis2=slip*numpy.sin(rake/180.*numpy.pi)
                dis3=0.
                #print 'dis1, dis2, dis3 = ', dis1, dis2, dis3
                
                my_source=numpy.array([x_cent, y_cent, d_cent,strike, dip, L, W, dis1, dis2, dis3])
                my_source=my_source.reshape((1,10))
                
                tsu_funct = okada_tsunami.earthquake_source(my_source, verbose=False)
                uz = tsu_funct(numpy.array([x_wanted]), numpy.array([y_wanted]))
              
                # Compute both relative and absolute versions of the error
                reltol = abs((uz - okada_values[i])/uz)
                abstol = abs(uz-okada_values[i])
                assert ((reltol<1.0e-03)|(abstol<1.0e-06))

    def test_reproduce_okadas_table_rotated(self):
        # This test is like 'reproduce_okadas_table' except we rotate the
        # earthquake first -- thus checking that rotation is done correctly in
        # the code   

        # Loop over a range of rotations
        for rotation in [0., 30., 90., 150., 210., 325.]:
            # Test cases 2 - 3 in Okada's table
            for j in range(2,4):
                #print 'Testing Case ', j
                if(j==2):
                    ##-------
                    # Case 2: Parameters from table 2
                    ##-------
                    okada_x_origin=0.
                    okada_y_origin=0
                    d = 4.
                    strike=rotation
                    dip = 70.
                    L = 3.
                    W = 2.
                    slip=1.0

                    # The desired x and y output location must be rotated too
                    x, y = okada_tsunami.rotate_coordinates([-3., 2.], -strike)

                    # uz for [strike, normal] faults
                    okada_values = [-2.747E-3, -3.564E-2]
                elif(j==3):
                    ##-------
                    # Case 3: Parameters from table 2
                    ##-------
                    okada_x_origin=0.
                    okada_y_origin=0
                    x = 0.
                    y = 0.
                    d = 4.
                    strike=rotation
                    dip = 90.
                    L = 3.
                    W = 2.
                    slip=1.0

                    # uz for [strike, normal] faults
                    okada_values = [0., 0.]
                 
                #Convert to the notation in tsunami_okada
                x_cent, y_cent, d_cent = okada_tsunami.okada_origin_2_slip_centroid([okada_x_origin*1000., okada_y_origin*1000.],\
                                                    d, L, W, strike, dip)
                
                x_wanted, y_wanted =[x*1000., y*1000.]
                

                # Set up earthquakes, and check that the vertical deformation is the same as in okada
                for i, rake in enumerate([0., 90.]):
                    dis1=slip*numpy.cos(rake/180.*numpy.pi)
                    dis2=slip*numpy.sin(rake/180.*numpy.pi)
                    dis3=0.
                    
                    my_source=numpy.array([x_cent, y_cent, d_cent,strike, dip, L, W, dis1, dis2, dis3])
                    my_source=my_source.reshape((1,10))
                    
                    tsu_funct = okada_tsunami.earthquake_source(my_source, verbose=False)
                    uz = tsu_funct(numpy.array([x_wanted]), numpy.array([y_wanted]))
                  
                    # Compute both relative and absolute versions of the error
                    reltol = abs((uz - okada_values[i])/uz)
                    abstol = abs(uz-okada_values[i])
                    assert ((reltol<1.0e-03)|(abstol<1.0e-06))

    def test_reproduce_okadas_table_rotated_and_offset(self):
        # This is like 'reproduce_okadas_table_rotated', except we also loop
        # over a range of earthquake origins and strikes. Hence, we check that
        # the translation/rotation is done correctly
        #
        # Actually this test subsumes the first two. However, it is useful to 
        # have separately, since separate tests will help isolate any problems

        # Define sets of origins for okada function (in physical space)
        okada_origins=numpy.array([[0., 0.], [1., 2.], [2., -1], [-3., -10.]])

        # Loop over all origins
        for k in range(okada_origins.shape[0]):
            okada_x_origin, okada_y_origin=okada_origins[k,:]
            # Loop over a range of rotations
            for rotation in [0., 30., 90., 150., 210., 325.]:
                # Loop over okada's test cases 2 and 3
                for j in range(2,4):
                    #print 'Testing Case ', j, ' Rotation = ', rotation, ' origin = ', [okada_x_origin, okada_y_origin]

                    if(j==2):
                        ##-------
                        # Case 2:
                        ##-------
                        #okada_x_origin=-1.
                        #okada_y_origin=10.
                        d = 4.
                        strike=rotation
                        dip = 70.
                        L = 3.
                        W = 2.
                        slip=1.0

                        # The desired x and y output location must be rotated too
                        x, y = okada_tsunami.rotate_coordinates([-3., 2.], -strike)
                        x = x+okada_x_origin
                        y = y+okada_y_origin

                        # uz for [strike, normal] faults
                        okada_values = [-2.747E-3, -3.564E-2]
                    elif(j==3):
                        ##-------
                        # Case 3:
                        ##-------
                        #okada_x_origin=0.
                        #okada_y_origin=0
                        d = 4.
                        strike=rotation
                        dip = 90.
                        L = 3.
                        W = 2.
                        slip=1.0
                        
                        x = 0. + okada_x_origin
                        y = 0. + okada_y_origin

                        # uz for [strike, normal] faults
                        okada_values = [0., 0.]
                     
                    # Compute new centroid coordinates
                    x_cent, y_cent, d_cent = okada_tsunami.okada_origin_2_slip_centroid([okada_x_origin*1000., okada_y_origin*1000.],\
                                                    d, L, W, strike, dip)
                    #print '##', x_cent, y_cent, d_cent
                    #assert 0==1

                    x_wanted, y_wanted =[x*1000., y*1000.]
                    

                    # Set up earthquakes, and check that the vertical deformation is the same as in okada
                    for i, rake in enumerate([0., 90.]):
                        dis1=slip*numpy.cos(rake/180.*numpy.pi)
                        dis2=slip*numpy.sin(rake/180.*numpy.pi)
                        dis3=0.
                        
                        my_source=numpy.array([x_cent, y_cent, d_cent,strike, dip, L, W, dis1, dis2, dis3])
                        my_source=my_source.reshape((1,10))
                        
                        tsu_funct = okada_tsunami.earthquake_source(my_source, verbose=False)
                        uz = tsu_funct(numpy.array([x_wanted]), numpy.array([y_wanted]))
                      
                        # Compute both relative and absolute versions of the error
                        reltol = abs((uz - okada_values[i])/uz)
                        abstol = abs(uz-okada_values[i])
                        assert ((reltol<1.0e-03)|(abstol<1.0e-06)), 'Okada_tsunami error for eq source: ' + str(my_source)
                        #print 'PASS'

#
    def test_superimpose_two_subfaults(self):
        # Here we generate a slip surface with 2 subfaults, and check that
        # their result adds to the sum of the original 2. This uses values from
        # Okada's table 2.

        # Point at which we calculate output
        x = -3. # x = 2 in Okada's coordinates -- rotate anticlockwise by 90 deg for physical space
        y = 2.  # y = 3 in Okada's coordinates

        okada_values = [-2.747E-3, -3.564E-2] # Values of slip for strike, normal fault respectively

        # Compute functions for Okada's subfaults 2 and 3 (from table 2 in the paper)
        # Translate case 3 to origin 2,3. Then, we know the deformation at
        # (2,3) should be the sum of the values provided by Okada.
        for j in range(2,4):
            if(j==2):
                ##-------
                # Case 2:
                ##-------
                okada_x_origin=0.
                okada_y_origin=0
                d = 4.
                strike=0.
                dip = 70.
                L = 3.
                W = 2.
                slip=1.0
                rake=0.0
                # uz for [strike, normal] faults at x,y
            elif(j==3):
                ##-------
                # Case 3:
                ##-------
                okada_x_origin=x
                okada_y_origin=y
                d = 4.
                strike=0.
                dip = 90.
                L = 3.
                W = 2.
                slip=1.0
                rake=0.0

                # uz for [strike, normal] faults
                #okada_values = [0., 0.] # at x,y
            
            #Convert to the notation in tsunami_okada
            #d_cent = d - (W/2.)*numpy.sin(dip/180.*numpy.pi) # Centroid depth of fault plain
            #x_cent=(okada_x_origin+ L/2.)*1000. # Centroid x in m
            #y_cent=(okada_y_origin + (W/2.)*numpy.cos(dip/180.*numpy.pi))*1000. # Centroid y in m
            x_cent, y_cent, d_cent = okada_tsunami.okada_origin_2_slip_centroid([okada_x_origin*1000., okada_y_origin*1000.],\
                                                    d, L, W, strike, dip)
            
            x_wanted=x*1000. # Desired values of x, y in m
            y_wanted=y*1000.

            # Set up earthquakes, and check that the vertical deformation is the same as in okada
            dis1=slip*numpy.cos(rake/180.*numpy.pi)
            dis2=slip*numpy.sin(rake/180.*numpy.pi)
            dis3=0.
            # Set up the 'my_source' array, or append to it
            if(j==2):
                my_source=numpy.array([x_cent, y_cent, d_cent,strike, dip, L, W, dis1, dis2, dis3])
                my_source=my_source.reshape((1,10))
            elif(j==3):
                my_source2=numpy.array([x_cent, y_cent, d_cent,strike, dip, L, W, dis1, dis2, dis3]) 
                my_source2=my_source2.reshape((1,10))

                mysource=numpy.concatenate((my_source, my_source2))
            else:
                raise Exception('j is ' + str(j) + ', it should not take this value')
            
        # Tsunami function with     
        tsu_funct = okada_tsunami.earthquake_source(my_source, verbose=False)
        uz = tsu_funct(numpy.array([x_wanted]), numpy.array([y_wanted]))
      
        # Compute both relative and absolute versions of the error
        reltol = abs((uz - okada_values[0])/uz)
        abstol = abs(uz-okada_values[0])
        assert ((reltol<1.0e-03)|(abstol<1.0e-06))
        #print 'PASS'

    def test_against_octave_code(self):
        # This test runs a test case from 
        # DENYS DUTYKH , DIMITRIOS MITSOTAKIS, XAVIER GARDEIL, AND FREDERIC DIAS 
        # ON THE USE OF THE FINITE FAULT SOLUTION FOR TSUNAMI
        # GENERATION PROBLEMS 
        # The copy I have seems to be a pre-print, not sure of the journal
        # We compare the results of our code with another code written purely in matlab
        # This code is by 'Francois Beauducel', and was obtained from matlab central (July 2012)
        # See: http://www.mathworks.com/matlabcentral/fileexchange/25982
        slip=2.5
        rake=95.
        dis1=slip*numpy.cos(rake/180.*numpy.pi)
        dis2=slip*numpy.sin(rake/180.*numpy.pi)
        strike=289.
        dip=10.
        length=80.9
        width=40.0
        focal_depth=20.
        depth=focal_depth + width/2.0*numpy.sin(dip/180.*numpy.pi) # T
        
        source=numpy.array([0., 0., depth, strike,dip , length, width, dis1, dis2,0.0])

        tsunami_fun=okada_tsunami.earthquake_source(source=source, verbose=False)

        # Make a grid
        mygrid=numpy.lib.index_tricks.nd_grid()
        grid_width=400000
        n=101
        grid2=mygrid[0:grid_width:(n*1j), 0:grid_width:(n*1j)]
        
        # Centre grid
        x=(grid2[0,:,:].reshape((n*n)) -grid_width/2.)
        y=(grid2[1,:,:].reshape((n*n)) -grid_width/2.)
        
        eq_source=tsunami_fun(x,y)
        
        ## Now read the same event from an octave code, which is completely
        ## independent of this one (i.e. they don't call okada's fortran)
        path=get_pathname_from_package('anuga.tsunami_source')
        octave=numpy.genfromtxt(path+sep+'tests'+sep+'okada_tsunami_octave_95.txt')
        octave_asvec=numpy.transpose(octave).reshape((1,101*101))

        # Estimate the differences between the 2 codes
        assert (abs(octave_asvec-eq_source)).max() < 1.0e-04
        
        ## Code used to run the scenario in octave inserted below for reference
        ## 
        ##
        ## % Here we implement the scenario described in 
        ## % DENYS DUTYKH , DIMITRIOS MITSOTAKIS, XAVIER GARDEIL, AND FREDERIC DIAS 
        ## % ON THE USE OF THE FINITE FAULT SOLUTION FOR TSUNAMI
        ## % GENERATION PROBLEMS 
        ## % The copy I have seems to be a pre-print, not sure of the journal
        ## % However, we get exactly the same answers
        ## 
        ## 
        ## strike=289.
        ## dip=10.
        ## length=80.9
        ## width=40.
        ## depth_base=20. + (width)*sin(dip/180*pi)
        ## rakes=[95.] %[0., 30., 90.]
        ## slip=2.5
        ## open=0.
        ## 
        ## n_pts=101
        ## 
        ## %Points to output
        ## x=2.
        ## y=3.
        ## 
        ## 
        ## % Compute depth as used in this code -- which are smaller than in Okada
        ## depth=depth_base-(width/2.0)*sin(dip/180.*pi)
        ## 
        ## % Compute x,y relative to coordinate system used in this code (which is centroid based)
        ## local_x=x-(length/2.0)
        ## local_y=y-(width/2.0)*cos(dip/180.*pi)
        ## 
        ## % 43-44, 49-50
        ## [u1,u2,u3] = okada85(local_x,local_y,depth,strike,dip,length,width,rakes(1),slip,open);
        ## 
        ## % Compute surface
        ## [E,N] = meshgrid(linspace(-200,200,n_pts), linspace(-200,200,n_pts));
        ## %	[uE,uN,uZ,uZE,uZN,uNN,uNE,uEN,uEE] = OKADA85(E,N,...
        ## %	   DEPTH,STRIKE,DIP,LENGTH,WIDTH,RAKE,SLIP,OPEN,NU)
        ## for i = 1:1
        ##     % call okada function
        ##     [uE,uN,uZ] = okada85(E,N,depth,strike,dip,length,width,rakes(i),slip,open);
        ##     %surf(E, N, uZ)
        ##     pic_filename=strcat('okada_octave_', num2str(rakes(i)), '.txt')
        ##     save('-ascii', file=pic_filename, 'uZ')
        ## end

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_okada_tsunami, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
