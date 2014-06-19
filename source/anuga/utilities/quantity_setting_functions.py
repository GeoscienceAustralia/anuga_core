"""

Function which can be useful when setting quantities

"""
import copy

def make_nearestNeighbour_quantity_function(quantity_xyValueIn, domain, threshold_distance = 9.0e+100, background_value = 9.0e+100):
    """
    Function which makes another function, which can be used in set_quantity 

    Idea: For every point x,y in the domain, we want to set a quantity based on
          the 'nearest-neighbours' from quantity_xyValue (a 3 column array with
          x,y,quantity-value), UNLESS the distance from x,y to the nearest-neighbour is > threshold_distance.
          In the latter case, we want to set the quantity value to 'background_value'
            
          We need a function f(x,y) to do that. This routine makes the
          function, with the desired quantity_xyValue points, threshold_distance, and
          background_value

    INPUTS:
        quantity_xyValueIn -- A 3 column array with:
                              x,y, Value 
                            defining the points used to set the new quantity values
        domain -- The ANUGA domain
        threshold_distance -- Points greater than this distance from their nearest quantity_xyValue point are set to background_value
        background_value -- see 'threshold_distance'

    OUTPUTS: 
        A function f which can be passed to domain.set_quantity('myQuantity', f)
    """
    import scipy
    import scipy.interpolate


    if(len(quantity_xyValueIn.shape)>1):
        quantity_xyValue=copy.copy(quantity_xyValueIn) # Pointer, no copy
    else:
        # Treat the single-point case
        quantity_xyValue=copy.copy(quantity_xyValueIn.reshape((1,3)))
    # Make a function which gives us the ROW-INDEX of the nearest xy point in quantity_xyValue
    quantity_xy_interpolator=scipy.interpolate.NearestNDInterpolator(quantity_xyValue[:,0:2], scipy.arange(len(quantity_xyValue[:,2])))

    #
    # Make a function of x,y which we can pass to domain.set_quantity
    def quant_NN_fun(x,y):
        """
        Function to assign quantity from the nearest point in quantity_xyValue,
        UNLESS the point is more than 'threshold_distance' away from the nearest point,
        in which case the background friction value is used
       
        """
        import scipy
        import scipy.interpolate
        # Since ANUGA stores x,y internally in non-georeferenced coordinates, we
        # adjust them here
        xll=domain.geo_reference.xllcorner
        yll=domain.geo_reference.yllcorner
        z=scipy.zeros(shape=(len(x), 2))
        z[:,0]=x+xll
        z[:,1]=y+yll
        # This will hold the quantity values
        quantity_output=x*0. +background_value
        # Compute the index of the nearest-neighbour in quantity_xyValue
        q_index=quantity_xy_interpolator(z)
        # Next find indices with distance < threshold_distance
        dist_lt_thresh=( (z[:,0]-quantity_xyValue[q_index,0])**2 + (z[:,1]-quantity_xyValue[q_index,1])**2 < threshold_distance**2)
        dist_lt_thresh=dist_lt_thresh.nonzero()[0]
        quantity_output[dist_lt_thresh] = quantity_xyValue[q_index[dist_lt_thresh],2] 
        return quantity_output
    # Return the quantity function
    return quant_NN_fun

###################################################################################################

def composite_quantity_setting_function(fun_poly_pairs, domain):
    """
        Make a 'composite function' to set quantities -- applies different functions inside different polygon regions.
             
        fun_poly_pairs = [ [f0, p0], [f1, p1], ...] where fi is a function and pi is a polygon, or None, or 'All'
              
        IMPORTANT: When polygons overlap, the first elements of the list are given priority. 
                   The approach is:
                       First f0 is applied to all points in p0, and we record that these points have been 'set'
                       Next f1 is applied to all points in p1 which have not been 'set', and then we record those points as being 'set'
                       Next f2 is applied to all points in p2 which have not been 'set', and then we record those points as being 'set'
                       ... etc

        INPUT: 
              fun_poly_pairs = [ [f0, p0], [f1, p1], ...]

                  where fi(x,y) are functions returning quantity values at points, and 
                  pi are polygons where we want to use fi inside
                
                  If any pi = 'All', then we assume that ALL unset points are set
                    using the function. This CAN ONLY happen in the last [fi,pi] pair where pi is
                    not None (since fi will be applied to all remaining points -- so anything else is probably an input mistake)
    
                  If any pi = None, then that [fi,pi] pair is skipped

              domain = ANUGA domain object


        OUTPUT: A function F(x,y) which can be used e.g. to set the quantity
                domain.set_quantity('elevation', F)
    """
    import scipy
    from matplotlib import path

    def F(x,y):
        """
            This is the function we return
        """
        isSet=scipy.zeros(len(x)) # Record if each point has been set
        quantityVal=x*0 # F value
        lfp=len(fun_poly_pairs)
        if(lfp<=0):
            raise Exception, 'Must have at least 1 fun-poly-pair'

        # Make an array of 'transformed' spatial coordinates, for checking
        # polygon inclusion
        xll=domain.geo_reference.xllcorner
        yll=domain.geo_reference.yllcorner
        xy_array_trans=scipy.vstack([x+xll,y+yll]).transpose()

        # Test that none of the pi polygons [except perhaps the last] is 'All'
        for i in range(lfp-1):
            if (fun_poly_pairs[i][1]=='All'):
                # This is only ok if all the othe fun_poly_pairs are None
                remaining_fun_poly_pairs_are_None=[ fun_poly_pairs[j][1] is None for j in range(i+1,lfp)]
                if(not all(remaining_fun_poly_pairs_are_None)):
                    raise Exception, 'Can only have the last polygon = All'

        # Apply the fi inside the pi
        for i in range(lfp):
            fi = fun_poly_pairs[i][0] # The function
            pi = fun_poly_pairs[i][1] # The polygon
            # Get a 'true/false' flag for whether the point is in pi, and not set
            if(pi is 'All'):
                fInside=(1-isSet)
            elif(pi is None):
                continue
            else:
                # Try matplotlib -- make the path a matplotlib path object
                pi_path=path.Path(scipy.array(pi))
                fInside=xy_array_trans[:,0]*0.
                for j in range(len(fInside)):
                    fInside[j] = pi_path.contains_point(xy_array_trans[j,:])
                
                fInside=fInside*(1-isSet)
            #
            # Find indices of points we need to adjust 
            fInds=(fInside==1.).nonzero()[0]
            if(len(fInds)>0):
                quantityVal[fInds] = fi(x[fInds], y[fInds])
                isSet[fInds]=1

        if(not min(isSet)==1):
            raise Exception, 'Some points were not inside any polygon'

        return quantityVal

    return F

##############################################################################

def quantityRasterFun(domain, rasterFile):
    """
    Make a function whick takes x,y in ANUGA coordinates, and returns the values
    on a raster rasterFile
    
    This can be used to set a quantity, and takes care of the manual conversion from 
    ANUGA coordinates to spatial coordinates.

    INPUTS: domain = ANUGA domain
            rasterFile = Filename of the raster to extract point values from
   
    OUTPUT: Function which takes x,y in ANUGA coordinates, and outputs their
            corresponding raster values 
    """
    import scipy
    from anuga.utilities.spatialInputUtil import rasterValuesAtPoints
    def QFun(x,y):
        xll=domain.geo_reference.xllcorner
        yll=domain.geo_reference.yllcorner
        inDat=scipy.vstack([x+xll,y+yll]).transpose()
        return rasterValuesAtPoints(xy=inDat,rasterFile=rasterFile)
    #
    return QFun

#################################################################################

def elevation_from_Pt_Pol_Data_and_Raster(Pt_Pol_Data, elevation_raster, domain):
    """
        Function to make a function F(x,y) which returns the elevation on elevation_raster,
            except if x,y is inside the polygon associated with any element of Pt_Pol_Data,
            in which case a Pt_Pol_-specific nearest neighbour interpolator is used.

        INPUT:
            Pt_Pol_Data = a list with [ [ Polygon_0, Pt_XYZ_0], 
                                        [ Polygon_1, Pt_XYZ_1], 
                                        ... ]
            elevation_raster = A GDAL-compatible elevation raster
            domain = ANUGA domain
    """

    # Function to set elevation from raster 
    elevFun1=quantityRasterFun(domain, rasterFile=elevation_raster)

    # List of function/polygon pairs defining the Pt_Pol_ elevation data
    elevFunChanList=[]
    for i in range(len(Pt_Pol_Data)):
        elevFunChanList.append([
             make_nearestNeighbour_quantity_function(Pt_Pol_Data[i][1], domain), 
             Pt_Pol_Data[i][0]
                              ])

    #
    elevFun=composite_quantity_setting_function(elevFunChanList+[[elevFun1, 'All']], domain)

    return elevFun
