"""

This little test script is used to test creating polygons at the end of a culvert inlet and outlet.
It needs to be called once only during an ANUGA run to create the required polygons for enquiring
the total energy at both ends of the culvert and transferring flow
May need to include ability to have a function that controls the blockage level of the culvert

"""

from anuga.geometry.polygon import plot_polygons
from culvert_polygons import create_culvert_polygons

# Culvert location
x0 = 5.0; y0 = 5.0  # One end
x1 = 15.0; y1 = 20.0  # Other end

# Culvert Size
culvert_width=2.4
culvert_height=1.2


# Call function
P = create_culvert_polygons(end_point0=[x0, y0],
                            end_point1=[x1, y1],
                            width=culvert_width,
                            height=culvert_height)



culv_polygon = [[x0,y0], [x1,y1]]
plot_polygons([culv_polygon,
               P['exchange_polygon0'],
               P['exchange_polygon1'],
               P['enquiry_polygon0'],
               P['enquiry_polygon1']],
              figname='culvert_polygon_example')


