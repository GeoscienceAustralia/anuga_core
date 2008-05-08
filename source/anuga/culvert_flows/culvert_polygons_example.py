"""

This little test script is used to test creating polygons at the end of a culvert inlet and outlet.
It needs to be called once only during an ANUGA run to create the required polygons for enquiring
the total energy at both ends of the culvert and transferring flow
May need to include ability to have a function that controls the blockage level of the culvert

"""

from anuga.utilities.polygon import read_polygon, plot_polygons
from culvert_polygons import create_culvert_polygons

# Read these values or passed from some where
culvert_description= '2.4mx1.2m Box Culvert'

# Culvert location
x1= 5.0; y1= 5.0  # One end
x2=15.0; y2=10.0  # Other end

# Culvert Size
culvert_width=2.4
culvert_height=1.2


# Call function
P = create_culvert_polygons(center1=[x1, y1],
                            center2=[x2, y2],
                            width=culvert_width,
                            height=culvert_height)



culv_polygon = [[x1,y1], [x2,y2]]
plot_polygons([culv_polygon,
               P['exchange_polygon1'],
               P['exchange_polygon2'],
               P['enquiry_polygon1'],
               P['enquiry_polygon2']],
              figname='Culvert Check')


