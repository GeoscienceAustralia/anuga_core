"""Common filenames and locations for topographic data, meshes and outputs.
"""
from os import sep, environ, getenv, getcwd
from os.path import expanduser
import sys
from time import localtime, strftime, gmtime
import anuga.utilities
from os.path import join

verbose = True
use_cache = False


###############################
# Domain definitions
###############################


# bounding polygon for study area
bounding_polygon = anuga.read_polygon('extent.csv')
print('Area of bounding polygon (km^2)', anuga.polygon_area(bounding_polygon)/1.0e+06)



houses_dir = 'houses'
poly_house000=anuga.read_polygon(join(houses_dir,'house000.csv'))
poly_house001=anuga.read_polygon(join(houses_dir,'house001.csv'))
poly_house002=anuga.read_polygon(join(houses_dir,'house002.csv'))
poly_house003=anuga.read_polygon(join(houses_dir,'house003.csv'))
poly_house004=anuga.read_polygon(join(houses_dir,'house004.csv'))
poly_house005=anuga.read_polygon(join(houses_dir,'house005.csv'))
poly_house006=anuga.read_polygon(join(houses_dir,'house006.csv'))
poly_house007=anuga.read_polygon(join(houses_dir,'house007.csv'))
poly_house008=anuga.read_polygon(join(houses_dir,'house008.csv'))
poly_house009=anuga.read_polygon(join(houses_dir,'house009.csv'))
poly_house010=anuga.read_polygon(join(houses_dir,'house010.csv'))
poly_house011=anuga.read_polygon(join(houses_dir,'house011.csv'))
poly_house012=anuga.read_polygon(join(houses_dir,'house012.csv'))
poly_house013=anuga.read_polygon(join(houses_dir,'house013.csv'))
poly_house014=anuga.read_polygon(join(houses_dir,'house014.csv'))
poly_house015=anuga.read_polygon(join(houses_dir,'house015.csv'))
poly_house016=anuga.read_polygon(join(houses_dir,'house016.csv'))
poly_house017=anuga.read_polygon(join(houses_dir,'house017.csv'))
poly_house018=anuga.read_polygon(join(houses_dir,'house018.csv'))
poly_house019=anuga.read_polygon(join(houses_dir,'house019.csv'))
poly_house020=anuga.read_polygon(join(houses_dir,'house020.csv'))
poly_house021=anuga.read_polygon(join(houses_dir,'house021.csv'))
poly_house022=anuga.read_polygon(join(houses_dir,'house022.csv'))
poly_house023=anuga.read_polygon(join(houses_dir,'house023.csv'))
poly_house024=anuga.read_polygon(join(houses_dir,'house024.csv'))
poly_house025=anuga.read_polygon(join(houses_dir,'house025.csv'))
poly_house026=anuga.read_polygon(join(houses_dir,'house026.csv'))
poly_house027=anuga.read_polygon(join(houses_dir,'house027.csv'))
poly_house028=anuga.read_polygon(join(houses_dir,'house028.csv'))
poly_house029=anuga.read_polygon(join(houses_dir,'house029.csv'))
#poly_house030=anuga.read_polygon(join(houses_dir,'house030.csv'))
poly_house031=anuga.read_polygon(join(houses_dir,'house031.csv'))
poly_house032_033=anuga.read_polygon(join(houses_dir,'house032_033.csv'))
poly_house032=anuga.read_polygon(join(houses_dir,'house032.csv'))
poly_house033=anuga.read_polygon(join(houses_dir,'house033.csv'))
poly_house034=anuga.read_polygon(join(houses_dir,'house034.csv'))
poly_house035=anuga.read_polygon(join(houses_dir,'house035.csv'))
poly_house036=anuga.read_polygon(join(houses_dir,'house036.csv'))
poly_house037=anuga.read_polygon(join(houses_dir,'house037.csv'))
poly_house038=anuga.read_polygon(join(houses_dir,'house038.csv'))
poly_house039=anuga.read_polygon(join(houses_dir,'house039.csv'))
poly_house040=anuga.read_polygon(join(houses_dir,'house040.csv'))
poly_house041=anuga.read_polygon(join(houses_dir,'house041.csv'))
poly_house042=anuga.read_polygon(join(houses_dir,'house042.csv'))
poly_house043=anuga.read_polygon(join(houses_dir,'house043.csv'))
poly_house044=anuga.read_polygon(join(houses_dir,'house044.csv'))
poly_house045=anuga.read_polygon(join(houses_dir,'house045.csv'))
poly_house046=anuga.read_polygon(join(houses_dir,'house046.csv'))
poly_house047=anuga.read_polygon(join(houses_dir,'house047.csv'))
poly_house048=anuga.read_polygon(join(houses_dir,'house048.csv'))
poly_house049=anuga.read_polygon(join(houses_dir,'house049.csv'))
poly_house050=anuga.read_polygon(join(houses_dir,'house050.csv'))
poly_house051=anuga.read_polygon(join(houses_dir,'house051.csv'))
poly_house052=anuga.read_polygon(join(houses_dir,'house052.csv'))
poly_house053=anuga.read_polygon(join(houses_dir,'house053.csv'))
poly_house054=anuga.read_polygon(join(houses_dir,'house054.csv'))
poly_house055=anuga.read_polygon(join(houses_dir,'house055.csv'))
poly_house056=anuga.read_polygon(join(houses_dir,'house056.csv'))
poly_house057=anuga.read_polygon(join(houses_dir,'house057.csv'))
poly_house058=anuga.read_polygon(join(houses_dir,'house058.csv'))

holes=[
poly_house000,
poly_house001,
poly_house002,
poly_house003,
poly_house004,
poly_house005,
poly_house006,
poly_house007,
poly_house008,
poly_house009,
poly_house010,
poly_house011,
poly_house012,
poly_house013,
poly_house014,
poly_house015,
poly_house016,
poly_house017,
poly_house018,
poly_house019,
poly_house020,
poly_house021,
poly_house022,
poly_house023,
poly_house024,
poly_house025,
poly_house026,
poly_house027,
poly_house028,
poly_house029,
poly_house031,
poly_house032_033,
#poly_house032,
#poly_house033, # house032 and house033 have a common edge.
poly_house034,
poly_house035,
poly_house036,
poly_house037,
poly_house038,
poly_house039,
poly_house040,
poly_house041,
poly_house042,
poly_house043,
poly_house044,
poly_house045,
poly_house046,
poly_house047,
poly_house048,
poly_house049,
poly_house050,
poly_house051,
poly_house052,
poly_house053,
poly_house054,
poly_house055,
poly_house056,
poly_house057,
poly_house058
]


###############################
# Interior region definitions
###############################

# interior polygons
poly_merewether = anuga.read_polygon('merewether.csv')

gauge_filename = 'gauges.csv'

# exporting asc grid
eastingmin = 382249.79174463
eastingmax = 382570.7715
northingmin = 6354265.4322858
northingmax = 6354681.40599876

# zone
zone = '56S'

