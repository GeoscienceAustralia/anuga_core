"""
Download data from anuga_case_studies_data repository at the given URL
location
"""

import os

URL = 'https://anuga.anu.edu.au/svn/anuga/trunk/anuga_case_studies_data/towradgi/DEM'

CMD = 'wget -r  --cut-dirs=5 --no-parent --reject "index.html" -nH %s'% URL

print CMD

os.system(CMD)

#wget -r  --cut-dirs=5 --no-parent --reject "index.html" -nH https://anuga.anu.edu.au/svn/anuga/trunk/anuga_case_studies_data/towradgi/DEM
