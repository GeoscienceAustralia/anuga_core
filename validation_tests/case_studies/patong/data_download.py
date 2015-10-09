"""
Download data from anuga_case_studies_data repository at the given URL
location
"""

import os

#URL = 'https://anuga.anu.edu.au/svn/anuga/trunk/anuga_case_studies_data/patong/boundaries'

#CMD = 'wget -r  --cut-dirs=5 --no-parent --reject "index.html" -nH %s'% URL

URL =  'https://sourceforge.net/projects/anuga/files/validation_data/patong-1.0/data.tgz'

CMD = 'wget %s'% URL

print CMD

try:
    import wget
    wget.download(URL)
except:
    print 'wget failed. Perhaps need to install wget via "pip install wget"'
    import sys
    sys.exit()

CMD = 'tar vzxf data.tgz'

print CMD
os.system(CMD)

CMD = 'mv data/thailand/patong_tsunami_scenario/anuga/* .'

print CMD
os.system(CMD)

CMD = 'rm -r data.tgz data'

print CMD
os.system(CMD)
