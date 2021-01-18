"""
Download data from anuga_case_studies_data repository at the given URL
location
"""

import os

URL =  'https://sourceforge.net/projects/anuga/files/validation_data/patong-1.0/data.tgz'

CMD = 'wget %s'% URL
print (CMD)

try:
    import wget
    wget.download(URL, )
except:
    print ('wget failed. Perhaps need to install wget via "pip install wget"')
    import sys
    sys.exit()


CMD = 'tar zxf data.tgz'
print ()
print (CMD)

import tarfile
tar = tarfile.open('data.tgz')
tar.extractall()
tar.close()


CMD = 'mv data/thailand/patong_tsunami_scenario/anuga/* .'
print (CMD)

import shutil
shutil.move('data/thailand/patong_tsunami_scenario/anuga/boundaries','.')
shutil.move('data/thailand/patong_tsunami_scenario/anuga/gauges','.')
shutil.move('data/thailand/patong_tsunami_scenario/anuga/meshes','.')
shutil.move('data/thailand/patong_tsunami_scenario/anuga/outputs','.')
shutil.move('data/thailand/patong_tsunami_scenario/anuga/polygons','.')
shutil.move('data/thailand/patong_tsunami_scenario/anuga/topographies','.')


CMD = 'rm -r data.tgz data'
print (CMD)

shutil.rmtree('data')
os.remove('data.tgz')






