"""
Download data from anuga_case_studies_data repository at the given URL
location
"""

import os

URL =  'https://sourceforge.net/projects/anuga/files/validation_data/towradgi-1.0/data.tgz'


CMD = 'wget %s'% URL
print(CMD)

try:
    import wget
    wget.download(URL)
except:
    print('wget failed. Perhaps you need to install wget for python via "pip install wget"')
    import sys
    sys.exit()

CMD = 'tar zxf data.tgz'
print()
print(CMD)

import tarfile
tar = tarfile.open('data.tgz')
tar.extractall()
tar.close()


CMD = 'rm data.tgz'
print(CMD)
os.remove('data.tgz')
