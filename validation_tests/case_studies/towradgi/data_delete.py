"""
delete the files that were downloaded from the
anuga_case_studies_data repository
"""

import shutil

for dir in ['DEM_bridges',  'Forcing',  'Model',  'Validation']:
    shutil.rmtree(dir)


