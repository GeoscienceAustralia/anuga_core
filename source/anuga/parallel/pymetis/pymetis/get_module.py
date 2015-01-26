
import os
import os.path
import shutil
from os.path import join

os.chdir('build')


for file in os.listdir('.'):
    if os.path.isdir(file) and file.startswith('lib'):
        os.chdir(file)
        locallist = os.listdir('.')
        for lfile in locallist:
            if  lfile=='metis.so' or lfile=='metis.pyd':
                print 'Found %s in %s\n' % (lfile, file)
                shutil.copy2(lfile, join('..','..','..'))
            else:
                pass
        os.chdir('..')

os.chdir('..')

   
        
                

    
   
   
   
   
    
"""    
for X in $(ls -d lib.*)
  do
  cd ${X}
  if [ -f metis.so ]
      then
      cp metis.so ../../../ -pf
      echo 'Found metis.so in' ${X}
  elif [ -f metis.pyd ]
      then
      cp metis.pyd ../../../ -pf
      echo 'Found metis.pyd in' ${X}
  fi
  cd ..
done
"""