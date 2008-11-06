import os

buildroot = os.getcwd()

os.chdir('source')

os.chdir('anuga')

#Complete horrible hack to decide which branch to take (Ole)
#try:
#    os.stat('inundation-numpy-branch')
#except:
#    os.chdir('inundation')
#else:
#    os.chdir('inundation-numpy-branch')    

print 'Changing to', os.getcwd()        

#entries = listdir('.')


#Attempt to compile all extensions

os.chdir('utilities')
execfile('compile.py')

os.chdir('..')
os.chdir('abstract_2d_finite_volumes')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

os.chdir('..')
os.chdir('advection')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')


os.chdir('..')
os.chdir('shallow_water')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

os.chdir('..')
os.chdir('mesh_engine')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

os.chdir(buildroot)    
#execfile('test_all.py')
    
