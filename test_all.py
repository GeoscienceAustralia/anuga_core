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

print 'Changing to', os.getcwd() #This is now different from buildroot   

execfile('test_all.py')

    
