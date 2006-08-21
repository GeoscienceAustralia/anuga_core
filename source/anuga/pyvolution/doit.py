import os

#Clean up

for file in os.listdir('.'):
    if file[-1] == '~':
    	os.remove(file)


#os.system('python compile.py')
os.system('python test_all.py')
os.system('python run_profile.py')
os.system('python show_balanced_limiters.py')
