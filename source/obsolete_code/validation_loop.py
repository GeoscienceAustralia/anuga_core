
import os
revisions = [3929, 3975,4000,4035 ]
for revision in revisions:
    s = 'svn up -r %s' %(revision)
    print s
    os.system(s)
    
    s = 'python set_pp_comp_test_all.py'
    print s
    os.system(s)

    root = os.sep    
    av_dir = os.path.join(root,'d','cit','1','dgray','validate_inundation','ga','anuga_validation','automated_validation_tests')

    os.chdir(av_dir)
    s = 'python set_pythonpath_validate.py'
    print s
    os.system(s)

    
