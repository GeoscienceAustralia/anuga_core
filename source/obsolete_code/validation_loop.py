
import os
revisions = [4204, 4106, 4182, 4175, 4174, 4160, 4130 ]
for revision in revisions:
    s = 'svn up -r %s source' %(revision)
    print s
    os.system(s)
    
    s = 'python set_pp_comp_test_all.py'
    print s
    os.system(s)

    #root = os.sep    
    #av_dir = os.path.join(root,'d','cit','1','dgray','validate_inundation','ga','anuga_validation','automated_validation_tests')
    av_dir = os.path.join('..','anuga_validation','automated_validation_tests','okushiri_tank_validation')
    os.chdir(av_dir)

    ## NOT SETTING PYTHONPATH
    s = 'python loading_pts.py'
    print s
    os.system(s)

    av_dir = os.path.join('..','..','..','anuga_core')
    os.chdir(av_dir)
