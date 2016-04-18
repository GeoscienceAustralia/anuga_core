import sys
import shutil
import argparse
import os

parser = argparse.ArgumentParser(description = 'Install the ANUGA sediment transport and vegetation operators into an existing ANUGA package.')

parser.add_argument('--link',
                    action = 'store_true',
                    help = 'Creates symbolic links within the installed ANUGA package instead of stand-alone files to preserve version control.')

args = parser.parse_args()

try:
    import anuga
except:
    print '\n'
    print '#' * 60
    print 'CANNOT INSTALL THE OPERATORS!'
    print 'ANUGA must be installed before running this script.'
    print 'Find it at https://github.com/GeoscienceAustralia/anuga_core'
    print '#' * 60
    print '\n'
    sys.exit()
    

path_to_anuga = anuga.__path__[0]

files = ['operators' + os.sep + 'sed_transport_operator.py',
         'operators' + os.sep + 'vegetation_operator.py']

full_path = os.path.realpath(__file__)
path, filename = os.path.split(full_path)
dir_name = full_path.split(os.sep)[-3]

def make_symlink():

    os.symlink(src, dst)
    print '\n' + src + ' --> ' + dst
    
    return counter + 1


if args.link:
    
    if dir_name == 'Downloads':
    
        print '\n'
        print '#' * 60
        print 'WARNING!'
        print 'The --link flag creates symbolic links from the current position'
        print 'of the operator scripts to the ANUGA package.'
        print "Make sure that you've moved the anugaSed folder to its final location"
        print "before running this script (ie. get it out of the Downloads folder!)."
        print '#' * 60
        print 'No symlinks were created.'
        print '\n'
        
    else:

        print '\n'
        print '#' * 60
        print 'Creating symbolic links:'
        print '#' * 60
    
        counter = 0

        for f in files:
    
            src = path + os.sep + f
            dst = path_to_anuga + os.sep + f
        
            try:
                counter = make_symlink()
                print counter
            except:
                print '\nThe file ' + f.split(os.sep)[-1] + ' already exists within the anuga package.'
                ovwrt = raw_input('Overwrite? (y/n): ')
            
                if ovwrt == 'y':
                
                    os.remove(dst)
                    counter = make_symlink()
                
                else:
                
                    print '\nNot linking ' + src

        skips = str(len(files) - counter)

        print '\n'
        print '#' * 60
        print 'Created ' + str(counter) + ' symlinks. Skipped ' + skips + ' files.\n'
    

else: # if not using the --link flag

    for f in files:
    
        src = path + os.sep + f
        dst = path_to_anuga + os.sep + f
        shutil.copy2(src, dst)

    print 'Success!'
    print 'Installed', str(len(files)), 'files to', path_to_anuga
    
    
    
    
    
    
    
    