"""
Script to run all the produce_results scripts in the
validation_tests/xxx/xxx/ directories
"""

import os
import time

import anuga
from anuga import indent

# Larger plot text and 200 dpi for every case's figures (see matplotlibrc)
os.environ['MATPLOTLIBRC'] = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'matplotlibrc')
#from anuga.validation_utilities.parameters import alg
#from anuga.validation_utilities.parameters import cfl


args = anuga.get_args()
alg = args.alg
np = args.np
verbose = args.verbose
no_run = getattr(args, 'no_run', False)
long_cases = getattr(args, 'long', False)

#---------------------------------
# Get the current svn revision
#---------------------------------
timestamp = time.asctime()
version = anuga.get_version()
version_split = version.split('.',3)

major_version = version_split[0]
minor_version = version_split[1]


#----------------------------------
# Now it is ok to create the latex
# macro file with run parameters
#
# FIXME: THis is a little dangerous as
# this is changed before all the tests
# are run.
#----------------------------------

f = open('saved_parameters.tex', 'w')
#f.write('\\newcommand{\\cfl}{\\UScore{%s}}\n' % str(cfl))
f.write('\\newcommand{\\alg}{\\UScore{%s}}\n' % str(alg))
f.write('\\newcommand{\\majorR}{\\UScore{%s}}\n' % str(major_version))
f.write('\\newcommand{\\minorR}{\\UScore{%s}}\n' % str(minor_version))
f.write('\\newcommand{\\timeR}{{%s}}\n' % str(timestamp))

f.close()

#---------------------------------
# Run the tests
#---------------------------------
os.chdir('..')
buildroot = os.getcwd()

Upper_dirs = os.listdir('.')
dir = '.'
Upper_dirs = [name for name in os.listdir(dir) if os.path.isdir(os.path.join(dir, name))]

try:
    Upper_dirs.remove('.svn')
except ValueError:
    pass

try:
    Upper_dirs.remove('reports')
except ValueError:
    pass

try:
    Upper_dirs.remove('case_studies')
except ValueError:
    pass

#print Upper_dirs
#os.chdir('./Tests')

#print 'Tests'
print(Upper_dirs)

time_total = 0.0
test_number = 1
failed = []
status = []      # (case, 'ok' | 'failed', seconds), for the combined report
for dir in Upper_dirs:
    upper_dir = dir

    os.chdir(dir)

    print(72 * '=')
    print('Directory: ' + dir)
    print(72 * '=')

    #print 'Changing to', os.getcwd()
    dir = '.'
    Lower_dirs =  [name for name in os.listdir(dir) if os.path.isdir(os.path.join(dir, name))]
    try:
        Lower_dirs.remove('.svn')
    except ValueError:
        pass
    #print Lower_dirs




    for l_dir in Lower_dirs:
        if not os.path.isfile(os.path.join(l_dir, 'produce_results.py')):
            # __pycache__ and the like
            continue

        os.chdir(l_dir)
        #print os.getcwd()
        print(60 * '=')
        print('Subdirectory %g: '% (test_number)  + l_dir)
        test_number += 1
        print(60 * '=')
        try:
            t0 = time.time()
            cmd = 'python produce_results.py -alg %s -np %s '% (str(alg),str(np))
            if verbose:
                cmd += '-v '
            if no_run:
                cmd += '-nr '
            if long_cases:
                cmd += '-l '
            print(2 * indent + 'Running: ' + cmd)
            # A case whose produce_results fails is reported at the end and
            # makes this script exit non-zero (the release workflow gates on
            # it), but the remaining cases still run and the report is still
            # typeset from whatever was produced.
            ok = os.system(cmd) == 0
            if not ok:
                failed.append(os.path.join(upper_dir, l_dir))
            t1 = time.time() - t0
            time_total += t1
            status.append((upper_dir + '/' + l_dir, 'ok' if ok else 'failed', t1))
            print(2 * indent + 'That took ' + str(t1) + ' secs')
        except Exception:
            print(2 * indent + 'Failed running produce_results in ' + os.getcwd())
            failed.append(os.path.join(upper_dir, l_dir))
            status.append((upper_dir + '/' + l_dir, 'failed', time.time() - t0))

        os.chdir('..')
        #print 'Changing to', os.getcwd()

    os.chdir('..')
    #print 'Changing to', os.getcwd()

os.chdir(buildroot)

print(72 * '=')
print('That took ' + str(time_total) + ' secs')
print(72 * '=')


# go back to reports directory to typeset report
os.chdir('reports')

# Per-case outcome and run time, read by validations_combine_algs.py
with open('case_status_alg_%s.csv' % str(alg), 'w') as f:
    f.write('case,status,seconds\n')
    for case, outcome, secs in status:
        f.write('%s,%s,%.1f\n' % (case, outcome, secs))


# A stale PDF from an earlier run must not pass for this one
if os.path.isfile('validations_report.pdf'):
    os.remove('validations_report.pdf')
os.system('python validations_typeset_report.py')

import sys

report = 'validations_report_alg_%s.pdf' % str(alg)
if os.path.isfile('validations_report.pdf'):
    os.replace('validations_report.pdf', report)
    print('Wrote ' + report)
else:
    print('validations_report.pdf was not produced (see validations_report.log)')
    failed.append('reports (typesetting)')

if failed:
    print(72 * '=')
    print('%d case(s) failed to produce results:' % len(failed))
    for name in failed:
        print(indent + name)
    print(72 * '=')
    sys.exit(1)
