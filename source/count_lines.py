"""Count total number of lines of code in ANUGA 
"""

import os

cmd_string = 'find . -type f -name "*.py" -print | xargs wc -l'
print cmd_string
os.system(cmd_string)


cmd_string = 'find . -type f -name "*.c" -print | xargs wc -l'
print cmd_string
os.system(cmd_string)


