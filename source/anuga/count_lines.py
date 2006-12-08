"""Count total number of lines of code in ANUGA 
"""

import os

s = 'find . -type f -name "*.py" -print | xargs wc -l'
print s
os.system(s)


s = 'find . -type f -name "*.c" -print | xargs wc -l'
print s
os.system(s)


