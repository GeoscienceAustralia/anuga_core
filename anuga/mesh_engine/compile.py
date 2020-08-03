"""compile.py - compile Python C-extension

   Commandline usage: 
     python compile.py <filename>

"""

from past.builtins import execfile
import os
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')
