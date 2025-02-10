project('anuga', 'c', 'cpp', 'cython', version: '1.0')

# Import python and cython modules
py3 = import('python').find_installation(pure: false)
dep_py = py3.dependency()

# Find NumPy
dep_numpy = dependency('numpy')

dependencies = [dep_py, dep_numpy]

# incdir_numpy = run_command(py3,
#     [
#       '-c',
#       '''import os
# os.chdir(os.path.join("..", "tools"))
# import numpy as np
# try:
#   incdir = os.path.relpath(np.get_include())
# except Exception:
#   incdir = np.get_include()
# print(incdir)
# ''',
# ]).stdout().strip()


#cc = meson.get_compiler('c')
#cpp = meson.get_compiler('cpp')
#cy = meson.get_compiler('cython')

#cython = find_program(cy.cmd_array()[0])


# Add subdirectories which contains python sources
subdir('anuga')



