project('anuga', 'c', 'cpp', 'cython', version: '1.0')

# Import python and cython modules
py3 = import('python').find_installation(pure: false)
dep_py = py3.dependency()

# Find NumPy
dep_numpy = dependency('numpy')

dependencies = [dep_py, dep_numpy]

# Add subdirectories which contains python sources
subdir('anuga')



