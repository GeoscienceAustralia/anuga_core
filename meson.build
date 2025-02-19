project('anuga', 'c', 'cpp', 'cython', version: '3.2.dev0')

# Import python
py3 = import('python').find_installation(pure: false)

message('Python found at ' + py3.get_variable('prefix'))

dep_py = py3.dependency()

# Find NumPy
dep_numpy = dependency('numpy', method : 'config-tool', required : true)

dependencies = [dep_py, dep_numpy]

# Add subdirectories which contains python sources
subdir('anuga')





