project('anuga', 'c', 'cpp', 'cython', version: '1.0')

# Import python and cython modules
py3 = import('python').find_installation(pure: false)
dep_py = py3.dependency()


#cython = import('cython')

# Find NumPy
dep_numpy = dependency('numpy')

dependencies = [dep_py, dep_numpy]


cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')
cy = meson.get_compiler('cython')

cython = find_program(cy.cmd_array()[0])

# Add subdirectories which contains python sources
subdir('anuga')



