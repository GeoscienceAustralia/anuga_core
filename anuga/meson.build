# project('anuga', 'cython', version: '1.0')

# # Import python and cython modules
# py3 = import('python').find_installation()
# dep_py = py3.dependency()


# #cython = import('cython')

# # Find NumPy
# dep_numpy = dependency('numpy')

# dependencies = [dep_py, dep_numpy]




python_sources = [
    'anuga_exceptions.py',
    '__config__.py',
    'config.py',
    'error_api.py',
    'extras.py',
    '__init__.py',
    '__metadata__.py',
    'old_setup.py',
    'revision.py',
]



py3.install_sources(
  python_sources,
  subdir: 'anuga'
)


# Add subdirectories for each module
subdir('abstract_2d_finite_volumes')
# subdir('advection')
# subdir('alpha_shape')
# subdir('caching')
# subdir('coordinate_transforms')
# subdir('culvert_flows')
# subdir('damage_modelling')
# subdir('file')
# subdir('file_conversion')
# subdir('fit_interpolate')
# subdir('geometry')
# subdir('geospatial_data')
# subdir('lib')
# subdir('load_mesh')
# subdir('mesh_engine')
# subdir('operators')
# subdir('parallel')
# subdir('pmesh')
# subdir('rain')
# subdir('shallow_water')
# subdir('simulation')
# subdir('structures')
# subdir('tsunami_source')
# subdir('utilities')
# subdir('validation_utilities')
