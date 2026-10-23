

python_sources = [
    'anuga_exceptions.py',
    'config.py',
    'error_api.py',
    'extras.py',
    '__init__.py',
    '__metadata__.py',
]

py3.install_sources(
  python_sources,
  subdir: 'anuga'
)

# Generate _version.py from git metadata at build time
_version_conf = configuration_data()
_version_conf.set('VERSION', meson.project_version())
_version_conf.set('GIT_SHA',  _git_sha)
_version_conf.set('GIT_DATE', _git_date)

configure_file(
  input:         '_version.py.in',
  output:        '_version.py',
  configuration: _version_conf,
  install:       true,
  install_dir:   py3.get_install_dir() / 'anuga',
)


# Add subdirectories for each module
subdir('abstract_2d_finite_volumes')
subdir('advection')
subdir('alpha_shape')
subdir('caching')
subdir('coordinate_transforms')
subdir('damage_modelling')
subdir('file')
subdir('file_conversion')
subdir('fit_interpolate')
subdir('geometry')
subdir('geospatial_data')
subdir('lib')
subdir('load_mesh')
subdir('mesh_engine')
subdir('operators')
subdir('parallel')
subdir('pmesh')
subdir('rain')
subdir('scenario')
subdir('shallow_water')
subdir('simulation')
subdir('structures')
subdir('tsunami_source')
subdir('utilities')
subdir('validation_utilities')
