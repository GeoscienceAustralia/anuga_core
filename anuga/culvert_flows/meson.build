# Title: Meson build file for culvert_flows

python_sources = [
  'culvert_class.py',
  'culvert_polygons.py',
  'culvert_routines.py',
  '__init__.py',
  'new_culvert_class.py',
]

py3.install_sources(
  python_sources,
  subdir: 'anuga/culvert_flows'
)

subdir('tests')