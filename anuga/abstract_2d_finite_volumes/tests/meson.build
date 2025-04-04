# project('anuga', 'cython', version: '1.0')

# # Import python and cython modules
# py3 = import('python').find_installation()

python_sources = [
    'test_ermapper.py',
    'test_gauge.py',
    'test_general_mesh.py',
    'test_generic_boundary_conditions.py',
    'test_generic_domain.py',
    'test_ghost.py',
    'test_neighbour_mesh.py',
    'test_pmesh2domain.py',
    'test_quantity.py',
    'test_region.py',
    'test_tag_region.py',
    'test_util.py',
    'show_balanced_limiters.py',
    '__init__.py',
]

py3.install_sources(
  python_sources,
  subdir: 'anuga/abstract_2d_finite_volumes/tests'
)