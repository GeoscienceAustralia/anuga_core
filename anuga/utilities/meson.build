

# Compile the Cython-generated C code and additional C code
py3.extension_module('cg_ext',
  sources: ['cg_ext.pyx'],
  dependencies: dependencies,
  subdir: 'anuga/utilities',
  install: true,
)

py3.extension_module('quad_tree_ext',
  sources: ['quad_tree_ext.pyx'],
  dependencies: dependencies,
  subdir: 'anuga/utilities',
  install: true,
)

py3.extension_module('sparse_ext',
  sources: ['sparse_ext.pyx'],
  dependencies: dependencies,
  subdir: 'anuga/utilities',
  install: true,
)

py3.extension_module('sparse_matrix_ext',
  sources: ['sparse_matrix_ext.pyx'],
  dependencies: dependencies,
  subdir: 'anuga/utilities',
  install: true,
)

py3.extension_module('util_ext',
  sources: ['util_ext_c.pyx'],
  dependencies: dependencies,
  subdir: 'anuga/utilities',
  install: true,
)



python_sources = [
'animate.py',
'argparsing.py',
'cg_solve.py',
'csv_tools.py',
'data_audit.py',
'data_audit_wrapper.py',
'file_utils.py',
'function_utils.py',
'__init__.py',
'interp.py',
'log_analyser.py',
'log.py',
'mem_time_equation.py',
'model_tools.py',
'norms.py',
'numerical_tools.py',
'parallel_abstraction.py',
'parse.py',
'parse_time.py',
'plot_utils.py',
'quantity_setting_functions.py',
'quickPlots.py',
'run_anuga_script.py',
'sparse.py',
'spatialInputUtil.py',
'sww_merge.py',
'system_tools.py',
'terminal_width.py',
'where_close.py',
'xml_tools.py',
]


py3.install_sources(
  python_sources,
  subdir: 'anuga/utilities'
)

subdir('tests')
