import os
import sys
import py_compile

# Function to build a .pyc from a .py
def build_pyc(target, source, env):
    python = sys.executable
    command = "import py_compile; py_compile.compile('%s', doraise=True)"
    if sys.platform == 'win32':
        command %= str(source[0]).replace('\\', '\\\\')
    else:
        command %= str(source[0])
    rv = os.system('%s -c "%s"' % (python, command))
    if not rv == 0:
        raise SyntaxError, "Could not compile %s" % str(source[0])
    return None

# Function to build a .pyo from a .py
def build_pyo(target, source, env):
    python = sys.executable
    command = "import py_compile; py_compile.compile('%s', doraise=True)"
    if sys.platform == 'win32':
        command %= str(source[0]).replace('\\', '\\\\')
    else:
        command %= str(source[0])
    options = '-' + 'O' * env['OPTIMISATION_LEVEL']
    rv = os.system('%s %s -c "%s"' % (python, options, command))
    if not rv == 0:
        raise SyntaxError, "Could not compile %s" % str(source[0])
    return None

# Create the builders
pyc_builder = Builder(action = build_pyc,
                      suffix = '.pyc',
                      src_suffix = '.py')
pyo_builder = Builder(action = build_pyo,
                      suffix = '.pyo',
                      src_suffix = '.py')

# Read in build options
opts = Options('build_options.py')
opts.Add('GCCFLAGS', 'Flags passed to GCC')
opts.Add('MSVCFLAGS', 'Flags passed to the MSVC compiler')
opts.Add('METIS_DIR', 'Location of the metis directory relative to the pymetis directory')
opts.Add('OPTIMISATION_LEVEL', 'Optimisation level for the building of .pyo files. Must be 1 or 2')

# Windows' site packages are in a different location to those on Linux
if sys.platform == 'win32':
    opts.Add(PathOption('PREFIX',
                        'Location to install compiled sources',
                        os.path.join(sys.exec_prefix, 'lib', 'site-packages')))
else:
    opts.Add(PathOption('PREFIX',
                        'Location to install compiled sources',
                        os.path.join(sys.exec_prefix, 'lib', 'python' + sys.version[:3], 'site-packages')))

opts.Add(BoolOption('INSTALL_PYTHON_SOURCE',
                    'Install the .py files as well as the .pyc/.pyo files',
                    0))

env = Environment(options = opts)
env.Append(BUILDERS={'Pyc' : pyc_builder,
                     'Pyo' : pyo_builder})

Help(opts.GenerateHelpText(env))

if not (env['OPTIMISATION_LEVEL'] == 1 or env['OPTIMISATION_LEVEL'] == 2):
    raise ValueError, "OPTIMISATION_LEVEL must be between 1 and 2 inclusive"

# Where to find the Python.h
if sys.platform == 'win32':
    # Prefer MinGW over MSVC
    Tool('mingw')(env)
    
    python_include = os.path.join(sys.exec_prefix, 'include')
    # Windows installs need to be told the lib path and the python library name
    # else it won't work right.
    python_libdir = os.path.join(sys.exec_prefix, 'libs')
    env.Append(LIBPATH=[python_libdir])
    python_libname = 'python%d%d' % (sys.version_info[0:2])
    env.Append(LIBS=[python_libname])
else:
    python_include = os.path.join(sys.exec_prefix, 'include', 'python' + sys.version[:3])
    
# Check existence of Python.h
headerfile = python_include + os.sep + 'Python.h'
try:
    open(headerfile, 'r')
except:
    raise """Did not find Python header file %s.
    Make sure files for Python C-extensions are installed.
    In debian linux, for example, you need to install a
    package called something like python2.3-dev""" %headerfile

env.Append(CPPPATH=[python_include])

# Set appropriate CCFLAGS for the compiler.
if env['CC'] == 'gcc':
    env.Append(CCFLAGS=['${GCCFLAGS}'])
elif env['CC'] == 'cl':
    env.Append(CCFLAGS=['${MSVCFLAGS}'])

# Suppress the "lib" prefix on built files
env['SHLIBPREFIX'] = ""

anuga_root = os.path.join('source', 'anuga')
install_root = os.path.join(env['PREFIX'], 'anuga')

# Build .pyc and .pyo files of every .py in here and below.
files = []
#dirs = filter(os.path.isdir, os.listdir('.'))
#Only build the source/anuga directory for now
dirs = [anuga_root]
while(dirs != []):
    dirs += filter(os.path.isdir, map(lambda x : os.path.join(dirs[0], x), os.listdir(dirs[0])))
    files += filter(lambda x : x[-3:] == '.py', map(lambda x : os.path.join(dirs[0], x), os.listdir(dirs[0])))
    dirs = dirs[1:]
for x in files:
    env.Pyc(x + 'c', x)
    env.Pyo(x + 'o', x)
    # We have to cut the first character off the result of os.path.dirname(x).replace('anuga_core, '')
    # or else we will start taking paths relative to the root directory.
    instdir = os.path.join(env['PREFIX'], os.path.dirname(x).replace('source', '')[1:])
    env.Install(instdir, x + 'c')
    env.Install(instdir, x + 'o')
    if env['INSTALL_PYTHON_SOURCE']:
        env.Install(instdir, x)

# Define the install target
env.Alias('install', '${PREFIX}')

# Mesh_engine
mesh_env = env.Copy()
mesh_env.Append(CPPDEFINES=[('TRILIBRARY', 1), ('NO_TIMER', 1)])

mesh_dir = os.path.join(anuga_root, 'mesh_engine')
mesh_install_dir = os.path.join(install_root, 'mesh_engine')

triang = mesh_env.SharedLibrary(os.path.join(mesh_dir, 'triang'), map(lambda s: os.path.join(mesh_dir, s), ['triangle.c', 'triang.c']))
triangle = mesh_env.SharedLibrary(os.path.join(mesh_dir, 'triangle'), map(lambda s: os.path.join(mesh_dir, s), ['triangle.c']))
env.Install(mesh_install_dir, triang)
env.Install(mesh_install_dir, triangle)

# Utilities
util_dir = os.path.join(anuga_root, 'utilities')
util_install_dir = os.path.join(install_root, 'utilities')

env.Install(util_install_dir, env.SharedLibrary(os.path.join(util_dir, 'polygon_ext'),
                                                map(lambda s: os.path.join(util_dir, s), ['polygon_ext.c'])))
env.Install(util_install_dir, env.SharedLibrary(os.path.join(util_dir, 'util_ext'),
                                                map(lambda s: os.path.join(util_dir, s), ['util_ext.c'])))
env.Install(util_install_dir, env.SharedLibrary(os.path.join(util_dir, 'sparse_ext'),
                                                map(lambda s: os.path.join(util_dir, s), ['sparse_ext.c'])))

# Abstract_2d_finite_volumes
a2fv_env = env.Copy()
a2fv_env.Append(CPPPATH=[os.path.join('#', anuga_root, 'utilities')])

a2fv_dir = os.path.join(anuga_root, 'abstract_2d_finite_volumes')
a2fv_install_dir = os.path.join(install_root, 'abstract_2d_finite_volumes')

env.Install(a2fv_install_dir, a2fv_env.SharedLibrary(os.path.join(a2fv_dir, 'quantity_ext'),
                                                     map(lambda s: os.path.join(a2fv_dir, s), ['quantity_ext.c'])))
env.Install(a2fv_install_dir, a2fv_env.SharedLibrary(os.path.join(a2fv_dir, 'shallow_water_kinetic'),
                                                     map(lambda s: os.path.join(a2fv_dir, s), ['shallow_water_kinetic_ext.c'])))

# Shallow_water
sw_env = env.Copy()
sw_env.Append(CPPPATH=[os.path.join('#', anuga_root, 'utilities')])

sw_dir = os.path.join(anuga_root, 'shallow_water')
sw_install_dir = os.path.join(install_root, 'shallow_water')

env.Install(sw_install_dir, sw_env.SharedLibrary(os.path.join(sw_dir, 'shallow_water_ext'),
                                                 map(lambda s: os.path.join(sw_dir, s), ['shallow_water_ext.c'])))
