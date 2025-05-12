

import sys


def configuration(parent_package='',top_path=None):

    from numpy.distutils.misc_util import Configuration

    config = Configuration('anuga',parent_package,top_path)
    config.add_subpackage('abstract_2d_finite_volumes')
    config.add_subpackage('advection')
    config.add_subpackage('alpha_shape')
    config.add_subpackage('caching')
    config.add_subpackage('coordinate_transforms')
    config.add_subpackage('culvert_flows')
    config.add_subpackage('damage_modelling')
    config.add_subpackage('file')
    config.add_subpackage('file_conversion')
    config.add_subpackage('fit_interpolate')
    config.add_subpackage('geometry')
    config.add_subpackage('geospatial_data')
    config.add_subpackage('lib')
    config.add_subpackage('load_mesh')
    config.add_subpackage('mesh_engine')
    config.add_subpackage('operators')
    config.add_subpackage('parallel')
    config.add_subpackage('pmesh')
    config.add_subpackage('simulation')
    config.add_subpackage('shallow_water')
    config.add_subpackage('structures')
    config.add_subpackage('tsunami_source')
    config.add_subpackage('utilities')
    config.add_subpackage('validation_utilities')


    try:
        import vtk
        config.add_subpackage('visualiser')
    except:
        pass

    config.make_config_py()

    return config

if __name__ == '__main__':
    print('This is the wrong setup.py file to run')
    #from numpy.distutils.core import setup
    #setup(**configuration(top_path='').todict())
