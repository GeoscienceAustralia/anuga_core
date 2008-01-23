"""f2py script -

"""

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('',parent_package,top_path)

    config.add_extension('advection_ext',
                         sources = ['advection_ext.pyf','advection_ext.c'])
    return config

#-------------------------------------------------------------
if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)


