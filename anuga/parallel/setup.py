
def configuration(parent_package='',top_path=None):
    
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    
    
    config = Configuration('parallel', parent_package, top_path)

    try:
        import mpi4py
    except:
        pass
    else:
        config.add_data_dir('tests')
        config.add_data_dir('data')
    
    return config
    
if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)



