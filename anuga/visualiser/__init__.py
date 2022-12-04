

try:
    import vtk

    from .realtime import RealtimeVisualiser
    from .offline import OfflineVisualiser
    from .visualiser import Visualiser
except:
    # we do this so that nosetests doesn't create an error when searching
    # this sub_package
    pass
