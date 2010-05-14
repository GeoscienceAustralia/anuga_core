from anuga.shallow_water.domain import Domain


class Model:
    def __init__(self):
        self.domain = None   # Domain class will be removed as we refactor the API
        self.geometry_data = None
        self.quantities = []
    
    def set_geometry(self, *args, **kwargs):
        if self.geometry_data:
            raise InputError, 'geometry already defined in model.'
        else:
            self.geometry_data = (args, kwargs)
        
    def set_name(self, name):
        self.name = name
        
    def set_quantity(self, *args, **kwargs):
        """Set values for named quantity
        """

        self.quantities.append((args, kwargs))

        
    def build(self):
        from anuga.shallow_water import Reflective_boundary
        from anuga.shallow_water import Dirichlet_boundary

        args, kwargs = self.geometry_data

        self.domain = Domain(*args, **kwargs)  # Create domain
        self.domain.set_name('channel1')                  # Output name

        #------------------------------------------------------------------------------
        # Setup initial conditions
        #------------------------------------------------------------------------------
        def topography(x, y):
            return -x/10                             # linear bed slope

        self.domain.set_quantity('elevation', topography) # Use function for elevation
        self.domain.set_quantity('friction', 0.01)        # Constant friction 
        self.domain.set_quantity('stage',                 # Dry bed
                            expression='elevation')          
        return
        
        
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
        if not self.name:
            raise InputError, 'Model needs a name: use model.set_name().'
        
        if not self.geometry_data:
            raise InputError, 'Model needs geometry. Set mesh or geometry.'

        if len(self.geometry_data[0]) == 1:
            raise InputError, 'Load mesh not implemented.'
        else:
            args, kwargs = self.geometry_data


            points, vertices, boundary = rectangular_cross(10, 5,
                                                           len1=10.0, len2=5.0) # Mesh

            domain = Domain(points, vertices, boundary)  # Create domain

        #    domain = Domain(*args, **kwargs)
            domain.set_name(self.name)

            def topography(x, y):
                return -x/10    
            domain.set_quantity('elevation', topography) # Use function for elevation
            domain.set_quantity('friction', 0.01)        # Constant friction 
            domain.set_quantity('stage',                 # Dry bed
                                expression='elevation')  
            
            print self.quantities
            
            for args, kwargs in self.quantities:
                self.domain.set_quantity(*args, **kwargs)
        
        
