from anuga.shallow_water.shallow_water_domain import Domain

# The different states the model can be in
STATE_SETUP, STATE_BUILT, STATE_RUNNING = range(0, 3)

class Model:   
    def __init__(self, name):
        self.domain = None   # Domain class will be removed as we refactor the API
        self.geometry_data = None
        self.quantities = []
        self.boundaries = []
        self._state = STATE_SETUP
        self.name = name

    
    def set_geometry(self, *args, **kwargs):
        if self.geometry_data:
            raise RuntimeError, 'geometry already defined in model.'
        elif self._state != STATE_SETUP:
            raise RuntimeError, 'Model already built - call this before build().'
        else:
            self.geometry_data = (args, kwargs)

        
    def set_quantity(self, *args, **kwargs):
        """Set values for named quantity
        """
        if self._state != STATE_SETUP:
            raise RuntimeError, 'Model already built - call this before build().'

        self.quantities.append((args, kwargs))

    def get_quantity(self, name):
        if self._state != STATE_RUNNING:
            raise RuntimeError, 'Model needs to be running - call this after run().'
            
        return self.domain.get_quantity(name)


    def set_boundary(self, *args, **kwargs):
        """Set the boundary properties.
        """
        if self._state == STATE_SETUP:
            # dont apply boundaries until model is built
            self.boundaries.append((args, kwargs))
        else:
            # model is running, apply boundaries right away
            self.domain.set_boundary(*args, **kwargs)
   
    def get_normals(self):
        # FIXME: This is a wrapper to allow reflective boundary to work.
        # Should be refactored and removed.
        return domain.get_normals()

    def build(self):        
        if not self.geometry_data:
            raise RuntimeError, 'Model needs geometry. Set mesh or geometry.'

        if self._state != STATE_SETUP:
            raise RuntimeError, 'Model already built.'

        if len(self.geometry_data[0]) == 1:
            raise RuntimeError, 'Load mesh not implemented.'
        else:
            args, kwargs = self.geometry_data

            self.domain = Domain(*args, **kwargs)
            self.domain.set_name(self.name)

            for args, kwargs in self.quantities:
                self.domain.set_quantity(*args, **kwargs)

            for args, kwargs in self.boundaries:
                self.domain.set_boundary(*args, **kwargs)

            self._state = STATE_BUILT
            self.geometry_data = None
            self.quantities = None
        

    def run(self, yieldstep, finaltime, verbose = False, timestep_function = None):
        if self._state != STATE_BUILT:
            msg = 'Model is not built, cannot run yet.'
            msg += ' Call build() before run().'
            raise RuntimeError, msg
 
        self._state = STATE_RUNNING
        
        for t in self.domain.evolve(yieldstep, finaltime):
            if verbose == True:
                print self.domain.timestepping_statistics()
                
            # call the user's function at every timestep
            if timestep_function:
                timestep_function(self, t)
