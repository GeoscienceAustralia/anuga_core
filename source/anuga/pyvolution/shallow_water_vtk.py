"""
Shallow water domain with VTK viewer. We are just over

Ole Nielsen, Duncan Gray
Geoscience Australia, 2006

Stephen Roberts,  Jack Kelly
ANU 2006
"""

#Subversion keywords:
#
#$LastChangedDate: 2006-01-13 12:33:38 +1100 (Fri, 13 Jan 2006) $
#$LastChangedRevision: 2205 $
#$LastChangedBy: steve $


from shallow_water import *
Shallow_Water_Domain = Domain #Rename

#Shallow water domain with VTK viewer
class Domain(Shallow_Water_Domain):

    def __init__(self, coordinates, vertices, boundary = None,
                 tagged_elements = None, geo_reference = None,
                 use_inscribed_circle=False):

        Shallow_Water_Domain.__init__(self, coordinates, vertices, boundary,
                                tagged_elements, geo_reference, use_inscribed_circle)



    def initialise_visualiser(self,scale_z=1.0,rect=None):
        #Realtime visualisation
        if self.visualiser is None:
            from vtk_realtime_visualiser import Visualiser
            self.visualiser = Visualiser(self,scale_z,rect)
            self.visualiser.coloring['stage'] = False
            self.visualiser.coloring['elevation'] = False
            self.visualiser.setup['elevation']=True
            self.visualiser.updating['stage']=True
            self.visualiser.qcolor['stage'] = (0.1,0.4,0.99)
        self.visualise = True
        if self.visualise_color_stage == True:
            self.visualiser.coloring['stage'] = True
            



    def evolve(self,
               yieldstep = None,
               finaltime = None,
               duration = None,
               skip_initial_step = False):
        """Specialisation of basic evolve method from parent class
        """

        import time
        first = True
        #Call check integrity here rather than from user scripts
        #self.check_integrity()

        msg = 'Parameter beta_h must be in the interval [0, 1['
        assert 0 <= self.beta_h < 1.0, msg
        msg = 'Parameter beta_w must be in the interval [0, 1['
        assert 0 <= self.beta_w < 1.0, msg

        self.distribute_to_vertices_and_edges()

        #Initialise real time viz if requested
        if self.visualise is True and self.time == 0.0 and self.visualiser is None:
            self.initialise_visualiser()
            #self.visualiser.update_timer()
            print "Warning: Enabling the visualiser with domain.visualise is not"
            print "recommended. Controlling the visualiser manually allows for much better"
            print "control over visualiser parameters."

        if self.visualise is True:
            self.visualiser.start()
            # Things go haywire if we start evolving before the vis is ready
            self.visualiser.idle.wait()
            self.visualiser.idle.clear()

        #Store model data, e.g. for visualisation
        if self.store is True and self.time == 0.0:
            self.initialise_storage()
            #print 'Storing results in ' + self.writer.filename
        else:
            pass

        #Call basic machinery from parent class
        for t in Generic_Domain.evolve(self,
                                       yieldstep=yieldstep,
                                       finaltime = finaltime,
                                       duration = duration,
                                       skip_initial_step = skip_initial_step):
            #Real time viz
            if self.visualise is True:
                self.visualiser.redraw_ready.set()
                self.visualiser.idle.wait()
                self.visualiser.idle.clear()
                self.visualiser.unpaused.wait()

            #Store model data, e.g. for subsequent visualisation
            if self.store is True:
                self.store_timestep(self.quantities_to_be_stored)

            #FIXME: Could maybe be taken from specified list
            #of 'store every step' quantities

            #Pass control on to outer loop for more specific actions
            yield(t)

#==================================== End of Updated Shallow Water VTK domain =============
if __name__ == "__main__":
    pass

# Profiling stuff
import profile
profiler = profile.Profile()
