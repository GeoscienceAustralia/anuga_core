
import anuga
from os.path import join


# We need to define this function and importantly the wrapping function so that
# pickle can find definitions of the functions. (Pickle cannot import the definitions otherwise)
# func = anuga.file_function(join('Forcing','Tide','Pioneer.tms'), quantities='rainfall')

# And pickle doesn't like lambda functions so we define an explicit function with name
#def wrapped_file_function(t):
#    func = anuga.file_function(join('Forcing','Tide','Pioneer.tms'), quantities='rainfall')
#    return [func(t), 0.0, 0.0]


def setup_boundaries(simulation):
    """
    Setup boundary conditions
    """

    domain  = simulation.domain

    Bd = anuga.Dirichlet_boundary([0,0,0])
    #Bw = anuga.Time_boundary(domain=domain, function=wrapped_file_function)

    func = anuga.file_function(join('Forcing','Tide','Pioneer.tms'), quantities='rainfall')
    Bw = anuga.Time_boundary(domain=domain, function=lambda t : [func(t), 0.0, 0.0] )


    domain.set_boundary({'west': Bd, 'south': Bd, 'north': Bd, 'east': Bw})
