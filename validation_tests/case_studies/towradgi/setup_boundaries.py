
import anuga
from os.path import join

def setup_boundaries(domain):
    """
    Setup boundary conditions
    """

    func = anuga.file_function(join('Forcing','Tide','Pioneer.tms'), quantities='rainfall')
    Bd = anuga.Dirichlet_boundary([0,0,0])
    Bw = anuga.Time_boundary(domain=domain, function=lambda t: [func(t)[0], 0.0, 0.0])
    #Bw = anuga.Time_boundary(domain=domain, function=func)

    domain.set_boundary({'west': Bd, 'south': Bd, 'north': Bd, 'east': Bw})

