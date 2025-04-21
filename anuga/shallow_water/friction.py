"""Friction forcing functions.

Constraints: See license in the user guide


"""


from warnings import warn
import numpy as num
from copy import copy

import anuga
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a


# --------------------------------------------------------------------------
# Friction Forcing Terms
# --------------------------------------------------------------------------


def manning_friction_implicit(domain):
    
    if domain.multiprocessor_mode in [0,1,2,3]:
        manning_friction_implicit_cpu(domain)
    elif domain.multiprocessor_mode == 4:
        manning_friction_implicit_gpu(domain)





def manning_friction_implicit_cpu(domain):
    """Apply (Manning) friction to water momentum
    Wrapper for c version.
    FIXME SR: Thi whole module should be replaced with a call to the C code
    in sw_domain_orig_ext.py
    """

    if domain.multiprocessor_mode == 2:
        from .sw_domain_openmp_ext import manning_friction_flat
        from .sw_domain_openmp_ext import manning_friction_sloped
    else:
        from .sw_domain_orig_ext import manning_friction_flat
        from .sw_domain_orig_ext import manning_friction_sloped

    
    xmom = domain.quantities['xmomentum']
    ymom = domain.quantities['ymomentum']

    # really only need this if using sloped mannings
    x = domain.get_vertex_coordinates()

    w = domain.quantities['stage'].centroid_values
    z = domain.quantities['elevation'].centroid_values
    zv = domain.quantities['elevation'].vertex_values

    uh = xmom.centroid_values
    vh = ymom.centroid_values
    eta = domain.quantities['friction'].centroid_values

    xmom_update = xmom.semi_implicit_update
    ymom_update = ymom.semi_implicit_update

    eps = domain.minimum_allowed_height
    g = domain.g

    if domain.use_sloped_mannings:
        manning_friction_sloped(g, eps, x, w, uh, vh, zv, eta, xmom_update, \
                                ymom_update)
    else:
        manning_friction_flat(g, eps, w, uh, vh, z, eta, xmom_update, \
                                ymom_update)




def manning_friction_explicit(domain):
    
    if domain.multiprocessor_mode in [0,1,2,3]:
        manning_friction_explicit_cpu(domain)
    elif domain.multiprocessor_mode == 4:
        manning_friction_explicit_gpu(domain)




def manning_friction_explicit_cpu(domain):
    """Apply (Manning) friction to water momentum
    Wrapper for c version
    """

    if domain.multiprocessor_mode == 2:
        from .sw_domain_openmp_ext import manning_friction_flat
        from .sw_domain_openmp_ext import manning_friction_sloped
    else:
        from .sw_domain_orig_ext import manning_friction_flat
        from .sw_domain_orig_ext import manning_friction_sloped


    xmom = domain.quantities['xmomentum']
    ymom = domain.quantities['ymomentum']

    x = domain.get_vertex_coordinates()

    w = domain.quantities['stage'].centroid_values
    z = domain.quantities['elevation'].centroid_values
    zv = domain.quantities['elevation'].vertex_values

    uh = xmom.centroid_values
    vh = ymom.centroid_values
    eta = domain.quantities['friction'].centroid_values

    xmom_update = xmom.explicit_update
    ymom_update = ymom.explicit_update

    eps = domain.minimum_allowed_height

    if domain.use_sloped_mannings:
        manning_friction_sloped(domain.g, eps, x, w, uh, vh, zv, eta, xmom_update, \
                            ymom_update)
    else:
        manning_friction_flat(domain.g, eps, w, uh, vh, z, eta, xmom_update, \
                            ymom_update)



#GPU version of manning_friction_implicit that'll call the kernal written in sw_domain_cuda
def manning_friction_implicit_gpu(domain):
    """Apply (Manning) friction to water momentum
    Wrapper for c version
    """
    if domain.use_sloped_mannings:
        domain.gpu_interface.compute_forcing_terms_manning_friction_sloped()
    else:
        domain.gpu_interface.compute_forcing_terms_manning_friction_flat()


#GPU version of manning_friction_explicit that'll call the kernal written in sw_domain_cuda
def manning_friction_explicit_gpu(domain):
    """Apply (Manning) friction to water momentum
    Wrapper for c version
    """
    if domain.use_sloped_mannings:
        domain.gpu_interface.compute_forcing_terms_manning_friction_sloped()
    else:
        domain.gpu_interface.compute_forcing_terms_manning_friction_flat()


# FIXME (Ole): This was implemented for use with one of the analytical solutions
def linear_friction(domain):
    """Apply linear friction to water momentum

    Assumes quantity: 'linear_friction' to be present
    """

    w = domain.quantities['stage'].centroid_values
    z = domain.quantities['elevation'].centroid_values
    h = w-z

    uh = domain.quantities['xmomentum'].centroid_values
    vh = domain.quantities['ymomentum'].centroid_values
    tau = domain.quantities['linear_friction'].centroid_values

    xmom_update = domain.quantities['xmomentum'].semi_implicit_update
    ymom_update = domain.quantities['ymomentum'].semi_implicit_update

    num_tris = len(domain)
    eps = domain.minimum_allowed_height

    for k in range(num_tris):
        if tau[k] >= eps:
            if h[k] >= eps:
                S = -tau[k]/h[k]

                #Update momentum
                xmom_update[k] += S*uh[k]
                ymom_update[k] += S*vh[k]

def depth_dependent_friction(domain, default_friction,
                             surface_roughness_data,
                             verbose=False):
    """Returns an array of friction values for each wet element adjusted for
            depth.

    Inputs:
        domain - computational domain object
        default_friction - depth independent bottom friction
        surface_roughness_data - N x 5 array of n0, d1, n1, d2, n2 values
        for each friction region.

    Outputs:
        wet_friction - Array that can be used directly to update friction as
                        follows:
                       domain.set_quantity('friction', wet_friction)



    """

    default_n0 = 0  # James - this was missing, don't know what it should be

    # Create a temp array to store updated depth dependent
    # friction for wet elements
    # EHR this is outwardly inneficient but not obvious how to avoid
    # recreating each call??????

    wet_friction    = num.zeros(len(domain), float)
    wet_friction[:] = default_n0  # Initially assign default_n0 to all array so
                                  # sure have no zeros values

    # create depth instance for this timestep
    depth = domain.create_quantity_from_expression('stage - elevation')
    # Recompute depth as vector
    d_vals = depth.get_values(location='centroids')

    # rebuild the 'friction' values adjusted for depth at this instant
    # loop for each wet element in domain

    for i in domain.get_wet_elements():
        # Get roughness data for each element
        d1 = float(surface_roughness_data[i, 1])
        n1 = float(surface_roughness_data[i, 2])
        d2 = float(surface_roughness_data[i, 3])
        n2 = float(surface_roughness_data[i, 4])


        # Recompute friction values from depth for this element

        if d_vals[i] <= d1:
            ddf = n1
        elif d_vals[i] >= d2:
            ddf = n2
        else:
            ddf = n1 + ((n2 - n1) / (d2 - d1)) * (d_vals[i] - d1)

        # Check sanity of result
        if ddf < 0.010 or ddf > 9999.0:
            log.critical('>>>> WARNING: computed depth_dependent friction '
                         'out of range, ddf%f, n1=%f, n2=%f'
                         % (ddf, n1, n2))

        # update depth dependent friction  for that wet element
        wet_friction[i] = ddf

    # EHR add code to show range of 'friction across domain at this instant as
    # sanity check?????????

    if verbose:
        # return array of domain nvals
        nvals = domain.get_quantity('friction').get_values(location='centroids')
        n_min = min(nvals)
        n_max = max(nvals)

        log.critical('         ++++ calculate_depth_dependent_friction - '
                     'Updated friction - range  %7.3f to %7.3f'
                     % (n_min, n_max))

    return wet_friction


