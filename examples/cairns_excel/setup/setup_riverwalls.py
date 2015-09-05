"""

Add riverwalls to the domain

Gareth Davies, Geoscience Australia 2014+

"""

def setup_riverwalls(domain, project):
    
    # #########################################################################
    #
    # Add Riverwalls [ must happen after distribute(domain) in parallel ]
    #
    # #########################################################################

    if not project.riverwalls == {}:
        domain.riverwallData.create_riverwalls(project.riverwalls,
                                               project.riverwall_par)
        domain.riverwallData.export_riverwalls_to_text(
            output_dir=project.output_dir + '/' +
            project.spatial_text_output_dir)

    return
