#!/usr/bin/python
"""
Setup box and pipe culverts from project data.

Supports Boyd_box_operator (type='boyd_box') and Boyd_pipe_operator
(type='boyd_pipe').  Both operator types accept either exchange-line files
or end-point coordinates to locate the culvert on the mesh.

Gareth Davies / ANUGA team
"""

from anuga.utilities import spatialInputUtil as su


def setup_culverts(domain, project):
    """Instantiate culvert operators from ``project.culvert_data``.

    ``project.culvert_data`` is a list of dicts produced by
    ``ProjectDataTOML._parse_culverts``.  Each dict contains all keyword
    arguments needed by Boyd_box_operator or Boyd_pipe_operator.

    Parameters
    ----------
    domain : anuga.Domain
        The simulation domain.
    project : ProjectData or ProjectDataTOML
        Scenario configuration object.  If the attribute ``culvert_data``
        is absent (e.g. an older Excel-based project file that predates
        culvert support) the function returns without error.
    """
    from anuga.parallel import numprocs
    if numprocs > 1:
        from anuga.parallel.parallel_operator_factory import (
            Boyd_box_operator, Boyd_pipe_operator)
    else:
        from anuga.structures.boyd_box_operator import Boyd_box_operator
        from anuga.structures.boyd_pipe_operator import Boyd_pipe_operator

    culvert_data = getattr(project, 'culvert_data', [])

    for cd in culvert_data:
        # Resolve geometry: exchange lines (file paths) or end points
        if cd['exchange_line_0'] is not None:
            exchange_lines = [
                su.read_polygon(cd['exchange_line_0']),
                su.read_polygon(cd['exchange_line_1']),
            ]
            end_points = None
        else:
            exchange_lines = None
            end_points = [cd['end_point_0'], cd['end_point_1']]

        common = dict(
            domain=domain,
            losses=cd['losses'],
            barrels=cd['barrels'],
            blockage=cd['blockage'],
            exchange_lines=exchange_lines,
            end_points=end_points,
            invert_elevations=cd['invert_elevations'],
            apron=cd['apron'],
            manning=cd['manning'],
            enquiry_gap=cd['enquiry_gap'],
            smoothing_timescale=cd['smoothing_timescale'],
            use_momentum_jet=cd['use_momentum_jet'],
            use_velocity_head=cd['use_velocity_head'],
            label=cd['label'],
            logging=True,
            verbose=True,
        )

        if cd['type'] == 'boyd_box':
            Boyd_box_operator(width=cd['width'], height=cd['height'],
                              z1=cd['z1'], z2=cd['z2'], **common)
        else:
            Boyd_pipe_operator(diameter=cd['diameter'], **common)
