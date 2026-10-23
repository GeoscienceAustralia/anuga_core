#!/usr/bin/python
"""
Setup weir / orifice structures from project data.

Uses Weir_orifice_trapezoid_operator, which computes discharge via combined
weir and orifice flow formulae transitioning smoothly between flow regimes.

Gareth Davies / ANUGA team
"""

from anuga.utilities import spatialInputUtil as su


def setup_weirs(domain, project):
    """Instantiate weir operators from ``project.weir_data``.

    ``project.weir_data`` is a list of dicts produced by
    ``ProjectDataTOML._parse_weirs``.  Each dict contains all keyword
    arguments needed by Weir_orifice_trapezoid_operator.

    Parameters
    ----------
    domain : anuga.Domain
        The simulation domain.
    project : ProjectData or ProjectDataTOML
        Scenario configuration object.  If the attribute ``weir_data``
        is absent (e.g. an older Excel-based project file that predates
        weir support) the function returns without error.
    """
    from anuga.parallel import numprocs
    if numprocs > 1:
        from anuga.parallel.parallel_operator_factory import (
            Weir_orifice_trapezoid_operator)
    else:
        from anuga.structures.weir_orifice_trapezoid_operator import (
            Weir_orifice_trapezoid_operator)

    weir_data = getattr(project, 'weir_data', [])

    for wd in weir_data:
        # Resolve geometry: exchange lines (file paths) or end points
        if wd['exchange_line_0'] is not None:
            exchange_lines = [
                su.read_polygon(wd['exchange_line_0']),
                su.read_polygon(wd['exchange_line_1']),
            ]
            end_points = None
        else:
            exchange_lines = None
            end_points = [wd['end_point_0'], wd['end_point_1']]

        Weir_orifice_trapezoid_operator(
            domain=domain,
            losses=wd['losses'],
            width=wd['width'],
            height=wd['height'],
            barrels=wd['barrels'],
            blockage=wd['blockage'],
            z1=wd['z1'],
            z2=wd['z2'],
            exchange_lines=exchange_lines,
            end_points=end_points,
            invert_elevations=wd['invert_elevations'],
            apron=wd['apron'],
            manning=wd['manning'],
            enquiry_gap=wd['enquiry_gap'],
            smoothing_timescale=wd['smoothing_timescale'],
            use_momentum_jet=wd['use_momentum_jet'],
            use_velocity_head=wd['use_velocity_head'],
            label=wd['label'],
            logging=True,
            verbose=True,
        )
