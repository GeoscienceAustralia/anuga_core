#!/usr/bin/python
"""

Set up erosion / scour operators

"""

import anuga
from anuga.operators.erosion_operators import (
    Erosion_operator,
    Bed_shear_erosion_operator,
    Flat_slice_erosion_operator,
    Flat_fill_slice_erosion_operator,
)
from anuga.operators.sanddune_erosion_operator import Sanddune_erosion_operator


# Five behaviours, not seven classes: Circular_erosion_operator and
# Polygonal_erosion_operator are thin region-specification wrappers over the
# base Erosion_operator and add no erosion physics, so 'simple' with a
# center/radius or a polygon reaches them without needing names of their own.
_OPERATORS = {
    'simple':     Erosion_operator,
    'bed_shear':  Bed_shear_erosion_operator,
    'flat_slice': Flat_slice_erosion_operator,
    'flat_fill':  Flat_fill_slice_erosion_operator,
    'sand_dune':  Sanddune_erosion_operator,
}

# Parameters each type accepts beyond the common set. Kept here as well as in
# the parser because this module can also be driven by the Excel interface,
# which does no such validation.
_TYPE_PARAMS = {
    'bed_shear':  ('shear_factor',),
    'flat_slice': ('elevation',),
    'flat_fill':  ('elevation',),
    'sand_dune':  ('Ra',),
}


def setup_erosion(domain, project):
    """
    Function to add erosion / scour operators to the domain
    """
    erosion_data = getattr(project, 'erosion_data', [])

    operators = []
    for e in erosion_data:
        etype = e['type']
        cls = _OPERATORS[etype]

        kwargs = {
            'threshold': e.get('threshold', 0.0),
            'base': e.get('base', 0.0),
        }

        # Region: polygon points were resolved by prepare_data; center/radius
        # pass through as given. Exactly one form is present — the parser
        # rejects both and neither.
        if e.get('polygon_points') is not None:
            kwargs['polygon'] = e['polygon_points']
        elif e.get('center') is not None:
            kwargs['center'] = e['center']
            kwargs['radius'] = e['radius']

        for prm in _TYPE_PARAMS.get(etype, ()):
            if e.get(prm) is not None:
                kwargs[prm] = e[prm]

        for prm in ('description', 'label', 'logging'):
            if e.get(prm) is not None:
                kwargs[prm] = e[prm]

        operators.append(cls(domain, **kwargs))

    return operators
