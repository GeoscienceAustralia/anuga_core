"""
Set up sediment transport from the [sediment] table of a TOML scenario.

The table maps onto ``Domain.initialize_sediment_operator`` and its setters,
each ``[[sediment.fractions]]`` onto ``Domain.add_sediment_fraction``. Only
keys the scenario gives are passed on, so an absent key leaves the domain's
own default; see the sediment page of the documentation for what those are.
"""


def setup_sediment(domain, project):
    """Configure sediment transport on *domain* from *project*.

    Returns the sediment operator, or None when the scenario has no
    [sediment] table.
    """
    sed = getattr(project, 'sediment_data', None)
    if not sed:
        return None

    def given(*keys):
        return {k: sed[k] for k in keys if k in sed}

    operator = domain.initialize_sediment_operator(
        **given('porosity', 'c_max', 'c_pack', 'bed_evolution', 'rho_w'))

    if 'shear_closure' in sed:
        domain.set_shear_closure(sed['shear_closure'],
                                 max_slope=sed.get('max_slope'),
                                 freeze_slope=sed.get('freeze_slope', False))

    if 'friction_mode' in sed:
        domain.set_sediment_friction(
            mode=sed['friction_mode'],
            **given('k_s', 'r_d', 'r_br', 'sigma_br', 'bed', 'grain_size'))

    if any(k in sed for k in ('bed_material', 'tau_crit', 'K_e', 'entrainment')):
        domain.set_bed_material(
            material=sed.get('bed_material', 'noncohesive'),
            **given('tau_crit', 'K_e', 'entrainment', 'de_leeuw_fit',
                    'skin_roughness'))

    if any(k in sed for k in ('deposition_law', 'tau_d', 'near_bed',
                              'reference_height_floor')):
        domain.set_deposition(
            law=sed.get('deposition_law', 'd_star'),
            **given('tau_d', 'near_bed', 'reference_height_floor'))

    if 'bedload' in sed:
        domain.set_bedload(
            formula=sed['bedload'],
            K=sed.get('bedload_K'), m=sed.get('bedload_m'),
            tau_c_star=sed.get('bedload_tau_c_star'),
            open_boundaries=sed.get('bedload_open_boundaries'))

    if 'angle_of_repose' in sed:
        domain.set_angle_of_repose(
            sed['angle_of_repose'],
            relax=sed.get('repose_relax', 1.0),
            max_sweeps=sed.get('repose_max_sweeps', 50))

    if 'erodible_base_elevation' in sed:
        domain.set_erodible_base(elevation=sed['erodible_base_elevation'])
    elif 'erodible_base_depth' in sed:
        domain.set_erodible_base(depth=sed['erodible_base_depth'])

    # Regions: polygon points were resolved by prepare_data; center/radius
    # pass through, as for the erosion operators.
    for r in sed.get('erodible_regions', []):
        kwargs = {'erodible': r.get('erodible', True)}
        if r.get('polygon_points') is not None:
            kwargs['polygon'] = r['polygon_points']
        elif r.get('center') is not None:
            kwargs['center'] = r['center']
            kwargs['radius'] = r['radius']
        domain.set_erodible_region(**kwargs)

    for fr in sed['fractions']:
        kwargs = {k: fr[k] for k in ('rho_s', 'tau_c_star', 'd_star',
                                     'initial_concentration',
                                     'reference_height', 'nu', 'C1', 'C2')
                  if k in fr}
        domain.add_sediment_fraction(fr['name'], fr['diameter'], **kwargs)
        for tag, value in fr.get('boundary', {}).items():
            domain.set_tracer_boundary(fr['name'], tag, value)

    return operator
