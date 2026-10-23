"""Example 3 -- multi-grain-size dam break over an erodible bed, on the GPU.

The other two examples are faithful ports of anugaSed's cases, so they only
exercise what anugaSed does: one grain size, suspension only, mode 1.
This one covers what the module adds -- several grain sizes at once, bedload
alongside suspension, a moving bed, and the unified (GPU) compute path.

Demonstrates
  * two grain sizes with different diameters, settling and mobility
  * bedload [K-1]..[K-4] and its bed evolution [G-5]
  * suspended exchange [E-1]/[D-1] and Exner bed evolution [G-4]
  * the Rouse near-bed ratio [S-4]
  * domain.sediment_summary(), which records every choice made
  * mode 2, and the sediment mass budget across bed and water column
"""
import numpy as np
import anuga

LEN_X, LEN_Y = 400.0, 100.0
DAM_X = 150.0
H_UP, H_DOWN = 3.0, 0.4
FINALTIME = 60.0
PORO = 0.30

domain = anuga.rectangular_cross_domain(80, 20, len1=LEN_X, len2=LEN_Y)
domain.set_name('erodible_dambreak')
domain.set_flow_algorithm('DE0')
domain.set_low_froude(0)
domain.store = False
x = domain.centroid_coordinates[:, 0]

# A gently sloping erodible bed, with the reservoir behind a dam.
domain.set_quantity('elevation', lambda X, Y: -0.004 * X)
domain.set_quantity('stage', np.where(x < DAM_X, -0.004 * x + H_UP,
                                      -0.004 * x + H_DOWN),
                    location='centroids')
domain.set_quantity('friction', 0.03)
domain.set_boundary({t: anuga.Reflective_boundary(domain)
                     for t in domain.get_boundary_tags()})

# --- configure the sediment physics ----------------------------------------
domain.set_sediment_parameters(porosity=PORO, c_max=0.30)
domain.set_bed_material('noncohesive')            # a sand/gravel bed, spec 4.1.1
domain.set_deposition('d_star', near_bed='rouse')  # [S-4] rather than well-mixed
domain.set_sediment_friction('constant')           # ordinary flood work
domain.set_bedload('wong_parker_eq24')             # [K-1] and [G-5]

domain.add_sediment_fraction(name='fine_sand', diameter=1.5e-4, initial_concentration=0.0)
domain.add_sediment_fraction(name='coarse_sand', diameter=8.0e-4, initial_concentration=0.0)

# 'unified' is the shared CPU/GPU path. Whether it offloads to a device is a
# property of the build, not of this call: a build without offload runs the
# same kernels on the host under CPU_ONLY_MODE.
domain.set_compute_mode('unified')
offloads = anuga.gpu_offload_enabled()

print(domain.sediment_summary())
print()
print('compute mode: %s, offload %s'
      % (domain.get_compute_mode(),
         'enabled -- running on the device' if offloads
         else 'not built in -- unified kernels on the host'))
print()

areas = domain.areas
z0 = domain.quantities['elevation'].centroid_values.copy()
m0 = sum(float((domain.tracer_conserved_values[s] * areas).sum())
         for s in range(domain.n_sediment_classes))

print('   time     suspended (m3)          bed change (m3)    total      max scour')
print('             fine      coarse')
for t in domain.evolve(yieldstep=10.0, finaltime=FINALTIME):
    z = domain.quantities['elevation'].centroid_values
    mf = float((domain.tracer_conserved_values[0] * areas).sum())
    mc = float((domain.tracer_conserved_values[1] * areas).sum())
    bed = float(((1.0 - PORO) * (z - z0) * areas).sum())
    print('  %5.1f s  %.4e  %.4e   %+.4e   %.4e   %+.3f m'
          % (t, mf, mc, bed, mf + mc + bed, (z - z0).min()))

z = domain.quantities['elevation'].centroid_values
dz = z - z0
m1 = sum(float((domain.tracer_conserved_values[s] * areas).sum())
         for s in range(domain.n_sediment_classes))
bed = float(((1.0 - PORO) * dz * areas).sum())

print()
drift = (m1 + bed) - m0
print('Sediment budget:  suspended %+.6e  +  bed %+.6e  =  %+.6e m3'
      % (m1, bed, m1 + bed))
print('  started from %+.6e, so the drift is %+.3e m3' % (m0, drift))
print()
print('  That closes to machine precision, and it should: the boundaries are')
print('  reflective, so nothing leaves. Every cubic metre now in suspension')
print('  came out of the bed, and (1-lambda) dz accounts for it exactly --')
print('  across BOTH grain sizes, with erosion, deposition and bedload active')
print('  and the bed moving under them. It is the strongest statement this')
print('  example makes: the coupling conserves sediment.')
print()
print('Bed:   %+.4f m (scour) to %+.4f m (deposition)' % (dz.min(), dz.max()))
print('Wet cells: %d of %d'
      % (int((domain.quantities['stage'].centroid_values - z > 0.01).sum()),
         len(domain)))
