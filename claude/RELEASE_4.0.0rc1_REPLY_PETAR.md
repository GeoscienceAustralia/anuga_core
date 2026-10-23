<!--
Reply to Petar during the 4.0.0rc1 soak, 2026-08-23. He asked what mesh size
was used for the Towradgi catchment and its refinement areas, so he could
reproduce the #229 numbers exactly.

Kept because the reproduction recipe is the thing anyone re-checking the #229
Towradgi figures will need: the mesh parameters, the exact commands, and the
two traps (the 0.32 m figure is old-vs-new write-back, NOT CPU-vs-GPU; and
concurrent runs corrupt the shared SWW, issue #232).

Numbers verified against run_small_towradgi.py and the SWW at the time of
writing: 256,688 triangles / 128,539 points at --scale 1.0.
See claude/archive/RELEASE_4.0.0rc1_ANNOUNCEMENT.md for the announcement it follows.
-->

Hi Petar,

Everything was the script defaults — `validation_tests/case_studies/towradgi/run_small_towradgi.py`
with `--scale 1.0`, so:

  bounding polygon              maximum_triangle_area = 1000 m^2
  Model/Bdy/Catchment.csv       scale * 100.0  =  100 m^2
  Model/Bdy/FineCatchment.csv   scale * 36.0   =   36 m^2
  Model/Bdy/CreekBanks.csv                        8 m^2   (fixed, not scaled)

That gives 256,688 triangles / 128,539 points, which matches the script's own
help text for scale=1.

One thing to check before you start: you asked about the two refinement areas,
and the version I ran has three — CreekBanks is the extra one, and note it is
hardcoded at 8 m^2 rather than multiplied by scale, so it stays put at coarse or
super-fine settings. If your copy only has Catchment and FineCatchment then our
meshes differ regardless of the other numbers, and that is worth sorting out
first.

The exact commands were:

    python run_small_towradgi.py -ft 3600 -ys 600 -mpm 1
    python run_small_towradgi.py -ft 3600 -ys 600 -mpm 2

so 1 hour of simulated time, 600 s yieldstep, DE1, everything else default.

Two things that will save you time:

1. The 0.32 m figure is NOT a CPU-vs-GPU comparison. It is the new structure
   write-back against the OLD one, both in mode 1 — I patched Inlet.set_average_depth
   back to the old uniform-depth behaviour in-process to get the "before" arm.
   Comparing mode 1 against mode 2 on the same code gives something much smaller:
   3.5e-02 m max, 99th percentile 0.54 mm, total water agreeing to 3.4e-06
   relative. So if you just run current develop twice you will not see 0.32 m —
   you need the old write-back to compare against.

2. Give each run its own output directory. domain_name is hardcoded in the
   script, so two runs both writing MODEL_OUTPUTS/Towradgi_historic_flood.sww
   interleave their frames into one file and produce a non-monotonic time axis.
   It cost me an entire invalid result before I noticed — the numbers looked
   plausible. It is issue #232. Worth asserting the time array is increasing
   before you difference anything.

For what it's worth, the reason Towradgi shows the largest effect we have seen
is the terrain: bed elevation spread across the 44 inlets of its 22 structures
has a median of 0.77 m, and the worst is 4.95 m across 12 triangles at Collins
St. The old code wrote a uniform depth across each inlet, which tilts a level
water surface onto the bed by roughly half that range, every timestep, even at
zero discharge. Flat-bed inlets are unaffected — three of our four flat-bed
validation cases are bit-identical before and after.

I think the new behaviour is right, but you know that catchment far better than
I do, so if the changed result looks wrong to you please say so — that would
reopen the decision rather than just the wording.

Cheers,
Steve
