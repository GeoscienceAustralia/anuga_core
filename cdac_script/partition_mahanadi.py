#!/usr/bin/env python3
"""
Partition Mahanadi Delta mesh for parallel execution

Usage:
    python partition_mahanadi.py --mesh-size 300 --nprocs 16
    python partition_mahanadi.py --mesh-size 100 --nprocs 32
    python partition_mahanadi.py --mesh-size 900 --nprocs 8
"""

import argparse
import time
import numpy as np
import linecache
import sys

def main():
    parser = argparse.ArgumentParser(description='Partition Mahanadi Delta mesh')
    parser.add_argument('--mesh-size', type=int, required=True, choices=[100, 300, 900],
                       help='Mesh resolution: 100, 300, or 900 sqm')
    parser.add_argument('--nprocs', type=int, required=True,
                       help='Number of partitions to create')
    parser.add_argument('--output', default='./partitions',
                       help='Output directory for partition files')
    parser.add_argument('--chilka-stage', type=float, default=0.15,
                       help='Initial stage inside Chilka polygon (default: 0.15)')
    parser.add_argument('--verbose', action='store_true',
                       help='Verbose output')

    args = parser.parse_args()

    # Construct filenames based on mesh size
    #base_name = f'delta11372sqkm_uniform_mesh_{args.mesh_size}sqm_Chilka{args.mesh_size}sqm'
    base_name = f'{args.mesh_size}sqm'
    mesh_filename = f'mesh_file/{base_name}.msh'  # Prefer .msh over .tsh for speed

    # Try .msh first, fall back to .tsh
    import os
    if not os.path.exists(mesh_filename):
        mesh_filename = f'mesh_file/{base_name}.tsh'
        if not os.path.exists(mesh_filename):
            print(f"Error: Mesh file not found: {mesh_filename}")
            print(f"Looking for: mesh_file/{base_name}.{{msh,tsh}}")
            sys.exit(1)

    # Data files (assuming same naming convention)
    friction_filename = 'friction_data/LULC_Apr2020_Norm.asc'
    elevation_filename = f'elevation_data/Elev_delta11372sqkm_uniform_mesh_300sqm_Chilka_300sqm_LIDAR1m_ALOS30m.csv'
    chilka_poly_file = 'polygons/Chilka_Modified_Poly.csv'

    # Check files exist
    for f in [friction_filename, elevation_filename]:
        if not os.path.exists(f):
            print(f"Warning: File not found: {f}")
            print(f"Will skip this quantity during partitioning")

    print("=" * 70)
    print(f"Partitioning Mahanadi Delta Mesh - {args.mesh_size} sqm resolution")
    print("=" * 70)
    print(f"Mesh file:       {mesh_filename}")
    print(f"Elevation file:  {elevation_filename}")
    print(f"Friction file:   {friction_filename}")
    print(f"Partitions:      {args.nprocs}")
    print(f"Output dir:      {args.output}")
    print("=" * 70)

    # Import ANUGA
    print("\nImporting ANUGA...")
    t_start = time.time()
    from mesh_cache import create_domain_from_file
    from anuga.parallel.sequential_distribute import Sequential_distribute
    import anuga
    print(f"  Done in {time.time() - t_start:.2f}s")

    # Step 1: Create domain from mesh
    print(f"\nStep 1: Loading mesh from {mesh_filename}...")
    t0 = time.time()
    domain = create_domain_from_file(mesh_filename)
    t_load = time.time() - t0

    n_triangles = len(domain)
    print(f"  Loaded {n_triangles:,} triangles in {t_load:.2f}s")

    # Step 2: Set domain parameters
    print("\nStep 2: Setting domain parameters...")
    domain.set_name(f'mahanadi_delta_{args.mesh_size}sqm')
    domain.set_flow_algorithm('DE1')
    domain.set_CFL(2.0)
    domain.set_minimum_allowed_height(0.008)
    domain.set_maximum_allowed_speed(1.0)

    # Step 3: Set friction from raster
    if os.path.exists(friction_filename):
        print(f"\nStep 3: Setting friction from {friction_filename}...")
        t0 = time.time()
        domain.set_quantity('friction',
                          filename=friction_filename,
                          location='centroids')
        print(f"  Done in {time.time() - t0:.2f}s")
    else:
        print(f"\nStep 3: Skipping friction (file not found)")

    # Step 4: Set elevation from CSV (custom per-vertex loading)
    if os.path.exists(elevation_filename):
        print(f"\nStep 4: Loading elevation from {elevation_filename}...")
        print(f"  Reading {n_triangles:,} triangles × 3 vertices = {n_triangles*3:,} values")
        t0 = time.time()

        elev_per_node = np.loadtxt(elevation_filename)
        triangle_index_elvation_main = elev_per_node[domain.triangles].tolist()
        triangle_index = list(range(n_triangles))

        t_read = time.time() - t0
        print(f"  Read elevation data in {t_read:.2f}s")

        # Set quantity
        t0 = time.time()
        domain.set_quantity('elevation',
                          numeric=triangle_index_elvation_main,
                          use_cache=False,
                          verbose=args.verbose,
                          alpha=0.1,
                          indices=triangle_index,
                          location='vertices')
        t_set = time.time() - t0

        print(f"  Set elevation quantity in {t_set:.2f}s")
        print(f"  Total elevation time: {t_read + t_set:.2f}s")

        # Clean up
        triangle_index = []
        triangle_index_elvation_main = []
        linecache.clearcache()
    else:
        print(f"\nStep 4: Skipping elevation (file not found)")

    # Step 5: Set initial stage
    print(f"\nStep 5: Setting initial stage...")
    t0 = time.time()

    # Set stage = elevation everywhere (dry bed)
    domain.set_quantity('stage', expression='elevation')

    # Add water inside Chilka polygon if polygon file exists
    if os.path.exists(chilka_poly_file):
        poly = anuga.read_polygon(chilka_poly_file)
        domain.add_quantity('stage', numeric=args.chilka_stage, polygon=poly)
        print(f"  Set Chilka stage = {args.chilka_stage} m")
    else:
        print(f"  Warning: Chilka polygon not found: {chilka_poly_file}")

    print(f"  Done in {time.time() - t0:.2f}s")

    # Step 6: Partition mesh
    print(f"\nStep 6: Partitioning mesh into {args.nprocs} parts...")
    print(f"  Using METIS graph partitioner...")

    t0 = time.time()
    partition = Sequential_distribute(domain, verbose=args.verbose, debug=False)
    partition.distribute(args.nprocs)
    t_partition = time.time() - t0

    print(f"\n  Partitioning complete!")
    print(f"    Time: {t_partition:.2f}s ({n_triangles/t_partition/1e6:.2f} M triangles/sec)")

    # Show partition distribution
    triangles_per_proc = partition.triangles_per_proc
    print(f"\n  Partition distribution:")
    for p in range(min(args.nprocs, 10)):
        pct = 100.0 * triangles_per_proc[p] / n_triangles
        print(f"    Rank {p:3d}: {triangles_per_proc[p]:8,} triangles ({pct:5.2f}%)")
    if args.nprocs > 10:
        print(f"    ... ({args.nprocs - 10} more ranks)")

    # Check balance
    min_tri = triangles_per_proc.min()
    max_tri = triangles_per_proc.max()
    avg_tri = n_triangles / args.nprocs
    imbalance = (max_tri - min_tri) / avg_tri * 100
    print(f"\n  Load balance:")
    print(f"    Min: {min_tri:,} triangles ({100*min_tri/avg_tri:.1f}% of average)")
    print(f"    Max: {max_tri:,} triangles ({100*max_tri/avg_tri:.1f}% of average)")
    print(f"    Imbalance: {imbalance:.2f}%")

    # Step 7: Save partitions
    print(f"\nStep 7: Saving partition files to {args.output}...")
    t0 = time.time()

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    # Save partitions manually (same approach as benchmark_partition.py)
    import pickle
    from os.path import join

    simulation_name = f'mahanadi_delta_{args.mesh_size}sqm'

    for p in range(args.nprocs):
        if args.verbose or args.nprocs <= 20 or p % max(1, args.nprocs // 10) == 0:
            print(f"  Saving partition {p}/{args.nprocs-1}...", end='', flush=True)

        tostore = partition.extract_submesh(p)

        pickle_name = simulation_name + f'_P{args.nprocs}_{p}.pickle'
        pickle_name = join(args.output, pickle_name)

        with open(pickle_name, 'wb') as f:
            lst = list(tostore)

            # Write points and triangles to their own files
            np.save(pickle_name + ".np1", tostore[1])
            lst[1] = pickle_name + ".np1.npy"
            np.save(pickle_name + ".np2", tostore[2])
            lst[2] = pickle_name + ".np2.npy"

            # Write each quantity to its own file
            for k in tostore[4]:
                np.save(pickle_name + ".np4." + k, np.array(tostore[4][k]))
                lst[4][k] = pickle_name + ".np4." + k + ".npy"

            pickle.dump(tuple(lst), f, protocol=pickle.HIGHEST_PROTOCOL)

        if args.verbose or args.nprocs <= 20 or p % max(1, args.nprocs // 10) == 0:
            print(" done")

    t_save = time.time() - t0
    print(f"\n  All partitions saved in {t_save:.2f}s ({t_save/args.nprocs:.2f}s per partition)")

    # Summary
    t_total = time.time() - t_start
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Mesh loading:       {t_load:8.2f}s")
    print(f"Partitioning:       {t_partition:8.2f}s")
    print(f"Saving partitions:  {t_save:8.2f}s")
    print(f"-" * 70)
    print(f"TOTAL TIME:         {t_total:8.2f}s ({t_total/60:6.2f} min)")
    print("=" * 70)

    print(f"\nPartition files created for {args.nprocs} MPI ranks")
    print(f"\nTo use in run_model_3.py, replace domain creation with:")
    print(f"  domain = anuga.sequential_distribute_load('{simulation_name}', partition_dir='{args.output}')")
    print()

if __name__ == '__main__':
    main()
