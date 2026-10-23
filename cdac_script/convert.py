#!/usr/bin/env python
"""
Convert .tsh (ASCII) mesh files to .msh (NetCDF binary) format

Binary .msh files load 10-20x faster than ASCII .tsh files.

Usage:
    python convert_tsh_to_msh.py input.tsh output.msh
    python convert_tsh_to_msh.py input.tsh  # Creates input.msh
"""

import sys
import os
import time

def convert_tsh_to_msh(tsh_file, msh_file=None):
    """Convert .tsh to .msh format"""

    if msh_file is None:
        # Auto-generate .msh filename
        base = os.path.splitext(tsh_file)[0]
        msh_file = base + '.msh'

    print("=" * 70)
    print("Converting .tsh (ASCII) to .msh (NetCDF binary)")
    print("=" * 70)
    print(f"Input:  {tsh_file}")
    print(f"Output: {msh_file}")

    # Check input exists
    if not os.path.exists(tsh_file):
        print(f"Error: Input file not found: {tsh_file}")
        sys.exit(1)

    # Import ANUGA
    from anuga.load_mesh.loadASCII import import_mesh_file, export_mesh_file

    # Get input file size
    input_size = os.path.getsize(tsh_file)
    print(f"\nInput file size: {input_size / 1e9:.2f} GB")

    # Load .tsh file
    print(f"\nLoading {tsh_file}...")
    t0 = time.time()
    mesh_dict = import_mesh_file(tsh_file)
    t_load = time.time() - t0

    n_triangles = len(mesh_dict['triangles'])
    n_vertices = len(mesh_dict['vertices'])

    print(f"  Loaded in {t_load:.2f}s")
    print(f"  Triangles: {n_triangles:,}")
    print(f"  Vertices:  {n_vertices:,}")
    print(f"  Loading rate: {n_triangles/t_load/1e6:.2f} M triangles/sec")

    # Export to .msh
    print(f"\nWriting {msh_file}...")
    t0 = time.time()
    export_mesh_file(msh_file, mesh_dict)
    t_write = time.time() - t0

    output_size = os.path.getsize(msh_file)

    print(f"  Written in {t_write:.2f}s")
    print(f"  Output file size: {output_size / 1e9:.2f} GB")
    print(f"  Compression: {100*(1 - output_size/input_size):.1f}%")
    print(f"  Writing rate: {n_triangles/t_write/1e6:.2f} M triangles/sec")

    # Verify by loading .msh
    print(f"\nVerifying {msh_file}...")
    t0 = time.time()
    test_dict = import_mesh_file(msh_file)
    t_verify = time.time() - t0

    print(f"  Loaded .msh in {t_verify:.2f}s")
    print(f"  Loading rate: {n_triangles/t_verify/1e6:.2f} M triangles/sec")
    print(f"  Speedup vs .tsh: {t_load/t_verify:.1f}x faster!")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f".tsh load time:  {t_load:8.2f}s ({t_load/60:6.2f} min)")
    print(f".msh load time:  {t_verify:8.2f}s ({t_verify/60:6.2f} min)")
    print(f"Speedup:         {t_load/t_verify:8.1f}x")
    print(f"File size:       {input_size/1e9:8.2f} GB → {output_size/1e9:.2f} GB ({100*(1-output_size/input_size):.1f}% smaller)")
    print("=" * 70)
    print(f"\nConversion successful!")
    print(f"You can now use '{msh_file}' instead of '{tsh_file}'")
    print()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    tsh_file = sys.argv[1]
    msh_file = sys.argv[2] if len(sys.argv) > 2 else None

    convert_tsh_to_msh(tsh_file, msh_file)
