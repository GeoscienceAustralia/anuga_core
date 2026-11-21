"""Put the NVIDIA HPC SDK run-time libraries into an nvc-built ANUGA wheel.

Usage: bundle_runtime.py RAW_WHEEL REDIST_LIB_DIR SDK_LICENSE NOTICE OUT_DIR

The wheel that nvc builds needs libnvomp, libacchost and the rest of the SDK's
OpenMP run time. auditwheel would copy them in, but it renames each one and
rewrites its ELF header, and the SDK licence gives no right to modify its
run-time files. This copies them byte for byte instead, into
anuga/_nvhpc_runtime/, and changes only ANUGA's own extension modules, which
get an RPATH pointing there.

Only files under a REDIST directory of the SDK may be distributed, so a
library found anywhere else stops the build.
"""

import hashlib
import pathlib
import shutil
import subprocess
import sys
import tempfile

RUNTIME_DIR = "_nvhpc_runtime"
# glibc and the gcc run time: on every Linux, never bundled.
SYSTEM = (
    "libc.so",
    "libm.so",
    "libdl.so",
    "libpthread.so",
    "librt.so",
    "libgcc_s.so",
    "libstdc++.so",
    "libgomp.so",
    "ld-linux",
    "libutil.so",
    "libatomic.so",
)


def run(*cmd):
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout


def needed(path):
    return run("patchelf", "--print-needed", str(path)).split()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(raw_wheel, redist_dir, sdk_license, notice, out_dir):
    redist = pathlib.Path(redist_dir)
    if "REDIST" not in redist.parts:
        sys.exit("%s is not under a REDIST directory of the SDK" % redist)
    out = pathlib.Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    work = pathlib.Path(tempfile.mkdtemp())
    run(sys.executable, "-m", "wheel", "unpack", "-d", str(work), raw_wheel)
    (tree,) = work.iterdir()
    package = tree / "anuga"
    runtime = package / RUNTIME_DIR
    runtime.mkdir()

    # The SDK libraries the extension modules need, and those they need.
    extensions = sorted(package.rglob("*.so"))
    wanted, queue = set(), []
    for ext in extensions:
        queue += needed(ext)
    while queue:
        name = queue.pop()
        if name in wanted or name.startswith(SYSTEM):
            continue
        source = redist / name
        if not source.exists():
            sys.exit("%s is needed but is not in %s: it may not be distributed" % (name, redist))
        wanted.add(name)
        queue += needed(source)

    for name in sorted(wanted):
        source = (redist / name).resolve()
        shutil.copyfile(source, runtime / name)
        if sha256(runtime / name) != sha256(source):
            sys.exit("%s changed while copying" % name)
        print("bundled %-28s %5.1f MB" % (name, source.stat().st_size / 1e6))

    shutil.copyfile(sdk_license, runtime / "LICENSE_NVIDIA_HPC_SDK.txt")
    shutil.copyfile(notice, runtime / "NOTICE_NVIDIA_RUNTIME.txt")

    # DT_RPATH, not DT_RUNPATH: it also applies to what the SDK libraries
    # load themselves (libnvomp needs libnvc, libnvc needs libnvcpumath).
    for ext in extensions:
        up = "/".join([".."] * (len(ext.relative_to(package).parts) - 1))
        rpath = "$ORIGIN/" + (up + "/" if up else "") + RUNTIME_DIR
        run("patchelf", "--force-rpath", "--set-rpath", rpath, str(ext))
    print("set RPATH on %d extension modules" % len(extensions))

    run(sys.executable, "-m", "wheel", "pack", "-d", str(out), str(tree))
    (wheel,) = out.glob("*.whl")
    print(
        run(
            sys.executable, "-m", "wheel", "tags", "--remove", "--platform-tag", "manylinux_2_28_x86_64", str(wheel)
        ).strip()
    )


if __name__ == "__main__":
    main(*sys.argv[1:])
