"""GMP-toolchain discovery, shared by `setup.py` (the Cython binding) and the
tests (the standalone oracle binary), so the build and the tests can't drift.

GMP is often off the compiler's default search path.  We locate it with
`pkg-config` (the standard, cross-platform mechanism; conda-forge and Homebrew
both ship a `gmp.pc`).  Activating a conda env does not by itself put
`$CONDA_PREFIX/lib/pkgconfig` on `PKG_CONFIG_PATH`, so we add it before
querying.  If `pkg-config` or the `.pc` is unavailable we fall back to probing
Homebrew's prefix on macOS and `$CONDA_PREFIX` directly.

Two views of the same discovery:
  * gmp_dirs()   -> (include_dirs, library_dirs) lists for setuptools'
                    Extension(...); empty when GMP is on the default path.
  * gmp_cflags() -> gcc flags (-I/-L plus an rpath) for command-line builds:
                    the test oracle and a raw `gcc ... unitri/na-query.c`.

Run as a script it prints `gmp_cflags()`, so a Makefile or shell build can
splice the flags in:

    gcc -O2 $(python3 _gmp.py) -DGMP -o na-query unitri/na-query.c -lgmp
"""
import os
import subprocess


def _pkg_config_env():
    """Environment for `pkg-config` with the active conda env's pkgconfig dir
    prepended to `PKG_CONFIG_PATH` (conda activation does not set this)."""
    env = os.environ.copy()
    conda = env.get("CONDA_PREFIX")
    if conda:
        pcdir = os.path.join(conda, "lib", "pkgconfig")
        env["PKG_CONFIG_PATH"] = os.pathsep.join(
            p for p in (pcdir, env.get("PKG_CONFIG_PATH", "")) if p)
    return env


def _pkg_config_dirs():
    """(include_dirs, library_dirs) from `pkg-config`, or None if it cannot
    locate GMP (or reports no special dirs -- i.e. GMP is on the default
    path, which the prefix probes below also resolve to empty)."""
    try:
        out = subprocess.check_output(
            ["pkg-config", "--cflags-only-I", "--libs-only-L", "gmp"],
            stderr=subprocess.DEVNULL, env=_pkg_config_env(),
        ).decode().split()
    except (OSError, subprocess.CalledProcessError):
        return None
    inc = [f[2:] for f in out if f.startswith("-I")]
    lib = [f[2:] for f in out if f.startswith("-L")]
    if not (inc or lib):
        return None
    return inc, lib


def _brew_gmp_prefix():
    """Homebrew's GMP prefix on macOS, or None."""
    try:
        prefix = subprocess.check_output(
            ["brew", "--prefix", "gmp"], stderr=subprocess.DEVNULL
        ).decode().strip()
    except (OSError, subprocess.CalledProcessError):
        return None
    return prefix or None


def _conda_gmp_prefix():
    """`$CONDA_PREFIX` if an active conda env actually has GMP, else None."""
    prefix = os.environ.get("CONDA_PREFIX")
    if prefix and os.path.exists(os.path.join(prefix, "include", "gmp.h")):
        return prefix
    return None


def gmp_dirs():
    """(include_dirs, library_dirs) locating GMP -- each a list, possibly empty
    when GMP is already on the default toolchain path."""
    found = _pkg_config_dirs()
    if found is not None:
        return found
    prefix = _brew_gmp_prefix() or _conda_gmp_prefix()
    if prefix:
        return [os.path.join(prefix, "include")], [os.path.join(prefix, "lib")]
    return [], []


def gmp_cflags():
    """gcc flags locating GMP, to be placed *before* the source on the command
    line (`-lgmp` is added separately by the caller).  Includes an rpath so the
    binary finds libgmp at run time without `LD_LIBRARY_PATH`.  Empty when GMP
    is on the default path."""
    inc, lib = gmp_dirs()
    flags = [f"-I{d}" for d in inc] + [f"-L{d}" for d in lib]
    flags += [f"-Wl,-rpath,{d}" for d in lib]
    return flags


if __name__ == "__main__":
    print(" ".join(gmp_cflags()))
