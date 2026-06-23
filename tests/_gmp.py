"""Shared GMP-toolchain helper for the test scripts.

On macOS (Homebrew) GMP is not on the default compiler search path, so the
test builds need its `-I`/`-L` from `brew --prefix gmp`. On Linux/conda GMP is
on the default path, so no extra flags are needed.
"""
import subprocess


def gmp_cflags():
    """`-I`/`-L` flags for Homebrew's GMP on macOS; an empty list elsewhere."""
    try:
        prefix = subprocess.check_output(
            ["brew", "--prefix", "gmp"], stderr=subprocess.DEVNULL
        ).decode().strip()
    except (OSError, subprocess.CalledProcessError):
        return []
    return [f"-I{prefix}/include", f"-L{prefix}/lib"]
