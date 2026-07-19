"""Shared fixtures for the unitri test suite."""
import os
import sys
import subprocess

import pytest

# the GMP-discovery helper lives at the repo root (shared with setup.py)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from _gmp import gmp_cflags

NA_QUERY_C = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "unitri", "na-query.c")


@pytest.fixture(scope="session")
def na_query_bin(tmp_path_factory):
    """Compile na-query.c with the exact big-integer (GMP) backend, once per
    session, and return the path to the binary.  The counting and symmetry
    suites query this directly so they pin exact counts rather than the
    default mod-prime backend."""
    out = str(tmp_path_factory.mktemp("unitri") / "na-query")
    subprocess.check_call(
        ["gcc", "-O2", *gmp_cflags(), "-DGMP", "-o", out, NA_QUERY_C, "-lgmp"])
    return out


@pytest.fixture(scope="session")
def na_query_mod_bin(tmp_path_factory):
    """Compile na-query.c with the default mod-prime backend (no -DGMP, hence no
    libgmp), once per session, and return the path to the binary.  The mod-prime
    suite queries this and reconstructs exact counts via crt_combine -- covering
    the default (most-compiled) build path that the GMP suite leaves untested."""
    out = str(tmp_path_factory.mktemp("unitri") / "na-query-mod")
    subprocess.check_call(["gcc", "-O2", "-o", out, NA_QUERY_C])
    return out
