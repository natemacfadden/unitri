# =============================================================================
#    Copyright (C) 2026  Nate MacFadden for the Liam McAllister Group
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
# =============================================================================
"""Shared fixtures for the unitri test suite."""
import os
import subprocess

import pytest

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
