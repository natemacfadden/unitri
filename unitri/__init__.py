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
"""unitri -- counting the fine (unimodular) triangulations of a lattice polygon.

    count_triangulations(points)          exact count (compiled GMP extension)
    count_triangulations_parallel(points) exact count, GMP-free: parallel
                                          mod-prime runs + CRT across cores

Both take a lattice point set.  na_query(m, n, upper, lower) is the low-level
profile API (GMP extension).  The package imports without the GMP extension, so
the mod-prime/parallel path works even when libgmp is absent.
"""
from .crt_parallel import count_parallel, count_triangulations_parallel
from .profiles import count_triangulations, points_to_profiles

try:
    # Binds na_query to the FUNCTION as unitri.na_query (this rebinds the
    # same-named compiled submodule).  The extension needs libgmp; without it
    # the import fails, `import unitri` still succeeds, and __getattr__ below
    # gives a pointed message if na_query is then used.
    from .na_query import na_query
except ImportError:
    pass

__all__ = [
    "count_triangulations",
    "count_triangulations_parallel",
    "count_parallel",
    "na_query",
    "points_to_profiles",
]


def __getattr__(name):
    if name == "na_query":     # only reached when the GMP extension didn't import
        raise ImportError(
            "unitri.na_query needs the compiled GMP extension (build with "
            "`pip install -e .`, which needs libgmp).  Without GMP, use "
            "unitri.count_triangulations_parallel (mod-prime + CRT)."
        )
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
