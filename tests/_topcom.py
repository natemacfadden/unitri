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
"""Shared TOPCOM (via CYTools) enumeration for the cross-check tests.

Both cross-checks -- check_topcom.py (fixed convex regions vs the GMP CLI) and
test_topcom_convex.py (random convex point sets vs count_triangulations) --
enumerate the same fine/primitive triangulations of a CYTools Polytope.  Keeping
that one enumeration (and its flag set) here means the two checks can't drift.
"""


def count_fine_triangulations(P, cap):
    """Number of fine triangulations of CYTools Polytope ``P`` via TOPCOM, or
    None if the count exceeds ``cap``.  The flags select exactly the fine
    (primitive) triangulations that unitri counts."""
    c = 0
    for _ in P.all_triangulations(only_fine=True, only_regular=False,
                                  only_star=False,
                                  include_points_interior_to_facets=True):
        c += 1
        if c > cap:
            return None
    return c
