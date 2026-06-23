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
"""Unimodular-invariance unit test (needs only the built unitri extension).

Each named point set is counted under a series of GL(2, Z)+translation images;
all images must give the same count, and it must equal the known value.

The first case is the regression that motivated the reflection fix: the region
under upper profile "1 . 2 2" over a flat floor (m=3, n=2).  In box coordinates
its profile has height_0 = 1 < height_m = 2, which lands in the half of the
table that the flat-floor recurrence does not store (it recovers that half by
the x<->m-x reflection).  Before the fix na_query reported "not_found" for this
orientation; count_triangulations now returns 134 in every orientation, matching
TOPCOM.
"""
import pytest

from transforms import invariant_count

# name, points, expected count (TOPCOM-verified)
CASES = [
    # region under "1 . 2 2" over flat floor: hull edge (0,1)->(2,2) skips a row
    ("'1 . 2 2' over flat floor (reflection regression)",
     [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2),
      (3, 0), (3, 1), (3, 2)], 134),
    # 3x2 rectangle, all 12 lattice points (check_topcom: 852)
    ("3x2 rectangle",
     [(x, y) for x in range(4) for y in range(3)], 852),
    # 2x2 square, all 9 lattice points
    ("2x2 square",
     [(x, y) for x in range(3) for y in range(3)], 64),
]


@pytest.mark.parametrize("name,pts,expected", CASES, ids=[c[0] for c in CASES])
def test_unimodular_invariance(name, pts, expected):
    # invariant_count asserts every GL(2, Z) image agrees and equals expected
    invariant_count(pts, expected)
