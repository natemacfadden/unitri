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
"""Unimodular-invariance helpers shared by the test suites.

The number of (fine) triangulations of a lattice point set is invariant under
the affine unimodular group: any map p -> M @ p + t with M in GL(2, Z)
(det M = +-1) and integer t carries the lattice to itself and triangulations to
triangulations, so the count must not change.  count_triangulations chooses its
own orientation internally, so feeding it several unimodular images of the same
set is a strong check that the orientation handling (and the x<->m-x reflection
fix) is correct.

Each M below has det +-1; verify with `assert _det(M) in (1, -1)`.
"""

# (M = (a, b, c, d) acting as [[a, b], [c, d]], translation (tx, ty)); det = +-1
UNIMODULAR = [
    ((1, 0, 0, 1),   (0, 0)),     # identity
    ((1, 1, 0, 1),   (0, 0)),     # shear in x
    ((1, 0, 1, 1),   (0, 0)),     # shear in y
    ((0, -1, 1, 0),  (3, -2)),    # 90-degree rotation (det +1) + shift
    ((-1, 0, 0, 1),  (5, 0)),     # reflection (det -1) + shift
    ((2, 1, 1, 1),   (-1, 4)),    # det +1, non-trivial
    ((1, -2, 0, 1),  (0, 7)),     # shear + shift
    ((-1, -1, 0, -1), (2, 2)),    # det +1
]

# The dihedral symmetries of the square (D4): det +-1 maps that PRESERVE the
# axis-aligned bounding-box dimensions, so the perpendicular extent n stays
# bounded.  Use these for random/large point sets where a general shear could
# blow up n (and hence na_query's ~ (n+2)^(m-1) cost); the full UNIMODULAR
# series above is for small fixed examples.
COMPACT = [
    ((1, 0, 0, 1),   (0, 0)),     # identity
    ((0, -1, 1, 0),  (0, 0)),     # rotate 90
    ((-1, 0, 0, -1), (0, 0)),     # rotate 180
    ((0, 1, -1, 0),  (0, 0)),     # rotate 270
    ((-1, 0, 0, 1),  (0, 0)),     # flip x
    ((1, 0, 0, -1),  (0, 0)),     # flip y
    ((0, 1, 1, 0),   (0, 0)),     # transpose
    ((0, -1, -1, 0), (0, 0)),     # anti-transpose
]


def _det(M):
    a, b, c, d = M
    return a * d - b * c


def apply(points, M, t):
    """Map every (x, y) by p -> M @ p + t (M = (a,b,c,d), t = (tx,ty))."""
    a, b, c, d = M
    tx, ty = t
    return [(a * x + b * y + tx, c * x + d * y + ty) for x, y in points]


def invariant_count(points, expected=None, transforms=UNIMODULAR):
    """Count `points` under every unimodular image in `transforms` and assert the
    counts all agree (and equal `expected`, if given).  Returns the common count.
    Raises AssertionError on any disagreement -- i.e. a count that depends on the
    orientation, which would be a correctness bug.
    """
    from unitri import count_triangulations
    counts = {(M, t): count_triangulations(apply(points, M, t))
              for M, t in transforms}
    vals = set(counts.values())
    assert len(vals) == 1, f"count not unimodular-invariant: {counts}"
    got = vals.pop()
    if expected is not None:
        assert got == expected, f"got {got}, expected {expected}"
    return got
