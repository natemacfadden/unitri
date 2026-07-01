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
"""Exact counts for deliberate, named *general* convex polygons.

The randomized fuzzer (test_topcom_convex.py) covers random convex point sets;
this pins specific non-rectangular shapes -- trapezoids, pentagons (one steep),
hexagons -- including two whose hull-tracing produces absent boundary vertices
(so the absent-mask path in na_query_compute is exercised, not just clean
rectangles). Convex polygons are the counter's least-tested regime.

Each expected count was verified against TOPCOM (CYTools); when CYTools is
installed, test_convex_polygon_matches_topcom re-verifies it live.
"""
import pytest

import unitri

try:
    from cytools import Polytope
    HAS_CYTOOLS = True
except ImportError:
    HAS_CYTOOLS = False

# (name, CCW vertices, TOPCOM-verified fine-triangulation count)
POLYGONS = [
    ("trapezoid",              [(0, 0), (4, 0), (3, 2), (1, 2)],                 140),
    ("hexagon (absent verts)", [(1, 0), (2, 0), (3, 1), (2, 2), (1, 2), (0, 1)], 24),
    ("wide trapezoid",         [(0, 0), (5, 0), (4, 2), (1, 2)],                 2046),
    ("steep pentagon",         [(0, 0), (2, 0), (3, 2), (2, 4), (0, 2)],         1392),
    ("hexagon",                [(0, 0), (3, 0), (4, 1), (4, 2), (2, 3), (0, 2)], 12135),
    ("asym pentagon (absent)", [(0, 0), (4, 0), (3, 3), (1, 3), (0, 1)],         10843),
]


def _lattice_points(verts):
    """Integer points inside or on the convex polygon `verts` (either winding).
    Matches CYTools' Polytope(verts).points() for these shapes, but needs no
    CYTools, so the count assertions run in CI."""
    xs = [x for x, _ in verts]
    ys = [y for _, y in verts]
    k = len(verts)

    def inside(px, py):
        signs = set()
        for i in range(k):
            x0, y0 = verts[i]
            x1, y1 = verts[(i + 1) % k]
            cr = (x1 - x0) * (py - y0) - (y1 - y0) * (px - x0)
            if cr:
                signs.add(cr > 0)
        return len(signs) <= 1

    return [(x, y) for x in range(min(xs), max(xs) + 1)
            for y in range(min(ys), max(ys) + 1) if inside(x, y)]


@pytest.mark.parametrize("name,verts,expected", POLYGONS, ids=[p[0] for p in POLYGONS])
def test_convex_polygon_count(name, verts, expected):
    assert unitri.count_triangulations(_lattice_points(verts)) == expected


@pytest.mark.skipif(not HAS_CYTOOLS, reason="cytools not installed")
@pytest.mark.parametrize("name,verts,expected", POLYGONS, ids=[p[0] for p in POLYGONS])
def test_convex_polygon_matches_topcom(name, verts, expected):
    # keep the pinned counts honest: re-enumerate with TOPCOM when it's available
    from _topcom import count_fine_triangulations
    assert count_fine_triangulations(Polytope(verts), cap=50000) == expected
