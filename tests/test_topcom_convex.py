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
"""Cross-check convex-region counts against TOPCOM (CYTools) -- convex polygons
are the counter's least-tested regime.  Two complementary checks:

  * curated, named convex polygons (trapezoids, pentagons, hexagons -- two with
    absent boundary vertices, exercising the absent-mask kernel path), with
    TOPCOM-verified counts; and
  * a randomized fuzzer over small convex point sets.

Both feed `count_triangulations` (which picks a hull-tracing orientation and
calls na_query) and require it to either match TOPCOM exactly or raise (sets
na_query can't represent without undercounting -- "long-diagonal" cases).  Each
decided set is also counted under the dihedral (D4) symmetries: the count must
be orientation-independent and equal TOPCOM.

The curated counts are asserted even without CYTools; the live TOPCOM
re-enumeration and the fuzzer are skipped when CYTools is absent.
"""
import random

import pytest

import unitri

try:
    from cytools import Polytope
    HAS_CYTOOLS = True
except ImportError:
    HAS_CYTOOLS = False

from _topcom import count_fine_triangulations
from transforms import COMPACT, invariant_count

# --- curated, named convex polygons -----------------------------------------
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


@pytest.mark.parametrize("name,verts,expected", POLYGONS, ids=[p[0] for p in POLYGONS])
def test_convex_polygon_count_parallel(name, verts, expected):
    # GMP-free path (parallel mod-prime + CRT) must agree with the exact count
    assert unitri.count_triangulations_parallel(_lattice_points(verts)) == expected


@pytest.mark.skipif(not HAS_CYTOOLS, reason="cytools not installed")
@pytest.mark.parametrize("name,verts,expected", POLYGONS, ids=[p[0] for p in POLYGONS])
def test_convex_polygon_matches_topcom(name, verts, expected):
    # keep the pinned counts honest: re-enumerate with TOPCOM when it's available
    assert count_fine_triangulations(Polytope(verts), cap=50000) == expected


# --- randomized fuzz over small convex point sets ---------------------------
SEED = 20260617
TRIALS = 6000
CHECK_LIMIT = 180    # stop after this many decided (matched/uncountable) cases
CAP = 8000           # skip regions TOPCOM can't enumerate quickly


@pytest.mark.skipif(not HAS_CYTOOLS, reason="cytools not installed")
def test_topcom_cross_check():
    rng = random.Random(SEED)
    matched = uncountable = 0
    problems = []
    for _ in range(TRIALS):
        coord = rng.choice([3, 4, 5, 6])
        cloud = {(rng.randint(0, coord), rng.randint(0, coord))
                 for _ in range(rng.randint(4, 10))}
        if len(cloud) < 3:
            continue
        try:
            P = Polytope(list(cloud))
        except Exception:
            continue
        if P.dim() != 2:
            continue
        pts = [tuple(int(v) for v in p) for p in P.points()]
        if len(pts) > 15:                 # keep TOPCOM enumeration cheap
            continue
        tc = count_fine_triangulations(P, CAP)
        if tc is None:
            continue
        try:
            # one count per dihedral image -- must all agree (orientation
            # independence) and equal TOPCOM
            dp = invariant_count(pts, transforms=COMPACT)
        except AssertionError as e:
            problems.append(f"NON-INVARIANT pts={pts}: {e}")
            continue
        except (ValueError, RuntimeError):
            uncountable += 1              # honest refusal -- fine
        else:
            if dp == tc:
                matched += 1
            else:
                problems.append(f"MISMATCH pts={pts} count_triangulations={dp} topcom={tc}")
        if matched + len(problems) + uncountable >= CHECK_LIMIT:
            break

    assert not problems, (
        f"{len(problems)} bad case(s) (matched={matched}, "
        f"raised={uncountable}):\n  " + "\n  ".join(problems))
    assert matched > 0, "no cases were decided (cytools/TOPCOM may be misconfigured)"
