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
"""na-query.c (big-integer GMP backend) against known exact counts.

Each case names a region (a profile of upper heights U over a floor L within an
m x n bounding box) and an expected exact count.  The query is run at the right
(m, n) and compared; the x<->m-x reflected profile (a unimodular map) is also
requeried and the count must be unchanged.

Sources for the expected values:
  1,3: https://arxiv.org/pdf/math/0211268
  2,4: https://arxiv.org/pdf/2602.16909
  6,7: TOPCOM
  5  : no independent check (expected value is conjectured)

Memory for this code's index skeleton grows like (n+2)^(m-1); cases needing
more than MEM_CAP are skipped (they require a big-memory machine).

Also checks the flat-rectangle f(m,k) table that na-query.c prints when given no
query on stdin -- the only path that exercises the table emitter.
"""
import re
import subprocess

import pytest

import unitri
from _cli import run_query

MEM_CAP = 8e9   # bytes; skip cases whose pointer skeleton exceeds this

# name, m, n, upper U, floor L (None = flat 0), expected count
CASES = [
    ("4x4 square",
     4, 4, [4, 4, 4, 4, 4], None, "736983568"),
    ("triangle (0,0),(3,0),(0,84)",          # arXiv 2602.16909, ~7.61e65
     3, 84, [84, 56, 28, 0], None,
     "761342982944289349099618507228200078481281500600912757801568059775"),
    ("4x10 square",
     4, 10, [10, 10, 10, 10, 10], None, "584455230176565718869688"),
    ("triangle (0,0),(7,0),(0,28)",
     7, 28, [28, 24, 20, 16, 12, 8, 4, 0], None,
     "153405520601065827395233041403916434565495619834253255337377215694753592286"),
    ("polygon L/U, box height 26",
     4, 26, [26, 22, 18, 13, 13], [0, 4, 8, 13, 13], "1189443740608860381225"),
    ("polygon L/U (TOPCOM 840021)",
     4, 4, [3, 4, 4, 4, 3], [3, 1, 0, 0, 1], "840021"),
    ("polygon L/U (TOPCOM 10653)",
     4, 3, [2, 3, 3, 3, 2], [2, 1, 0, 0, 1], "10653"),
]


def skeleton_bytes(m, n):
    n2 = n + 2
    return (m + 1) * (n2 ** (m - 3)) * (n2 + n2 * n2) * 8


def _reflect(profile):
    # x <-> m-x is a unimodular map; reversing the profile must not change the
    # count.  (This is exactly the symmetry the flat-floor recurrence exploits,
    # so it directly exercises the height_0 < height_m reflection handling.)
    return None if profile is None else profile[::-1]


@pytest.mark.parametrize("name,m,n,U,L,expected", CASES, ids=[c[0] for c in CASES])
def test_count(na_query_bin, name, m, n, U, L, expected):
    mem = skeleton_bytes(m, n)
    if mem > MEM_CAP:
        pytest.skip(f"needs ~{mem/1e9:.0f} GB (m={m}, n={n}); big-memory machine")

    got = run_query(na_query_bin, m, n, U, L)
    assert got is not None, "no query_value (the run may have failed)"
    if expected.startswith("~"):                 # approximate expectation
        ratio = int(got) / float(expected[1:])
        assert 0.5 < ratio < 2.0, f"got {got}, ~{ratio:.3g}x expected {expected}"
    else:
        assert got == expected, f"got {got}, expected {expected}"

    # unimodular (reflection) invariance: the reversed profile must agree
    rev = run_query(na_query_bin, m, n, _reflect(U), _reflect(L))
    assert rev == got, f"reflected profile gave {rev}, expected {got}"


# flat-rectangle f(m,k) table (na-query.c with no query on stdin).  Exact and
# cross-consistent -- f(m,k)=f(k,m), e.g. f(3,4)=f(4,3)=2822648 -- and matches
# the known rectangle counts (f(3,2)=852, f(4,4)=736983568).
TABLES = [
    (3, 5, [20, 852, 46456, 2822648, 182881520]),
    (4, 4, [70, 12170, 2822648, 736983568]),
]


def _ftable(binary, m, n):
    """Run in table mode (empty stdin) and return the list of f(m,k) values."""
    out = subprocess.run([binary, str(m), str(n)],
                         input="", capture_output=True, text=True).stdout
    return [int(v) for v in re.findall(r"f\(\d+,\d+\) = \*\) (\d+)", out)]


@pytest.mark.parametrize("m,n,expected", TABLES, ids=[f"{m}x{n}" for m, n, _ in TABLES])
def test_ftable_mode(na_query_bin, m, n, expected):
    assert _ftable(na_query_bin, m, n) == expected


# The README Performance-table regions as convex-hull vertices, so the point-set
# API is exercised on the big counts too -- both the exact GMP path and the
# parallel mod-prime + CRT path must reproduce them.  (These vertices hull to the
# same regions as the profiles above.)
VERTEX_CASES = [
    ("3x2 rectangle",  [(0, 0), (0, 2), (3, 0), (3, 2)],                          852),
    ("polygon 10653",  [(0, 2), (1, 3), (2, 0), (3, 0), (3, 3), (4, 1), (4, 2)],  10653),
    ("polygon 840021", [(0, 3), (1, 1), (1, 4), (2, 0), (3, 0), (3, 4), (4, 1), (4, 3)], 840021),
    ("4x4 square",     [(0, 0), (0, 4), (4, 0), (4, 4)],                          736983568),
    ("4x10 square",    [(0, 0), (0, 10), (4, 0), (4, 10)],                        584455230176565718869688),
    ("triangle h84",   [(0, 0), (0, 84), (3, 0)],
     761342982944289349099618507228200078481281500600912757801568059775),
]


@pytest.mark.parametrize("name,verts,expected", VERTEX_CASES, ids=[c[0] for c in VERTEX_CASES])
def test_benchmark_region_pointset(name, verts, expected):
    assert unitri.count_triangulations(verts) == expected
    assert unitri.count_triangulations_parallel(verts) == expected
