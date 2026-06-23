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
"""Tough cases for the x<->m-x reflection used by the flat-floor recurrence.

The recurrence stores only profiles with height_0 >= height_m and recovers the
rest as mirror images.  A query whose profile has height_0 < height_m therefore
lands in the un-stored half; before the fix that returned "not_found".  Every
case below is chosen to exercise that half hard:

  * asymmetric endpoints with the small height on the left (height_0 < height_m),
  * absent vertices ('.') placed asymmetrically (long hull edges skipping rows),
  * nontrivial upper AND lower profiles (the need_full_table path; reversing the
    query reverses the floor too, so reflection is checked there as well).

For each case three things are checked against na-query.c (GMP build):
  1. the query is found (not "not_found") and equals the expected count,
  2. the reversed profile -- a unimodular reflection -- gives the same count,
  3. where the region is small enough (< 17 lattice points), TOPCOM agrees.

Counts marked TOPCOM-verified were cross-checked against CYTools; the larger
ones are regression pins (TOPCOM cannot enumerate >= ~17 points in reasonable
time) and are still checked for reflection invariance.

    pytest tests/test_symmetry.py
"""
import subprocess

import pytest

TOPCOM_MAX_POINTS = 17     # CYTools enumeration is infeasible at/above this

# name, m, n, upper U, floor L (None = flat 0), expected count
# (tokens: integers or '.' for an absent vertex)
CASES = [
    # --- TOPCOM-verified (< 17 lattice points) ----------------------------
    ("absent vertex, ends 1<2",      3, 2, "1 . 2 2",     None, 134),
    ("absent vertex, ends 0<3",      3, 3, "0 . 3 3",     None, 186),
    ("absent vertex, ends 0<4",      3, 4, "0 . 3 4",     None, 522),
    ("absent vertex, ends 0<5",      3, 5, "0 . 3 5",     None, 717),
    ("absent vertex, ends 1<2 (m4)", 4, 2, "1 . 2 2 2",   None, 2042),
    ("two absent, ends 0<4",         4, 4, "0 . 3 . 4",   None, 15546),
    ("absent, rising 1<4",           4, 4, "1 . 2 3 4",   None, 17863),
    # --- nontrivial U AND L, TOPCOM-verified (full-table path; the reversed
    #     query reverses the floor too, so this checks reflection there) ------
    ("nontrivial U/L, symmetric",    3, 3, "2 3 3 2", "1 0 0 1", 1736),
    ("nontrivial U/L, asymmetric",   3, 3, "3 3 3 2", "1 0 0 1", 4768),
    ("nontrivial U/L, asymmetric 2", 3, 4, "1 3 4 3", "0 0 1 2", 1564),
    # --- regression pins (>= 17 points; TOPCOM infeasible) ----------------
    ("absent, ends 0<5 (m4)",        4, 5, "0 . 3 4 5",   None, 237625),
    ("two absent, ends 0<6",         4, 6, "0 . 3 . 6",   None, 408826),
    ("three absent, ends 0<5 (m5)",  5, 5, "0 . 3 . . 5", None, 7799342),
    ("nontrivial U/L, asymmetric",   4, 6, "6 5 4 3 3", "0 1 2 3 3", 14303),
    ("nontrivial U/L, reversed",     4, 6, "3 3 4 5 6", "3 3 2 1 0", 14303),
]


def query(binary, m, n, U, L):
    """Return the count string, or None on not_found / an invalid profile."""
    inp = " ".join(U) + "\n" + (" ".join(L) + "\n" if L else "")
    out = subprocess.run([binary, str(m), str(n)],
                         input=inp, capture_output=True, text=True).stdout
    for line in out.splitlines():
        if line.startswith("query_value"):
            tok = line.split()[1]
            return None if tok == "not_found" else tok
    return None


def _absent(tok, n):
    return tok == "." or tok == str(n + 1)


def _height_at(prof, x, n):
    if not _absent(prof[x], n):
        return int(prof[x]), 1
    lo, hi = x, x
    while _absent(prof[lo], n):
        lo -= 1
    while _absent(prof[hi], n):
        hi += 1
    return int(prof[lo]) * (hi - x) + int(prof[hi]) * (x - lo), hi - lo


def region_points(m, n, U, L):
    pts = []
    for x in range(m + 1):
        un, ud = _height_at(U, x, n)
        top = un // ud                       # floor of the upper boundary
        bot = 0 if L is None else -((-_height_at(L, x, n)[0]) // _height_at(L, x, n)[1])
        pts += [(x, y) for y in range(bot, top + 1)]
    return pts


def topcom_count(m, n, U, L):
    """('ok', count) | ('skip', reason).  Uses CYTools; convex regions only."""
    pts = sorted(set(region_points(m, n, U, L)))
    if len(pts) >= TOPCOM_MAX_POINTS:
        return ("skip", f"{len(pts)} points (TOPCOM infeasible)")
    try:
        from cytools import Polytope
    except ImportError:
        return ("skip", "cytools not installed")
    P = Polytope(pts)
    if P.dim() != 2 or len(P.points()) != len(pts):
        return ("skip", "non-convex region")
    c = 0
    for _ in P.all_triangulations(only_fine=True, only_regular=False,
                                  only_star=False,
                                  include_points_interior_to_facets=True):
        c += 1
    return ("ok", c)


@pytest.mark.parametrize("name,m,n,U,L,expected", CASES, ids=[c[0] for c in CASES])
def test_symmetry(na_query_bin, name, m, n, U, L, expected):
    Us, Ls = U.split(), (L.split() if L else None)

    got = query(na_query_bin, m, n, Us, Ls)
    assert got is not None, "not_found / invalid profile"
    assert got == str(expected), f"count {got} != expected {expected}"

    # the reversed profile is a unimodular reflection -> same count
    rev = query(na_query_bin, m, n, Us[::-1], Ls[::-1] if Ls else None)
    assert rev == got, f"not reflection-invariant (reversed={rev})"

    # where the region is small enough, TOPCOM must agree
    status, tc = topcom_count(m, n, Us, Ls)
    if status == "ok":
        assert str(tc) == got, f"TOPCOM disagrees ({tc})"
