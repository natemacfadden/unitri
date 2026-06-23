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
"""Convert a set of 2-D lattice points into na_query (m, n, upper, lower).

na_query counts triangulations of the region *under the column profiles*.  For
that to equal the number of triangulations of a point set's convex hull, the
profiles must trace the hull boundary.  We do this exactly: in a chosen
orientation we read the hull's upper/lower boundary at each integer column, and
where the boundary crosses *between* lattice rows (a hull edge spanning more than
one column) we mark that column ABSENT with na_query's n+1 sentinel ("."), so the
profiles describe the true hull rather than a per-column-max polyline that dips
inside it.

Triangulation counts are unimodular-invariant, so we pick a minimal-width
orientation (width >= 3, na_query's lower bound) to keep na_query's
~ (n+2)^(m-1) cost down.  The result is cross-checked against TOPCOM.  If every
orientation is too thin (width < 3) we raise rather than guess.
"""
from math import gcd


def _ext_gcd(a, b):
    """Return (g, x, y) with a*x + b*y = g = gcd(a, b)."""
    if b == 0:
        return (abs(a), (1 if a >= 0 else -1) if a else 0, 0)
    g, x, y = _ext_gcd(b, a % b)
    return (g, y, x - (a // b) * y)


def _width(pts, a, b):
    vals = [a * x + b * y for x, y in pts]
    return max(vals) - min(vals)


def _convex_hull(pts):
    """CCW convex-hull vertices of `pts` (Andrew's monotone chain)."""
    pts = sorted(set(pts))

    def cross(o, a, b):
        return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])

    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def _column_bounds(hull, x):
    """Hull's [lower, upper] y-range at integer column `x`, each as a reduced
    (num, den) with den > 0 (a lattice point exactly when num % den == 0).
    Robust to vertical edges (they only occur at the extreme columns)."""
    lo = hi = None
    V = len(hull)
    for i in range(V):
        x0, y0 = hull[i]
        x1, y1 = hull[(i + 1) % V]
        if x0 == x1:                              # vertical edge (extreme column)
            if x0 == x:
                for yy in (y0, y1):
                    if hi is None or yy * hi[1] > hi[0]:
                        hi = (yy, 1)
                    if lo is None or yy * lo[1] < lo[0]:
                        lo = (yy, 1)
            continue
        if min(x0, x1) <= x <= max(x0, x1):
            num = y0 * (x1 - x0) + (y1 - y0) * (x - x0)
            den = x1 - x0
            if den < 0:
                num, den = -num, -den
            if hi is None or num * hi[1] > hi[0] * den:
                hi = (num, den)
            if lo is None or num * lo[1] < lo[0] * den:
                lo = (num, den)
    return lo, hi


def _orientations(pts):
    """Primitive orientations (a, b) with lattice width >= 3, smallest width
    first (smaller width => smaller m => cheaper na_query)."""
    xs = [x for x, _ in pts]
    ys = [y for _, y in pts]
    K = max(max(xs) - min(xs), max(ys) - min(ys), 1)
    cands = []
    for a in range(0, K + 1):
        for b in range(-K, K + 1):
            if a == 0 and b <= 0:
                continue
            if gcd(abs(a), abs(b)) != 1:
                continue
            w = _width(pts, a, b)
            if w >= 3:
                cands.append((w, (a, b)))
    cands.sort()
    return [ab for _, ab in cands]


def _build_profile(pts, a, b):
    """(m, n, upper, lower) reading the convex-hull boundary in the
    orientation u = (a, b).  A column the boundary crosses at a non-integer
    height (passing *between* lattice rows) is marked ABSENT with the n+1
    sentinel -- na_query's "." vertex -- so the profiles describe the true hull,
    not a dipped per-column-max polyline."""
    _, x0, y0 = _ext_gcd(a, b)                  # a*x0 + b*y0 = 1
    c, d = -y0, x0
    tpts = [(a * x + b * y, c * x + d * y) for x, y in pts]
    sx = min(X for X, _ in tpts)
    sy = min(Y for _, Y in tpts)
    tpts = [(X - sx, Y - sy) for X, Y in tpts]
    hull = _convex_hull(tpts)
    m = max(X for X, _ in tpts)
    present_U, raw_U, raw_L = [], [], []
    for X in range(m + 1):
        (ln, ld), (un, ud) = _column_bounds(hull, X)
        u = un // ud if un % ud == 0 else None
        lo = ln // ld if ln % ld == 0 else None
        raw_U.append(u)
        raw_L.append(lo)
        if u is not None:
            present_U.append(u)
    n = max(present_U)
    absent = n + 1
    upper = [absent if u is None else u for u in raw_U]
    lower = [absent if lo is None else lo for lo in raw_L]
    return m, n, upper, lower


def points_to_profiles(points):
    """Convert 2-D lattice points to na_query's (m, n, upper, lower), reading
    the convex-hull boundary in a minimal-width orientation (absent columns get
    the n+1 sentinel).
    """
    pts = sorted({(int(x), int(y)) for x, y in points})
    if len(pts) < 3:
        raise ValueError("need at least 3 distinct points")
    ors = _orientations(pts)
    if not ors:
        raise ValueError("point set is too thin: width < 3 in every orientation "
                         "(na_query needs m >= 3)")
    a, b = ors[0]
    m, n, upper, lower = _build_profile(pts, a, b)
    return m, n, upper, lower


def count_triangulations(points):
    """Count the fine triangulations of a 2-D lattice point set -> Python int.

    Uses a hull-tracing orientation, smallest lattice width first (na_query needs
    m >= 3 and costs ~ (n+2)^(m-1)).  The count is orientation-invariant and the
    hull+absent profile is exact, so the smallest-width orientation suffices; we
    fall through to the next only if na_query refuses one (a rare profile it can't
    represent), since a different orientation may avoid it.  Requires the built
    extension.
    """
    from .na_query import na_query
    pts = sorted({(int(x), int(y)) for x, y in points})
    if len(pts) < 3:
        raise ValueError("need at least 3 distinct points")
    ors = _orientations(pts)
    if not ors:
        raise ValueError("point set is too thin: width < 3 in every orientation "
                         "(na_query needs m >= 3)")
    last = None
    for a, b in ors[:24]:                       # bound the retries
        m, n, upper, lower = _build_profile(pts, a, b)
        try:
            return na_query(m, n, upper, lower)
        except RuntimeError as e:               # na_query rejected this orientation
            last = e
    raise ValueError(f"na_query rejected all tried orientations (last: {last})")
