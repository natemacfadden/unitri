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

To count a point set's triangulations, the profiles must trace its convex hull.
In a chosen orientation we read the hull boundary at each integer column; where
the boundary height is non-integer (the hull edge skips a row) we mark that
column with the n+1 sentinel ("."), an absent vertex.

Counts are unimodular-invariant, so we pick a minimal-width orientation (width
>= 3, na_query's lower bound) to keep the ~ (n+2)^(m-1) cost down, and cross-check
against TOPCOM. If every orientation has width < 3 we raise.
"""
from __future__ import annotations

from collections.abc import Iterable
from math import gcd

Point = tuple[int, int]


def _ext_gcd(a: int, b: int) -> tuple[int, int, int]:
    """Return (g, x, y) with a*x + b*y = g = gcd(a, b)."""
    if b == 0:
        return (abs(a), (1 if a >= 0 else -1) if a else 0, 0)
    g, x, y = _ext_gcd(b, a % b)
    return (g, y, x - (a // b) * y)


def _width(pts: list[Point], a: int, b: int) -> int:
    vals = [a * x + b * y for x, y in pts]
    return max(vals) - min(vals)


def _convex_hull(pts: list[Point]) -> list[Point]:
    """CCW convex-hull vertices of `pts` (Andrew's monotone chain)."""
    pts = sorted(set(pts))

    def cross(o: Point, a: Point, b: Point) -> int:
        return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])

    lower: list[Point] = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper: list[Point] = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def _column_bounds(hull: list[Point], x: int
                   ) -> tuple[tuple[int, int] | None, tuple[int, int] | None]:
    """Hull's [lower, upper] y-range at integer column `x`, each as a reduced
    (num, den) with den > 0 (a lattice point exactly when num % den == 0).
    Robust to vertical edges (they only occur at the extreme columns)."""
    lo: tuple[int, int] | None = None
    hi: tuple[int, int] | None = None
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


def _orientations(pts: list[Point]) -> list[tuple[int, int]]:
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


def _build_profile(pts: list[Point], a: int, b: int
                   ) -> tuple[int, int, list[int], list[int]]:
    """(m, n, upper, lower) reading the convex-hull boundary in orientation
    u = (a, b). A column whose boundary height is non-integer gets the n+1
    sentinel ("."), an absent vertex."""
    _, x0, y0 = _ext_gcd(a, b)                  # a*x0 + b*y0 = 1
    c, d = -y0, x0
    tpts = [(a * x + b * y, c * x + d * y) for x, y in pts]
    sx = min(X for X, _ in tpts)
    sy = min(Y for _, Y in tpts)
    tpts = [(X - sx, Y - sy) for X, Y in tpts]
    hull = _convex_hull(tpts)
    m = max(X for X, _ in tpts)
    present_U: list[int] = []
    raw_U: list[int | None] = []
    raw_L: list[int | None] = []
    for X in range(m + 1):
        col_lo, col_hi = _column_bounds(hull, X)
        # columns 0..m are inside the hull's x-range, so both bounds exist
        assert col_lo is not None and col_hi is not None
        (ln, ld), (un, ud) = col_lo, col_hi
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


def points_to_profiles(points: Iterable[tuple[int, int]]
                       ) -> tuple[int, int, list[int], list[int]]:
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


def count_triangulations(points: Iterable[tuple[int, int]]) -> int:
    """Count the fine triangulations of a 2-D lattice point set -> Python int.

    Uses a minimal-width hull-tracing orientation (na_query needs m >= 3 and
    costs ~ (n+2)^(m-1)); the count is orientation-invariant. Falls through to
    the next orientation only if na_query rejects one. Requires the built
    extension.
    """
    try:
        from .na_query import na_query
    except ImportError as exc:
        raise ImportError(
            "count_triangulations needs the compiled GMP extension (build with "
            "`pip install -e .`, which needs libgmp).  Without GMP, use "
            "count_triangulations_parallel (mod-prime + CRT)."
        ) from exc
    pts = sorted({(int(x), int(y)) for x, y in points})
    if len(pts) < 3:
        raise ValueError("need at least 3 distinct points")
    ors = _orientations(pts)
    if not ors:
        raise ValueError("point set is too thin: width < 3 in every orientation "
                         "(na_query needs m >= 3)")
    last: RuntimeError | None = None
    for a, b in ors[:24]:                       # bound the retries
        m, n, upper, lower = _build_profile(pts, a, b)
        try:
            return na_query(m, n, upper, lower)
        except RuntimeError as e:               # na_query rejected this orientation
            last = e
    raise ValueError(f"na_query rejected all tried orientations (last: {last})")
