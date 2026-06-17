#!/usr/bin/env python3
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
r"""Uniform sampler for fine (unimodular) triangulations of a lattice polygon.

Method (provably uniform): fix a boundary edge e=(a,b).  In every triangulation
e belongs to exactly one triangle (a,b,c), so
    f(S) = sum_c f(S \ triangle(a,b,c))
is a sign-free partition (no inclusion-exclusion, no signs).  Removing the
triangle either inserts c into the boundary (c interior) or splits S into two
sub-polygons at c (c on the boundary), where f then multiplies.  Exact counts
use Python bignums; sampling picks c with probability f(child)/f(S) and recurses.

This is the clean general-position version; it reuses the same fine-triangulation
count we validated against TOPCOM and the Orevkov DP.
"""
import sys, random
from math import gcd
from functools import lru_cache

def area2(poly):
    s = 0
    for i in range(len(poly)):
        x1, y1 = poly[i]; x2, y2 = poly[(i+1) % len(poly)]
        s += x1*y2 - x2*y1
    return s

def cross(o, a, b):
    return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])

def boundary_points(poly):
    """All lattice points on the boundary, in CCW order (every edge primitive)."""
    pts = []
    n = len(poly)
    for i in range(n):
        a = poly[i]; b = poly[(i+1) % n]
        dx = b[0]-a[0]; dy = b[1]-a[1]; g = gcd(abs(dx), abs(dy))
        for k in range(g):
            pts.append((a[0]+dx*k//g, a[1]+dy*k//g))
    return pts

def point_in_poly(p, poly):
    """True if p is inside or on boundary of simple polygon poly."""
    x, y = p; n = len(poly); inside = False
    for i in range(n):
        ax, ay = poly[i]; bx, by = poly[(i+1) % n]
        # on-segment test
        if min(ax, bx) <= x <= max(ax, bx) and min(ay, by) <= y <= max(ay, by):
            if (bx-ax)*(y-ay) - (by-ay)*(x-ax) == 0:
                return True
        if (ay > y) != (by > y):
            xint = ax + (bx-ax)*(y-ay)/(by-ay)
            if x < xint:
                inside = not inside
    return inside

def seg_proper_cross(a, b, c, d):
    d1 = cross(c, d, a); d2 = cross(c, d, b)
    d3 = cross(a, b, c); d4 = cross(a, b, d)
    if d1 == 0 or d2 == 0 or d3 == 0 or d4 == 0:
        return False  # endpoint-touching / collinear is not a proper crossing
    return ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0))

def orient(poly):
    """Reduce collinear vertices, force CCW, rotate to start at the
    lexicographically smallest vertex -- but DO NOT translate, so coordinates
    stay absolute.  Both count() and sample() run on orient(poly), guaranteeing
    they choose the same first edge (poly[0],poly[1])."""
    if area2(poly) < 0:
        poly = poly[::-1]
    v = []
    n = len(poly)
    for i in range(n):
        if cross(poly[i-1], poly[i], poly[(i+1) % n]) != 0:
            v.append(poly[i])
    if len(v) < 3:
        v = list(poly)
    i0 = min(range(len(v)), key=lambda i: v[i])
    return tuple(v[i0:] + v[:i0])

def canon(poly):
    """Translation-normalized memo key (counts are translation-invariant)."""
    v = orient(poly)
    minx = min(p[0] for p in v); miny = min(p[1] for p in v)
    return tuple((p[0]-minx, p[1]-miny) for p in v)

def children(poly):
    """Yield (apex c, [child polygons]) for removing the triangle on edge
    (poly[0],poly[1]); child polygons multiply."""
    bpts = boundary_points(poly)
    n = len(bpts)
    a = bpts[0]; b = bpts[1]
    inside_pts = [p for p in lattice_points(poly)]
    res = []
    for c in inside_pts:
        if c == a or c == b:
            continue
        if cross(a, b, c) != 1:           # empty unimodular tri on interior side
            continue
        # diagonals a-c, b-c must not properly cross the boundary
        ok = True
        for i in range(n):
            e1 = bpts[i]; e2 = bpts[(i+1) % n]
            if seg_proper_cross(a, c, e1, e2) or seg_proper_cross(b, c, e1, e2):
                ok = False; break
        if not ok:
            continue
        cent = ((a[0]+b[0]+c[0])/3.0, (a[1]+b[1]+c[1])/3.0)
        if not point_in_poly(cent, poly):
            continue
        if c in bpts:                     # split at boundary vertex c
            k = bpts.index(c)
            polyA = bpts[1:k+1]           # b .. c
            polyB = [bpts[0]] + bpts[k:]  # a, c .. last
            parts = [p for p in (polyA, polyB) if len(p) >= 3]
        else:                             # c interior: insert into boundary
            parts = [[bpts[0], c] + bpts[1:]]
        res.append((c, parts))
    return res

# region lattice points (cached per canonical polygon)
_lp_cache = {}
def lattice_points(poly):
    key = orient(poly)            # absolute (untranslated): points stay correct
    if key in _lp_cache:
        return _lp_cache[key]
    xs = [p[0] for p in poly]; ys = [p[1] for p in poly]
    pts = [(x, y) for x in range(min(xs), max(xs)+1)
                  for y in range(min(ys), max(ys)+1)
                  if point_in_poly((x, y), poly)]
    _lp_cache[key] = pts
    return pts

_count_cache = {}
def count(poly):
    key = canon(poly)             # translation-merged memo key
    if key in _count_cache:
        return _count_cache[key]
    poly = list(orient(poly))     # same first edge as sample(), abs coords
    if len(lattice_points(poly)) <= 3:
        _count_cache[key] = 1
        return 1
    total = 0
    for c, parts in children(poly):
        prod = 1
        for part in parts:
            prod *= count(part)
        total += prod
    _count_cache[key] = total
    return total

def sample(poly, rng):
    """Return a list of triangles (each a frozenset of 3 points), in the
    polygon's ABSOLUTE coordinates.  Runs on orient(poly) -- same first edge as
    count() -- so the branch weights are consistent."""
    poly = list(orient(poly))
    if len(lattice_points(poly)) <= 3:
        pts = lattice_points(poly)
        return [frozenset(pts)] if len(pts) == 3 else []
    bpts = boundary_points(poly)
    a, b = bpts[0], bpts[1]
    opts = children(poly)
    weights = []
    for c, parts in opts:
        w = 1
        for part in parts:
            w *= count(part)
        weights.append(w)
    total = sum(weights)
    r = rng.randrange(total)
    acc = 0
    for (c, parts), w in zip(opts, weights):
        acc += w
        if r < acc:
            tris = [frozenset([a, b, c])]
            for part in parts:
                tris += sample(part, rng)
            return tris
    raise RuntimeError("sampling failed")

def all_triangulations(poly):
    """Yield every fine triangulation (each a list of triangle frozensets),
    exactly once.  Same canonical-edge partition as count(): the triangle on the
    first edge is unique per triangulation, so iterating all apex choices (and,
    for split children, the product over the pieces) enumerates with no
    duplicates.  Cost is O(count(poly)) -- small polygons only."""
    import itertools
    poly = list(orient(poly))
    pts = lattice_points(poly)
    if len(pts) <= 3:
        yield [frozenset(pts)] if len(pts) == 3 else []
        return
    bpts = boundary_points(poly)
    a, b = bpts[0], bpts[1]
    for c, parts in children(poly):
        for combo in itertools.product(*(all_triangulations(p) for p in parts)):
            tris = [frozenset([a, b, c])]
            for part_tris in combo:
                tris += part_tris
            yield tris

if __name__ == "__main__":
    # default demo polygon: [0,4]^2
    import json
    if len(sys.argv) > 1:
        poly = tuple(map(tuple, json.loads(sys.argv[1])))
    else:
        poly = ((0, 0), (4, 0), (4, 4), (0, 4))
    print("polygon:", poly)
    print("fine triangulations:", count(poly))
