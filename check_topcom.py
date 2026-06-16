#!/usr/bin/env python3
"""
Independent cross-check of the floor-aware counter (na-query.c, big-integer GMP
back-end) against TOPCOM, via CYTools' Polytope.all_triangulations.

TOPCOM is a separate codebase and algorithm, so agreement is real evidence the
recurrence -- including the non-flat-floor (need_full_table) generalization and
the small widths m=3 / m>=4 -- counts fine (primitive) triangulations correctly.

Scope/limits:
  - Only CONVEX regions are comparable: TOPCOM triangulates a point
    configuration (its convex hull), which equals the region between floor L
    and upper U only when that region is convex.  Non-convex (L,U) are skipped.
  - Enumeration is one triangulation at a time, so this is feasible only for
    SMALL regions (few thousand triangulations), not f(5,6)-scale counts.

This validates on small cases; it is not a proof for all m,n.

Usage:  python3 check_topcom.py     (needs cytools + a GMP toolchain)
"""

import subprocess
from cytools import Polytope

MEM_CAP = 8e9   # bytes; skip configs whose pointer skeleton exceeds this


def skeleton_bytes(m, n):
    n2 = n + 2
    return (m + 1) * (n2 ** (m - 3)) * (n2 + n2 * n2) * 8


def build_dp():
    # m, n are runtime arguments now, so one binary serves every region
    gmp = subprocess.check_output(["brew", "--prefix", "gmp"]).decode().strip()
    out = "/tmp/na_qg_check"
    subprocess.check_call(
        ["gcc", "-O2",
         f"-I{gmp}/include", f"-L{gmp}/lib",
         "-DGMP", "-o", out, "na-query.c", "-lgmp"])
    return out


def dp_query(binary, m, n, U, L=None):
    inp = " ".join(map(str, U)) + "\n"
    if L is not None:
        inp += " ".join(map(str, L)) + "\n"
    out = subprocess.run([binary, str(m), str(n)],
                         input=inp, capture_output=True, text=True).stdout
    for line in out.splitlines():
        if line.startswith("query_value"):
            return int(line.split()[1])
    return None


def region_points(U, L):
    # lattice points of {(x,y): 0<=x<=m, L[x] <= y <= U[x]} for present profiles
    return [[x, y] for x, (lo, hi) in enumerate(zip(L, U)) for y in range(lo, hi + 1)]


def topcom_count(U, L, cap=500000):
    pts = region_points(U, L)
    p = Polytope(pts)
    if len(p.points()) != len(pts):
        return ("non-convex", None)          # region != its convex hull
    c = 0
    for _ in p.all_triangulations(only_fine=True, only_regular=False,
                                  only_star=False,
                                  include_points_interior_to_facets=True):
        c += 1
        if c > cap:
            return ("too-many", None)
    return ("ok", c)


# (m, n, [(name, upper U, floor L), ...]); every U,L has m+1 heights in [0,n]
CONFIGS = [
    (4, 6, [
        ("4x1 rectangle (floor 0)",      [1, 1, 1, 1, 1], [0, 0, 0, 0, 0]),
        ("4x2 rectangle (floor 0)",      [2, 2, 2, 2, 2], [0, 0, 0, 0, 0]),
        ("flat top, valley floor",       [2, 2, 2, 2, 2], [2, 1, 0, 1, 2]),
        ("linear floor (parallelogram)", [1, 2, 3, 4, 5], [0, 1, 2, 3, 4]),
        ("flat floor h=1",               [3, 3, 3, 3, 3], [1, 1, 1, 1, 1]),
        ("concave top, floor 0",         [0, 1, 2, 1, 0], [0, 0, 0, 0, 0]),
    ]),
    (3, 4, [
        ("3x1 rectangle (floor 0)",      [1, 1, 1, 1],    [0, 0, 0, 0]),
        ("3x2 rectangle (floor 0)",      [2, 2, 2, 2],    [0, 0, 0, 0]),
        ("flat top, valley floor",       [2, 2, 2, 2],    [1, 0, 0, 1]),
        ("linear floor (parallelogram)", [1, 2, 3, 4],    [0, 1, 2, 3]),
        ("flat floor h=1",               [3, 3, 3, 3],    [1, 1, 1, 1]),
        ("concave top, floor 0",         [0, 2, 2, 0],    [0, 0, 0, 0]),
    ]),
]


def main():
    dp = build_dp()
    fails = 0
    for m, n, cases in CONFIGS:
        mem = skeleton_bytes(m, n)
        if mem > MEM_CAP:
            print(f"== m={m} n={n} == [SKIP] needs ~{mem/1e9:.0f} GB")
            continue
        print(f"== m={m} n={n} ==")
        for name, U, L in cases:
            status, tc = topcom_count(U, L)
            q = dp_query(dp, m, n, U, L)
            if status != "ok":
                print(f"  [skip] {name:32s} ({status})")
                continue
            ok = (tc == q)
            fails += not ok
            print(f"  [{'OK ' if ok else 'BAD'}] {name:32s} topcom={tc:<8} dp={q}")
    raise SystemExit(1 if fails else 0)


if __name__ == "__main__":
    main()
