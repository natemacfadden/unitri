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
"""Timing benchmark: na_query (transfer-matrix counting) vs TOPCOM (enumeration).

na_query counts triangulations without enumerating them, so its cost depends on
the bounding box (m, n), not the number of triangulations; TOPCOM enumerates one
at a time and cannot reach large regions. Needs a GMP toolchain; the TOPCOM
column needs cytools (skipped, but na_query is still timed, if it is absent).

    python benchmarks/benchmark.py
"""
import os
import subprocess
import time

try:
    from cytools import Polytope
    HAS_CYTOOLS = True
except ImportError:
    HAS_CYTOOLS = False

NA_QUERY_C = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "unitri", "na-query.c")
TRIALS = 5   # na_query is fast; report the min over a few runs

# (name, m, n, upper U, floor L (None = flat 0), topcom-feasible)
CASES = [
    ("3x2 rectangle",       3, 4,  [2, 2, 2, 2],          [0, 0, 0, 0],     True),
    ("polygon (U/L)",       4, 3,  [2, 3, 3, 3, 2],       [2, 1, 0, 0, 1],  True),
    ("polygon (U/L)",       4, 4,  [3, 4, 4, 4, 3],       [3, 1, 0, 0, 1],  True),
    ("4x4 square",          4, 4,  [4, 4, 4, 4, 4],       None,             False),
    ("4x10 square",         4, 10, [10, 10, 10, 10, 10],  None,             False),
    ("triangle, height 84", 3, 84, [84, 56, 28, 0],       None,             False),
]


def gmp_cflags():
    # GMP is on the default path on Linux/conda; Homebrew on macOS puts it elsewhere
    try:
        prefix = subprocess.check_output(["brew", "--prefix", "gmp"],
                                         text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return []
    return [f"-I{prefix}/include", f"-L{prefix}/lib"]


def build_gmp():
    out = "/tmp/na_query_bench"
    subprocess.check_call(
        ["gcc", "-O2", *gmp_cflags(), "-DGMP", "-o", out, NA_QUERY_C, "-lgmp"])
    return out


def na_query(binary, m, n, U, L):
    inp = " ".join(map(str, U)) + "\n"
    if L is not None:
        inp += " ".join(map(str, L)) + "\n"
    out = subprocess.run([binary, str(m), str(n)],
                         input=inp, capture_output=True, text=True).stdout
    for line in out.splitlines():
        if line.startswith("query_value"):
            return int(line.split()[1])
    return None


def topcom_count(U, L, cap=2_000_000):
    """Enumerate fine triangulations with TOPCOM; None if non-convex or > cap."""
    pts = [[x, y] for x, (lo, hi) in enumerate(zip(L, U)) for y in range(lo, hi + 1)]
    P = Polytope(pts)
    if len(P.points()) != len(pts):
        return None                          # region is not its own convex hull
    c = 0
    for _ in P.all_triangulations(only_fine=True, only_regular=False,
                                  only_star=False,
                                  include_points_interior_to_facets=True):
        c += 1
        if c > cap:
            return None
    return c


def time_na_query(binary, m, n, U, L):
    best, count = float("inf"), None
    for _ in range(TRIALS):
        t0 = time.perf_counter()
        count = na_query(binary, m, n, U, L)
        best = min(best, time.perf_counter() - t0)
    return count, best


def main():
    binary = build_gmp()
    print(f"{'region':22} {'triangulations':>18}  {'na_query':>10}  {'TOPCOM':>12}")
    print("-" * 68)
    for name, m, n, U, L, feasible in CASES:
        count, na = time_na_query(binary, m, n, U, L)
        shown = str(count) if count < 10**9 else f"~{count:.1e}"
        if not feasible:
            topcom = "infeasible"
        elif not HAS_CYTOOLS:
            topcom = "no cytools"
        else:
            floor = L if L is not None else [0] * (m + 1)
            t0 = time.perf_counter()
            tc = topcom_count(U, floor)
            topcom = f"{time.perf_counter() - t0:.2f} s" if tc is not None else "too many"
        print(f"{name:22} {shown:>18}  {na * 1000:7.1f} ms  {topcom:>12}")


if __name__ == "__main__":
    main()
