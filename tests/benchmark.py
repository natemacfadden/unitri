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
at a time and cannot reach large regions. Needs cytools + a GMP toolchain.

    python tests/benchmark.py
"""
import time

from check_topcom import build_dp, dp_query, topcom_count

# (name, m, n, upper U, floor L (None = flat 0), topcom-feasible)
CASES = [
    ("3x2 rectangle",       3, 4,  [2, 2, 2, 2],          [0, 0, 0, 0],     True),
    ("polygon (U/L)",       4, 3,  [2, 3, 3, 3, 2],       [2, 1, 0, 0, 1],  True),
    ("polygon (U/L)",       4, 4,  [3, 4, 4, 4, 3],       [3, 1, 0, 0, 1],  True),
    ("4x4 square",          4, 4,  [4, 4, 4, 4, 4],       None,             False),
    ("4x10 square",         4, 10, [10, 10, 10, 10, 10],  None,             False),
    ("triangle, height 84", 3, 84, [84, 56, 28, 0],       None,             False),
]
TRIALS = 5   # na_query is fast; report the min over a few runs


def time_na_query(binary, m, n, U, L):
    best, count = float("inf"), None
    for _ in range(TRIALS):
        t0 = time.perf_counter()
        count = dp_query(binary, m, n, U, L)
        best = min(best, time.perf_counter() - t0)
    return count, best


def main():
    binary = build_dp()
    print(f"{'region':22} {'triangulations':>18}  {'na_query':>10}  {'TOPCOM':>12}")
    print("-" * 68)
    for name, m, n, U, L, feasible in CASES:
        count, na = time_na_query(binary, m, n, U, L)
        shown = str(count) if count < 10**9 else f"~{count:.1e}"
        if feasible:
            floor = L if L is not None else [0] * (m + 1)
            t0 = time.perf_counter()
            status, _ = topcom_count(U, floor, cap=2_000_000)
            topcom = f"{time.perf_counter() - t0:.2f} s" if status == "ok" else status
        else:
            topcom = "infeasible"
        print(f"{name:22} {shown:>18}  {na * 1000:7.1f} ms  {topcom:>12}")


if __name__ == "__main__":
    main()
