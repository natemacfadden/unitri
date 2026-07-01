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
"""Timing benchmark: na_query (dynamic-programming counting) vs TOPCOM (enumeration).

na_query counts triangulations without enumerating them, so its cost depends on
the bounding box (m, n), not the number of triangulations; TOPCOM enumerates one
at a time and cannot reach large regions.

Both columns are timed the same way: warm up once (so one-time costs -- caches /
first-touch allocation for na_query, the cytools import + numba JIT for TOPCOM --
are excluded), then measure several batches and report the per-call mean +/-
stdev.  na_query is timed *in process* through the compiled `unitri.na_query`
extension (no subprocess spawn), with auto-scaled inner repetitions so the
microsecond cases are not perf_counter-resolution bound.  TOPCOM is one call per
batch; because it is slow, its repeats stop early once their cumulative time
exceeds TOPCOM_BUDGET, so a multi-minute case does not balloon.

    pip install -e .          # build the extension first
    python benchmarks/benchmark.py

The TOPCOM column needs cytools (the na_query column is still produced if it is
absent).
"""
import statistics
import time

try:
    from cytools import Polytope
    HAS_CYTOOLS = True
except ImportError:
    HAS_CYTOOLS = False

try:
    from unitri.na_query import na_query
except ImportError as e:                     # pragma: no cover - setup guidance
    raise SystemExit(
        "could not import unitri.na_query -- build the extension first with "
        "`pip install -e .` (see the project README).\n"
        f"(import error: {e})")

WARMUP = 1            # untimed calls before measuring (one-time cost exclusion)
MIN_BATCH = 0.05      # grow inner reps until a timed batch lasts >= this (seconds)
REPEAT = 7            # number of timed batches; mean +/- stdev is over these
TOPCOM_BUDGET = 60.0  # seconds; cap repeated runs of the slow TOPCOM baseline

# (name, m, n, upper U, floor L (None = flat 0), topcom-feasible)
CASES = [
    ("3x2 rectangle",       3, 4,  [2, 2, 2, 2],          [0, 0, 0, 0],     True),
    ("polygon",             4, 3,  [2, 3, 3, 3, 2],       [2, 1, 0, 0, 1],  True),
    ("polygon",             4, 4,  [3, 4, 4, 4, 3],       [3, 1, 0, 0, 1],  True),
    ("4x4 square",          4, 4,  [4, 4, 4, 4, 4],       None,             False),
    ("4x10 square",         4, 10, [10, 10, 10, 10, 10],  None,             False),
    ("triangle, height 84", 3, 84, [84, 56, 28, 0],       None,             False),
]


def measure(fn):
    """Per-call (mean, stdev) seconds for a fast routine, timeit-style: warm up,
    auto-scale the inner repetition count so each batch runs >= MIN_BATCH, then
    take REPEAT batches so the spread (not just the best case) is visible."""
    for _ in range(WARMUP):
        fn()

    reps = 1
    while True:
        t0 = time.perf_counter()
        for _ in range(reps):
            fn()
        dt = time.perf_counter() - t0
        if dt >= MIN_BATCH:
            break
        # scale up toward MIN_BATCH (with headroom); guard against dt == 0
        reps = max(reps + 1, int(reps * MIN_BATCH / dt * 1.2)) if dt > 0 else reps * 2

    per_call = []
    for _ in range(REPEAT):
        t0 = time.perf_counter()
        for _ in range(reps):
            fn()
        per_call.append((time.perf_counter() - t0) / reps)
    mean = statistics.fmean(per_call)
    sd = statistics.stdev(per_call) if len(per_call) > 1 else 0.0
    return mean, sd


def measure_slow(fn, repeat=REPEAT, budget=TOPCOM_BUDGET):
    """Per-call (mean, stdev, n_batches) seconds for a slow, expensive baseline
    (one call per batch).  Warm up once, then time up to `repeat` calls, stopping
    early once cumulative measured time exceeds `budget`.  Returns None if the
    warmup call itself returns None (e.g. TOPCOM gave up: non-convex or > cap)."""
    if fn() is None:                          # warmup; also a validity probe
        return None
    times = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn()
        times.append(time.perf_counter() - t0)
        if sum(times) >= budget:
            break
    mean = statistics.fmean(times)
    sd = statistics.stdev(times) if len(times) > 1 else 0.0
    return mean, sd, len(times)


def fmt_profile(U, L):
    """Compact upper/lower boundary heights that define the region's geometry;
    a flat floor (lower is None) shows as 'flat'."""
    upper = ",".join(map(str, U))
    lower = "flat" if L is None else ",".join(map(str, L))
    return f"{upper} / {lower}"


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


def main():
    print(f"{'region':20} {'upper / lower':21} {'triangulations':>16}  "
          f"{'na_query (mean +/- sd)':>24}  {'TOPCOM (mean +/- sd)':>22}")
    print("-" * 111)
    for name, m, n, U, L, feasible in CASES:
        count = na_query(m, n, U, L)
        mean, sd = measure(lambda: na_query(m, n, U, L))
        shown = str(count) if count < 10**9 else f"~{count:.1e}"
        na = f"{mean * 1e3:.3f} +/- {sd * 1e3:.3f} ms"

        if not feasible:
            topcom = "infeasible"
        elif not HAS_CYTOOLS:
            topcom = "no cytools"
        else:
            floor = L if L is not None else [0] * (m + 1)
            res = measure_slow(lambda: topcom_count(U, floor))
            topcom = "too many" if res is None else f"{res[0]:.2f} +/- {res[1]:.2f} s"
        print(f"{name:20} {fmt_profile(U, L):21} {shown:>16}  {na:>24}  {topcom:>22}")


if __name__ == "__main__":
    main()
