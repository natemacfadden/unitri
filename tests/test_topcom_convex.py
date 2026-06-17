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
"""Cross-check unitri.count_triangulations against TOPCOM on convex point sets.

Generates random small convex lattice point sets and, for each, compares
count_triangulations (which picks a hull-tracing orientation and feeds na_query)
against TOPCOM (CYTools).  The contract under test: count_triangulations either
returns a count that EXACTLY matches TOPCOM, or it raises (point sets na_query
cannot represent without undercounting -- those needing "long-diagonal"
triangulations).  A returned-but-wrong count is a failure.

    python3 tests/test_topcom_convex.py
"""
import random
import sys

from cytools import Polytope

import unitri

SEED = 20260617
TRIALS = 4000
CHECK_LIMIT = 140    # stop after this many decided (matched/uncountable) cases
CAP = 5000           # skip regions TOPCOM can't enumerate quickly


def topcom_count(P, cap):
    c = 0
    for _ in P.all_triangulations(only_fine=True, only_regular=False,
                                  only_star=False,
                                  include_points_interior_to_facets=True):
        c += 1
        if c > cap:
            return None
    return c


def main():
    rng = random.Random(SEED)
    matched = mismatched = uncountable = 0
    for _ in range(TRIALS):
        coord = rng.choice([3, 4, 5])
        cloud = {(rng.randint(0, coord), rng.randint(0, coord))
                 for _ in range(rng.randint(4, 9))}
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
        tc = topcom_count(P, CAP)
        if tc is None:
            continue
        try:
            dp = unitri.count_triangulations(pts)
        except (ValueError, RuntimeError):
            uncountable += 1              # honest refusal -- fine
        else:
            if dp == tc:
                matched += 1
            else:
                mismatched += 1
                print(f"  MISMATCH pts={pts} count_triangulations={dp} topcom={tc}")
        if matched + mismatched + uncountable >= CHECK_LIMIT:
            break

    print(f"matched TOPCOM            = {matched}")
    print(f"raised (na_query can't)   = {uncountable}")
    print(f"MISMATCHED (must be 0)    = {mismatched}")
    print("OK" if mismatched == 0 else "FAILED")
    sys.exit(1 if mismatched else 0)


if __name__ == "__main__":
    main()
