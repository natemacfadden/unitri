"""Shared TOPCOM (via CYTools) enumeration for the cross-check tests.

Both cross-checks -- check_topcom.py (fixed convex regions vs the GMP CLI) and
test_topcom_convex.py (random convex point sets vs count_triangulations) --
enumerate the same fine/primitive triangulations of a CYTools Polytope.  Keeping
that one enumeration (and its flag set) here means the two checks can't drift.
"""


def count_fine_triangulations(P, cap):
    """Number of fine triangulations of CYTools Polytope ``P`` via TOPCOM, or
    None if the count exceeds ``cap``.  The flags select exactly the fine
    (primitive) triangulations that unitri counts."""
    c = 0
    for _ in P.all_triangulations(only_fine=True, only_regular=False,
                                  only_star=False,
                                  include_points_interior_to_facets=True):
        c += 1
        if c > cap:
            return None
    return c
