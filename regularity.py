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
"""Regularity check for a fine lattice triangulation.

Ported (verbatim, minus the surrounding package) from dualGNN:
src/dualgnn/geometry.py :: is_regular.  Uses `regfans` to build a fan and test
whether the triangulation is regular (i.e. induced by a lifting / lower convex
hull).  A triangulation of a 2D lattice polygon is regular iff such a lifting
exists; not all fine triangulations are (e.g. P_{4,4} has 1,553,020 irregular
ones).

Requires `regfans` and `numpy`.
"""
import ctypes
import signal
import warnings

import numpy as np
from regfans import VectorConfiguration

# ppl messes up rounding types... we have to fix them periodically :(
_libc = ctypes.CDLL(None)
_libc.fesetround(0)


def is_regular(pts: np.ndarray, simps: np.ndarray) -> bool:
    """
    Check if a triangulation is regular via `regfans`. This requires
        1) homogenizing pts,
        2) building a Fan in regfans, and then
        3) checking the regularity of said Fan.
    A 60s SIGALRM timeout guards the call (regfans can hang on degenerate
    inputs); regfans warnings are promoted to errors. Both timeout and
    promoted-warning are reported and return False; any other exception
    propagates.

    Parameters
    ----------
    pts : ndarray
        `(Npts, 2)` int. Lattice points of the polygon.
    simps : ndarray
        `(Nsimps, 3)` int. Each row is a triple of indices into `pts`. Together
        they should cover the polygon (i.e., form a triangulation); the function
        does not check that.

    Returns
    -------
    regular : bool
        True iff the triangulation is regular. If regfans hangs >60s or emits a
        warning, this is set to False and logged.
    """
    if len(simps) == 1:
        return True

    # homogenize
    ones = np.ones((pts.shape[0], 1), dtype=pts.dtype)
    vecs = np.hstack([pts, ones])

    # make the fan
    vc = VectorConfiguration(vecs, labels=list(range(len(pts))))
    fan = vc.triangulate(cells=simps)

    # check regularity, with a SIGALRM watchdog for regfans hangs.
    # Note: SIGALRM only works on the main thread of the main interpreter,
    # so calling is_regular from a worker thread will raise ValueError -- if
    # you ever need threaded harvesting, replace this with a thread-safe
    # watchdog (e.g. threading.Timer + _thread.interrupt_main).
    def raise_timeout(signum, frame):
        raise TimeoutError("is_regular hung")

    old = signal.signal(signal.SIGALRM, raise_timeout)
    signal.alarm(60)
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("error", module="regfans")
            return fan.is_regular()
    except (TimeoutError, Warning) as e:
        print(f"[is_regular] {type(e).__name__}: {e} "
              f"(Npts={len(pts)}, Nsimps={len(simps)})", flush=True)
        return False
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)

        # fix ppl's rounding bug
        _libc.fesetround(0)


def is_regular_from_triangles(triangles) -> bool:
    """Convenience wrapper for the local samplers' output format: `triangles`
    is an iterable of 3-point sets/tuples (each point an (x, y) pair), as
    returned by sample_triangulation.sample().  Builds the (pts, simps) index
    form and calls is_regular."""
    index, pts, simps = {}, [], []
    for tri in triangles:
        idx = []
        for p in tri:
            p = tuple(p)
            if p not in index:
                index[p] = len(pts)
                pts.append(p)
            idx.append(index[p])
        simps.append(idx)
    return is_regular(np.asarray(pts, dtype=int), np.asarray(simps, dtype=int))
