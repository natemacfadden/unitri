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
"""GMP-free exact counting: parallel mod-prime runs + CRT reconstruction.

Each prime's count is an independent run of the *default* (mod-prime) na-query
build, so they run across cores; crt_combine then rebuilds the exact integer.
This needs only a C compiler -- no libgmp -- so it is both the parallel path for
large counts and the way to get exact counts when you don't have GMP.

``count_triangulations_parallel(points)`` is the entry ~everyone uses;
``count_parallel(m, n, upper, lower)`` is the (m, n, profile) form underneath.
"""
import os
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor

from .crt_combine import PRIMES, combine
from .profiles import _build_profile, _orientations

_NA_QUERY_C = os.path.join(os.path.dirname(os.path.abspath(__file__)), "na-query.c")
_mod_binary_path = None


def _mod_binary():
    """Compile the default mod-prime na-query binary once (no -DGMP, no libgmp)
    and cache its path."""
    global _mod_binary_path
    if _mod_binary_path and os.path.exists(_mod_binary_path):
        return _mod_binary_path
    out = os.path.join(tempfile.gettempdir(), "unitri-na-query-mod")
    subprocess.check_call([os.environ.get("CC", "cc"), "-O2", "-o", out, _NA_QUERY_C])
    _mod_binary_path = out
    return out


def _residue(binary, m, n, upper, lower, prime_index):
    """Count modulo prime[prime_index] via the mod-prime CLI; None if not found."""
    inp = " ".join(map(str, upper)) + "\n"
    if lower is not None:
        inp += " ".join(map(str, lower)) + "\n"
    out = subprocess.run([binary, str(m), str(n), str(prime_index)],
                         input=inp, capture_output=True, text=True).stdout
    for line in out.splitlines():
        if line.startswith("query_value"):
            tok = line.split()[1]
            return None if tok == "not_found" else int(tok)
    return None


def count_parallel(m, n, upper, lower=None, workers=None):
    """Exact count of the (m, n, upper, lower) region via parallel mod-prime runs
    + CRT -- no libgmp.  Runs primes across ``workers`` threads (each blocks on
    its own na-query process, so the C runs are genuinely parallel), adding
    primes until the reconstructed value stabilises (their product exceeds the
    count).  Raises ValueError if the region is not representable."""
    binary = _mod_binary()
    workers = workers or os.cpu_count() or 1
    residues, value, i = [], None, 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        while i < len(PRIMES):
            batch = list(range(i, min(i + workers, len(PRIMES))))
            got = list(pool.map(
                lambda j: _residue(binary, m, n, upper, lower, j), batch))
            if any(r is None for r in got):
                raise ValueError("region not representable (na-query: not_found)")
            residues.extend(got)
            i += len(batch)
            new_value = combine([(r, PRIMES[k]) for k, r in enumerate(residues)])[0]
            if new_value == value:
                return new_value          # stable: prime product exceeds the count
            value = new_value
    raise RuntimeError("exhausted tabulated primes before the count stabilised")


def count_triangulations_parallel(points, workers=None):
    """Exact fine-triangulation count of a lattice point set, GMP-free -- the
    parallel mod-prime + CRT counterpart of count_triangulations.  Traces the
    convex hull in a minimal-width orientation (retrying orientations na_query
    rejects) and counts across cores."""
    pts = sorted({(int(x), int(y)) for x, y in points})
    if len(pts) < 3:
        raise ValueError("need at least 3 distinct points")
    ors = _orientations(pts)
    if not ors:
        raise ValueError("point set is too thin: width < 3 in every orientation "
                         "(na_query needs m >= 3)")
    last = None
    for a, b in ors[:24]:
        m, n, upper, lower = _build_profile(pts, a, b)
        try:
            return count_parallel(m, n, upper, lower, workers=workers)
        except ValueError as e:          # this orientation unrepresentable; retry
            last = e
    raise ValueError(f"na_query rejected all tried orientations (last: {last})")
