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
"""End-to-end test of the default mod-prime backend (na-query.c built without
-DGMP).  Every other shipped test compiles -DGMP, so the default -- and most
commonly built -- arithmetic path is otherwise untested.  Here we collect the
count modulo successive primes and reconstruct the exact integer via crt_combine,
checking it against a known value."""
from unitri.crt_combine import combine, PRIMES

from _cli import run_query

# 4x4 square: known exact count (cf. test_counts.py and crt_combine's docstring)
M, N, U, L = 4, 4, [4, 4, 4, 4, 4], None
EXACT = 736983568


def test_single_residue_matches(na_query_mod_bin):
    # the default backend should compute the count reduced modulo prime[0]
    tok = run_query(na_query_mod_bin, M, N, U, L, prime_index=0)
    assert tok is not None, "no query_value (the run may have failed)"
    assert int(tok) == EXACT % PRIMES[0]


def test_residues_crt_to_exact(na_query_mod_bin):
    # gather residues over enough primes that their product exceeds the count,
    # then CRT-reconstruct the exact integer
    k, prod = 0, 1
    while prod <= EXACT:
        prod *= PRIMES[k]
        k += 1
    congruences = []
    for i in range(k):
        tok = run_query(na_query_mod_bin, M, N, U, L, prime_index=i)
        assert tok is not None, f"no query_value for prime_index {i}"
        congruences.append((int(tok), PRIMES[i]))
    value, _ = combine(congruences)
    assert value == EXACT
