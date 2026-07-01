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
"""The default mod-prime backend and its CRT reconstruction.

Two layers, both on the default (no -DGMP) path that every other test skips:
  * unit tests for crt_combine.py -- the CRT that rebuilds an exact count from
    per-prime residues; and
  * end to end: build na-query.c without -DGMP, collect the count modulo
    successive primes, and CRT-reconstruct the exact integer -- checked against
    a known value.
"""
import pytest

from unitri.crt_combine import PRIMES, combine, combine_pair, parse_residues
from _cli import run_query

# 4x4 square: known exact count (cf. test_counts.py and crt_combine's docstring)
M, N, U, L = 4, 4, [4, 4, 4, 4, 4], None
EXACT = 736983568


# --- crt_combine unit tests -------------------------------------------------
def test_combine_pair_basic():
    # x = 2 (mod 3), x = 3 (mod 5)  ->  x = 8 (mod 15)
    assert combine_pair(2, 3, 3, 5) == (8, 15)


def test_combine_three_coprime():
    # x = 2 (mod 3), 3 (mod 5), 2 (mod 7)  ->  23 (mod 105)
    assert combine([(2, 3), (3, 5), (2, 7)]) == (23, 105)


def test_combine_inconsistent_raises():
    # x = 0 (mod 4) and x = 1 (mod 6): gcd 2 does not divide (1 - 0) -> impossible
    with pytest.raises(ValueError):
        combine_pair(0, 4, 1, 6)


def test_parse_residues_strips_comments_and_blanks():
    lines = ["25278   # mod prime[0]", "", "# a full-line comment", "10602"]
    assert parse_residues(lines) == [25278, 10602]


def test_documented_4x4_example():
    # crt_combine.py's docstring: these two residues reconstruct the 4x4 count
    value, modulus = combine([(25278, PRIMES[0]), (10602, PRIMES[1])])
    assert value == EXACT
    assert modulus == PRIMES[0] * PRIMES[1]


def test_roundtrip_known_count():
    # reduce the known count modulo the first few primes, then reconstruct it
    k, prod = 0, 1
    while prod <= EXACT:
        prod *= PRIMES[k]
        k += 1
    value, _ = combine([(EXACT % PRIMES[i], PRIMES[i]) for i in range(k)])
    assert value == EXACT


# --- default backend end to end (residues -> CRT -> exact count) ------------
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
