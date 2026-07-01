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
"""Unit tests for crt_combine.py -- the CRT step that reconstructs an exact
count from the default mod-prime backend's per-prime residues."""
import pytest

from unitri.crt_combine import combine, combine_pair, parse_residues, PRIMES


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
    # crt_combine.py's own docstring: these two residues reconstruct the 4x4
    # square's exact count, 736983568.
    value, modulus = combine([(25278, PRIMES[0]), (10602, PRIMES[1])])
    assert value == 736983568
    assert modulus == PRIMES[0] * PRIMES[1]


def test_roundtrip_known_count():
    # reduce a known exact count modulo the first few primes, then reconstruct
    # it -- enough primes that their product exceeds the count.
    n = 736983568
    k, prod = 0, 1
    while prod <= n:
        prod *= PRIMES[k]
        k += 1
    value, _ = combine([(n % PRIMES[i], PRIMES[i]) for i in range(k)])
    assert value == n
