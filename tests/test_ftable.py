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
"""Flat-rectangle f(m,k) table mode (na-query.c run with no query on stdin).

This is the only path that prints the f(m,k) table; every other test sends a
query, so without this the table emitter is uncovered.  The values are exact
(GMP build) and cross-consistent -- f(m,k)=f(k,m), e.g. f(3,4)=f(4,3)=2822648 --
and match the known rectangle counts (f(3,2)=852, f(4,4)=736983568).
"""
import re
import subprocess

import pytest

# (m, n, [f(m,1), f(m,2), ..., f(m,n)])
TABLES = [
    (3, 5, [20, 852, 46456, 2822648, 182881520]),
    (4, 4, [70, 12170, 2822648, 736983568]),
]


def _ftable(binary, m, n):
    """Run in table mode (empty stdin) and return the list of f(m,k) values."""
    out = subprocess.run([binary, str(m), str(n)],
                         input="", capture_output=True, text=True).stdout
    return [int(v) for v in re.findall(r"f\(\d+,\d+\) = \*\) (\d+)", out)]


@pytest.mark.parametrize("m,n,expected", TABLES, ids=[f"{m}x{n}" for m, n, _ in TABLES])
def test_ftable_mode(na_query_bin, m, n, expected):
    assert _ftable(na_query_bin, m, n) == expected
