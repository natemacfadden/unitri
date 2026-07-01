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
"""Shared helper to drive the na-query CLI in the cross-check tests.

Every cross-check reinvented "run the binary, read the `query_value` line"; this
is the one copy.  Returns the value as a string (callers ``int()`` it when they
want a number), so both the exact-integer and mod-prime-residue cases work.
"""
import subprocess


def run_query(binary, m, n, upper=None, lower=None, prime_index=None, stdin=None):
    """Run the na-query CLI once and return its ``query_value`` token as a
    string, ``None`` if the query is not found / no value is printed.

    Provide either ``upper`` (and optional ``lower``) as height sequences --
    ints or ``.``-style tokens both work -- or a raw ``stdin`` string for inputs
    the sequences can't express.  Set ``prime_index`` to choose the modulus in a
    default (mod-prime) build.
    """
    if stdin is None:
        stdin = " ".join(map(str, upper)) + "\n"
        if lower is not None:
            stdin += " ".join(map(str, lower)) + "\n"
    args = [binary, str(m), str(n)]
    if prime_index is not None:
        args.append(str(prime_index))
    out = subprocess.run(args, input=stdin, capture_output=True, text=True).stdout
    for line in out.splitlines():
        if line.startswith("query_value"):
            tok = line.split()[1]
            return None if tok == "not_found" else tok
    return None
