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
"""
Test suite for na-query.c (big-integer GMP back-end) against known counts.

Each case names a region (a profile of upper heights U over a floor L within an
m x n bounding box) and an expected exact count.  The script builds na-query.c
with -DGMP at the right (m, n), runs the query, and compares.

Sources for the expected values:
  1,3: https://arxiv.org/pdf/math/0211268
  2,4: https://arxiv.org/pdf/2602.16909
  6,7: TOPCOM
  5  : no independent check (expected value is conjectured)

Memory for this code's index skeleton grows like (n+2)^(m-1); cases needing
more than MEM_CAP are skipped (they require a big-memory machine).
"""

import os
import subprocess

GMP = subprocess.check_output(["brew", "--prefix", "gmp"]).decode().strip()
NA_QUERY_C = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "unitri", "na-query.c")
MEM_CAP = 8e9   # bytes; skip cases whose pointer skeleton exceeds this

# name, m, n, upper U, floor L (None = flat 0), expected count
CASES = [
    ("4x4 square",
     4, 4, [4, 4, 4, 4, 4], None, "736983568"),
    ("triangle (0,0),(3,0),(0,84)",          # arXiv 2602.16909, ~7.61e65
     3, 84, [84, 56, 28, 0], None,
     "761342982944289349099618507228200078481281500600912757801568059775"),
    ("4x10 square",
     4, 10, [10, 10, 10, 10, 10], None, "584455230176565718869688"),
    ("triangle (0,0),(7,0),(0,28)",
     7, 28, [28, 24, 20, 16, 12, 8, 4, 0], None,
     "153405520601065827395233041403916434565495619834253255337377215694753592286"),
    ("polygon L/U, box height 26",
     4, 26, [26, 22, 18, 13, 13], [0, 4, 8, 13, 13], "1189443740608860381225"),
    ("polygon L/U (TOPCOM 840021)",
     4, 4, [3, 4, 4, 4, 3], [3, 1, 0, 0, 1], "840021"),
    ("polygon L/U (TOPCOM 10653)",
     4, 3, [2, 3, 3, 3, 2], [2, 1, 0, 0, 1], "10653"),
]


def skeleton_bytes(m, n):
    n2 = n + 2
    return (m + 1) * (n2 ** (m - 3)) * (n2 + n2 * n2) * 8


def build():
    # m, n are runtime arguments now, so one binary serves every case
    out = "/tmp/na_test"
    subprocess.check_call(
        ["gcc", "-O2",
         f"-I{GMP}/include", f"-L{GMP}/lib",
         "-DGMP", "-o", out, NA_QUERY_C, "-lgmp"])
    return out


def run_query(binary, m, n, U, L):
    inp = " ".join(map(str, U)) + "\n"
    if L is not None:
        inp += " ".join(map(str, L)) + "\n"
    out = subprocess.run([binary, str(m), str(n)],
                         input=inp, capture_output=True, text=True).stdout
    for line in out.splitlines():
        if line.startswith("query_value"):
            return line.split()[1]
    return None


def verdict(expected, got):
    if got is None:
        return "ERROR (no query_value; run may have failed)"
    if expected.startswith("~"):                 # approximate expectation
        exp = float(expected[1:])
        ratio = int(got) / exp
        ok = 0.5 < ratio < 2.0
        return f"{'OK ' if ok else 'OFF'} (got {got}, ~{ratio:.3g}x the expected {expected})"
    return "OK" if got == expected else f"MISMATCH (got {got}, expected {expected})"


def main():
    binary = build()
    fails = 0
    for name, m, n, U, L, expected in CASES:
        mem = skeleton_bytes(m, n)
        if mem > MEM_CAP:
            print(f"[SKIP] {name}\n       needs ~{mem/1e9:.0f} GB (m={m}, n={n}); "
                  f"run on a big-memory machine")
            continue
        got = run_query(binary, m, n, U, L)
        v = verdict(expected, got)
        fails += not v.startswith("OK")
        print(f"[{'PASS' if v.startswith('OK') else 'FAIL'}] {name}\n       {v}")
    raise SystemExit(1 if fails else 0)


if __name__ == "__main__":
    main()
