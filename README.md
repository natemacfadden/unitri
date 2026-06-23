# Counting unimodular triangulations of a lattice polygon

See https://arxiv.org/abs/2602.16909.

Exactly counts the unimodular (fine) triangulations of a lattice polygon -- given by an `m x n` bounding box plus upper/lower boundary profiles, so general regions and not just rectangles.

Based on the original program by **Stepan Orevkov**
(<http://picard.ups-tlse.fr/~orevkov>), reworked and generalized by Nate
MacFadden (with Claude Opus 4.8); an earlier minor cleanup was by Michael
Stepniczka and Nate MacFadden.

> **If there are any bugs/issues in this code, assume they are due to Nate
> MacFadden's rework and *not* Stepan Orevkov's original code.**

## Build

The counting code is the single-header `na_query.h` (stb-style); `na-query.c`
is a thin CLI that pulls in its implementation. Width `m` and height `n` are
**runtime arguments**, so one binary handles every size at runtime -- no
recompiling. The default build computes **modulo a prime**:

```sh
gcc -O2 -o na-query unitri/na-query.c
```

Add `-DGMP` for the **big-integer** back-end (one run gives the whole
count, no CRT needed); this links libgmp:

```sh
gcc -O2 -DGMP -o na-query unitri/na-query.c -lgmp
```

If GMP isn't on the default search path (e.g. Homebrew on macOS), point the
compiler at it with `-I` (headers) and `-L` (library). Apple Silicon paths are
shown; on Intel macOS use `/usr/local/...`:

```sh
gcc -O2 -DGMP -I/opt/homebrew/opt/gmp/include -L/opt/homebrew/opt/gmp/lib \
    -o na-query unitri/na-query.c -lgmp
```

## Querying a region (upper / lower boundaries)

Run `./na-query <m> <n> [prime_index]` (width `m >= 3`; `prime_index` selects
the modulus in the default build, ignored under `-DGMP`). The program reads an
optional target region from **stdin**, one profile per line, and reports the
count for that region:

- **Line 1 -- upper boundary** (the query): `m+1` heights `h_0 h_1 ... h_m`,
  each an integer in `[0, n]`, or `.` for an *absent* vertex (the boundary
  passes between lattice points there). The endpoints `h_0` and `h_m` must be
  present.
- **Line 2 -- lower boundary / floor** (optional): same format. A blank line or
  Ctrl-D leaves the floor flat at 0. The upper profile must lie on or above the
  floor.

With no input at all (an empty pipe / immediate Ctrl-D), the program skips the
query and just prints the flat-square `f(m,k)` table for `k = 1..n`.

Examples (big-integer build, `-DGMP`).

The full 4x4 square (floor defaults to flat 0) -- prints `query_value 736983568`:

```sh
echo "4 4 4 4 4" | ./na-query 4 4
```

A region over a non-flat floor -- prints `query_value 14032211`:

```sh
printf '4 4 4 4 4\n0 1 0 1 0\n' | ./na-query 4 4
```

An absent vertex (`.`) on the upper boundary -- prints `query_value 35`:

```sh
echo "0 . 3 . 0" | ./na-query 4 4
```

The result prints as `query_value <count>` -- the whole integer with `-DGMP`, or a
residue mod the chosen prime in the default build (combine several primes with
`unitri/crt_combine.py` to recover the exact count).

## Performance

`na_query` counts triangulations with a transfer-matrix recurrence -- it never
enumerates them -- so its cost depends only on the bounding box `(m, n)`, not on
the (often astronomically large) number of triangulations. TOPCOM, by contrast,
enumerates one triangulation at a time, so its cost scales with the count and
cannot reach large regions at all.

Exact (GMP) build vs TOPCOM (via CYTools) on an Intel i5-10600K, Ubuntu 24.04,
gcc 13.3 (`na_query` times are the min of 5 runs):

| region | triangulations | `na_query` | TOPCOM |
|---|---|---|---|
| 3x2 rectangle | 852 | 0.5 ms | 0.10 s |
| polygon (upper/lower) | 10,653 | 0.9 ms | 1.6 s |
| polygon (upper/lower) | 840,021 | 1.4 ms | 146 s |
| 4x4 square | 736,983,568 | 1.1 ms | infeasible |
| 4x10 square | ~5.8e23 | 22 ms | infeasible |
| triangle, height 84 | ~7.6e65 | 2.1 s | infeasible |

Counts agree exactly with TOPCOM wherever TOPCOM can finish. Reproduce with
`python benchmarks/benchmark.py`.

## Python (Cython binding)

`unitri/na_query.pyx` wraps the in-process counting API so you can count from
Python directly -- no subprocess, no stdout parsing. Build the extension (needs
Cython and libgmp):

```sh
pip install -e .                      # or: python3 setup.py build_ext --inplace
```

Then:

```python
import unitri
unitri.na_query(4, 4, [4, 4, 4, 4, 4])           # 736983568  (the 4x4 square)
unitri.na_query(3, 12, [12, 8, 4, 0])            # 668517487  (a base-3 triangle)
unitri.na_query(4, 4, [4,4,4,4,4], [0,1,0,1,0])  # 14032211   (over a non-flat floor)
```

`unitri.na_query(m, n, upper, lower=None)` returns the exact count as a Python
int (arbitrary precision). `upper`/`lower` are the `m+1` boundary heights, as in
the CLI; omit `lower` for a flat floor at 0.

To count an arbitrary lattice **point set** (not just a profile), use
`unitri.count_triangulations(points)`: it traces the convex hull in a
minimal-width orientation, builds the profile, and calls `na_query`.
`unitri.points_to_profiles(points)` returns that `(m, n, upper, lower)` without
counting. There is also a CLI -- `python -m unitri [points_file]` (reads stdin
if no file; accepts a pasted numpy array, `[x, y]` lists, or `x y` per line).

Run the tests with `pip install -e .[test]` (adds pytest, plus cytools for the
TOPCOM cross-checks) then `pytest tests/`.

## Organization

```
unitri/
├── unitri/
│   ├── na_query.h      # the counter: stb-style single header, mod-prime (default) or GMP (-DGMP)
│   ├── na-query.c      # thin CLI wrapper around na_query.h
│   ├── na_query.pyx    # Cython binding: in-process na_query(m, n, upper, lower)
│   ├── profiles.py     # lattice point set -> (m, n, upper, lower); count_triangulations
│   ├── crt_combine.py  # combine the default build's per-prime residues into the exact count
│   ├── __main__.py     # CLI: python -m unitri (count a lattice point set)
│   └── __init__.py
├── tests/
│   ├── test_counts.py        # exact counts vs literature / TOPCOM
│   ├── test_symmetry.py      # hard x<->m-x reflection cases
│   ├── test_unimodular.py    # unimodular invariance
│   ├── test_topcom_convex.py # randomized TOPCOM cross-check
│   ├── check_topcom.py       # standalone TOPCOM cross-check (small convex regions)
│   ├── conftest.py           # shared fixtures (builds the GMP binary)
│   ├── transforms.py         # GL(2,Z) invariance helpers
│   └── _gmp.py               # locate libgmp (Homebrew on macOS)
├── benchmarks/
│   ├── benchmark.py          # na_query vs TOPCOM timing (the Performance table)
│   └── profile.sh            # wall-time + peak-memory profiler for any command
├── pyproject.toml
├── setup.py
├── MANIFEST.in
├── CITATION.cff
└── LICENSE
```
