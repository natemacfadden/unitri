# unitri

**Paper:** [Further Bounding the Kreuzer-Skarke Landscape](https://arxiv.org/abs/2602.16909) (arXiv:2602.16909)

Exactly counts the fine (unimodular) triangulations of a lattice polygon. For a
convex polygon, hand it the polygon's lattice points and get the exact count;
for more general regions it also takes an explicit upper/lower boundary
description. It counts *without enumerating*, so it reaches regions with
astronomically many triangulations that enumeration tools cannot.

Based on the original program by **Stepan Orevkov**
(<http://picard.ups-tlse.fr/~orevkov>), reworked and generalized by Nate
MacFadden (with Claude Opus 4.8); an earlier minor cleanup was by Michael
Stepniczka and Nate MacFadden.

> **If there are any bugs/issues in this code, assume they are due to Nate
> MacFadden's rework and *not* Stepan Orevkov's original code.**

## Install

```sh
pip install -e .          # builds the Cython extension; needs a C toolchain + libgmp
```

libgmp is **optional**: it powers the fast single-run exact counter
(`count_triangulations`, `na_query`). It comes from `apt install libgmp-dev`
(Debian/Ubuntu), `brew install gmp` (macOS), or `conda install -c conda-forge
gmp`, and the build finds it automatically (`pkg-config`, falling back to
Homebrew/conda via the bundled `_gmp.py`). **Without GMP**, `pip install` simply
skips that extension and you count via `count_triangulations_parallel`
(mod-prime + CRT, needs only a C compiler) -- see below.

## Counting a polygon's unimodular triangulations

Give `count_triangulations` the polygon's lattice points and it returns the
exact number of fine (unimodular) triangulations:

```python
import unitri

unitri.count_triangulations([(0,0), (4,0), (3,2), (1,2)])          # 140    (a trapezoid)
unitri.count_triangulations([(0,0), (4,0), (3,3), (1,3), (0,1)])   # 10843  (a pentagon)
```

There is also a command-line front end. It reads a point set from a file or
stdin in almost any bracket/comma/whitespace format -- a pasted numpy array,
`[x, y]` lists, or `x y` per line:

```sh
echo "[[0,0],[4,0],[3,3],[1,3],[0,1]]" | python -m unitri     # -> 10843
python -m unitri points.txt
```

The count is exact (arbitrary precision). The point-set path is for convex
regions; for non-convex regions -- valleys, concave tops -- describe them
directly with boundary profiles, below.

For large counts, or when you don't have GMP, count the same thing **in parallel
across cores** via the mod-prime + CRT path (needs only a C compiler, no libgmp):

```python
unitri.count_triangulations_parallel([(0,0), (4,0), (3,2), (1,2)])   # 140, in parallel
```

It runs the mod-prime counter for successive primes across cores and CRT-combines
them into the exact integer -- the GMP-free way to get exact counts, and often
faster than the single GMP run for very large ones.

## The counting core (C CLI and boundary profiles)

Under the point-set convenience is a single-header C counter, `na_query.h`
(stb-style), driven by an `m x n` bounding box and upper/lower boundary
*profiles*. This is the lower-level, more general interface -- it handles any
region between a lower and an upper boundary, including non-convex ones -- and
it needs no Python. `na-query.c` is a thin CLI that pulls in the implementation;
width `m` and height `n` are runtime arguments, so one binary handles every size:

```sh
gcc -O2 -o na-query unitri/na-query.c                 # default: counts modulo a prime
gcc -O2 -DGMP -o na-query unitri/na-query.c -lgmp     # big-integer: the whole count
```

If GMP isn't on the compiler's default path (Homebrew, or a conda env), splice
in the bundled locator's flags: `gcc -O2 $(python3 _gmp.py) -DGMP -o na-query unitri/na-query.c -lgmp`.
Or just use the `Makefile`: `make na-query-mod` (mod-prime), `make na-query`
(GMP), or `make both`.

Run `./na-query <m> <n> [prime_index]` and pipe the region to stdin, one profile
per line:

- **Line 1 -- upper boundary** (the query): `m+1` heights `h_0 h_1 ... h_m`, each
  an integer in `[0, n]`, or `.` for an *absent* vertex (the boundary passes
  between lattice points there). The endpoints `h_0` and `h_m` must be present.
- **Line 2 -- lower boundary / floor** (optional): same format; a blank line or
  Ctrl-D leaves the floor flat at 0. The upper profile must lie on or above it.

```sh
echo "4 4 4 4 4" | ./na-query 4 4                  # query_value 736983568  (full 4x4 square)
printf '4 4 4 4 4\n0 1 0 1 0\n' | ./na-query 4 4    # query_value 14032211   (over a non-flat floor)
echo "0 . 3 . 0" | ./na-query 4 4                  # query_value 35          (an absent vertex)
```

The result prints as `query_value <count>` -- the whole integer under `-DGMP`, or
a residue mod the chosen prime in the default build (combine several primes with
`unitri/crt_combine.py` to recover the exact count). With no input at all,
`na-query <m> <n>` prints the flat-rectangle `f(m,k)` table for `k = 1..n`.

The same profile interface is available in-process from Python as
`unitri.na_query(m, n, upper, lower=None)` -- what `count_triangulations` calls
under the hood; `upper`/`lower` are the `m+1` boundary heights, omit `lower` for
a flat floor at 0:

```python
unitri.na_query(4, 4, [4, 4, 4, 4, 4])              # 736983568  (the 4x4 square)
unitri.na_query(3, 12, [12, 8, 4, 0])               # 668517487  (a base-3 triangle)
unitri.na_query(4, 4, [4,4,4,4,4], [0,1,0,1,0])     # 14032211   (over a non-flat floor)
```

Run the tests with `pip install -e .[test]` (adds pytest, plus cytools for the
TOPCOM cross-checks) then `pytest tests/`.

## Performance

`na_query` counts triangulations with a dynamic-programming recurrence -- it never
enumerates them -- so its cost depends only on the bounding box `(m, n)`, not on
the (often astronomically large) number of triangulations. TOPCOM, by contrast,
enumerates one triangulation at a time, so its cost scales with the count and
cannot reach large regions at all.

Exact (GMP) build vs TOPCOM (via CYTools) on an Intel Core Ultra 7 270K Plus,
Ubuntu 26.04, gcc 15.2. Both columns are warmed up once, then reported as the
per-call mean ± stdev: `na_query` in process through the compiled extension (not
a subprocess), with auto-scaled repetitions over 7 batches; TOPCOM one run per
batch, its repeats stopping once their cumulative time exceeds a 60 s budget.
The `vertices` column gives each region's convex-hull corners; pass them to
`count_triangulations` (or `count_triangulations_parallel`) to reproduce.

| region | vertices | triangulations | `na_query` | TOPCOM |
|---|---|---|---|---|
| 3x2 rectangle | (0,0),(0,2),(3,0),(3,2) | 852 | 0.020 ± 0.000 ms | 0.04 ± 0.01 s |
| polygon | (0,2),(1,3),(2,0),(3,0),(3,3),(4,1),(4,2) | 10,653 | 0.190 ± 0.001 ms | 0.56 ± 0.00 s |
| polygon | (0,3),(1,1),(1,4),(2,0),(3,0),(3,4),(4,1),(4,3) | 840,021 | 0.517 ± 0.001 ms | 55.0 ± 0.3 s |
| 4x4 square | (0,0),(0,4),(4,0),(4,4) | 736,983,568 | 0.330 ± 0.004 ms | infeasible |
| 4x10 square | (0,0),(0,10),(4,0),(4,10) | ~5.8e23 | 13.10 ± 0.04 ms | infeasible |
| triangle, height 84 | (0,0),(0,84),(3,0) | ~7.6e65 | 1.26 ± 0.003 s | infeasible |

Counts agree exactly with TOPCOM wherever TOPCOM can finish. Reproduce with
`pip install -e . && python benchmarks/benchmark.py`.

## Organization

```
unitri/
├── unitri/
│   ├── profiles.py     # point set -> (m,n,upper,lower); count_triangulations (the main entry)
│   ├── na_query.pyx    # Cython binding: in-process na_query(m, n, upper, lower)
│   ├── na_query.h      # the counter: stb-style single header, mod-prime (default) or GMP (-DGMP)
│   ├── na-query.c      # thin CLI wrapper around na_query.h
│   ├── crt_combine.py  # combine the default build's per-prime residues into the exact count
│   ├── crt_parallel.py # GMP-free exact counts: parallel mod-prime runs + CRT
│   ├── __main__.py     # CLI: python -m unitri (count a lattice point set)
│   └── __init__.py
├── tests/
│   ├── test_topcom_convex.py    # convex polygons vs TOPCOM: curated + randomized (the main use-case)
│   ├── test_counts.py           # exact counts vs literature / TOPCOM; f(m,k) table mode
│   ├── test_mod_prime.py        # default mod-prime backend + crt_combine (CRT reconstruction)
│   ├── test_symmetry.py         # hard x<->m-x reflection cases
│   ├── test_unimodular.py       # unimodular invariance
│   ├── test_readme_examples.py  # the examples in this README, asserted
│   ├── check_topcom.py          # standalone TOPCOM cross-check (fixed convex regions)
│   ├── conftest.py              # shared fixtures (builds the GMP/mod-prime binaries)
│   ├── _cli.py                  # shared "run the na-query CLI" helper
│   ├── _topcom.py               # shared TOPCOM enumeration helper
│   └── transforms.py            # GL(2,Z) invariance helpers
├── benchmarks/
│   ├── benchmark.py             # na_query vs TOPCOM timing (the Performance table)
│   └── profile.sh               # wall-time + peak-memory profiler for any command
├── Makefile                     # build the na-query CLI (make both / na-query / na-query-mod)
├── pyproject.toml
├── setup.py
├── environment.yml              # conda env with the GMP backend (recommended)
├── environment-nogmp.yml        # conda env without GMP (mod-prime + CRT path)
├── _gmp.py                      # locate GMP (pkg-config -> Homebrew/conda); shared by setup.py + tests
├── MANIFEST.in
├── CITATION.cff
└── LICENSE
```

## License

[GPLv3](LICENSE). Copyright (c) 2026 Nate MacFadden.
