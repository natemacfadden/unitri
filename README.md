# Counting primitive lattice triangulations of an m x n rectangle

Based on the original program by **Stepan Orevkov**
(<http://picard.ups-tlse.fr/~orevkov>), reworked and generalized by Nate
MacFadden (with Claude Opus 4.8); an earlier minor cleanup was by Michael
Stepniczka and Nate MacFadden.

> **If there are any bugs/issues in this code, assume they are due to Nate
> MacFadden's rework and *not* Stepan Orevkov's original code.**

## Files

| file | what it does |
|------|--------------|
| `orig.c` | Orevkov's original counter (unmodified); computes **modulo a prime**. |
| `na-query.c` | Reworked/generalized counter. Two compile-time back-ends: **modulo a prime** (default) or **exact** big integers (`-DGMP`, links libgmp). Optionally reads an upper (and lower/floor) boundary profile from stdin and reports that cell's count. |
| `crt_combine.py` | Combines per-prime residues from `na-query.c` (default build) into the exact count via the Chinese Remainder Theorem. |
| `check_topcom.py` | Independent cross-check of the floor logic against TOPCOM (via CYTools), on small convex regions. |
| `run_tests.py` | Test suite: checks several regions against known counts (literature / TOPCOM). |
| `profile.sh` | Reports wall time (min/mean over `ITERS` runs) and peak memory of a command. |
| `sample_triangulation.py` | Self-contained count / uniform-sample / enumerate of fine triangulations of a lattice polygon (general polygons; small only). |
| `regularity.py` | Checks whether a triangulation is regular, via `regfans` (`pip install regfans`). |
| `baseline.txt` | Reference outputs used to check the rework. |

## Build

Default build computes **modulo a prime** (`m, n` default to 5, 6):

```sh
gcc -O2 -o na-query na-query.c
```

Override the rectangle size with `-D` (width `m` must be >= 3; height `n` is
unconstrained):

```sh
gcc -O2 -Dm=4 -Dn=4 -o na-query na-query.c
```

Add `-DGMP` for the **exact** big-integer back-end (one run gives the true
count, no CRT needed); this links libgmp:

```sh
gcc -O2 -DGMP -o na-query na-query.c -lgmp
```

If GMP isn't on the default search path (e.g. Homebrew on macOS), point the
compiler at it with `-I` (headers) and `-L` (library). Apple Silicon paths are
shown; on Intel macOS use `/usr/local/...`:

```sh
gcc -O2 -DGMP -I/opt/homebrew/opt/gmp/include -L/opt/homebrew/opt/gmp/lib \
    -o na-query na-query.c -lgmp
```

## Querying a region (upper / lower boundaries)

`na-query` (either back-end) reads an optional target region from **stdin**, one
profile per line, and reports the count for that region:

- **Line 1 -- upper boundary** (the query): `m+1` heights `h_0 h_1 ... h_m`,
  each an integer in `[0, n]`, or `.` for an *absent* vertex (the boundary
  passes between lattice points there). The endpoints `h_0` and `h_m` must be
  present.
- **Line 2 -- lower boundary / floor** (optional): same format. A blank line or
  Ctrl-D leaves the floor flat at 0. The upper profile must lie on or above the
  floor.

With no input at all (an empty pipe / immediate Ctrl-D), the program skips the
query and just prints the flat-square `f(m,k)` table.

Examples (exact build, `-DGMP`, with `m=4, n=4`).

The full 4x4 square (floor defaults to flat 0) -- prints `query_value 736983568`:

```sh
echo "4 4 4 4 4" | ./na-query
```

A region over a non-flat floor -- prints `query_value 14032211`:

```sh
printf '4 4 4 4 4\n0 1 0 1 0\n' | ./na-query
```

An absent vertex (`.`) on the upper boundary -- prints `query_value 35`:

```sh
echo "0 . 3 . 0" | ./na-query
```

The result prints as `query_value <count>` -- exact with `-DGMP`, or the
residue mod the chosen prime in the default build (combine several primes with
`crt_combine.py` to recover the exact count).

