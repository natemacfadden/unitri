# Counting primitive lattice triangulations of an m × n rectangle

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
| `na-query.c` | Reworked/generalized counter, **modulo a prime**. Optionally reads an upper (and lower/floor) boundary profile from stdin and reports that cell's count. |
| `na-query-gmp.c` | Same as `na-query.c` but **exact** (GMP big integers); one run gives the true count, no CRT needed. |
| `crt_combine.py` | Combines per-prime residues from `na-query.c` into the exact count via the Chinese Remainder Theorem. |
| `baseline.txt` | Reference outputs used to check the rework. |

## Build

```sh
# m,n default to 5,6
gcc -O2 -o na-query na-query.c

# override the rectangle size
gcc -O2 -Dm=4 -Dn=4 -o na-query na-query.c

# needs the GMP library
gcc -O2 -o na-query-gmp na-query-gmp.c -lgmp
```

If GMP isn't on the default search path (e.g. Homebrew on macOS), point the
compiler at it with `-I` (headers) and `-L` (library):

```sh
# Apple Silicon paths shown; on Intel macOS use /usr/local/... instead
gcc -O2 -I/opt/homebrew/opt/gmp/include -L/opt/homebrew/opt/gmp/lib \
    -o na-query-gmp na-query-gmp.c -lgmp
```

