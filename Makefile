# Build the na-query CLI counter.  Two backends:
#
#   na-query-mod  default: counts modulo a prime.  Needs only a C compiler.
#                 Combine several primes with unitri/crt_combine.py (or, in
#                 parallel, unitri.count_triangulations_parallel) for the exact
#                 count -- this is the GMP-free path.
#   na-query      exact big-integer counts in one run.  Needs libgmp
#                 (apt install libgmp-dev / brew install gmp / conda-forge gmp).
#
# `make both` is recommended.  libgmp is a light but extra dependency; if you
# don't have it, `make na-query-mod` alone still gives exact counts via the
# mod-prime + CRT path.

CC     ?= gcc
CFLAGS ?= -O2

.PHONY: both clean
both: na-query na-query-mod

na-query-mod: unitri/na-query.c          # mod-prime build (no libgmp)
	$(CC) $(CFLAGS) -o $@ $<

na-query: unitri/na-query.c              # exact big-integer build (needs libgmp)
	$(CC) $(CFLAGS) $$(python3 _gmp.py) -DGMP -o $@ $< -lgmp

clean:
	rm -f na-query na-query-mod
