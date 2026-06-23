# na_query.pyx -- Cython binding for the lattice-triangulation counter.
#
# Wraps na_query.h's in-process API (na_query_count) so Python can count fine
# triangulations directly -- no subprocess, no stdout parsing.  Built against
# the big-integer (GMP) back-end, so the returned count is the exact integer.
#
# starting scaffold: the C API may still change.  If na_query_count's signature
# in na_query.h changes, update the `cdef extern` block below to match.

from libc.stdlib cimport malloc, free

# minimal GMP declarations: just enough to read the mpz_t result back out
cdef extern from "gmp.h":
    ctypedef struct __mpz_struct:
        pass
    ctypedef __mpz_struct mpz_t[1]
    void mpz_init(mpz_t)
    void mpz_clear(mpz_t)
    char *mpz_get_str(char *, int, const mpz_t)
    size_t mpz_sizeinbase(const mpz_t, int)

cdef extern from "na_query.h":
    int na_query_count(int m, int n, const int *upper, const int *lower,
                       mpz_t *out_count)


def na_query(int m, int n, upper, lower=None):
    """Count fine triangulations of the region between ``upper`` and an optional
    ``lower`` floor in an ``m x n`` box.  Returns the exact count as a Python int.

    Parameters
    ----------
    m, n : int
        Bounding-box width (``m >= 3``) and height.
    upper : sequence of int, length m+1
        Upper-boundary heights ``h_0 .. h_m``.
    lower : sequence of int, length m+1, optional
        Floor heights; omit (``None``) for a flat floor at 0.
    """
    if m < 3:
        raise ValueError(f"m must be >= 3, got {m}")
    if len(upper) != m + 1:
        raise ValueError(f"upper must have m+1={m + 1} heights, got {len(upper)}")
    if lower is not None and len(lower) != m + 1:
        raise ValueError(f"lower must have m+1={m + 1} heights, got {len(lower)}")

    cdef int i
    cdef int *U = <int*>malloc((m + 1) * sizeof(int))
    cdef int *L = NULL
    cdef mpz_t c
    cdef int st
    cdef size_t ndig
    cdef char *buf
    cdef bytes digits

    if U == NULL:
        raise MemoryError()
    try:
        for i in range(m + 1):
            U[i] = upper[i]
        if lower is not None:
            L = <int*>malloc((m + 1) * sizeof(int))
            if L == NULL:
                raise MemoryError()
            for i in range(m + 1):
                L[i] = lower[i]

        mpz_init(c)
        st = na_query_count(m, n, U, L, &c)
        if st != 0:
            mpz_clear(c)
            raise RuntimeError(f"na_query_count failed (status {st})")

        ndig = mpz_sizeinbase(c, 10) + 2          # digits + sign + NUL
        buf = <char*>malloc(ndig)
        if buf == NULL:
            mpz_clear(c)
            raise MemoryError()
        mpz_get_str(buf, 10, c)
        digits = buf                              # copy C string into Python bytes
        free(buf)
        mpz_clear(c)
        return int(digits.decode("ascii"))
    finally:
        free(U)
        if L != NULL:
            free(L)
