"""Build the unitri package.

The exact-count Cython extension (unitri.na_query) is built against libgmp.  If
libgmp isn't available the extension is skipped and the package still installs
-- the GMP-free mod-prime + CRT path (unitri.count_triangulations_parallel)
needs only a C compiler.
"""
import os
import subprocess
import sys
import tempfile

from setuptools import Extension, setup

# locate GMP via the shared helper (pkg-config, falling back to Homebrew/conda)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _gmp import gmp_dirs


def _have_gmp(inc, lib):
    """True if a trivial gmp.h program compiles and links.  When false the GMP
    extension is skipped -- install still succeeds and count_triangulations_parallel
    (mod-prime + CRT) works without GMP."""
    src = ("#include <gmp.h>\n"
           "int main(void){ mpz_t z; mpz_init(z); mpz_clear(z); return 0; }\n")
    cc = os.environ.get("CC", "cc")
    with tempfile.TemporaryDirectory() as d:
        csrc = os.path.join(d, "gmp_probe.c")
        with open(csrc, "w") as f:
            f.write(src)
        cmd = ([cc] + [f"-I{p}" for p in inc] + [f"-L{p}" for p in lib]
               + [csrc, "-lgmp", "-o", os.path.join(d, "gmp_probe")])
        try:
            subprocess.check_call(cmd, stdout=subprocess.DEVNULL,
                                  stderr=subprocess.DEVNULL)
            return True
        except (OSError, subprocess.CalledProcessError):
            return False


gmp_inc, gmp_lib = gmp_dirs()
ext_modules = []
if _have_gmp(gmp_inc, gmp_lib):
    from Cython.Build import cythonize
    ext_modules = cythonize(
        [Extension(
            "unitri.na_query",
            sources=["unitri/na_query.pyx"],
            include_dirs=["unitri"] + gmp_inc,
            library_dirs=gmp_lib,
            libraries=["gmp"],
            # rpath so the built extension finds libgmp without LD_LIBRARY_PATH
            extra_link_args=[f"-Wl,-rpath,{d}" for d in gmp_lib],
            # pull na_query.h's implementation in, big-integer back-end
            define_macros=[("NA_QUERY_IMPLEMENTATION", None), ("GMP", None)],
            extra_compile_args=["-O3"],
            language="c",
        )],
        compiler_directives={"language_level": "3"},
    )
else:
    sys.stderr.write(
        "unitri: libgmp not found -- installing WITHOUT the exact GMP extension; "
        "use unitri.count_triangulations_parallel (mod-prime + CRT, no GMP).\n")

setup(ext_modules=ext_modules)
