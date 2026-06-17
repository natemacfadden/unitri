"""Build the unitri Cython extension (the lattice-triangulation counter).

    python3 setup.py build_ext --inplace     # or:  pip install -e .

produces an importable `unitri.na_query` extension wrapping unitri/na_query.h
(built with the big-integer GMP back-end, so results are exact).  Requires
Cython and libgmp.
"""
import os
import subprocess

from setuptools import Extension, setup
from Cython.Build import cythonize

# locate GMP (Homebrew on macOS; otherwise assume it's on the default path)
gmp_inc, gmp_lib = [], []
try:
    prefix = subprocess.check_output(["brew", "--prefix", "gmp"]).decode().strip()
    gmp_inc, gmp_lib = [os.path.join(prefix, "include")], [os.path.join(prefix, "lib")]
except Exception:
    pass

setup(
    ext_modules=cythonize(
        [Extension(
            "unitri.na_query",
            sources=["unitri/na_query.pyx"],
            include_dirs=["unitri"] + gmp_inc,
            library_dirs=gmp_lib,
            libraries=["gmp"],
            # pull na_query.h's implementation into the extension, big-int back-end
            define_macros=[("NA_QUERY_IMPLEMENTATION", None), ("GMP", None)],
            extra_compile_args=["-O3"],
            language="c",
        )],
        compiler_directives={"language_level": "3"},
    )
)
