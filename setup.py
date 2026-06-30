"""Build the unitri Cython extension (the lattice-triangulation counter).

    python3 setup.py build_ext --inplace     # or:  pip install -e .

produces an importable `unitri.na_query` extension wrapping unitri/na_query.h
(built with the big-integer GMP back-end, so results are exact).  Requires
Cython and libgmp.
"""
import os
import sys

from setuptools import Extension, setup
from Cython.Build import cythonize

# locate GMP via the shared helper (pkg-config, falling back to Homebrew/conda)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _gmp import gmp_dirs

gmp_inc, gmp_lib = gmp_dirs()

setup(
    ext_modules=cythonize(
        [Extension(
            "unitri.na_query",
            sources=["unitri/na_query.pyx"],
            include_dirs=["unitri"] + gmp_inc,
            library_dirs=gmp_lib,
            libraries=["gmp"],
            # rpath so the built extension finds libgmp without LD_LIBRARY_PATH
            extra_link_args=[f"-Wl,-rpath,{d}" for d in gmp_lib],
            # pull na_query.h's implementation into the extension, big-int back-end
            define_macros=[("NA_QUERY_IMPLEMENTATION", None), ("GMP", None)],
            extra_compile_args=["-O3"],
            language="c",
        )],
        compiler_directives={"language_level": "3"},
    )
)
