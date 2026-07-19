"""Executable checks of the worked examples in README.md, so the documented
inputs keep producing the documented outputs (doc drift then fails CI).

The expected values mirror the README examples verbatim -- if you change an
example there, update it here too.
"""
import pytest

import unitri
from _cli import run_query

# README "Python (Cython binding)" section: (name, m, n, upper, lower, count)
PY_EXAMPLES = [
    ("4x4 square",            4, 4,  [4, 4, 4, 4, 4], None,            736983568),
    ("base-3 triangle",       3, 12, [12, 8, 4, 0],   None,            668517487),
    ("4x4 over a non-flat floor", 4, 4, [4, 4, 4, 4, 4], [0, 1, 0, 1, 0], 14032211),
]

# README "Querying a region" section: (name, m, n, stdin, count).  These go
# through the CLI so the absent-vertex ('.') input path is covered too.
CLI_EXAMPLES = [
    ("full 4x4 square",        4, 4, "4 4 4 4 4\n",             736983568),
    ("over a non-flat floor",  4, 4, "4 4 4 4 4\n0 1 0 1 0\n",  14032211),
    ("absent vertex",          4, 4, "0 . 3 . 0\n",            35),
]


@pytest.mark.parametrize("name,m,n,upper,lower,expected", PY_EXAMPLES,
                         ids=[c[0] for c in PY_EXAMPLES])
def test_readme_python_example(name, m, n, upper, lower, expected):
    assert unitri.na_query(m, n, upper, lower) == expected


@pytest.mark.parametrize("name,m,n,stdin,expected", CLI_EXAMPLES,
                         ids=[c[0] for c in CLI_EXAMPLES])
def test_readme_cli_example(na_query_bin, name, m, n, stdin, expected):
    got = run_query(na_query_bin, m, n, stdin=stdin)
    assert got is not None, "no query_value (the run may have failed)"
    assert int(got) == expected
