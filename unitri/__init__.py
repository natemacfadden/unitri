# =============================================================================
#    Copyright (C) 2026  Nate MacFadden for the Liam McAllister Group
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
# =============================================================================
"""unitri -- counting unimodular triangulations of a lattice polygon.

The polygon is given by an m x n bounding box plus upper/lower boundary
profiles. Exposes the in-process counter from the Cython extension
(na_query.pyx, which wraps the C single-header na_query.h).  See `na_query`
for the entry point.
"""
from .na_query import na_query
from .profiles import (
    count_triangulations,
    points_to_profiles,
)

__all__ = [
    "count_triangulations",
    "na_query",
    "points_to_profiles",
]
