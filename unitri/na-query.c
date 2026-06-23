/* =============================================================================
 *    Copyright (C) 2026  Nate MacFadden for the Liam McAllister Group
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 * ============================================================================= */
// CLI driver for the lattice-triangulation counter.  The counting code lives
// in na_query.h (STB single-header style); this translation unit pulls in its
// implementation and forwards to the entry point.  m and n are runtime
// arguments; the only build-time option is the arithmetic backend (-DGMP).
#define NA_QUERY_IMPLEMENTATION
#include "na_query.h"

int main(int argc, char **argv){
    return na_query_run(argc, argv);
}
