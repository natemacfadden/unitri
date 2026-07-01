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
#ifndef NA_QUERY_H
#define NA_QUERY_H

/*
   Computation of the number of primitive lattice triangulations
   of a rectangle m x n (m is small, n is big)

   Author: Stepan Orevkov  http://picard.ups-tlse.fr/~orevkov

   Last update:  Feb 10, 2017

*/
/*
   Provenance and modifications
   ----------------------------
   The header above and the algorithm it describes are the original work of
   Stepan Orevkov, which is the base of all the other methods here.

   Michael Stepniczka and Nate MacFadden performed a minor cleanup of that
   original code.  Subsequently, Nate MacFadden, together
   with Claude Opus 4.8, performed a more systematic cleanup and generalization.

   This file is the query counter, with two arithmetic back-ends selected at
   compile time:

     default (no -DGMP) : evaluates the area-graded recurrence modulo a prime.
                          Combine several primes via the Chinese Remainder
                          Theorem (unitri/crt_combine.py) to recover the exact count.

     -DGMP              : arbitrary-precision big-integer arithmetic (GMP / mpz_t);
                          a single run yields the true integer, no CRT needed
                          (link with -lgmp).

   Either way, given a target upper-boundary profile (optionally over a
   lower/floor profile) on stdin, it reports that cell's count.
*/

// HEADER
// ======

// Value back-end (select with -DGMP before including): a count is a big-integer
// mpz_t (-DGMP) or an int residue mod a prime (default).  VAL is that type and
// VAL_INIT initializes one (mpz_init under GMP, a no-op otherwise) -- a caller of
// na_query_count declares `VAL c; VAL_INIT(c);` and passes &c.
#ifdef GMP
#include <gmp.h>
typedef mpz_t VAL;
#define VAL_INIT(x)  mpz_init(x)
#else
typedef int VAL;
#define VAL_INIT(x)  ((void)0)
#endif

/*
CLI driver.  Parses ``<m> <n> [prime_index]`` and an optional target profile on
stdin, runs the recurrence, and prints either the flat-square f-table (no query)
or the queried cell's count.  The actual counting is na_query_count below -- the
in-process API that this CLI and any binding both call.

Compile-time configuration (define before the implementation include):
    GMP : big-integer mpz_t arithmetic (link -lgmp); otherwise the recurrence is
          evaluated modulo a prime.

Parameters
----------
argc, argv : int, char**
    Command line.  argv[1] = m (width, >= 3), argv[2] = n (height); optional
    argv[3] = prime_index, which selects the modulus in the default build and is
    ignored under -DGMP.  An optional target profile is read from stdin.

Returns
-------
int
    Process exit status: 0 on success, 1 on a usage error, 2 if another call is
    already in progress.
*/
int na_query_run(int argc, char **argv);

// status codes returned by na_query_count (also documented in its Returns)
enum {
  NA_OK          =  0,
  NA_ERR_BAD_M   = -1,   // m not in [3, NA_MAX_M]
  NA_ERR_PROFILE = -2,   // invalid query (e.g. upper profile below the floor)
  NA_ERR_BUSY    = -3,   // a call is already in progress (not reentrant)
};

// largest supported width: the file-scope work arrays are sized [100] and
// accessed through +1-shifted pointers, so the maximum valid index is m = 98
#define NA_MAX_M 98

/*
Count fine (primitive) triangulations of the region between an upper profile and
an optional floor, inside an ``m x n`` box.  This is the in-process counting API
(the CLI na_query_run and the Cython binding both call it): it writes the count
through ``out_count`` and returns a status code.  Not reentrant -- a concurrent
or recursive call returns NA_ERR_BUSY (parallelize with separate processes).

Parameters
----------
m : int
    Bounding-box width; must be in [3, NA_MAX_M] (= [3, 98]).
n : int
    Bounding-box height; profile heights range over [0, n].
upper : const int* restrict
    Upper boundary: m+1 heights h_0 .. h_m.  Use n+1 as the sentinel for an
    absent vertex (the boundary passes between lattice points there).
lower : const int* restrict
    Floor: m+1 heights, or NULL for a flat floor at 0.  The upper profile must
    lie on or above the floor.
out_count : VAL* restrict
    Written with the count.  A VAL is an mpz_t (-DGMP) or an int residue mod a
    prime (default); declare ``VAL c; VAL_INIT(c);`` and pass ``&c``.

Returns
-------
int
    Status code:
        NA_OK          ( 0) : success
        NA_ERR_BAD_M   (-1) : m not in [3, NA_MAX_M]
        NA_ERR_PROFILE (-2) : invalid query (e.g. upper profile below the floor)
        NA_ERR_BUSY    (-3) : a call is already in progress (not reentrant)
*/
int na_query_count(int m, int n,
                   const int * restrict upper, const int * restrict lower,
                   VAL * restrict out_count);

// IMPLEMENTATION
// ==============
#ifdef NA_QUERY_IMPLEMENTATION

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdatomic.h> // single-call-at-a-time reentrancy guard (see na_query_run)

// internal back-end helpers.  VAL (the cell type) and VAL_INIT are public and
// live in the HEADER above; here we add ACC -- the per-cell recurrence
// accumulator, a wider long in the mod-p build to sum signed products before
// reducing -- and the macros used only inside the implementation.
#ifdef GMP
typedef mpz_t ACC;
#define VAL_SET_UI(x,u)        mpz_set_ui((x), (u))
#define ACC_ZERO(a)            mpz_set_ui((a), 0)
#define ACC_FMA(a,sign,cell)   do { if ((sign) > 0) mpz_add((a),(a),(cell)); \
                                    else             mpz_sub((a),(a),(cell)); } while (0)
#define STORE_RESULT(cell,a)   recurrence_value_with_lower_base((cell), (a))
#define PRINT_VAL(cell)        gmp_printf("%Zd", (cell))
#define VAL_OUT_SET(outp,src)  mpz_set(*(outp), (src))   // *outp is mpz_t -> ptr
#else
typedef long ACC;
#define VAL_SET_UI(x,u)        ((x) = (u))
#define ACC_ZERO(a)            ((a) = 0)
#define ACC_FMA(a,sign,cell)   ((a) += (sign) * (cell))
#define STORE_RESULT(cell,a)   ((cell) = recurrence_value_with_lower_base(a))
#define PRINT_VAL(cell)        printf("%d", (cell))
#define VAL_OUT_SET(outp,src)  (*(outp) = (src))
#endif


// problem config
// --------------
// Bounding-box width m (assumed WLOG <= height, for efficiency) and height n
// are set at run time -- see na_query_run, which reads them and derives
// n1 = n+1, n2 = n+2.  m >= 3 is required (checked there): the indexing uses
// height[m-4] and a shape-index table over height[1..m-3]; m >= 4 is the
// general case and m = 3 degenerates cleanly (no interior heights, one shape
// slot).
int m, n;    // bounding-box width / height
int n1, n2;  // n+1, n+2


// useful macros
// -------------
// min/max of two numbers
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))


// debugging flags
// ---------------
// debug flag
//#define DEBUG


#ifdef GMP
// big-integer arithmetic
// ----------------------
// this back-end keeps full arbitrary-precision counts: no modulus, no CRT --
// one run gives the whole integer
#else
// memory optimization
// -------------------
// program does calculations mod primes, for memory considerations. By doing
// calculations for N properly selected primes, you can then use CRT to get
// the full answer

// note: a single run uses one prime; it does not combine primes itself

int modulus;      // do calculations modulo modulus
int n_primes=200; // number of pre-computed primes to select from
int prime[200]=
   {29443, 29453, 29473, 29483, 29501, 29527, 29531, 29537, 29567, 29569,
    29573, 29581, 29587, 29599, 29611, 29629, 29633, 29641, 29663, 29669,
    29671, 29683, 29717, 29723, 29741, 29753, 29759, 29761, 29789, 29803,
    29819, 29833, 29837, 29851, 29863, 29867, 29873, 29879, 29881, 29917,
    29921, 29927, 29947, 29959, 29983, 29989, 30011, 30013, 30029, 30047,
    30059, 30071, 30089, 30091, 30097, 30103, 30109, 30113, 30119, 30133,
    30137, 30139, 30161, 30169, 30181, 30187, 30197, 30203, 30211, 30223,
    30241, 30253, 30259, 30269, 30271, 30293, 30307, 30313, 30319, 30323,
    30341, 30347, 30367, 30389, 30391, 30403, 30427, 30431, 30449, 30467,
    30469, 30491, 30493, 30497, 30509, 30517, 30529, 30539, 30553, 30557,
    30559, 30577, 30593, 30631, 30637, 30643, 30649, 30661, 30671, 30677,
    30689, 30697, 30703, 30707, 30713, 30727, 30757, 30763, 30773, 30781,
    30803, 30809, 30817, 30829, 30839, 30841, 30851, 30853, 30859, 30869,
    30871, 30881, 30893, 30911, 30931, 30937, 30941, 30949, 30971, 30977,
    30983, 31013, 31019, 31033, 31039, 31051, 31063, 31069, 31079, 31081,
    31091, 31121, 31123, 31139, 31147, 31151, 31153, 31159, 31177, 31181,
    31183, 31189, 31193, 31219, 31223, 31231, 31237, 31247, 31249, 31253,
    31259, 31267, 31271, 31277, 31307, 31319, 31321, 31327, 31333, 31337,
    31357, 31379, 31387, 31391, 31393, 31397, 31469, 31477, 31481, 31489,
    31511, 31513, 31517, 31531, 31541, 31543, 31547, 31567, 31573, 31583};
#endif

// global state
// ------------
// amount of allocated memory, used only for progress and summary output
long alloc_pointers=0;    // for pointers
long alloc_total=0;       // total instantaneous
long max_alloc_total=0;   // maximum instantaneous
long *alloc_memH;         // [2*m*n+1] byte size of each H area-layer value page

// table of subshape triangulation counts, graded by twice the shape area.
// Each value is a VAL (a big-integer mpz_t, or an int residue mod a prime);
// H[v] points into the value page memH[v], and old H pointer layers are reused
// once their volume is no longer needed.
//
// Index order: H[twice_area][shape_index][height_m2][height_m1][height_m]
// (cf. the [area][coord][m2][m1] note at H init). twice_area is the area grade
// (2*m*n+1 layers); shape_index is the flattened index over the interior boundary
// heights height[1..m-3] (built in shape_index_prefix); the last three levels are
// the heights of the three rightmost columns m-2, m-1, m. A cell is the count of
// fine triangulations of the shape with that area and profile, filled bottom-up
// in area by the one-step inclusion-exclusion recurrence below.
VAL *****H;       // [2*m*n+1] area layers; allocated once m,n are known
VAL **memH;       // [2*m*n+1] value pages
VAL memH0;
#ifdef GMP
long *memH_cells; // [2*m*n+1] number of mpz_t cells in each page (for clearing)
#endif

// size and next free offset for the contiguous H[twice_area] value page
long mem_page_size, allo;

// recurrence used below (inclusion-exclusion): f(S) = sum (-1)^(k-1) f(S'),
// where S' ranges over proper subshapes of S obtained by removing a set of k
// pairwise-disjoint S-maximal primitive tiles (k>=1; k=1 terms add, k=2
// subtract, k=3 add, ...). one-step only: this relation is not transitive
//
// work arrays use shifted pointers so that subshape_code[-1] and
// inclusion_sign[-1] are valid
// sentinel values in the lexicographic one-step subshape iterators
int inclusion_sign_array[100] = {-1};
int *inclusion_sign = &(inclusion_sign_array[1]);

int subshape_code_array[100] = {0};
int *subshape_code = &(subshape_code_array[1]);

int subshape_code_max[100];


// height[coord] gives the upper boundary height at x=coord
// twice_area_prefix stores twice the accumulated area contribution
// shape_index_prefix stores the flattened H index for height[1],...,height[m-3]
int twice_area_prefix[100];
int shape_index_prefix[100];
int subshape_index_delta[100];
int interpolated_height[100];


int subshape_right_vertex_array[100] = {-1};
int *subshape_right_vertex = &(subshape_right_vertex_array[1]);
int height_options_0[3], height_options_m1[3], height_options_m2[3];
int removed_area_increment[100];
int index_increment[100];
int subshape_index_increment[100];
int subshape_code_index_delta[100][3];
int segment_length[100];
int removed_twice_area_array[100] = {0};
int *removed_twice_area = &(removed_twice_area_array[1]);

// nearest present boundary vertices to the left and right; a deleted vertex is
// represented by height n1 and belongs to the long segment joining these
int left_vertex[100], right_vertex[100];
int left_idx,right_idx;
int left_height,right_height;


int twice_area;

// heights of the upper boundary at integer x-coordinates
int height[100];
int height_m;  // height[m]
int height_m1; // height[m-1]
int height_m2; // height[m-2]
int height_m_min, height_m_max;
int height_m1_min, height_m1_max;
int height_m2_min, height_m2_max;

int cursor, height_equation_value, height_0, shape_index;

// target profile query read from stdin; disabled when query_enabled is 0
int query_enabled = 0;
int query_found = 0;

// suppress all diagnostic output (the stderr progress line); set by the
// in-process API (na_query_count) so a library call is silent, left 0 by the CLI
int quiet = 0;
int query_area = -1;
int query_height[100];
int lower_enabled = 0;
int lower_is_flat = 1;
int lower_height[100];
VAL query_value;          // count at the query profile
#ifdef GMP
// reusable global mpz_t accumulator for the per-cell recurrence: a local would
// force an mpz_init/mpz_clear on every cell.  (The mod-p build instead declares
// its accumulator locally -- a plain long stays in a register through the hot
// loop, which a global would not.)
ACC recurrence_sum;
#endif

// non-flat floor breaks the x<->m-x reflection we normally use to store only
// the height_0 >= height_m half of each layer; set this to store the full
// layer and read cells directly
int need_full_table = 0;

// helper functions
// ----------------
#ifndef GMP
int positive_mod(long value);
#endif

// euclidean GCD algorithm
int gcd( int x, int y ){
  // ensure inputs have the right signs
  if (x<0) x=-x;
  if (y<0) y=-y;

  // edge case: if y==0, return x
  if (y==0) {
    return x;
  }

  int remainder = x%y;
  while( remainder > 0 ){
    x = y;
    y = remainder;
    remainder = x%y;
  }
  
  return y;
}


static int parse_query_height_token(const char *token, int *height){
  char *end = NULL;
  long value;

  if (strcmp(token, ".") == 0 || strcmp(token, "*") == 0) {
    *height = n1;
    return 1;
  }

  value = strtol(token, &end, 10);
  if (*token == '\0' || *end != '\0') return 0;
  if (value < 0 || value > n) return 0;

  *height = (int)value;
  return 1;
}

static void read_profile_token(
    const char *token,
    int *profile,
    int coord){
  if (!parse_query_height_token(token, &profile[coord])) {
    fprintf(stderr, "profile height must be an integer in [0,%d] or .\n", n);
    exit(1);
  }
}

static void validate_profile(const int *profile, const char *name){
  if (profile[0] == n1 || profile[m] == n1) {
    fprintf(stderr, "%s endpoints must be present heights\n", name);
    exit(1);
  }

  for (int coord=0; coord<=m; coord++) {
    int left = coord;
    int right;

    if (profile[coord] != n1) continue;

    while (left >= 0 && profile[left] == n1) left--;
    right = coord;
    while (right <= m && profile[right] == n1) right++;

    if (left < 0 || right > m) {
      fprintf(stderr, "%s absent vertices need present neighbors\n", name);
      exit(1);
    }

    for (; coord<right; coord++) {
      int numerator = profile[left]*(right-coord)
          + profile[right]*(coord-left);
      int denominator = right - left;

      if (numerator % denominator == 0) {
        fprintf(stderr,
            "%s absent vertex at x=%d is a lattice point on segment\n",
            name,
            coord);
        exit(1);
      }
    }
    coord = right - 1;
  }
}

static int profile_area(const int *profile){
  int area = 0;
  int left = 0;

  while (left < m) {
    int right = left + 1;

    while (right <= m && profile[right] == n1) right++;
    if (right > m) {
      fprintf(stderr, "invalid profile\n");
      exit(1);
    }

    area += (right - left) * (profile[left] + profile[right]);
    left = right;
  }

  return area;
}

static void profile_height_at(
    const int *profile,
    int coord,
    int *numerator,
    int *denominator){
  int left = coord;
  int right = coord;

  if (profile[coord] != n1) {
    *numerator = profile[coord];
    *denominator = 1;
    return;
  }

  while (left >= 0 && profile[left] == n1) left--;
  while (right <= m && profile[right] == n1) right++;

  if (left < 0 || right > m) {
    fprintf(stderr, "invalid absent vertex in profile comparison\n");
    exit(1);
  }

  *numerator = profile[left]*(right-coord)
      + profile[right]*(coord-left);
  *denominator = right - left;
}


static int profile_geq(const int *upper, const int *lower){
  for (int coord=0; coord<=m; coord++) {
    int upper_num;
    int upper_den;
    int lower_num;
    int lower_den;

    profile_height_at(upper, coord, &upper_num, &upper_den);
    profile_height_at(lower, coord, &lower_num, &lower_den);

    if ((long long)upper_num*lower_den < (long long)lower_num*upper_den) {
      return 0;
    }
  }

  return 1;
}

static int current_profile_matches(const int *profile){
  for (int coord=0; coord<=m; coord++) {
    if (height[coord] != profile[coord]) return 0;
  }

  return 1;
}

static int current_profile_above_lower(void){
  if (!lower_enabled) return 1;
  return profile_geq(height, lower_height);
}

static int lower_profile_is_zero(void){
  if (!lower_enabled) return 1;

  for (int coord=0; coord<=m; coord++) {
    if (lower_height[coord] != 0) return 0;
  }

  return 1;
}

#ifdef GMP
// Write the cell's count (into dst) for the current upper profile over the
// fixed floor lower_height[]:
#else
// Cell value for the current upper profile over the fixed floor lower_height[]:
#endif
//   upper == floor          -> 1   (empty region)
//   upper dips below floor  -> 0
//   otherwise               -> the accumulated recurrence sum
// A non-flat floor needs the full table (need_full_table): the original code
// kept only the height_0 >= height_m half via the x<->m-x reflection, which is
// wrong for an asymmetric floor.
//
// na_query counts triangulations of the region under the given profile. To
// count a point set's hull instead, the profile must trace the hull; that is
// what unitri.points_to_profiles builds (cross-checked against TOPCOM).
#ifdef GMP
static void recurrence_value_with_lower_base(mpz_t dst,
                                             const mpz_t recurrence_sum){
  if (lower_enabled && current_profile_matches(lower_height)) {
    mpz_set_ui(dst, 1);
    return;
  }
  if (!current_profile_above_lower()) {
    mpz_set_ui(dst, 0);
    return;
  }
  mpz_set(dst, recurrence_sum);
#else
static int recurrence_value_with_lower_base(long recurrence_sum){
  if (lower_enabled && current_profile_matches(lower_height)) return 1;
  if (!current_profile_above_lower()) return 0;

  return positive_mod(recurrence_sum);
#endif
}

static int profile_is_flat_present(const int *profile){
  int value = profile[0];

  if (value == n1) return 0;
  for (int coord=1; coord<=m; coord++) {
    if (profile[coord] != value) return 0;
  }

  return 1;
}

// read one line holding exactly m+1 height tokens into profile[].  returns 1 on
// success, or 0 on EOF / a blank line (meaning "none given").  a line with the
// wrong number of tokens is a hard error.  reading one whole line at a time is
// what makes the floor profile genuinely optional: a blank line (just Enter)
// ends input, whereas scanf("%s") would skip the newline and keep blocking.
static int read_profile_line(int *profile, const char *name){
  char line[4096];

  if (!fgets(line, sizeof line, stdin)) return 0;   // EOF -> none

  int count = 0;
  for (char *tok = strtok(line, " \t\r\n"); tok; tok = strtok(NULL, " \t\r\n")) {
    if (count >= m+1) {
      fprintf(stderr, "%s profile needs exactly %d heights\n", name, m+1);
      exit(1);
    }
    read_profile_token(tok, profile, count++);
  }

  if (count == 0) return 0;                          // blank line -> none
  if (count != m+1) {
    fprintf(stderr, "%s profile needs exactly %d heights\n", name, m+1);
    exit(1);
  }
  return 1;
}

// With a flat floor the recurrence exploits the x<->m-x reflection symmetry and
// stores only shapes with height_0 >= height_m, recovering the other half by
// reflection; a query whose profile has height_0 < height_m therefore lives in
// the unstored half and is never matched (it would report "not_found").  The
// triangulation count is reflection-invariant, so reflect such a query (and its
// flat, hence symmetric, floor) into the stored half.  A non-flat floor forces
// the full table (need_full_table) and no reflection applies.
static void reflect_query_into_stored_half(void){
  if (need_full_table || query_height[0] >= query_height[m]) return;
  for (int i=0, j=m; i<j; i++, j--) {
    int t = query_height[i]; query_height[i] = query_height[j]; query_height[j] = t;
    if (lower_enabled) {
      t = lower_height[i]; lower_height[i] = lower_height[j]; lower_height[j] = t;
    }
  }
}

// Finish configuring a query once query_height[] (and lower_height[] when
// lower_enabled) are populated: compute the area, derive the floor flags,
// require the query to lie on/above the floor, and fold into the stored half.
// Returns NA_OK, or NA_ERR_PROFILE if the query dips below the floor.  Shared
// by query_from_stdin and query_from_arrays so the two paths cannot drift.
// (Structural validation is NOT here: query_from_stdin's validate_profile
// exit()s, which the library path must avoid -- see query_from_arrays.)
static int finalize_query(void){
  query_area = profile_area(query_height);
  if (lower_enabled) {
    lower_is_flat   = profile_is_flat_present(lower_height);
    need_full_table = !lower_is_flat;
    if (!profile_geq(query_height, lower_height)) return NA_ERR_PROFILE;
  } else {
    lower_is_flat   = 1;
    need_full_table = 0;
  }
  reflect_query_into_stored_half();
  return NA_OK;
}

void query_from_stdin(void){
  // upper (query) profile; a blank line or EOF here means "no query at all".
  // (Input format is documented in the README; see na_query_run's usage line.)
  if (!read_profile_line(query_height, "query")) {
    query_enabled = 0;
    return;
  }
  query_enabled = 1;
  validate_profile(query_height, "query");

  // optional floor profile; a blank line or EOF means "no floor"
  if (read_profile_line(lower_height, "lower")) {
    lower_enabled = 1;
    validate_profile(lower_height, "lower");
  } else {
    lower_enabled = 0;
  }

  // shared finalization; its only failure is the query dipping below the floor
  if (finalize_query() != NA_OK) {
    fprintf(stderr, "query profile must lie on or above lower profile\n");
    exit(1);
  }
}

int current_profile_is_query(void){
  if (!query_enabled || twice_area != query_area) return 0;

  for (int coord=0; coord<=m; coord++) {
    if (height[coord] != query_height[coord]) return 0;
  }

  return 1;
}

#ifdef GMP
void maybe_record_query(const mpz_t value){
#else
void maybe_record_query(int value){
#endif
  if (!current_profile_is_query()) return;

  query_found = 1;
#ifdef GMP
  mpz_set(query_value, value);
#else
  query_value = value;
#endif
}

// print a boundary profile as space-separated heights ('.' = absent vertex)
void print_profile(const int *profile){
  for (int coord=0; coord<=m; coord++) {
    if (coord) putchar(' ');
    if (profile[coord] == n1) putchar('.');
    else printf("%d", profile[coord]);
  }
#ifndef GMP
}

int positive_mod(long value){
  return (value % modulus + modulus) % modulus;
#endif
}

int corner_subshape_code_max(int left, int mid, int right){
  if (mid <= 0) return 0;
  return 2*mid == left + right + 1 ? 2 : 1;
}

int has_nonprimitive_edge(int from_coord, int to_coord, int to_height){
  int delta_x = to_coord - from_coord;

  if (from_coord <= 0 || delta_x <= 1) return 0;
  return gcd(to_height - height[from_coord], delta_x) > 1;
}

void reset_shape_state(void){
  for (int coord=0; coord<=m; coord++) {
    height[coord] = 0;
    twice_area_prefix[coord] = 0;
    shape_index_prefix[coord] = 0;
  }
}

void reset_shape_and_subshape_state(void){
  for (int coord=0; coord<=m; coord++) {
    height[coord] = 0;
    twice_area_prefix[coord] = 0;
    subshape_code_max[coord] = 0;
    shape_index_prefix[coord] = 0;
  }
}

void set_boundary_neighbors(int absent_mask){
  left_vertex[0] = 0;
  left_vertex[1] = 0;
  for (int coord=2; coord<=m; coord++) {
    if (absent_mask & (1 << (coord-1))) {
      left_vertex[coord] = left_vertex[coord-1];
    } else {
      left_vertex[coord] = coord-1;
    }
  }

  right_vertex[m] = m;
  right_vertex[m-1] = m;
  for (int coord=m-2; coord>=0; coord--) {
    if (absent_mask & (1 << (coord+1))) {
      right_vertex[coord] = right_vertex[coord+1];
    } else {
      right_vertex[coord] = coord+1;
    }
  }
}

// return the height of the coord-th vertex in the current subshape
int subshape_height_at(int coord){
  switch(subshape_code[coord]){
    case 0: return height[coord];
    case 1: if (height[coord] <= n) {
              if (segment_length[coord]==2) return height[coord]-1;
              else return n1;
            }
            else return interpolated_height[coord];
    case 2: return n1;
  }
  return n1;
}

#ifndef GMP
int parse_prime_index(int argc, char *argv[]){
  int prime_index = 0;

  if (argc < 4) {
    printf("Using default modulus=prime[%d]=%d...\n",
           prime_index, prime[prime_index]);
    return prime_index;
  }

  if ((1 != sscanf(argv[3], "%d", &prime_index)) ||
      (prime_index < 0) ||
      (prime_index >= n_primes)) {
    printf("Usage: %s <m> <n> [prime_index]\n", argv[0]);
    printf("prime_index is an integer in range [0,%d) "
           "for computations mod prime[prime_index]\n",
           n_primes);
    exit(0);
  }

  return prime_index;
}

#endif
// ================
// main calculation
// ================
// The counter keeps its working state in file-scope globals, so it is not
// reentrant.  Both public entry points (na_query_run, na_query_count) take the
// atomic flag below, which refuses a second concurrent or recursive entry
// (returning a busy status) rather than letting two calls silently corrupt the
// shared tables -- e.g. a threaded caller that released the GIL.  Parallelize
// with separate processes instead.  (Belt-and-suspenders: a Cython binding also
// holds the GIL by default, serializing calls; this guard backstops anyone who
// releases it.)
static atomic_int na_query_busy = 0;

// Run the area-graded recurrence over the configured region.  Caller must have
// set m, n, n1, n2, the modulus (mod-p build), and the query-config globals
// (via query_from_stdin or query_from_arrays).  Allocates and frees the tables, leaves
// Print one row of the flat-square f(m,k) table -- only in table mode (no
// query), and only when twice_area completes a k*(2m) block; a no-op otherwise.
static void emit_ftable(long *n2pow){
  if (!query_enabled && twice_area%(2*m) == 0) {
    shape_index = 0;
    for (int coord=1; coord<=m-3; coord++) shape_index += n2pow[m-3-coord];
    int k = twice_area/(2*m);
    printf("(* f(%d,%d) = *) ", m, k);
    PRINT_VAL(H[twice_area][k*shape_index][k][k][k]);
    if (k<n) puts(","); else puts("};");
  }
}

// the queried count in query_value (when query_enabled), and prints the f-table
// only in non-query mode.  Internal: does not take the reentrancy guard -- the
// public entry points do.
static int na_query_compute(void){
  alloc_pointers = alloc_total = max_alloc_total = 0;   // reset diagnostics

  // area-graded tables, sized now that m,n are known (freed at the end)
  long n_areas = 2*m*n + 1;
  H          = (VAL *****)malloc(n_areas * sizeof *H);
  memH       = (VAL **)   calloc(n_areas, sizeof *memH);   // NULL-init: cleanup frees non-NULL pages
  alloc_memH = (long *)   malloc(n_areas * sizeof *alloc_memH);
#ifdef GMP
  memH_cells = (long *)   malloc(n_areas * sizeof *memH_cells);
  // query_value / recurrence_sum are reusable globals: init once, ever
  static int mpz_inited = 0;
  if (!mpz_inited) { mpz_init(query_value); mpz_init(recurrence_sum); mpz_inited = 1; }
#endif

  // helper variables
  // ----------------
  // cache powers of n+2
  long n2pow[m+1];

  n2pow[0] = 1;
  for (int coord=0; coord<m; coord++) {
    n2pow[coord+1]=n2pow[coord]*n2;
  }

  // cache shape_index_prefix[coord] += (n+2)^{m-3-coord}
  for (int coord=1; coord<m-2; coord++) {
    index_increment[coord] = n2pow[m-3-coord];
  }

  // allocate memory for H (just pointers to memH)
  // ---------------------------------------------
  // for twice_area>m, use the pre-allocated memory at H[twice_area-m-1]
  // H[twice_area] does not depend on H[twice_area-m-1]
  long mem;

  for (twice_area=0; twice_area<=m; twice_area++) {
    mem = n2pow[m-3]*sizeof(VAL ***);
    alloc_pointers += mem;
    H[twice_area] = (VAL ****)malloc( mem );

    for (int coord=0; coord<n2pow[m-3]; coord++) {
      mem = n2*sizeof(VAL **);
      alloc_pointers += mem;
      H[twice_area][coord] = (VAL ***)malloc( mem );

      for (height_m2=0; height_m2<n2; height_m2++) {
        mem = n2*sizeof(VAL *);
        alloc_pointers += mem;
        H[twice_area][coord][height_m2] = (VAL **)malloc( mem );
      }
    }
  }

  // all memory allocated (so far) is in pointers
  alloc_total = alloc_pointers;

  // subshape_code value 2 means that the vertex is absent in the subshape
  // (the [m]/height_0 corner only ever has codes 0/1, so it needs no [2] slot)
  height_options_m1[2] = height_options_m2[2] = n+1;

  // base case: the empty shape has one (trivial) triangulation
  VAL_INIT(memH0);
  memH[0] = &memH0;
#ifdef GMP
  memH_cells[0] = 1;
#endif
  // false positive: gcc -Wmaybe-uninitialized can't see through the three
  // malloc'd pointer levels (the H[0] skeleton is allocated for twice_area
  // 0..m at the top of this function), so it flags this write
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
  H[0][0][0][0] = memH[0];                                  // [area][coord][m2][m1]
#pragma GCC diagnostic pop
  VAL_SET_UI(H[0][0][0][0][0], lower_profile_is_zero() ? 1 : 0);  // [..][m]

  // =========
  // main loop
  // =========
  for (twice_area=1; twice_area<=2*m*n; twice_area++) {

    // on a query the f(m,k) table is suppressed, so show progress instead
    // (stderr, one rewriting line); no query: the table is the progress.
    // quiet suppresses it entirely (library calls via na_query_count)
    if (query_enabled && !quiet) {
      fprintf(stderr, "\r  computing area %d / %d ...", twice_area, query_area);
    }

    // for H[twice_area>m], just use memory allocated for H[twice_area-m-1]
    // (said memory is already allocated; H[twice_area-m-1] is no longer needed)
    if (twice_area>m) {
      H[twice_area] = H[twice_area-m-1];
    }

    // free unnecessary value memory
#ifdef GMP
    if (twice_area>m+1) {
      long old = twice_area-m-1;
      for (long i=0; i<memH_cells[old]; i++) mpz_clear(memH[old][i]);
      free(memH[old]);
      memH[old] = NULL;
#else
    if (twice_area>m+1) {
      free(memH[twice_area-m-1]);
      memH[twice_area-m-1] = NULL;
#endif
      alloc_total -= alloc_memH[twice_area-m-1];
    }

    // compute page size
    // -----------------
    allo = 0;
    reset_shape_state();

    // loop over all (m-3)-tuples height[1],...,height[m-3]
    // for height[coord] <= n
    while (1) {
      for (height_m2=0; height_m2<=n; height_m2++) {
        twice_area_prefix[m-2]=twice_area_prefix[m-3]+2*height_m2;
        if (twice_area_prefix[m-2] > twice_area) break;

        for (height_m1=0; height_m1<=n; height_m1++) {
          twice_area_prefix[m-1]=twice_area_prefix[m-2]+2*height_m1;
          if (twice_area_prefix[m-1] > twice_area) break;

          // restrict to range [0,n]
          height_m_min = MAX(0, twice_area - twice_area_prefix[m-1] - n);
          height_m_max = MIN(n, need_full_table
              ? (twice_area - twice_area_prefix[m-1])
              : (twice_area - twice_area_prefix[m-1])/2);

          if (height_m_min <= height_m_max) {
            allo += height_m_max-height_m_min+1;
          }
        } /* end of loop by height[m-1] */
      } /* end of loop by height[m-2] */

      // pass to next (m-3)-tuple height[1],...,height[m-3]
      // for height[coord] <= n
      for (cursor=m-3; cursor>0; cursor--) {
        if (height[cursor]<n && twice_area_prefix[cursor]<twice_area - 1) break;
      }
      if (cursor<=0) break;

      height[cursor]++;
      twice_area_prefix[cursor] += 2;
      shape_index_prefix[cursor] += index_increment[cursor];

      while (++cursor < m-2) {
        height[cursor]=0;
        twice_area_prefix[cursor]=twice_area_prefix[cursor-1];
        shape_index_prefix[cursor]=shape_index_prefix[cursor-1];
      }
    }

    // set boundary heights from the page-size loop
    height[m-1] = height_m1-1;
    height[m-2] = height_m2-1;
    height_options_m1[0] = height_m1-1;
    height_options_m1[1] = height_m1-2;
    height_options_m2[0] = height_m2-1;
    height_options_m2[1] = height_m2-2;


    // the coord-th bit of absent_mask is 1 <==> height[coord]==n1
    for (int absent_mask=2; absent_mask<(1<<m); absent_mask+=2) {
      set_boundary_neighbors(absent_mask);

      segment_length[0] = right_vertex[0];
      segment_length[m] = m - left_vertex[m];
      twice_area_prefix[0]=shape_index_prefix[0]=0;
      for (int coord=1; coord<m; coord++) {
        shape_index_prefix[coord] = shape_index_prefix[coord-1];
        twice_area_prefix[coord] = twice_area_prefix[coord-1];
        left_idx = left_vertex[coord];
        segment_length[coord] = right_vertex[coord] - left_idx;
        if( (absent_mask&(1<<coord)) == 0 ){
          height[coord] = 0;
          if (0 < left_idx && left_idx < coord-1 && height[left_idx]==0) {
            height[coord] = 1;
            twice_area_prefix[coord] += segment_length[coord];
            shape_index_prefix[coord] += index_increment[coord];
          } 
        } else {
          height[coord]=n1;
          shape_index_prefix[coord] += n1*index_increment[coord];
        }
      }

      // loop over   height[1], ...., height[m-3]
      while(1){
        if( absent_mask&(1<<(m-2)) ) {
          height_m2_min = height_m2_max = n1;
        } else {
          height_m2_min = 0;
          height_m2_max = MIN(
            n,
            (twice_area - twice_area_prefix[m-3]) / segment_length[m-2]);
        }
        
        for (height_m2=height_m2_min; height_m2<=height_m2_max; height_m2++) {
          if(height_m2 <= n) {
            left_idx = left_vertex[m-2];
            if (has_nonprimitive_edge(left_idx, m-2, height_m2)) continue;
          }
          height_options_m2[0]=height[m-2]=height_m2;
          if (segment_length[m-2] == 2) height_options_m2[1] = height_m2-1;

          if (height_m2<=n) {
            twice_area_prefix[m-2] = twice_area_prefix[m-3]
              + height_m2*segment_length[m-2];
            if (twice_area_prefix[m-2] > twice_area) continue;
          } else {
            twice_area_prefix[m-2] = twice_area_prefix[m-3];
          }

          if (absent_mask&(1<<(m-1))) {
            height_m1_min = height_m1_max = n1;
          } else {
            height_m1_min = 0;
            height_m1_max = MIN(
                n,
                (twice_area - twice_area_prefix[m-2]) / segment_length[m-1]);
          }

          for (height_m1=height_m1_min;
                  height_m1<=height_m1_max;
                  height_m1++) {
            if(height_m1 <= n) {
              left_idx = left_vertex[m-1];
              if (has_nonprimitive_edge(left_idx, m-1, height_m1)) continue;
            }

            height_options_m1[0]=height[m-1]=height_m1;
            if( segment_length[m-1] == 2 )height_options_m1[1] = height_m1-1;
            
            if (height_m1<=n) {
              twice_area_prefix[m-1] = twice_area_prefix[m-2]
                + height_m1*segment_length[m-1];
              if (twice_area_prefix[m-1] > twice_area) continue;
            } else {
              twice_area_prefix[m-1] = twice_area_prefix[m-2];
            }

            height_m_min = MAX(
                0,
                (twice_area - twice_area_prefix[m-1]
                 - n*segment_length[0] + segment_length[m] - 1)
                / segment_length[m]);
            height_m_max = MIN(
                n,
                (twice_area - twice_area_prefix[m-1])
                / (need_full_table ? segment_length[m]
                                   : segment_length[0] + segment_length[m]));

            if (height_m_min <= height_m_max) {
              allo += height_m_max - height_m_min + 1;
            }
          } /* end of loop by height[m-1] */
        } /* end of loop by height[m-2] */

         /* pass to next compatible (m-3)-tuple */

        loop1:
        for (cursor = m-3; cursor>0; cursor--)
          if (height[cursor] < n
              && twice_area_prefix[cursor]
                  <= twice_area - segment_length[cursor]) {
            break;
          }
        if (cursor <= 0) break;

                height[cursor]++;
        twice_area_prefix[cursor] += segment_length[cursor];

        shape_index_prefix[cursor] += index_increment[cursor];
        left_idx = left_vertex[cursor];

        if (has_nonprimitive_edge(left_idx, cursor, height[cursor])) {
          if (height[cursor] < n
              && twice_area_prefix[cursor]
                  <= twice_area - segment_length[cursor]) {
            height[cursor]++;
            twice_area_prefix[cursor] += segment_length[cursor];
            shape_index_prefix[cursor] += index_increment[cursor];
          } else {
            goto loop1;
          }
        }


        while (++cursor< m-2) {
          twice_area_prefix[cursor] = twice_area_prefix[cursor-1];
          if (height[cursor]==n1) {
            shape_index_prefix[cursor] = shape_index_prefix[cursor-1]
              + n1*index_increment[cursor];
            continue;
          }
          
          shape_index_prefix[cursor] = shape_index_prefix[cursor-1];
          height[cursor]=0;
          left_idx = left_vertex[cursor];
          
          if (left_idx >= 0
              && cursor-left_idx > 1
              && (height[left_idx] == 0
                  || (height[left_idx] > 1
                      && gcd(height[left_idx], cursor-left_idx) > 1))) {
            height[cursor]=1;
            twice_area_prefix[cursor] += segment_length[cursor];
            shape_index_prefix[cursor] += index_increment[cursor];
          }
        }
      }
    } /* end of loop by absent_mask */

    // allocate the H[twice_area] value page
    mem_page_size = allo;

     mem = mem_page_size*sizeof(VAL);

     alloc_memH[twice_area] = mem;
     alloc_total += mem;

     if( alloc_total > max_alloc_total ) max_alloc_total = alloc_total;

#ifdef GMP
    memH[twice_area] = (mpz_t *)malloc( mem );
    memH_cells[twice_area] = mem_page_size;
    for (long i=0; i<mem_page_size; i++) mpz_init(memH[twice_area][i]);
#else
    memH[twice_area] = (int *)malloc( mem );
#endif

     allo = 0;

     reset_shape_and_subshape_state();

     // fill H[twice_area] for shapes with all vertices present

     while (1) {
      shape_index = shape_index_prefix[m-3];

      for (int coord=2; coord<m-3; coord++) {
        subshape_code_max[coord] = corner_subshape_code_max(
          height[coord-1], height[coord], height[coord+1]);
      }


       for (height_m2=0; height_m2<=n; height_m2++) {
         twice_area_prefix[m-2] = twice_area_prefix[m-3] + 2*height_m2;
         if (twice_area_prefix[m-2] > twice_area) break;
         height[m-2] = height_m2;
         if (m >= 4)
           subshape_code_max[m-3] = corner_subshape_code_max(
             height[m-4], height[m-3], height_m2);
         height_options_m2[1] = (height_options_m2[0] = height_m2)-1;

         for (height_m1=0; height_m1<=n; height_m1++) {
           twice_area_prefix[m-1] = twice_area_prefix[m-2] + 2*height_m1;
           if (twice_area_prefix[m-1] > twice_area) break;
           height[m-1] = height_m1;
           subshape_code_max[m-2] = corner_subshape_code_max(
             height[m-3], height_m2, height_m1);
           height_options_m1[1] = (height_options_m1[0] = height_m1)-1;

           height_m_min = MAX(0, twice_area - twice_area_prefix[m-1] - n);
           height_m_max = MIN(n, need_full_table
               ? (twice_area - twice_area_prefix[m-1])
               : (twice_area - twice_area_prefix[m-1])/2);

           if( height_m_min <= height_m_max ){
                H[twice_area][shape_index][height_m2][height_m1] =
                  memH[twice_area] + allo - height_m_min;
                allo += (height_m_max - height_m_min + 1);
           }
           for( height_m=height_m_min; height_m<=height_m_max; height_m++ ){
             height[m] = height_m;
             height_0 = height[0]
                 = twice_area - twice_area_prefix[m-1] - height_m;
             height_options_0[0] = height_m;
             height_options_0[1] = height_m-1;
             subshape_code_max[m-1] = corner_subshape_code_max(
               height_m2, height_m1, height_m);
             subshape_code_max[1] = corner_subshape_code_max(
               height_0, height[1], height[2]);
             subshape_code_max[0] = (height_0>0 ? 1 : 0);
             subshape_code_max[m] = (height_m>0 ? 1 : 0);
#ifdef GMP

            ACC_ZERO(recurrence_sum);     // reset the reusable global accumulator
#else

            ACC recurrence_sum = 0;        // local: stays in a register
#endif

             for (int coord=0; coord<=m; coord++) {
                 subshape_code[coord] = 0;
                 removed_twice_area[coord] = 0;
                 subshape_index_delta[coord]=0;
                 inclusion_sign[coord]=-1;
                 subshape_index_increment[coord] =
                     (n+2 - height[coord])*index_increment[coord];
             }

             while(1){      /* loop over one-step subshapes */

                            /* subshape_code encodes one-step subshape a_prime
                               obtained by removing a set of disjoint S-maximal
                               primitive tiles; do not reuse this as
                               a transitive subshape test
                               0: a_prime[coord] = height[coord];
                               1: a_prime[coord] = height[coord] - 1;
                               2: a_prime[coord] = n1, vertex absent */

                 for (cursor=m; cursor>=0; cursor--) {
                     if (subshape_code[cursor-1] == 0
                         && subshape_code[cursor] < subshape_code_max[cursor]) {
                       break;
                     }
                 }
                 if(cursor<0)break;
                 if(cursor>0){
                     if (subshape_code[cursor]++) {
                       subshape_index_delta[cursor] +=
                           subshape_index_increment[cursor];
                       removed_twice_area[cursor]--;
                     }
                     else {
                       subshape_index_delta[cursor] -= index_increment[cursor];
                       removed_twice_area[cursor] += (cursor==m ? 1 : 2);
                     }
                 }
                 else {
                   subshape_code[0]++;
                   subshape_index_delta[0] = 0;
                   removed_twice_area[cursor]++;
                 }
                 inclusion_sign[cursor]=-inclusion_sign[cursor-1];
                 while (++cursor <= m) {
                  subshape_code[cursor]=0;
                  removed_twice_area[cursor]=removed_twice_area[cursor-1];
                  inclusion_sign[cursor]=inclusion_sign[cursor-1];
                  subshape_index_delta[cursor]=subshape_index_delta[cursor-1];
                }

                 if (need_full_table
                     || height_0 > height_m
                     || subshape_code[0] <= subshape_code[m]) {
                     ACC_FMA(recurrence_sum, inclusion_sign[m], H[twice_area - removed_twice_area[m]]
                            [shape_index + subshape_index_delta[m-3]]
                            [height_options_m2[subshape_code[m-2]]]
                            [height_options_m1[subshape_code[m-1]]]
                            [height_options_0[subshape_code[m]]]);
                 }
                 else {
                  // set encoded_index
                  int encoded_index, coord;

                  for (encoded_index=0, coord=m-1; coord > 2; coord--) {
                    if (subshape_code[coord] == 2) {
                      encoded_index = encoded_index * n2 + n1;
                    } else {
                      encoded_index = encoded_index*n2
                          + (height[coord] - subshape_code[coord]);
                    }
                  }

                  // increment recurrence_sum
                  int index1, index2;

                  if (subshape_code[1] == 2) {
                      index1 = n1;
                  } else {
                      index1 = height[1]-subshape_code[1];
                  }

                  if (subshape_code[2] == 2) {
                      index2 = n1;
                  } else {
                      index2 = height[2]-subshape_code[2];
                  }

                  ACC_FMA(recurrence_sum, inclusion_sign[m], H[twice_area - removed_twice_area[m]]
                         [encoded_index][index2][index1][height_0-1]);
                 }
             }

             STORE_RESULT(H[twice_area][shape_index][height_m2][height_m1][height_m], recurrence_sum);
             maybe_record_query(
                H[twice_area][shape_index][height_m2][height_m1][height_m]);

           } /* end of loop by height[m] */
         } /* end of loop by height[m-1] */
       } /* end of loop by height[m-2] */


       /* pass to next (m-3)-tuple */
       for (cursor=m-3; cursor>0; cursor--) {
         if (height[cursor] < n
             && twice_area_prefix[cursor] < twice_area - 1) {
           break;
         }
       }
       if( cursor<=0 )break;
       height[cursor]++;
      twice_area_prefix[cursor] += 2;
       shape_index_prefix[cursor] += index_increment[cursor];
       while (++cursor < m-2) {
         height[cursor] = 0;
         twice_area_prefix[cursor] = twice_area_prefix[cursor-1];
         shape_index_prefix[cursor] = shape_index_prefix[cursor-1];
       }

     }

     emit_ftable(n2pow);
     if (twice_area==2*n*m) break;

     // fill H[twice_area] for shapes with one or more absent upper-boundary
     // absent_mask marks missing vertices by setting height[coord] = n1

     for (int absent_mask=2; absent_mask<(1<<m); absent_mask+=2){

       set_boundary_neighbors(absent_mask);
       segment_length[0] = right_vertex[0];
       segment_length[m] = m - left_vertex[m];
       twice_area_prefix[0]=shape_index_prefix[0]=0;
       for (int coord=1; coord<m; coord++) {
           subshape_code_max[coord] = 0;
           shape_index_prefix[coord] = shape_index_prefix[coord-1];
        twice_area_prefix[coord] = twice_area_prefix[coord-1];
           left_idx = left_vertex[coord];
           segment_length[coord] = right_vertex[coord] - left_idx;
           if( (absent_mask&(1<<coord)) == 0 ){
               height[coord] = 0;
               if (0 < left_idx && left_idx < coord-1
                   && height[left_idx] == 0) {
                    height[coord] = 1;
                    twice_area_prefix[coord] += segment_length[coord];
                    shape_index_prefix[coord] += index_increment[coord];
               } 
           }
           else {
             height[coord] = n1;
             shape_index_prefix[coord] += n1*index_increment[coord];
           }
       }
       shape_index = shape_index_prefix[m-3];
       height_options_m1[2] = height_options_m2[2] = n1;
       if( segment_length[m]   > 1 )height_options_0[1]=n1;
       if( segment_length[m-1] > 2 )height_options_m1[1]=n1;
       if( segment_length[m-2] > 2 )height_options_m2[1]=n1;

       while (1) {

           if( absent_mask&(1<<(m-2)) ) height_m2_min = height_m2_max = n1;
           else{
             height_m2_min = 0;
             height_m2_max = MIN(
            n,
            (twice_area - twice_area_prefix[m-3]) / segment_length[m-2]);
           }
           for (height_m2=height_m2_min; height_m2<=height_m2_max; height_m2++){
             if(height_m2 <= n) {
              left_idx = left_vertex[m-2];
              if (has_nonprimitive_edge(left_idx, m-2, height_m2)) continue;
              } 
             height_options_m2[0]=height[m-2]=height_m2;
             if( segment_length[m-2] == 2 )height_options_m2[1] = height_m2-1;

             if( height_m2<=n ){twice_area_prefix[m-2] = twice_area_prefix[m-3]
              + height_m2*segment_length[m-2];
            if (twice_area_prefix[m-2] > twice_area) continue;}
             else twice_area_prefix[m-2] = twice_area_prefix[m-3];

             if( absent_mask&(1<<(m-1)) )height_m1_min = height_m1_max = n1;
             else{
               height_m1_min = 0;
               height_m1_max = MIN(
                n,
                (twice_area - twice_area_prefix[m-2]) / segment_length[m-1]);
             }

             for (height_m1=height_m1_min;
                  height_m1<=height_m1_max;
                  height_m1++) {
               if (height_m1 <= n) {
                 left_idx = left_vertex[m-1];
                 if (has_nonprimitive_edge(left_idx, m-1, height_m1)) continue;
               }
               height_options_m1[0]=height[m-1]=height_m1;
               if( segment_length[m-1] == 2 )height_options_m1[1] = height_m1-1;

               if (height_m1 <= n) {
                 twice_area_prefix[m-1] = twice_area_prefix[m-2]
                     + height_m1*segment_length[m-1];
                 if (twice_area_prefix[m-1] > twice_area) continue;
               }
               else twice_area_prefix[m-1] = twice_area_prefix[m-2];

//   hm*sm + h0*s0 = twice_area - twice_area_prefix[m-1]
//   height[0] >= height[m] >= 0
//   h0 = (twice_area - prefix - hm*sm) / s0
//   hm <= h0 <= n
//   lower and upper h0 bounds imply the height_m range below


               height_m_min = MAX(
                   0,
                   (twice_area - twice_area_prefix[m-1]
                    - n*segment_length[0] + segment_length[m] - 1)
                   / segment_length[m]);
               height_m_max = MIN(
                   n,
                   (twice_area - twice_area_prefix[m-1])
                   / (need_full_table ? segment_length[m]
                                      : segment_length[0] + segment_length[m]));

               if( height_m_min <= height_m_max ){
                  H[twice_area][shape_index][height_m2][height_m1] =
                  memH[twice_area] + allo - height_m_min;
                  allo += (height_m_max - height_m_min + 1);
               }
               for( height_m=height_m_min; height_m<=height_m_max; height_m++ ){
                left_idx=left_vertex[m];
                 if (has_nonprimitive_edge(left_idx, m, height_m)) continue;

                 height[m] = height_m;
                 height_equation_value = twice_area - twice_area_prefix[m-1]
                     - height_m*segment_length[m];
                 if (height_equation_value % segment_length[0]) continue;
                 if( height_equation_value < 0 )continue;
                 height_0 = height[0] = height_equation_value/segment_length[0];
                 right_idx = right_vertex[0];
                 if (right_idx > 1
                     && gcd(height[right_idx] - height_0, right_idx) > 1) {
                   continue;
                 }

                 height_options_0[1] = (height_options_0[0] = height_m)-1;
#ifdef GMP

                 ACC_ZERO(recurrence_sum);     // reset the reusable global accumulator
#else

                 ACC recurrence_sum = 0;        // local: stays in a register
#endif
                 subshape_code[0] = 0;
                 subshape_code[m] = 0;
                 subshape_code_max[0] = 0;
                 subshape_code_max[m] = 0;
                 removed_twice_area[0] = 0;
                 removed_twice_area[m] = 0;
                 subshape_index_delta[0] = 0;

                 inclusion_sign[0] = inclusion_sign[m] = -1;
                 subshape_right_vertex[0] = -1;
                 if (right_vertex[0] == 1) {
                   subshape_code_max[0] = (height_0 == 0 ? 0 : 1);
                   removed_area_increment[0] = 1;
                 }
                 if (left_vertex[m] == m-1) {
                  subshape_code_max[m] = (height_m==0 ? 0 : 1);
                  removed_area_increment[m] = 1;
                }

                 for (int coord=1; coord<m; coord++) {
                     inclusion_sign[coord] = -1;
                     subshape_right_vertex[coord] = -1;
                     int lenseg = segment_length[coord];
                     subshape_code[coord] = 0;
                     subshape_code_max[coord] = 0;
                     removed_twice_area[coord] = 0;
                     subshape_index_delta[coord] = 0;
                     if ((absent_mask & (1 << coord)) == 0) {
                         if (lenseg == 2) {
                             if( height[coord] > 0 ){
                                 subshape_code_max[coord] = 1;
                                 if (coord <= m-3) {
                                   subshape_code_index_delta[coord][0] =
                                       -index_increment[coord];
                                 }
                                 removed_area_increment[coord] = 2;
                                 if (2*height[coord] - 1
                                     == height[coord-1] + height[coord+1]) {
                                     subshape_code_max[coord] = 2;
                                     if (coord <= m-3) {
                                       subshape_code_index_delta[coord][1] =
                                           (n+1 - height[coord])
                                           * n2pow[m-3-coord];
                                     }
                                 }
                             }
                         }
                         else{                 // there is a long adj segment
                            left_idx = left_vertex[coord];
                            right_idx = right_vertex[coord];
                            left_height = height[left_idx];
                            right_height = height[right_idx];
                          /* 
                              det( lenseg             coord  -   left_idx  )
                                 ( height[right_idx] - height[left_idx]
                                   height[coord] - height[left_idx] )
                                                                   */
                             if (lenseg*(height[coord] - left_height)
                                 - (coord-left_idx)
                                   *(right_height-left_height) == 1) {
                                 subshape_code_max[coord] = 1;
                                 if (coord <= m-3) {
                                   subshape_code_index_delta[coord][0] =
                                       (n+1 - height[coord])
                                       * n2pow[m-3-coord];
                                 }
                                 removed_area_increment[coord] = 1;
                             }
                         }
                     }
                     else{                 // no vertex at coord
                         left_idx =  left_vertex[coord];
                         right_idx = right_vertex[coord];
                         left_height = height[left_idx];
                         right_height = height[right_idx];
                         height_equation_value =
                             (left_height*(right_idx-coord)
                              + right_height*(coord-left_idx)) / lenseg;
                      /* 
                          det( coord - left_idx     lenseg      )
                             ( height_equation_value - height[left_idx]
                             height[right_idx] - height[left_idx] )
                                                         */ 
                         if ((coord-left_idx)*(right_height-left_height)
                             - (height_equation_value-left_height)*lenseg
                             == 1) {
                             subshape_code_max[coord] = 1;
                             interpolated_height[coord] = height_equation_value;

                             if (coord <= m-3) {
                               subshape_code_index_delta[coord][0] =
                                   (height_equation_value - n - 1)
                                   * n2pow[m-3-coord];
                             }
                             else if (coord == m-2) {
                               height_options_m2[1] = height_equation_value;
                             }
                             else if (coord == m-1) {
                               height_options_m1[1] = height_equation_value;
                             }
                             removed_area_increment[coord] = 1;
                         }
                     }
                 }

                 while(1){ /* loop over one-step subshapes */

                           /* subshape_code encodes one-step subshape a_prime
                              obtained by removing a set of disjoint S-maximal
                              primitive tiles; do not reuse this as
                              a transitive subshape test
                              0: a_prime[coord] = height[coord];
                              1, no long segment:
                                      a_prime[coord] = height[coord] - 1;
                              1, long segment:
                                      a_prime[coord] = n1;
                              1, height[coord] = n1: a_prime[coord] < n1;
                              2, no adjacent long segment:
                                      long segment; then a_prime[coord] = n1 */

                     for (cursor=m; cursor>=0; cursor--) {
                       if (subshape_right_vertex[cursor-1] < cursor
                           && subshape_code[cursor]
                              < subshape_code_max[cursor]) {
                         break;
                       }
                     }
                     if(cursor<0)break;
                     subshape_right_vertex[cursor] = right_vertex[cursor];
                     if (cursor <= m-3 && cursor > 0) {
                       subshape_index_delta[cursor] =
                           subshape_index_delta[cursor-1]
                           + subshape_code_index_delta[cursor]
                               [subshape_code[cursor]];
                     }
                     removed_twice_area[cursor] = removed_twice_area[cursor-1]
                         + (subshape_code[cursor]++
                            ? 1 : removed_area_increment[cursor]);
                     inclusion_sign[cursor]=-inclusion_sign[cursor-1];
                     while( ++cursor <= m ){
                         subshape_right_vertex[cursor] =
                             subshape_right_vertex[cursor-1];
                         subshape_code[cursor] = 0;
                         removed_twice_area[cursor] =
                             removed_twice_area[cursor-1];
                         inclusion_sign[cursor] = inclusion_sign[cursor-1];
                         subshape_index_delta[cursor] =
                             subshape_index_delta[cursor-1];
                     }
                     if (need_full_table
                         || height_0 > height_m
                         || subshape_code[0] <= subshape_code[m]) {
                         ACC_FMA(recurrence_sum, inclusion_sign[m], H[twice_area - removed_twice_area[m]]
                                [shape_index + subshape_index_delta[m-3]]
                                [height_options_m2[subshape_code[m-2]]]
                                [height_options_m1[subshape_code[m-1]]]
                                [height_options_0[subshape_code[m]]]);
                     }
                     else{
                      // set encoded_index
                      int encoded_index,coord;

                      for (encoded_index=0,coord=m-1; coord>2; coord--) {
                        encoded_index = encoded_index*n2
                            + subshape_height_at(coord);
                      }

                      // increment recurrence_sum
                      ACC_FMA(recurrence_sum, inclusion_sign[m], H[twice_area - removed_twice_area[m]]
                             [encoded_index]
                             [subshape_height_at(2)]
                             [subshape_height_at(1)]
                             [height_0-1]);
                     }
                 }     // end of loop over subshape_code

                 STORE_RESULT(H[twice_area][shape_index][height_m2][height_m1][height_m], recurrence_sum);
             maybe_record_query(
                H[twice_area][shape_index][height_m2][height_m1][height_m]);

               } /* end of loop by height[m] */
             } /* end of loop by height[m-1] */
           } /* end of loop by height[m-2] */

           /* pass to next compatible (m-3)-tuple */

    loop:
           for (cursor=m-3; cursor>0; cursor--) {
             if (height[cursor] < n
                 && twice_area_prefix[cursor]
                    <= twice_area - segment_length[cursor]) {
               break;
             }
           }
           if (cursor<=0) break;
           height[cursor]++;
           subshape_code_max[cursor]=1;
           twice_area_prefix[cursor] += segment_length[cursor];
           shape_index_prefix[cursor] += index_increment[cursor];
           left_idx = left_vertex[cursor];
           if (has_nonprimitive_edge(left_idx, cursor, height[cursor])) {
               if( height[cursor]<n && twice_area_prefix[cursor]
                  <= twice_area - segment_length[cursor] ){
                   height[cursor]++;
                   subshape_code_max[cursor] = 1;
                   twice_area_prefix[cursor] += segment_length[cursor];
                   shape_index_prefix[cursor] += index_increment[cursor];
               }
               else goto loop;
           }
           while( ++cursor < m-2 ){
               twice_area_prefix[cursor] = twice_area_prefix[cursor-1];
               if (height[cursor]==n1) {
                shape_index_prefix[cursor] = shape_index_prefix[cursor-1]
              + n1*index_increment[cursor];
                continue;
              }
              
              shape_index_prefix[cursor] = shape_index_prefix[cursor-1];
              height[cursor] = subshape_code_max[cursor]=0;
              left_idx = left_vertex[cursor];
              if (left_idx >= 0 && cursor-left_idx > 1) {
                if (height[left_idx]==0) {
                  height[cursor]=1;
                  twice_area_prefix[cursor] += segment_length[cursor];
                  shape_index_prefix[cursor] += index_increment[cursor];
                } else if (height[left_idx] > 1
                           && gcd(height[left_idx], cursor-left_idx) > 1) {
                  height[cursor]=1;
                  twice_area_prefix[cursor] += segment_length[cursor];
                  shape_index_prefix[cursor] += index_increment[cursor];
                }
              }
           }
           shape_index = shape_index_prefix[m-3];
       }
     }
     // a query cell lives at twice_area == query_area; nothing larger feeds it,
     // so stop here once it is computed (the lower-boundary base cases it needs
     // all sit at smaller area and are already done)
     if (query_enabled && twice_area == query_area) break;
   } /* end of loop by twice_area */

   // free the per-area H pointer skeleton: only areas 0..m are distinct (later
   // areas alias them via H[ta] = H[ta-m-1]), so free those m+1 skeletons once
   for (long ta=0; ta<=m; ta++) {
     for (long coord=0; coord<n2pow[m-3]; coord++) {
       for (long m2=0; m2<n2; m2++) free(H[ta][coord][m2]);
       free(H[ta][coord]);
     }
     free(H[ta]);
   }
   // free the value pages still live: aged-out pages were freed + NULLed in the
   // loop, unreached areas stay NULL (calloc), and memH[0] is the static cell
   for (long ta=1; ta<n_areas; ta++) {
     if (!memH[ta]) continue;
#ifdef GMP
     for (long i=0; i<memH_cells[ta]; i++) mpz_clear(memH[ta][i]);
#endif
     free(memH[ta]);
   }
#ifdef GMP
   mpz_clear(memH0);   // the base-case cell memH[0], init'd at the top
#endif

   free(H);
   free(memH);
   free(alloc_memH);
#ifdef GMP
   free(memH_cells);
#endif
   return NA_OK;
}

// set the query-config globals from in-memory profiles -- the array-based,
// no-I/O, no-exit counterpart of query_from_stdin, used by na_query_count.
// upper has m+1 heights; lower has m+1 heights, or is NULL for a flat floor.
static int query_from_arrays(const int *upper, const int *lower){
  query_found = 0;
  // Validate heights at the API boundary so this stays the "no-exit" path: an
  // out-of-range value (notably n+1, which equals the absent-vertex marker n1)
  // would otherwise reach exit() in profile_area/profile_height_at.  The CLI
  // does the equivalent check while parsing each height token.
  for (int c=0; c<=m; c++)
    if (upper[c] < 0 || upper[c] > n) return NA_ERR_PROFILE;
  if (lower)
    for (int c=0; c<=m; c++)
      if (lower[c] < 0 || lower[c] > n) return NA_ERR_PROFILE;

  for (int c=0; c<=m; c++) query_height[c] = upper[c];
  query_enabled = 1;
  if (lower) {
    for (int c=0; c<=m; c++) lower_height[c] = lower[c];
    lower_enabled = 1;
  } else {
    lower_enabled = 0;
  }
  return finalize_query();
}

// In-process counting API (box_enum style): count fine triangulations of the
// region between `upper` and an optional `lower` floor in an m x n box.  Writes
// the count through *out_count (a VAL: mpz_t under -DGMP, else int -- pass &c
// for `VAL c; VAL_INIT(c);`) and returns a status code.  Not reentrant; a
// concurrent/recursive call returns NA_ERR_BUSY.
int na_query_count(int m_arg, int n_arg,
                   const int * restrict upper, const int * restrict lower,
                   VAL * restrict out_count){
  if (atomic_exchange(&na_query_busy, 1)) return NA_ERR_BUSY;
  if (m_arg < 3 || m_arg > NA_MAX_M) { atomic_store(&na_query_busy, 0); return NA_ERR_BAD_M; }
  m = m_arg;  n = n_arg;  n1 = n + 1;  n2 = n + 2;
  quiet = 1;   // a library call is silent (no stderr progress)
#ifndef GMP
  if (modulus == 0) modulus = prime[0];   // default modulus if caller set none
#endif
  int st = query_from_arrays(upper, lower);
  if (st != NA_OK)  { atomic_store(&na_query_busy, 0); return st; }
  na_query_compute();
  if (!query_found) { atomic_store(&na_query_busy, 0); return NA_ERR_PROFILE; }
  VAL_OUT_SET(out_count, query_value);
  atomic_store(&na_query_busy, 0);
  return NA_OK;
}

// CLI entry: parse <m> <n> [prime_index] + an optional profile on stdin, run the
// recurrence, and print the f-table or the queried cell's count.
int na_query_run( int argc, char *argv[] ){
  if (atomic_exchange(&na_query_busy, 1)) {
    fprintf(stderr, "na_query: concurrent/reentrant call detected; this build "
                    "is single-call-at-a-time. Parallelize with separate "
                    "processes.\n");
    return 2;          // busy: do not clear the flag -- the in-progress call owns it
  }

  // runtime configuration:  <m> <n> [prime_index]   (profile on stdin)
  if (argc < 3) {
    fprintf(stderr,
            "Usage: %s <m> <n> [prime_index]   (m>=3; target profile on stdin)\n",
            argv[0]);
    atomic_store(&na_query_busy, 0);
    return 1;
  }
  m = atoi(argv[1]);
  n = atoi(argv[2]);
  if (m < 3 || m > NA_MAX_M) {
    fprintf(stderr,
            "m must be in [3, %d] (widths m < 3 need special-case code; larger m "
            "exceeds the fixed work arrays); got m=%d\n", NA_MAX_M, m);
    atomic_store(&na_query_busy, 0);
    return 1;
  }
  n1 = n + 1;
  n2 = n + 2;
  quiet = 0;   // CLI: show the stderr progress line
#ifndef GMP
  modulus = prime[parse_prime_index(argc, argv)];
#endif

  // read the optional target profile (blocks on stdin), then announce the run
  query_from_stdin();
  if (query_enabled)
    printf("(* fine triangulations of the queried region in a %dx%d box *)\n", m, n);
  else
    printf("(* fine triangulations of the %dx%d rectangle *)\n", m, n);
#ifdef GMP
  printf("(* arbitrary-precision (big-integer) counts *)\n\n");
#else
  printf("(* calculations mod %d *)\n\n", modulus);
#endif
  fflush(stdout);

  na_query_compute();

  printf("(* Maximal memory allocation: %6.2lfg\n"
         " (including %6.2lfg for pointers) *)\n",
         ((double)max_alloc_total)/1e9, ((double)alloc_pointers)/1e9);
  if (query_enabled) {
    fprintf(stderr, "\n");   // finish the progress line
    if (query_found) {
      printf("(* query: top [");
      print_profile(query_height);
      printf("]");
      if (lower_enabled) {
        printf(" over floor [");
        print_profile(lower_height);
        printf("]");
      }
#ifdef GMP
      printf("  (big integer) *)\n");
#else
      printf("  (mod %d) *)\n", modulus);
#endif
      printf("query_value ");
      PRINT_VAL(query_value);
      printf("\n");
    } else {
      printf("query_value not_found\n");
    }
  }

  atomic_store(&na_query_busy, 0);
  return 0;
}

// internal macros are scoped to the implementation; the public VAL / VAL_INIT
// (defined in the HEADER) intentionally remain
#undef VAL_SET_UI
#undef ACC_ZERO
#undef ACC_FMA
#undef STORE_RESULT
#undef PRINT_VAL
#undef VAL_OUT_SET

#endif // NA_QUERY_IMPLEMENTATION
#endif // NA_QUERY_H
