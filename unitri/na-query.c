// CLI driver for the lattice-triangulation counter.  The counting code lives
// in na_query.h (STB single-header style); this translation unit pulls in its
// implementation and forwards to the entry point.  m and n are runtime
// arguments; the only build-time option is the arithmetic backend (-DGMP).
#define NA_QUERY_IMPLEMENTATION
#include "na_query.h"

int main(int argc, char **argv){
    return na_query_run(argc, argv);
}
