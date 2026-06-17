// CLI driver for the lattice-triangulation counter.  The counting code lives
// in na_query.h (STB single-header style); this translation unit pulls in its
// implementation and forwards to the entry point.  Build-time configuration
// (m, n, GMP) is set with -D on the compiler command line.
#define NA_QUERY_IMPLEMENTATION
#include "na_query.h"

int main(int argc, char **argv){
    return na_query_run(argc, argv);
}
