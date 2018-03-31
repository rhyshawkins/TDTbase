#ifndef tetrahedron_h
#define tetrahedron_h

#include "manifold.h"

manifold_t *
tetrahedron_create(int degree);

/*
 * Counting functions
 */
int 
tetrahedron_nvertices(int depth);

int 
tetrahedron_nedges(int depth);

int 
tetrahedron_ntriangles(int depth);



#endif /* tetrahedron_h */
