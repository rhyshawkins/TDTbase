#ifndef icosahedron_h
#define icosahedron_h

#include "manifold.h"

manifold_t *
icosahedron_create(int degree);

/*
 * Counting functions
 */
int 
icosahedron_nvertices(int depth);

int 
icosahedron_nedges(int depth);

int 
icosahedron_ntriangles(int depth);

/*
 * Edge subtended angle -> resolution at depth
 */
double
icosahedron_angle(int depth);



#endif /* icosahedron_h */
