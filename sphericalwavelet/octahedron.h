#ifndef octahedron_h
#define octahedron_h

#include "manifold.h"

manifold_t *
octahedron_create(int degree);

void
octahedron_destroy(manifold_t *o);

int
octahedron_parent(manifold_t *o, int index);

int
octahedron_children(manifold_t *o, int index, int *child_indices);

/*
 * Printing functions
 */
void
octahedron_print_triangles(manifold_t *o, int depth);

/*
 * Counting functions
 */
int 
octahedron_nvertices(int depth);

int 
octahedron_nedges(int depth);

int 
octahedron_ntriangles(int depth);



#endif /* octahedron_h */
