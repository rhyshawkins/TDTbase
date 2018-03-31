
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "octahedron.h"

#include "slog.h"

manifold_t *
octahedron_create(int degree)
{
  manifold_t *o;
  int i;

  o = manifold_create(degree,
		      octahedron_nvertices,
		      octahedron_nedges,
		      octahedron_ntriangles);
  
  /*
   * Everything allocated, fill in data for depth 0 (6 vertices, 12 edges, 8 triangles)
   */

  if (manifold_set_vertex(o, 0,
			  0.0, 0.0, 1.0, 0, 0) < 0) {
    return NULL;
  }
  if (manifold_set_vertex(o, 1,
			  0.0, 0.0, -1.0, 0, 0) < 0) {
    return NULL;
  }
  if (manifold_set_vertex(o, 2,
			  1.0, 0.0, 0.0, 0, 0) < 0) {
    return NULL;
  }
  if (manifold_set_vertex(o, 3,
			  0.0, 1.0, 0.0, 0, 0) < 0) {
    return NULL;
  }
  if (manifold_set_vertex(o, 4,
			  -1.0, 0.0, 0.0, 0, 0) < 0) {
    return NULL;
  }
  if (manifold_set_vertex(o, 5,
			  0.0, -1.0, 0.0, 0, 0) < 0) {
    return NULL;
  }

  /* Edges Panel 0 */
  if (manifold_set_edge(o, 0, 0,
			2, 0) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 1,
			2, 3) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 2,
			2, 1) < 0) {
    return NULL;
  }

  /* Edges Panel 1 */
  if (manifold_set_edge(o, 0, 3,
			3, 0) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 4,
			3, 4) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 5,
			3, 1) < 0) {
    return NULL;
  }

  /* Edges Panel 2 */
  if (manifold_set_edge(o, 0, 6,
			4, 0) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 7,
			4, 5) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 8,
			4, 1) < 0) {
    return NULL;
  }

  /* Edges Panel 3 */
  if (manifold_set_edge(o, 0, 9,
			5, 0) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 10,
			5, 2) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 11,
			5, 1) < 0) {
    return NULL;
  }
  
  /* Triangles Panel 0 */
  if (manifold_set_triangle(o, 0, 0,
			    0, 3, 2,
			    3, 1, 0) < 0) {
    return NULL;
  }
  if (manifold_set_triangle(o, 0, 1,
			    1, 2, 3,
			    2, 1, 5) < 0) {
    return NULL;
  }

  /* Triangles Panel 1 */
  if (manifold_set_triangle(o, 0, 2,
			    0, 4, 3,
			    6, 4, 3) < 0) {
    return NULL;
  }

  if (manifold_set_triangle(o, 0, 3,
			    1, 3, 4,
			    5, 4, 8) < 0) {
    return NULL;
  }

  /* Triangles Panel 2 */
  if (manifold_set_triangle(o, 0, 4,
			    0, 5, 4,
			    9, 7, 6) < 0) {
    return NULL;
  }

  if (manifold_set_triangle(o, 0, 5,
			    1, 4, 5,
			    8, 7, 11) < 0) {
    return NULL;
  }

  /* Triangles Panel 3 */
  if (manifold_set_triangle(o, 0, 6,
			    0, 2, 5,
			    0, 10, 9) < 0) {
    return NULL;
  }
  
  if (manifold_set_triangle(o, 0, 7,
			    1, 5, 2,
			    11, 10, 2) < 0) {
    return NULL;
  }

  /*
   * For remaining levels, subdivide
   */
  for (i = 1; i <= degree; i ++) {
    if (manifold_subdivide(o, i) < 0) {
      ERROR("failed to subdivide level %d", i);
      return NULL;
    }
  }

  for (i = degree; i >= 0; i --) { 
    if (manifold_compute_areas(o, i) < 0) {
      ERROR("failed to compute area for level %d", i);
      return NULL;
    }
  }

  if (manifold_build_neighbors(o) < 0) {
    ERROR("failed to build neighbors");
    return NULL;
  }
  
  return o;
}

int
octahedron_parent(manifold_t *o, int vertex_index)
{
  return -1;
}

int
octahedron_children(manifold_t *o, int index, int *child_indices)
{
  return -1;
}

void
octahedron_print_triangles(manifold_t *o, int depth)
{
  int i;

  if (o == NULL ||
      depth < 0 || depth > o->degree) {
    return;
  }

  for (i = 0; i < o->ntriangles[depth]; i ++) {
    printf("%d: %d %d %d (%d %d %d)\n",
	   i,
	   o->triangles[depth][i].a,
	   o->triangles[depth][i].b,
	   o->triangles[depth][i].c,
	   o->triangles[depth][i].ab,
	   o->triangles[depth][i].bc,
	   o->triangles[depth][i].ca);
  }
}

int 
octahedron_nvertices(int depth)
{
  if (depth >= 0 && depth < 16) {
    return (1 << (2*depth + 2)) + 2;
  } 

  return -1;
}

int 
octahedron_nedges(int depth)
{
  if (depth >= 0 && depth < 16) {
    return (1 << (2*depth + 2)) + (1 << (2*depth + 3));
  } 

  return -1;
}

int 
octahedron_ntriangles(int depth)
{
  if (depth >= 0 && depth < 16) {
    return 1 << (2*depth + 3);
  } 

  return -1;
}


