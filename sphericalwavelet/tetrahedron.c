
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "tetrahedron.h"

#include "slog.h"

manifold_t *
tetrahedron_create(int degree)
{
  manifold_t *o;
  int i;
  double phi;
  
  o = manifold_create(degree,
		      tetrahedron_nvertices,
		      tetrahedron_nedges,
		      tetrahedron_ntriangles);

  /*
   * Everything allocated, fill in data for depth 0 (4 vertices, 6 edges, 4 triangles)
   */

  phi = 1.0/sqrt(2.0);

  if (manifold_set_vertex(o, 0,
			  -1.0, 0.0, phi, 0, -1) < 0) {
    return NULL;
  }
  if (manifold_set_vertex(o, 1,
			  0.0, -1.0, -phi, 0, -1) < 0) {
    return NULL;
  }
  if (manifold_set_vertex(o, 2,
			  1.0, 0.0, phi, 0, -1) < 0) {
    return NULL;
  }
  if (manifold_set_vertex(o, 3,
			  0.0, 1.0, -phi, 0, -1) < 0) {
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
			3, 2) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 5,
			3, 1) < 0) {
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
			    0, 2, 3,
			    0, 4, 3) < 0) {
    return NULL;
  }
  if (manifold_set_triangle(o, 0, 3,
			    1, 3, 2,
			    5, 4, 2) < 0) {
    return NULL;
  }

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
tetrahedron_nvertices(int depth)
{
  if (depth >= 0 && depth < 16) {
    return (1 << (2*depth + 1)) + 2;
  }

  return -1;
}

int 
tetrahedron_nedges(int depth)
{
  if (depth >= 0 && depth < 16) {
    return (1 << (2*depth + 1)) + (1 << (2*depth + 2));
  }

  return -1;
}

int 
tetrahedron_ntriangles(int depth)
{
  if (depth >= 0 && depth < 16) {
    return (1 << (2*depth + 2));
  }

  return -1;
}
