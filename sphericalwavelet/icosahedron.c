//
//    Spherical Subdivision/Wavelet library
//    
//    Copyright (C) 2014 - 2018 Rhys Hawkins
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "icosahedron.h"

#include "transform.h"

#include "spherical.h"

#include "slog.h"

manifold_t *
icosahedron_create(int degree)
{
  manifold_t *o;
  int i;
  double phi;
  double ax, ay, az, theta;
  double transform[16];
  
  o = manifold_create(degree,
		      icosahedron_nvertices,
		      icosahedron_nedges,
		      icosahedron_ntriangles);

  /*
   * Everything allocated, fill in data for depth 0 (12 vertices, 30 edges, 20 triangles)
   */
  phi = (1.0 + sqrt(5.0))/2.0;

  /* Poles */
  if (manifold_set_vertex(o, 0,
			  1.0, 0.0, phi, 0, MANIFOLD_NORMALIZE) < 0) {
    return NULL;
  }
  if (manifold_set_vertex(o, 1,
			  -1.0, 0.0, -phi, 0, MANIFOLD_NORMALIZE) < 0) {
    return NULL;
  }

  /* First ring */
  if (manifold_set_vertex(o, 2, 
			  phi, -1.0, 0.0, 0, MANIFOLD_NORMALIZE) < 0) {
    return NULL;
  }
  if (manifold_set_vertex(o, 3,
			  phi, 1.0, 0.0, 0, MANIFOLD_NORMALIZE) < 0) {
    return NULL;
  }
  if (manifold_set_vertex(o, 4,
			  0.0, phi, 1.0, 0, MANIFOLD_NORMALIZE) < 0) {
    return NULL;
  }
  if (manifold_set_vertex(o, 5,
			  -1.0, 0.0, phi, 0, MANIFOLD_NORMALIZE) < 0) {
    return NULL;
  }
  if (manifold_set_vertex(o, 6,
			  0.0, -phi, 1.0, 0, MANIFOLD_NORMALIZE) < 0) {
    return NULL;
  }

  /* Second ring */
  if (manifold_set_vertex(o, 7,
			  0.0, -phi, -1.0, 0, MANIFOLD_NORMALIZE) < 0) {
    return NULL;
  }
  if (manifold_set_vertex(o, 8,
			  1.0, 0.0, -phi, 0, MANIFOLD_NORMALIZE) < 0) {
    return NULL;
  }
  if (manifold_set_vertex(o, 9,
			  0.0, phi, -1.0, 0, MANIFOLD_NORMALIZE) < 0) {
    return NULL;
  }
  if (manifold_set_vertex(o, 10,
			  -phi, 1.0, 0.0, 0, MANIFOLD_NORMALIZE) < 0) {
    return NULL;
  }
  if (manifold_set_vertex(o, 11,
			  -phi, -1.0, 0.0, 0, MANIFOLD_NORMALIZE) < 0) {
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
			2, 8) < 0) {
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
			3, 9) < 0) {
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
			4, 10) < 0) {
    return NULL;
  }

  /* Edges Panel 3 */
  if (manifold_set_edge(o, 0, 9,
			5, 0) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 10,
			5, 6) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 11,
			5, 11) < 0) {
    return NULL;
  }

  /* Edges Panel 4 */
  if (manifold_set_edge(o, 0, 12,
			6, 0) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 13,
			6, 2) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 14,
			6, 7) < 0) {
    return NULL;
  }

  /* Edges Panel 5 */
  if (manifold_set_edge(o, 0, 15,
			7, 2) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 16,
			7, 8) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 17,
			7, 1) < 0) {
    return NULL;
  }
  
  /* Edges Panel 6 */
  if (manifold_set_edge(o, 0, 18,
			8, 3) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 19,
			8, 9) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 20,
			8, 1) < 0) {
    return NULL;
  }

  /* Edges Panel 7 */
  if (manifold_set_edge(o, 0, 21,
			9, 4) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 22,
			9, 10) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 23,
			9, 1) < 0) {
    return NULL;
  }

  /* Edges Panel 8 */
  if (manifold_set_edge(o, 0, 24,
			10, 5) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 25,
			10, 11) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 26,
			10, 1) < 0) {
    return NULL;
  }
  
  /* Edges Panel 9 */
  if (manifold_set_edge(o, 0, 27,
			11, 6) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 28,
			11, 7) < 0) {
    return NULL;
  }
  if (manifold_set_edge(o, 0, 29,
			11, 1) < 0) {
    return NULL;
  }
  
  /* Triangles Panel 0 */
  if (manifold_set_triangle(o, 0, 0,
			    0, 3, 2,
			    3, 1, 0) < 0) {
    return NULL;
  }
  if (manifold_set_triangle(o, 0, 1,
			    8, 2, 3,
			    2, 1, 18) < 0) {
    return NULL;
  }

  /* Triangles Panel 1 */
  if (manifold_set_triangle(o, 0, 2,
			    0, 4, 3,
			    6, 4, 3) < 0) {
    return NULL;
  }
  if (manifold_set_triangle(o, 0, 3,
			    9, 3, 4,
			    5, 4, 21) < 0) {
    return NULL;
  }

  /* Triangles Panel 2 */
  if (manifold_set_triangle(o, 0, 4,
			    0, 5, 4,
			    9, 7, 6) < 0) {
    return NULL;
  }
  if (manifold_set_triangle(o, 0, 5,
			    10, 4, 5,
			    8, 7, 24) < 0) {
    return NULL;
  }

  /* Triangles Panel 3 */
  if (manifold_set_triangle(o, 0, 6,
			    0, 6, 5,
			    12, 10, 9) < 0) {
    return NULL;
  }
  if (manifold_set_triangle(o, 0, 7,
			    11, 5, 6,
			    11, 10, 27) < 0) {
    return NULL;
  }

  /* Triangles Panel 4 */
  if (manifold_set_triangle(o, 0, 8,
			    0, 2, 6,
			    0, 13, 12) < 0) {
    return NULL;
  }
  if (manifold_set_triangle(o, 0, 9,
			    7, 6, 2,
			    14, 13, 15) < 0) {
    return NULL;
  }

  /* Triangles Panel 5 */
  if (manifold_set_triangle(o, 0, 10,
			    2, 8, 7,
			    2, 16, 15) < 0) {
    return NULL;
  }
  if (manifold_set_triangle(o, 0, 11,
			    1, 7, 8,
			    17, 16, 20) < 0) {
    return NULL;
  }

  /* Triangles Panel 6 */
  if (manifold_set_triangle(o, 0, 12,
			    3, 9, 8,
			    5, 19, 18) < 0) {
    return NULL;
  }
  if (manifold_set_triangle(o, 0, 13,
			    1, 8, 9,
			    20, 19, 23) < 0) {
    return NULL;
  }

  /* Triangles Panel 7 */
  if (manifold_set_triangle(o, 0, 14,
			    4, 10, 9,
			    8, 22, 21) < 0) {
    return NULL;
  }
  if (manifold_set_triangle(o, 0, 15,
			    1, 9, 10,
			    23, 22, 26) < 0) {
    return NULL;
  }

  /* Triangles Panel 8 */
  if (manifold_set_triangle(o, 0, 16,
			    5, 11, 10,
			    11, 25, 24) < 0) {
    return NULL;
  }
  if (manifold_set_triangle(o, 0, 17,
			    1, 10, 11,
			    26, 25, 29) < 0) {
    return NULL;
  }

  /* Triangles Panel 9 */
  if (manifold_set_triangle(o, 0, 18,
			    6, 7, 11,
			    14, 28, 27) < 0) {
    return NULL;
  }
  if (manifold_set_triangle(o, 0, 19,
			    1, 11, 7,
			    29, 28, 17) < 0) {
    return NULL;
  }


  /*
   * Generate a transform that rotates (1, 0, phi) to (0, 0, 1)
   */
  transform_compute_axis_angle(1.0, 0.0, phi,
			       0.0, 0.0, 1.0,
			       &ax, &ay, &az, &theta);
  transform_from_axis_angle(ax, ay, az, theta,
			    transform);
  

  /* 
   * Rotate all vertices to align poles 
   */
  for (i = 0; i < 12; i ++) {
    vertex3_transform(&(o->vertices[i]), transform);
  }

  /*
   * Check north and south poles
   */
  if (fabs(o->vertices[0].x) > 1.0e-9 ||
      fabs(o->vertices[0].y) > 1.0e-9 ||
      fabs(o->vertices[0].z - 1.0) > 1.0e-9) {
    ERROR("rotation failed to move vertex 0 to north pole (%g %g %g)",
	  o->vertices[0].x,
	  o->vertices[0].y,
	  o->vertices[0].z);
    return NULL;
  }
    
  if (fabs(o->vertices[1].x) > 1.0e-9 ||
      fabs(o->vertices[1].y) > 1.0e-9 ||
      fabs(o->vertices[1].z + 1.0) > 1.0e-9) {
    ERROR("rotation failed to move vertex 1 to south pole (%g %g %g)",
	  o->vertices[1].x,
	  o->vertices[1].y,
	  o->vertices[1].z);
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
icosahedron_nvertices(int depth)
{
  if (depth >= 0 && depth < 16) {
    return (10 * (1 << (2*depth))) + 2;
  }

  return -1;
}

int 
icosahedron_nedges(int depth)
{
  if (depth >= 0 && depth < 16) {
    return 30 * (1 << (2*depth));
  }

  return -1;
}

int 
icosahedron_ntriangles(int depth)
{
  if (depth >= 0 && depth < 16) {
    return 10 * (1 << (2*depth + 1));
  }

  return -1;
}

double
icosahedron_angle(int depth)
{
  manifold_t *m = icosahedron_create(depth);
  int v1, v2;

  double lon1, lat1;
  double lon2, lat2;
  
  if (m == NULL) {
    ERROR("failed to create icosahedron manifold");
    return -1.0;
  }

  /* First edge */
  v1 = m->edges[depth][0].a;
  v2 = m->edges[depth][0].b;

  
  vertex3_carttosph(m->vertices[v1].x, m->vertices[v1].y, m->vertices[v1].z,
		    &lon1, &lat1);
  
  vertex3_carttosph(m->vertices[v2].x, m->vertices[v2].y, m->vertices[v2].z,
		    &lon2, &lat2);

  manifold_destroy(m);

  return spherical_gcdist(lon1, lat1, lon2, lat2, 1.0);
}
