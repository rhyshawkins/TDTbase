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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>

#include <gsl/gsl_linalg.h>

#include "triangle.h"

void
triangle_init(triangle_t *t)
{
  t->a = -1;
  t->b = -1;
  t->c = -1;

  t->ab = -1;
  t->bc = -1;
  t->ca = -1;

  t->parent = -1;
  t->child_triangles[0] = -1;
  t->child_triangles[1] = -1;
  t->child_triangles[2] = -1;
  t->child_triangles[3] = -1;
}

double
triangle_area(triangle_t *t,
	      const vertex3_t *vertices)
{
  double nx, ny, nz;

  vertex3_cross(vertices[t->b].x - vertices[t->a].x,
		vertices[t->b].y - vertices[t->a].y,
		vertices[t->b].z - vertices[t->a].z,
		vertices[t->a].x - vertices[t->c].x,
		vertices[t->a].y - vertices[t->c].y,
		vertices[t->a].z - vertices[t->c].z,
		&nx, &ny, &nz);

  return 0.5 * sqrt(nx*nx + ny*ny + nz*nz);
}

int
triangle_centroid(triangle_t *t,
		  const vertex3_t *vertices,
		  double *x, double *y, double *z)
{
  double cx, cy, cz;
  double l;
  
  cx = (vertices[t->a].x + vertices[t->b].x + vertices[t->c].x)/3.0;
  cy = (vertices[t->a].y + vertices[t->b].y + vertices[t->c].y)/3.0;
  cz = (vertices[t->a].z + vertices[t->b].z + vertices[t->c].z)/3.0;

  l = sqrt(cx*cx + cy*cy + cz*cz);
  if (l < 1.0e-12) {
    return -1;
  }
  
  cx /= l;
  cy /= l;
  cz /= l;

  *x = cx;
  *y = cy;
  *z = cz;

  return 0;
}

struct triangle_point_in_triangle_workspace {
  gsl_matrix *U;
  gsl_matrix *V;
  gsl_vector *S;
  gsl_vector *W;

  gsl_vector *b;
  gsl_vector *x;
};

triangle_workspace_t *
triangle_point_in_triangle_create_workspace(void)
{
  triangle_workspace_t *w;

  w = malloc(sizeof(triangle_workspace_t));
  if (w == NULL) {
    return NULL;
  }

  w->U = gsl_matrix_alloc(3,3);
  w->V = gsl_matrix_alloc(3,3);
  w->S = gsl_vector_alloc(3);
  w->W = gsl_vector_alloc(3);

  w->b = gsl_vector_alloc(3);
  w->x = gsl_vector_alloc(3);

  if (w->U == NULL ||
      w->V == NULL ||
      w->S == NULL ||
      w->W == NULL ||
      w->b == NULL ||
      w->x == NULL) {
    return NULL;
  }

  return w;
}

void
triangle_point_in_triangle_free_workspace(triangle_workspace_t *w)
{
  if (w != NULL) {
    gsl_vector_free(w->x);
    gsl_vector_free(w->b);
    
    gsl_vector_free(w->W);
    gsl_vector_free(w->S);
    gsl_matrix_free(w->V);
    gsl_matrix_free(w->U);
    
    free(w);
  }
}

static double truncate_barycentre_coordinate(double bc, double epsilon)
{
  if (fabs(bc) < epsilon) {
    return 0.0;
  } else {
    return bc;
  }
}

const double DEFAULT_TRIANGLE_EPSILON = 1.0e-14;

int
triangle_point_in_triangle(triangle_workspace_t *workspace,
			   triangle_t *t,
			   const vertex3_t *vertices,
			   double x, double y, double z,
			   double *ba, double *bb, double *bc,
			   double epsilon)
{
  double den;
  vertex3_t p;

  p.x = x;
  p.y = y;
  p.z = z;

  den = vertex3_determinant(&(vertices[t->a]),
			    &(vertices[t->b]),
			    &(vertices[t->c]));

  *ba = truncate_barycentre_coordinate(vertex3_determinant(&p, 
							   &(vertices[t->b]),
							   &(vertices[t->c]))/den,
				       epsilon);
  
  *bb = truncate_barycentre_coordinate(vertex3_determinant(&(vertices[t->a]), 
							   &p,
							   &(vertices[t->c]))/den,
				       epsilon);

  *bc = truncate_barycentre_coordinate(vertex3_determinant(&(vertices[t->a]),
							   &(vertices[t->b]),
							   &p)/den,
				       epsilon);

  return ((*ba) >= 0.0 && ((*bb) >= 0.0) && ((*bc) >= 0.0));
}
