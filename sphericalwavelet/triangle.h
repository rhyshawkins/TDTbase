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

#ifndef triangle_h
#define triangle_h

#include "vertex3.h"

struct _triangle {
  /* Indices into vertex list */
  int a, b, c;
    
  /* Indices into edge list */
  int ab, bc, ca;

  int parent;
  int child_triangles[4];

  double area;
};
typedef struct _triangle triangle_t;

void
triangle_init(triangle_t *t);

double
triangle_area(triangle_t *t,
	      const vertex3_t *vertices);

int
triangle_centroid(triangle_t *t,
		  const vertex3_t *vertices,
		  double *x, double *y, double *z);


typedef struct triangle_point_in_triangle_workspace triangle_workspace_t;

triangle_workspace_t *
triangle_point_in_triangle_create_workspace(void);

void
triangle_point_in_triangle_free_workspace(triangle_workspace_t *w);

extern const double DEFAULT_TRIANGLE_EPSILON;

int
triangle_point_in_triangle(triangle_workspace_t *workspace,
			   triangle_t *t,
			   const vertex3_t *vertices,
			   double x, double y, double z,
			   double *ba, double *bb, double *bc,
			   double epsilon);

#endif /* triangle_h */
