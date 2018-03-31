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
