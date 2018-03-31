
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vertex3.h"

#include "slog.h"

void
vertex3_initialize(vertex3_t *v)
{
  v->x = 0.0;
  v->y = 0.0;
  v->z = 0.0;

  v->depth = -1;
  
  v->parent = -1;

  v->children[0] = -1;
  v->children[1] = -1;
  v->children[2] = -1;
  v->children[3] = -1;

  v->v[0] = -1;
  v->v[1] = -1;
  
  v->f[0] = -1;
  v->f[1] = -1;

  v->e[0] = -1;
  v->e[1] = -1;
  v->e[2] = -1;
  v->e[3] = -1;
  
  v->n[0] = -1;
  v->n[1] = -1;
  v->n[2] = -1;
  v->n[3] = -1;
  v->n[4] = -1;
  v->n[5] = -1;

  v->area = 0.0;
}

double
vertex3_dot(const vertex3_t *a,
	    const vertex3_t *b)
{
  return a->x*b->x + a->y*b->y + a->z*b->z;
}

void
vertex3_subtract(const vertex3_t *a,
		 const vertex3_t *b,
		 vertex3_t *m)
{
  m->x = a->x - b->x;
  m->y = a->y - b->y;
  m->z = a->z - b->z;
}

void
vertex3_midpoint(const vertex3_t *a,
		 const vertex3_t *b,
		 vertex3_t *m)
{
  m->x = (a->x + b->x)/2.0;
  m->y = (a->y + b->y)/2.0;
  m->z = (a->z + b->z)/2.0;
}

int
vertex3_normalize(vertex3_t *a)
{
  double l = sqrt(vertex3_dot(a, a));

  if (l > 0.0) {
    a->x /= l;
    a->y /= l;
    a->z /= l;

    return 0;
  }

  return -1;
}

void
vertex3_transform(vertex3_t *a,
		  double *matrix)
{
  double v[4];
  double w[4];
  int i;
  int j;

  v[0] = a->x;
  v[1] = a->y;
  v[2] = a->z;
  v[3] = 1.0;

  w[0] = 0.0;
  w[1] = 0.0;
  w[2] = 0.0;
  w[3] = 0.0;

  for (j = 0; j < 4; j ++) {
    for (i = 0; i < 4; i ++) {
      w[j] += matrix[4*j + i]*v[i];
    }
  }

  a->x = w[0]/w[3];
  a->y = w[1]/w[3];
  a->z = w[2]/w[3];
}

int
vertex3_add_child(vertex3_t *a,
		  int ci)
{
  int i;
  
  if (a == NULL || ci < 0) {
    ERROR("invalid parameters (%p %d)", a, ci);
    return -1;
  }

  for (i = 0; i < 4; i ++) {
    if (a->children[i] < 0) {
      a->children[i] = ci;
      return 0;
    }
  }

  ERROR("no empty slots");
  return -1;
}

int
vertex3_add_neighbor(vertex3_t *v,
		     int ni)
{
  int i;

  for (i = 0; i < 6; i ++) {

    if (v->n[i] < 0) {
      v->n[i] = ni;
      return 0;
    }

    if (v->n[i] == ni) {
      ERROR("neighbor already present");
      return -1;
    }
  }

  ERROR("neighbors full (adding %d)", ni);
  for (i = 0; i < 6; i ++ ) {
    ERROR("  %d", v->n[i]);
  }
  return -1;
}

void
vertex3_cross(double x1, double y1, double z1,
	      double x2, double y2, double z2,
	      double *nx, double *ny, double *nz)
{
  *nx = y1*z2 - z1*y2;
  *ny = z1*x2 - x1*z2;
  *nz = x1*y2 - y1*x2;
}

void
vertex3_sphtocart(double lon, double lat,
		  double *x, double *y, double *z)
{
  double theta;
  double phi;

  theta = lon * M_PI/180.0;
  phi = lat * M_PI/180.0;
  
  *x = cos(theta)*cos(phi);
  *y = sin(theta)*cos(phi);
  *z = sin(phi);
}

void
vertex3_carttosph(double x, double y, double z,
		  double *lon, double *lat)

{
  double r;

  r = sqrt(x*x + y*y + z*z);

  *lat = asin(z/r) * 180.0/M_PI;
  *lon = atan2(y, x) * 180.0/M_PI;
}

int
vertex3_midpoint_plane(const vertex3_t *a,
		       const vertex3_t *b,
		       double pa, double pb, double pc, double pd)
{
  vertex3_t t;
  vertex3_t n;

  vertex3_midpoint(a, b, &t);

  vertex3_subtract(a, b, &n);
  if (vertex3_normalize(&n) < 0) {
    return -1;
  }

  pa = n.x;
  pb = n.y;
  pc = n.z;

  pd = -vertex3_dot(&t, &n);
  
  return 0;
}

double
vertex3_determinant(const vertex3_t *a,
		    const vertex3_t *b,
		    const vertex3_t *c)
{
  return
    (a->x * b->y * c->z) +
    (b->x * c->y * a->z) +
    (c->x * a->y * b->z) -
    (c->x * b->y * a->z) -
    (b->x * a->y * c->z) -
    (a->x * c->y * b->z);
}
