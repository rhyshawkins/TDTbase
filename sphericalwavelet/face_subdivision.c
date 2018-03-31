
#include <stdio.h>

#include "face_subdivision.h"

/*
 * Full Transforms 
 */
int
face_subdivision_forward(manifold_t *m,
			 double *coeff)
{
  int depth;

  for (depth = m->degree; depth > 0; depth --) {
    if (face_subdivision_forward_step(m, coeff, depth) < 0) {
      return -1;
    }
  }

  return 0;
}

int
face_subdivision_inverse(manifold_t *m,
			 double *coeff)
{
  int depth;

  for (depth = 1; depth <= m->degree; depth ++) {
    if (face_subdivision_inverse_step(m, coeff, depth) < 0) {
      return -1;
    }
  }

  return 0;
}

/*
 * Single Steps
 */
int
face_subdivision_forward_step(manifold_t *m,
			      double *coeff,
			      int depth)
{
  int pt;
  int ct;
  double mean;

  int poffset;
  int coffset;

  int i;
  
  if (depth <= 0 || depth > m->degree) {
    return -1;
  }

  coffset = 0;
  for (i = 0; i < depth; i ++) {
    poffset = coffset;
    coffset += m->ntrianglesatdepth(i);
  }

  for (pt = 0; pt < m->ntriangles[depth - 1]; pt ++) {

    /* Compute the mean of the 4 child triangles */
    mean = 0.0;
    for (ct = 0; ct < 4; ct ++) {
      mean += coeff[coffset + m->triangles[depth - 1][pt].child_triangles[ct]];
    }
    mean /= 4.0;


    /* Set the mean of the children to the parent */
    coeff[poffset + pt] = mean;

    /* Remove the mean from the 4 children */
    for (ct = 0; ct < 4; ct ++) {
      coeff[coffset + m->triangles[depth - 1][pt].child_triangles[ct]] -= mean;
    }
  }

  return 0;
}

int
face_subdivision_inverse_step(manifold_t *m,
			      double *coeff,
			      int depth)
{
  int pt;
  int ct;

  int poffset;
  int coffset;
  
  int i;

  double mean;
  
  if (depth <= 0 || depth > m->degree) {
    return -1;
  }

  coffset = 0;
  for (i = 0; i < depth; i ++) {
    poffset = coffset;
    coffset += m->ntrianglesatdepth(i);
  }
  
  for (pt = 0; pt < m->ntriangles[depth - 1]; pt ++) {

    mean = coeff[poffset + pt];

    /* And the parent value to the 4 children */
    for (ct = 0; ct < 4; ct ++) {
      coeff[coffset + m->triangles[depth - 1][pt].child_triangles[ct]] += mean;
    }

  }

  return 0;
}
