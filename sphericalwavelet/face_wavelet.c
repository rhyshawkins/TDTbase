
#include <stdio.h>

#include "face_wavelet.h"

#include "slog.h"

/*
 * There is a bit of messiness here since we have a coefficient for each
 * triangle in the hierarchy and for the wavelet transform we have a 
 * duplication of the scaling coefficient from the parent triangle to its
 * central sub-division triangle. This is to maintain compatibility with
 * the simple face subdivision code which uses the whole hierarchy of
 * triangles whereas biohaar could be done inplace.
 */

/*
 * Full transforms
 */
int
face_wavelet_biohaar_forward(manifold_t *m,
			     double *coeff)
{
  int depth;

  for (depth = m->degree; depth > 0; depth --) {
    if (face_wavelet_biohaar_forward_step(m, coeff, depth) < 0) {
      return -1;
    }
  }

  return 0;
}

int
face_wavelet_biohaar_inverse(manifold_t *m,
			     double *coeff)
{
  int depth;

  for (depth = 1; depth <= m->degree; depth ++) {
    if (face_wavelet_biohaar_inverse_step(m, coeff, depth) < 0) {
      return -1;
    }
  }

  return 0;
}


/*
 * Individual steps
 */
int
face_wavelet_biohaar_forward_step(manifold_t *m,
				  double *coeff,
				  int depth)
{
  int pt;
  int ct;
  double scaling;

  int poffset;
  int coffset;

  int i;
  int ci;
  
  if (depth <= 0 || depth > m->degree) {
    return -1;
  }

  coffset = 0;
  for (i = 0; i < depth; i ++) {
    poffset = coffset;
    coffset += m->ntrianglesatdepth(i);
  }

  for (pt = 0; pt < m->ntriangles[depth - 1]; pt ++) {

    /* Compute the scaling term from the values at the 4 child triangles */
    scaling = 0.0;
    for (ct = 0; ct < 4; ct ++) {
      ci = m->triangles[depth - 1][pt].child_triangles[ct];

      scaling +=
	coeff[coffset + ci] *
	m->triangles[depth][ci].area;
    }
    scaling /= m->triangles[depth - 1][pt].area;

    /* Copy the scaling to the parent */
    coeff[poffset + pt] = scaling;

    /* 
     * Compute the resulting wavelet coefficients for the 3 children.
     * In manifold code, the last child is the central triangle which
     * stores the scaling coefficient, but we copied that to the parent
     * above.
     */
    
    for (ct = 0; ct < 3; ct ++) {
      ci = m->triangles[depth - 1][pt].child_triangles[ct];
      coeff[coffset + ci] = 0.5 * (coeff[coffset + ci] - scaling);
    }

    /* 
     * This isn't strictly necessary, but we set the central triangle 
     * coeff to 0.
     */
    ci = m->triangles[depth - 1][pt].child_triangles[3];
    coeff[coffset + ci] = 0.0;
  }

  return 0;
  
}

int
face_wavelet_biohaar_inverse_step(manifold_t *m,
				  double *coeff,
				  int depth)
{
  int pt;
  int ct;

  int poffset;
  int coffset;
  
  int i;
  int ci;

  double scaling;
  double diff;
  
  if (depth <= 0 || depth > m->degree) {
    return -1;
  }

  coffset = 0;
  for (i = 0; i < depth; i ++) {
    poffset = coffset;
    coffset += m->ntrianglesatdepth(i);
  }
  
  for (pt = 0; pt < m->ntriangles[depth - 1]; pt ++) {

    scaling = coeff[poffset + pt];

    /* And the parent value to the 3 children */
    diff = 0.0;
    for (ct = 0; ct < 3; ct ++) {
      ci = m->triangles[depth - 1][pt].child_triangles[ct];
      coeff[coffset + ci] = 2.0 * coeff[coffset + ci] + scaling;
      diff += coeff[coffset + ci] * m->triangles[depth][ci].area;
    }

    /* 
     * Lastly we compute the central triangle from difference between scaling and
     * other coefficients.
     */

    ci = m->triangles[depth - 1][pt].child_triangles[3];
    coeff[coffset + ci] = (m->triangles[depth - 1][pt].area * scaling - diff)/m->triangles[depth][ci].area;
  }

  return 0;
  
}

int
face_wavelet_biohaar_shell_forward(manifold_t *m,
				   double *coeff,
				   int ncoeff,
				   double *workspace,
				   face_wavelet_radial_forward_step_t radial_forward_step)
{
  int depth;
  int rowstride;
  int radial_size;
  int j;
  int vend;
  
  rowstride = m->ntotaltriangles;
  radial_size = 1 << m->degree;

  if (ncoeff != (rowstride * radial_size)) {
    ERROR("size mismatch");
    return -1;
  }
  
  for (depth = m->degree; depth > 0; depth --) {

    radial_size = 1 << depth;

    /*
     * Lateral transform
     */
    for (j = 0; j < radial_size; j ++) {
      if (face_wavelet_biohaar_forward_step(m, coeff + rowstride*j, depth) < 0) {
	ERROR("failed to do lateral transform");
	return -1;
      }
    }

    /*
     * Radial transform
     */
    vend = m->nverticesatdepth(depth);

    for (j = 0; j < vend; j ++) {
      if (radial_forward_step(coeff + j, radial_size, rowstride, workspace) < 0) {
	ERROR("failed to do radial transform");
	return -1;
      }
    }
  }

  return 0;
}

int
face_wavelet_biohaar_shell_inverse(manifold_t *m,
				   double *coeff,
				   int ncoeff,
				   double *workspace,
				   face_wavelet_radial_inverse_step_t radial_inverse_step)
{
  int depth;
  int radial_size;
  int vend;
  int j;
  int rowstride;
  
  rowstride = m->ntotaltriangles;
  radial_size = 1 << m->degree;

  if (ncoeff != (rowstride * radial_size)) {
    ERROR("size mismatch");
    return -1;
  }

  for (depth = 1; depth <= m->degree; depth ++) {

    radial_size = 1 << depth;

    /*
     * Inverse Radial transform
     */
    vend = m->nverticesatdepth(depth);

    for (j = 0; j < vend; j ++) {
      if (radial_inverse_step(coeff + j, radial_size, rowstride, workspace) < 0) {
	ERROR("failed to do radial transform");
	return -1;
      }
    }
    
    /*
     * Lateral transform
     */
    for (j = 0; j < radial_size; j ++) {
      if (face_wavelet_biohaar_inverse_step(m, coeff + rowstride*j, depth) < 0) {
	ERROR("failed to do lateral transform");
	return -1;
      }
    }
  }

  return 0;
}
