
#include <stdio.h>

#include "cdf97_lift.h"

#include "vertex_wavelet.h"

#include "slog.h"

int
vertex_wavelet_butterfly_forward(manifold_t *m,
				 double *coeff)
{
  int depth;
  
  for (depth = m->degree; depth > 0; depth --) {
    if (vertex_wavelet_butterfly_forward_step(m, coeff, depth) < 0) {
      return -1;
    }
  }

  return 0;
}

int
vertex_wavelet_butterfly_inverse(manifold_t *m,
				 double *coeff)
{
  int depth;
  
  for (depth = 1; depth <= m->degree; depth ++) {
    if (vertex_wavelet_butterfly_inverse_step(m, coeff, depth) < 0) {
      return -1;
    }
  }

  return 0;
}

int
vertex_wavelet_butterfly_forward_lifted(manifold_t *m,
					double *coeff)
{
  int depth;
  
  for (depth = m->degree; depth > 0; depth --) {
    if (vertex_wavelet_butterfly_forward_lifted_step(m, coeff, depth) < 0) {
      return -1;
    }
  }

  return 0;
}

int
vertex_wavelet_butterfly_inverse_lifted(manifold_t *m,
					double *coeff)
{
  int depth;
  
  for (depth = 1; depth <= m->degree; depth ++) {
    if (vertex_wavelet_butterfly_inverse_lifted_step(m, coeff, depth) < 0) {
      return -1;
    }
  }

  return 0;
}

int
vertex_wavelet_butterfly_forward_step(manifold_t *m,
				      double *coeff,
				      int depth)
{
  int i;
  
  int vstart;
  int vend;

  if (depth <= 0 || depth > m->degree) {
    return -1;
  }

  vstart = m->nverticesatdepth(depth - 1);
  vend = m->nverticesatdepth(depth);

  for (i = vstart; i < vend; i ++) {

    coeff[i] -= (0.5 * (coeff[m->vertices[i].v[0]] +
			coeff[m->vertices[i].v[1]])
		 +
		 0.125 * (coeff[m->vertices[i].f[0]] +
			  coeff[m->vertices[i].f[1]])
		 -
		 0.0625 * (coeff[m->vertices[i].e[0]] +
			   coeff[m->vertices[i].e[1]] +
			   coeff[m->vertices[i].e[2]] +
			   coeff[m->vertices[i].e[3]]));
    
  }
		      
  return 0;  
}

int
vertex_wavelet_butterfly_forward_lifted_step(manifold_t *m,
					     double *coeff,
					     int depth)
{
  int i;
  int vstart;
  int vend;
  double lift;
  
  if (vertex_wavelet_butterfly_forward_step(m, coeff, depth) < 0) {
    return -1;
  }

  vstart = m->nverticesatdepth(depth - 1);
  vend = m->nverticesatdepth(depth);

  for (i = vstart; i < vend; i ++) {

    lift = coeff[i] * m->vertices[i].area;

    coeff[m->vertices[i].v[0]] += lift;
    coeff[m->vertices[i].v[1]] += lift;

  }

  return 0;
}

int
vertex_wavelet_butterfly_inverse_step(manifold_t *m,
				      double *coeff,
				      int depth)
{
  int i;
  int vstart;
  int vend;

  if (depth <= 0 || depth > m->degree) {
    return -1;
  }

  vstart = m->nverticesatdepth(depth - 1);
  vend = m->nverticesatdepth(depth);

  for (i = vstart; i < vend; i ++) {
    coeff[i] += (0.5 * (coeff[m->vertices[i].v[0]] +
			coeff[m->vertices[i].v[1]])
		 +
		 0.125 * (coeff[m->vertices[i].f[0]] +
			  coeff[m->vertices[i].f[1]])
		 -
		 0.0625 * (coeff[m->vertices[i].e[0]] +
			   coeff[m->vertices[i].e[1]] +
			   coeff[m->vertices[i].e[2]] +
			   coeff[m->vertices[i].e[3]]));
  }
		      
  return 0;  
  
}

int
vertex_wavelet_butterfly_inverse_lifted_step(manifold_t *m,
					     double *coeff,
					     int depth)
{
  int i;
  int vstart;
  int vend;
  double lift;
  
  if (depth <= 0 || depth > m->degree) {
    return -1;
  }

  vstart = m->nverticesatdepth(depth - 1);
  vend = m->nverticesatdepth(depth);

  for (i = vstart; i < vend; i ++) {

    lift = coeff[i] * m->vertices[i].area;

    coeff[m->vertices[i].v[0]] -= lift;
    coeff[m->vertices[i].v[1]] -= lift;

  }

  if (vertex_wavelet_butterfly_inverse_step(m, coeff, depth) < 0) {
    return -1;
  }

  return 0;
}

/*
 * Shell transforms : radial degree == manifold degree
 */
int
vertex_wavelet_butterfly_shell_forward(manifold_t *m,
				       double *coeff,
				       int ncoeff,
				       double *workspace)
{
  int depth;
  int rowstride;
  int radial_size;
  int j;
  int vend;
  
  rowstride = m->nvertices;
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
      if (vertex_wavelet_butterfly_forward_step(m, coeff + rowstride*j, depth) < 0) {
	ERROR("failed to do lateral transform");
	return -1;
      }
    }

    /*
     * Radial transform
     */
    vend = m->nverticesatdepth(depth);

    for (j = 0; j < vend; j ++) {
      if (cdf97_lift_forward1d_cdf97_step(coeff + j, radial_size, rowstride, workspace) < 0) {
	ERROR("failed to do radial transform");
	return -1;
      }
    }
  }

  return 0;
}

int
vertex_wavelet_butterfly_shell_inverse(manifold_t *m,
				       double *coeff,
				       int ncoeff,
				       double *workspace)
{
  int depth;
  int radial_size;
  int vend;
  int j;
  int rowstride;
  
  rowstride = m->nvertices;
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
      if (cdf97_lift_inverse1d_cdf97_step(coeff + j, radial_size, rowstride, workspace) < 0) {
	ERROR("failed to do radial transform");
	return -1;
      }
    }
    
    /*
     * Lateral transform
     */
    for (j = 0; j < radial_size; j ++) {
      if (vertex_wavelet_butterfly_inverse_step(m, coeff + rowstride*j, depth) < 0) {
	ERROR("failed to do lateral transform");
	return -1;
      }
    }
  }

  return 0;
}

int
vertex_wavelet_butterfly_shell_forward_lifted(manifold_t *m,
					      double *coeff,
					      int ncoeff,
					      double *workspace)
{
  int depth;
  int rowstride;
  int radial_size;
  int j;
  int vend;
  
  rowstride = m->nvertices;
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
      if (vertex_wavelet_butterfly_forward_lifted_step(m, coeff + rowstride*j, depth) < 0) {
	ERROR("failed to do lateral transform");
	return -1;
      }
    }

    /*
     * Radial transform
     */
    vend = m->nverticesatdepth(depth);

    for (j = 0; j < vend; j ++) {
      if (cdf97_lift_forward1d_cdf97_step(coeff + j, radial_size, rowstride, workspace) < 0) {
	ERROR("failed to do radial transform");
	return -1;
      }
    }
  }

  return 0;
}

int
vertex_wavelet_butterfly_shell_inverse_lifted(manifold_t *m,
					      double *coeff,
					      int ncoeff,
					      double *workspace)
{
  int depth;
  int radial_size;
  int vend;
  int j;
  int rowstride;
  
  rowstride = m->nvertices;
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
      if (cdf97_lift_inverse1d_cdf97_step(coeff + j, radial_size, rowstride, workspace) < 0) {
	ERROR("failed to do radial transform");
	return -1;
      }
    }
    
    /*
     * Lateral transform
     */
    for (j = 0; j < radial_size; j ++) {
      if (vertex_wavelet_butterfly_inverse_lifted_step(m, coeff + rowstride*j, depth) < 0) {
	ERROR("failed to do lateral transform");
	return -1;
      }
    }
  }

  return 0;
}
