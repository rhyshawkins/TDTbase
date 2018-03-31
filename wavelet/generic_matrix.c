
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_blas.h>

#include "generic_matrix.h"

#include "slog.h"

int
generic_matrix_interlace_index(int i, int width)
{
  if (i < 0 || i >= width) {
    return -1;
  }

  if (i < width/2) {

    return 2*i;

  } else {

    i -= width/2;

    return 2*i + 1;
    
  }
}

int
generic_matrix_deinterlace_index(int i, int width)
{
  if (i < 0 || i >= width) {
    return -1;
  }

  if (i % 2 == 0) {
    return i/2;
  } else {
    return width/2 + (i - 1)/2;
  }
}

int
generic_matrix_fill_forward2d_col_A(int j_max,
				    int si,
				    int nscaling,
				    const double *scaling,
				    const int *soffset,
				    gsl_matrix *m,
				    int row_offset,
				    int col_offset,
				    generic_matrix_edge_func_t edge)
{
  int width;

  int i;
  int orow;
  int j;
  int ci;
  int ocol;

  width = 1 << j_max;

  if ((row_offset + width) > m->size1 ||
      (col_offset + width*width) > m->size2) {
    ERROR("matrix too small (%d %d) (%d %d)",
	  row_offset + width, (int)m->size1,
	  col_offset + width*width, (int)m->size2);
    return -1;
  }

  for (i = 0; i < width; i ++) {

    orow = i;

    for (j = 0; j < nscaling; j ++) {

      ci = edge(2*si + soffset[j], width);

      ocol = i + width*ci;

      gsl_matrix_set(m, row_offset + orow, col_offset + ocol,
		     gsl_matrix_get(m, row_offset + orow, col_offset + ocol) +
		     scaling[j]);
    }
  }
		     
  return 0;
}

int
generic_matrix_fill_forward2d_col_B(int j_max,
				    int wi,
				    int nwavelet,
				    const double *wavelet,
				    const int *woffset,
				    gsl_matrix *m,
				    int row_offset,
				    int col_offset,
				    generic_matrix_edge_func_t edge)
{
  int width;

  int i;
  int orow;
  int j;
  int ci;
  int ocol;

  width = 1 << j_max;

  if ((row_offset + width) > m->size1 ||
      (col_offset + width*width) > m->size2) {
    ERROR("matrix too small (%d %d) (%d %d)",
	  row_offset + width, (int)m->size1,
	  col_offset + width*width, (int)m->size2);
    return -1;
  }

  for (i = 0; i < width; i ++) {

    orow = i;

    for (j = 0; j < nwavelet; j ++) {

      ci = edge(2*wi + 1 + woffset[j], width);

      ocol = i + width*ci;

      gsl_matrix_set(m, row_offset + orow, col_offset + ocol,
		     gsl_matrix_get(m, row_offset + orow, col_offset + ocol) +
		     wavelet[j]);
    }
  }
		     
  return 0;
}

int
generic_matrix_fill_forward2d_col_step(int j_max,
				       int nscaling,
				       const double *scaling,
				       const int *soffset,
				       int nwavelet,
				       const double *wavelet,
				       const int *woffset,
				       gsl_matrix *m,
				       generic_matrix_edge_func_t edge)
{
  int rowstride;
  int size;
  
  int width;
  int hwidth;
  
  int i;

  
  if (j_max < 1) {
    ERROR("invalid j parameters (%d)", j_max);
    return -1;
  }
  
  rowstride = 1 << j_max;
  size = rowstride * rowstride;

  if (m->size1 < size ||
      m->size2 < size) {
    ERROR("matrix too small");
    return -1;
  }

  width = rowstride;
  hwidth = width/2;

  for (i = 0; i < hwidth; i ++ ) {

    if (generic_matrix_fill_forward2d_col_A(j_max,
					    i,
					    nscaling,
					    scaling,
					    soffset,
					    m,
					    width * i,
					    0,
					    edge) < 0) {
      return -1;
    }

    if (generic_matrix_fill_forward2d_col_B(j_max,
					    i,
					    nwavelet,
					    wavelet,
					    woffset,
					    m,
					    hwidth*width + width*i,
					    0,
					    edge) < 0) {
      return -1;
    }
  }

  return 0;
}

int
generic_matrix_create_forward2d_col_step(int j_max,
					 int nscaling,
					 const double *scaling,
					 const int *soffset,
					 int nwavelet,
					 const double *wavelet,
					 const int *woffset,
					 gsl_matrix **_m,
					 generic_matrix_edge_func_t edge)
{
  int rowstride;
  int size;
  gsl_matrix *m;
  
  rowstride = 1 << j_max;
  size = rowstride * rowstride;

  m = gsl_matrix_alloc(size, size);
  if (m == NULL) {
    ERROR("failed to allocate matrix");
    return -1;
  }

  gsl_matrix_set_all(m, 0.0);

  if (generic_matrix_fill_forward2d_col_step(j_max,
					     nscaling,
					     scaling,
					     soffset,
					     nwavelet,
					     wavelet,
					     woffset,
					     m,
					     edge) < 0) {
    ERROR("failed to fill matrix");
    return -1;
  }

  *_m = m;
						      
  return 0;
}

int
generic_matrix_create_forward2d_col_substep(int j_max,
					    int j0,
					    int nscaling,
					    const double *scaling,
					    const int *soffset,
					    int nwavelet,
					    const double *wavelet,
					    const int *woffset,
					    gsl_matrix **_m,
					    generic_matrix_edge_func_t edge)
{
  int rowstride;
  int size;
  int rowstride0;
  int size0;
  int i;
  
  gsl_matrix *m;
  
  rowstride = 1 << j_max;
  size = rowstride * rowstride;
  rowstride0 = 1 << j0;
  size0 = rowstride0 * rowstride0;

  m = gsl_matrix_alloc(size, size);
  if (m == NULL) {
    ERROR("failed to allocate matrix");
    return -1;
  }

  gsl_matrix_set_all(m, 0.0);

  if (generic_matrix_fill_forward2d_col_step(j0,
					     nscaling,
					     scaling,
					     soffset,
					     nwavelet,
					     wavelet,
					     woffset,
					     m,
					     edge) < 0) {
    ERROR("generic_matrix_create_forward2d_col_substep: failed to fill matrix\n");
    return -1;
  }

  for (i = size0; i < size; i ++) {
    gsl_matrix_set(m, i, i, 1.0);
  }

  *_m = m;
						      
  return 0;
}

int
generic_matrix_fill_forward2d_row_A(int j_max,
				    int nscaling,
				    const double *scaling,
				    const int *soffset,
				    gsl_matrix *m,
				    int row_offset,
				    int col_offset,
				    generic_matrix_edge_func_t edge)
{
  int width;
  int height;
  int row;
  int col;
  int i;

  width = 1 << j_max;
  height = width/2;

  if (row_offset + height > m->size1 ||
      col_offset + width > m->size2) {
    return -1;
  }

  for (row = 0; row < height; row ++) {

    for (i = 0; i < nscaling; i ++) {
      col = edge(2 * row + soffset[i], width);

      gsl_matrix_set(m, row_offset + row, col_offset + col,
		     gsl_matrix_get(m, row_offset + row, col_offset + col) + scaling[i]);

    }
  }

  return 0;
}

int
generic_matrix_fill_forward2d_row_B(int j_max,
				    int nwavelet,
				    const double *wavelet,
				    const int *woffset,
				    gsl_matrix *m,
				    int row_offset,
				    int col_offset,
				    generic_matrix_edge_func_t edge)
{
  int width;
  int height;
  int row;
  int col;
  int i;
  
  width = 1 << j_max;
  height = width/2;

  if (row_offset + height > m->size1 ||
      col_offset + width > m->size2) {
    return -1;
  }

  for (row = 0; row < height; row ++) {

    for (i = 0; i < nwavelet; i ++) {
      col = edge(2 * row + 1 + woffset[i], width);

      gsl_matrix_set(m, row_offset + row, col_offset + col,
		     gsl_matrix_get(m, row_offset + row, col_offset + col) + wavelet[i]);

    }
  }

  return 0;
}

int
generic_matrix_fill_forward2d_row_step(int j_max,
				       int nscaling,
				       const double *scaling,
				       const int *soffset,
				       int nwavelet,
				       const double *wavelet,
				       const int *woffset,
				       gsl_matrix *m,
				       generic_matrix_edge_func_t edge)
{
  int rowstride;
  int size;
  int hsize;
  
  int width;
  int hwidth;
  
  int i;
  int c;
  
  if (j_max < 1) {
    ERROR("invalid j parameters (%d)", j_max);
    return -1;
  }
  
  rowstride = 1 << j_max;
  size = rowstride * rowstride;

  if (m->size1 < size ||
      m->size2 < size) {
    ERROR("matrix too small");
    return -1;
  }

  width = rowstride;
  hwidth = width/2;
  hsize = size/2;

  for (c = 0; c < 2; c ++) {
    
    for (i = 0; i < hwidth; i ++) {
      
      if (generic_matrix_fill_forward2d_row_A(j_max,
					      nscaling,
					      scaling,
					      soffset,
					      m,
					      c*hsize + hwidth*i,
					      c*hsize + width*i,
					      edge) < 0) {
	return -1;
      }
      
      if (generic_matrix_fill_forward2d_row_B(j_max,
					      nwavelet,
					      wavelet,
					      woffset,
					      m,
					      c*hsize + hwidth*i + hwidth*hwidth,
					      c*hsize + width*i,
					      edge) < 0) {
	return -1;
      }
    }
  }
  
  return 0;
}

int
generic_matrix_create_forward2d_row_step(int j_max,
					 int nscaling,
					 const double *scaling,
					 const int *soffset,
					 int nwavelet,
					 const double *wavelet,
					 const int *woffset,
					 gsl_matrix **_m,
					 generic_matrix_edge_func_t edge)
{
  int rowstride;
  int size;
  gsl_matrix *m;
  
  rowstride = 1 << j_max;
  size = rowstride * rowstride;

  m = gsl_matrix_alloc(size, size);
  if (m == NULL) {
    ERROR("failed to allocate matrix");
    return -1;
  }

  gsl_matrix_set_all(m, 0.0);

  if (generic_matrix_fill_forward2d_row_step(j_max,
					     nscaling,
					     scaling,
					     soffset,
					     nwavelet,
					     wavelet,
					     woffset,
					     m,
					     edge) < 0) {
    ERROR("failed to fill matrix");
    return -1;
  }

  *_m = m;
						      
  return 0;
}  

int
generic_matrix_create_forward2d_row_substep(int j_max,
					    int j0,
					    int nscaling,
					    const double *scaling,
					    const int *soffset,
					    int nwavelet,
					    const double *wavelet,
					    const int *woffset,
					    gsl_matrix **_m,
					    generic_matrix_edge_func_t edge)
{
  int rowstride;
  int size;
  int rowstride0;
  int size0;
  int i;
  
  gsl_matrix *m;
  
  rowstride = 1 << j_max;
  size = rowstride * rowstride;
  rowstride0 = 1 << j0;
  size0 = rowstride0 * rowstride0;

  m = gsl_matrix_alloc(size, size);
  if (m == NULL) {
    ERROR("failed to allocate matrix");
    return -1;
  }

  gsl_matrix_set_all(m, 0.0);

  if (generic_matrix_fill_forward2d_row_step(j0,
					     nscaling,
					     scaling,
					     soffset,
					     nwavelet,
					     wavelet,
					     woffset,
					     m,
					     edge) < 0) {
    ERROR("failed to fill matrix");
    return -1;
  }

  for (i = size0; i < size; i ++) {
    gsl_matrix_set(m, i, i, 1.0);
  }

  *_m = m;
						      
  return 0;
}

int
generic_matrix_fill_forward2d_step(int jmax,
				   int nscaling,
				   const double *scaling,
				   const int *soffset,
				   int nwavelet,
				   const double *wavelet,
				   const int *woffset,
				   gsl_matrix *m,
				   generic_matrix_edge_func_t edge)
{
  gsl_matrix *r;
  gsl_matrix *c;
  gsl_matrix *combined;

  int i;
  int j;
  
  if (generic_matrix_create_forward2d_col_step(jmax,
					       nscaling,
					       scaling,
					       soffset,
					       nwavelet,
					       wavelet,
					       woffset,
					       &c,
					       edge) < 0) {
    return -1;
  }
  
  if (generic_matrix_create_forward2d_row_step(jmax,
					       nscaling,
					       scaling,
					       soffset,
					       nwavelet,
					       wavelet,
					       woffset,
					       &r,
					       edge) < 0) {
    return -1;
  }

  if (c->size1 != r->size1 ||
      c->size2 != r->size2) {
    return -1;
  }

  if (m->size1 < c->size1 ||
      m->size2 < c->size2) {
    return -1;
  }

  if (m->size1 == c->size1 &&
      m->size2 == c->size2) {

    /*
     * Multiply directly into output matrix in the case where the size matches
     */
    if (gsl_blas_dgemm(CblasNoTrans,
		       CblasNoTrans,
		       1.0,
		       r,
		       c,
		       0.0,
		       m) < 0) {
      return -1;
    }
    
  } else {
    
    combined = gsl_matrix_alloc(c->size1, c->size2);
    if (combined == NULL) {
      return -1;
    }
    
    if (gsl_blas_dgemm(CblasNoTrans,
		       CblasNoTrans,
		       1.0,
		       r,
		       c,
		       0.0,
		       combined) < 0) {
      return -1;
    }
    
    for (i = 0; i < combined->size1; i ++) {
      for (j = 0; j < combined->size2; j ++) {
	
	gsl_matrix_set(m, i, j, gsl_matrix_get(combined, i, j));
	
      }
    }
    
    gsl_matrix_free(combined);
  }
  
  gsl_matrix_free(r);
  gsl_matrix_free(c);

  return 0;
}

int
generic_matrix_create_forward2d_step(int j_max,
				     int nscaling,
				     const double *scaling,
				     const int *soffset,
				     int nwavelet,
				     const double *wavelet,
				     const int *woffset,
				     gsl_matrix **_m,
				     generic_matrix_edge_func_t edge)
{
  int rowstride;
  int size;
  
  gsl_matrix *m;
  
  rowstride = 1 << j_max;
  size = rowstride * rowstride;

  m = gsl_matrix_alloc(size, size);
  if (m == NULL) {
    ERROR("failed to allocate matrix");
    return -1;
  }

  gsl_matrix_set_all(m, 0.0);

  if (generic_matrix_fill_forward2d_step(j_max,
					 nscaling,
					 scaling,
					 soffset,
					 nwavelet,
					 wavelet,
					 woffset,
					 m,
					 edge) < 0) {
    ERROR("failed to fill matrix");
    return -1;
  }

  *_m = m;
						      
  return 0;
}

int
generic_matrix_create_forward2d_substep(int j_max,
					int j0,
					int nscaling,
					const double *scaling,
					const int *soffset,
					int nwavelet,
					const double *wavelet,
					const int *woffset,
					gsl_matrix **_m,
					generic_matrix_edge_func_t edge)
{
  int rowstride;
  int size;
  int rowstride0;
  int size0;
  int i;
  
  gsl_matrix *m;
  
  rowstride = 1 << j_max;
  size = rowstride * rowstride;
  rowstride0 = 1 << j0;
  size0 = rowstride0 * rowstride0;

  m = gsl_matrix_alloc(size, size);
  if (m == NULL) {
    ERROR("failed to allocate matrix");
    return -1;
  }

  gsl_matrix_set_all(m, 0.0);

  if (generic_matrix_fill_forward2d_step(j0,
					 nscaling,
					 scaling,
					 soffset,
					 nwavelet,
					 wavelet,
					 woffset,
					 m,
					 edge) < 0) {
    ERROR("failed to fill matrix");
    return -1;
  }

  for (i = size0; i < size; i ++) {
    gsl_matrix_set(m, i, i, 1.0);
  }

  *_m = m;
						      
  return 0;
}

int
generic_matrix_create_forward2d(int j_max,
				int nscaling,
				const double *scaling,
				const int *soffset,
				int nwavelet,
				const double *wavelet,
				const int *woffset,
				gsl_matrix **_m,
				generic_matrix_edge_func_t edge)
{
  gsl_matrix *m;
  gsl_matrix *a;
  gsl_matrix *b;

  int j;

  if (generic_matrix_create_forward2d_step(j_max,
					   nscaling,
					   scaling,
					   soffset,
					   nwavelet,
					   wavelet,
					   woffset,
					   &m,
					   edge) < 0) {
    return -1;
  }

  for (j = j_max - 1; j > 0; j --) {
    if (generic_matrix_create_forward2d_substep(j_max,
						j,
						nscaling,
						scaling,
						soffset,
						nwavelet,
						wavelet,
						woffset,
						&a,
						edge) < 0) {
      return -1;
    }

    b = m;

    m = gsl_matrix_alloc(b->size1, b->size2);
    
    if (gsl_blas_dgemm(CblasNoTrans,
		       CblasNoTrans,
		       1.0,
		       a,
		       b,
		       0.0,
		       m) < 0) {
    return -1;
    }
    
    gsl_matrix_free(a);
    gsl_matrix_free(b);
  }

  *_m = m;
  return 0;
}


int
generic_matrix_fill_inverse2d_col_A(int j_max,
				    int si,
				    int nscaling,
				    const double *scaling,
				    const int *soffset,
				    gsl_matrix *m,
				    int row_offset,
				    int col_offset,
				    generic_matrix_edge_func_t edge)
{
  int width;

  int i;
  int orow;
  int j;
  int ci;
  int ocol;

  width = 1 << j_max;

  if ((row_offset + width) > m->size1 ||
      (col_offset + width*width) > m->size2) {
    ERROR("matrix too small (%d %d) (%d %d)",
	  row_offset + width, (int)m->size1,
	  col_offset + width*width, (int)m->size2);
    return -1;
  }

  for (i = 0; i < width; i ++) {

    orow = i;

    for (j = 0; j < nscaling; j ++) {

      ci = edge(2*si + soffset[j], width);

      if (ci % 2 == 0) {
	ocol = i + width*ci/2;
      } else {
	ocol = width*width/2 + i + width*(ci - 1)/2;
      }

      gsl_matrix_set(m, row_offset + orow, col_offset + ocol,
		     gsl_matrix_get(m, row_offset + orow, col_offset + ocol) +
		     scaling[j]);
    }
  }
		     
  return 0;
}

int
generic_matrix_fill_inverse2d_col_B(int j_max,
				    int wi,
				    int nwavelet,
				    const double *wavelet,
				    const int *woffset,
				    gsl_matrix *m,
				    int row_offset,
				    int col_offset,
				    generic_matrix_edge_func_t edge)
{
  int width;

  int i;
  int orow;
  int j;
  int ci;
  int ocol;

  width = 1 << j_max;

  if ((row_offset + width) > m->size1 ||
      (col_offset + width*width) > m->size2) {
    ERROR("matrix too small (%d %d) (%d %d)",
	  row_offset + width, (int)m->size1,
	  col_offset + width*width, (int)m->size2);
    return -1;
  }

  for (i = 0; i < width; i ++) {

    orow = i;

    for (j = 0; j < nwavelet; j ++) {

      ci = edge(2*wi + 1 + woffset[j], width);

      if (ci % 2 == 0) {
	ocol = i + width*ci/2;
      } else {
	ocol = width*width/2 + i + width*(ci - 1)/2;
      }

      gsl_matrix_set(m, row_offset + orow, col_offset + ocol,
		     gsl_matrix_get(m, row_offset + orow, col_offset + ocol) +
		     wavelet[j]);
    }
  }
		     
  return 0;
}

int
generic_matrix_fill_inverse2d_col_step(int j_max,
				       int nscaling,
				       const double *scaling,
				       const int *soffset,
				       int nwavelet,
				       const double *wavelet,
				       const int *woffset,
				       gsl_matrix *m,
				       generic_matrix_edge_func_t edge)
{
  int rowstride;
  int size;
  
  int width;
  int hwidth;
  
  int i;
  
  if (j_max < 1) {
    ERROR("invalid j parameters (%d)", j_max);
    return -1;
  }
  
  rowstride = 1 << j_max;
  size = rowstride * rowstride;

  if (m->size1 < size ||
      m->size2 < size) {
    ERROR("matrix too small");
    return -1;
  }

  width = rowstride;
  hwidth = width/2;

  for (i = 0; i < hwidth; i ++) {

    if (generic_matrix_fill_inverse2d_col_A(j_max,
					    i,
					    nscaling,
					    scaling,
					    soffset,
					    m,
					    2*width*i,
					    0,
					    edge) < 0) {
      return -1;
    }

    if (generic_matrix_fill_inverse2d_col_B(j_max,
					    i,
					    nwavelet,
					    wavelet,
					    woffset,
					    m,
					    2*width*i + width,
					    0,
					    edge) < 0) {
      return -1;
    }
  }

  return 0;
}

int
generic_matrix_create_inverse2d_col_step(int j_max,
					 int nscaling,
					 const double *scaling,
					 const int *soffset,
					 int nwavelet,
					 const double *wavelet,
					 const int *woffset,
					 gsl_matrix **_m,
					 generic_matrix_edge_func_t edge)
{
  int rowstride;
  int size;
  gsl_matrix *m;
  
  rowstride = 1 << j_max;
  size = rowstride * rowstride;

  m = gsl_matrix_alloc(size, size);
  if (m == NULL) {
    ERROR("failed to allocate matrix");
    return -1;
  }

  gsl_matrix_set_all(m, 0.0);

  if (generic_matrix_fill_inverse2d_col_step(j_max,
					     nscaling,
					     scaling,
					     soffset,
					     nwavelet,
					     wavelet,
					     woffset,
					     m,
					     edge) < 0) {
    ERROR("failed to fill matrix");
    return -1;
  }

  *_m = m;
						      
  return 0;
}

int
generic_matrix_create_inverse2d_col_substep(int j_max,
					    int j0,
					    int nscaling,
					    const double *scaling,
					    const int *soffset,
					    int nwavelet,
					    const double *wavelet,
					    const int *woffset,
					    gsl_matrix **_m,
					    generic_matrix_edge_func_t edge)
{
  int rowstride;
  int size;
  int rowstride0;
  int size0;
  int i;
  
  gsl_matrix *m;
  
  rowstride = 1 << j_max;
  size = rowstride * rowstride;
  rowstride0 = 1 << j0;
  size0 = rowstride0 * rowstride0;

  m = gsl_matrix_alloc(size, size);
  if (m == NULL) {
    ERROR("failed to allocate matrix");
    return -1;
  }

  gsl_matrix_set_all(m, 0.0);

  if (generic_matrix_fill_inverse2d_col_step(j0,
					     nscaling,
					     scaling,
					     soffset,
					     nwavelet,
					     wavelet,
					     woffset,
					     m,
					     edge) < 0) {
    ERROR("failed to fill matrix");
    return -1;
  }

  for (i = size0; i < size; i ++) {
    gsl_matrix_set(m, i, i, 1.0);
  }

  *_m = m;
						      
  return 0;
}

int
generic_matrix_fill_inverse2d_row_A(int j_max,
				    int nscaling,
				    const double *scaling,
				    const int *soffset,
				    gsl_matrix *m,
				    int row_offset,
				    int col_offset,
				    generic_matrix_edge_func_t edge)
{
  int orow;
  int ocol;
  int width;
  int hwidth;
  int size;

  int i;
  int j;
  int k;

  width = 1 << j_max;
  hwidth = width/2;

  size = width * hwidth;
  if (m->size1 < row_offset + size ||
      m->size2 < col_offset + size) {
    ERROR("matrix too small (%d %d) (%d %d)",
	  (int)m->size1, row_offset + size,
	  (int)m->size2, col_offset + size);
    return -1;
  }
    
  for (k = 0; k < hwidth; k ++) {

    for (i = 0; i < hwidth; i ++) {
      
      orow = k*width + 2*i;
      
      for (j = 0; j < nscaling; j ++) {
	
	ocol = edge(2*i + soffset[j], width);
	
	if (ocol % 2 == 0) {
	  ocol /= 2;
	} else {
	  ocol = hwidth*hwidth + (ocol - 1)/2;
	}
	
	gsl_matrix_set(m,
		       row_offset + orow,
		       col_offset + ocol + k*hwidth,
		       gsl_matrix_get(m, row_offset + orow, col_offset + ocol + k*hwidth) + scaling[j]);
      }
    }
  }

  return 0;
}

int
generic_matrix_fill_inverse2d_row_B(int j_max,
				    int nwavelet,
				    const double *wavelet,
				    const int *woffset,
				    gsl_matrix *m,
				    int row_offset,
				    int col_offset,
				    generic_matrix_edge_func_t edge)
{
  int orow;
  int ocol;
  int width;
  int hwidth;
  int size;

  int i;
  int j;
  int k;

  width = 1 << j_max;
  hwidth = width/2;

  size = width * hwidth;
  if (m->size1 < row_offset + size ||
      m->size2 < col_offset + size) {
    ERROR("matrix too small (%d %d) (%d %d)",
	  (int)m->size1, row_offset + size,
	  (int)m->size2, col_offset + size);
    return -1;
  }
    
    
  for (k = 0; k < hwidth; k ++) {

    for (i = 0; i < hwidth; i ++) {

      orow = k*width + 2*i + 1;
      
      for (j = 0; j < nwavelet; j ++) {
	
	ocol = edge(2*i + 1 + woffset[j], width);
	
	if (ocol % 2 == 0) {
	  ocol /= 2;
	} else {
	  ocol = hwidth*hwidth + (ocol - 1)/2;
	}
	
	gsl_matrix_set(m,
		       row_offset + orow,
		       col_offset + ocol + k*hwidth,
		       gsl_matrix_get(m, row_offset + orow, col_offset + ocol + k*hwidth) + wavelet[j]);
      }
    }
  }

  return 0;
}


int
generic_matrix_fill_inverse2d_row_step(int j_max,
				       int nscaling,
				       const double *scaling,
				       const int *soffset,
				       int nwavelet,
				       const double *wavelet,
				       const int *woffset,
				       gsl_matrix *m,
				       generic_matrix_edge_func_t edge)
{
  int rowstride;
  int size;
  
  int i;
  
  if (j_max < 1) {
    ERROR("invalid j parameters (%d)", j_max);
    return -1;
  }
  
  rowstride = 1 << j_max;
  size = rowstride * rowstride;

  if (m->size1 < size ||
      m->size2 < size) {
    ERROR("matrix too small");
    return -1;
  }

  for (i = 0; i < 2; i ++) {
    if (generic_matrix_fill_inverse2d_row_A(j_max,
					    nscaling,
					    scaling,
					    soffset,
					    m,
					    size/2 * i,
					    size/2 * i,
					    edge) < 0) {
      return -1;
    }
    
    if (generic_matrix_fill_inverse2d_row_B(j_max,
    					    nwavelet,
    					    wavelet,
    					    woffset,
    					    m,
    					    size/2 * i,
    					    size/2 * i,
					    edge) < 0) {
      return -1;
    }
  }

  return 0;
}

int
generic_matrix_create_inverse2d_row_step(int j_max,
					 int nscaling,
					 const double *scaling,
					 const int *soffset,
					 int nwavelet,
					 const double *wavelet,
					 const int *woffset,
					 gsl_matrix **_m,
					 generic_matrix_edge_func_t edge)
{
  int rowstride;
  int size;
  gsl_matrix *m;
  
  rowstride = 1 << j_max;
  size = rowstride * rowstride;

  m = gsl_matrix_alloc(size, size);
  if (m == NULL) {
    ERROR("failed to allocate matrix");
    return -1;
  }

  gsl_matrix_set_all(m, 0.0);

  if (generic_matrix_fill_inverse2d_row_step(j_max,
					     nscaling,
					     scaling,
					     soffset,
					     nwavelet,
					     wavelet,
					     woffset,
					     m,
					     edge) < 0) {
    ERROR("failed to fill matrix");
    return -1;
  }

  *_m = m;
						      
  return 0;
}

int
generic_matrix_create_inverse2d_row_substep(int j_max,
					    int j0,
					    int nscaling,
					    const double *scaling,
					    const int *soffset,
					    int nwavelet,
					    const double *wavelet,
					    const int *woffset,
					    gsl_matrix **_m,
					    generic_matrix_edge_func_t edge)
{
  int rowstride;
  int size;
  int rowstride0;
  int size0;
  int i;
  
  gsl_matrix *m;
  
  rowstride = 1 << j_max;
  size = rowstride * rowstride;
  rowstride0 = 1 << j0;
  size0 = rowstride0 * rowstride0;

  m = gsl_matrix_alloc(size, size);
  if (m == NULL) {
    ERROR("failed to allocate matrix");
    return -1;
  }

  gsl_matrix_set_all(m, 0.0);

  if (generic_matrix_fill_inverse2d_row_step(j0,
					     nscaling,
					     scaling,
					     soffset,
					     nwavelet,
					     wavelet,
					     woffset,
					     m,
					     edge) < 0) {
    ERROR("failed to fill matrix");
    return -1;
  }

  for (i = size0; i < size; i ++) {
    gsl_matrix_set(m, i, i, 1.0);
  }

  *_m = m;
						      
  return 0;
}

int
generic_matrix_fill_inverse2d_step(int j_max,
				   int nscaling,
				   const double *scaling,
				   const int *soffset,
				   int nwavelet,
				   const double *wavelet,
				   const int *woffset,
				   gsl_matrix *m,
				   generic_matrix_edge_func_t edge)
{
  gsl_matrix *r;
  gsl_matrix *c;
  gsl_matrix *combined;

  int i;
  int j;
  
  if (generic_matrix_create_inverse2d_col_step(j_max,
					       nscaling,
					       scaling,
					       soffset,
					       nwavelet,
					       wavelet,
					       woffset,
					       &c,
					       edge) < 0) {
    return -1;
  }
  
  if (generic_matrix_create_inverse2d_row_step(j_max,
					       nscaling,
					       scaling,
					       soffset,
					       nwavelet,
					       wavelet,
					       woffset,
					       &r,
					       edge) < 0) {
    return -1;
  }

  if (c->size1 != r->size1 ||
      c->size2 != r->size2) {
    return -1;
  }

  if (m->size1 < c->size1 ||
      m->size2 < c->size2) {
    return -1;
  }

  if (m->size1 == c->size1 &&
      m->size2 == c->size2) {

    /*
     * Multiply directly into output matrix in the case where the size matches
     */
    if (gsl_blas_dgemm(CblasNoTrans,
		       CblasNoTrans,
		       1.0,
		       c,
		       r,
		       0.0,
		       m) < 0) {
      return -1;
    }
    
  } else {
    
    combined = gsl_matrix_alloc(c->size1, c->size2);
    if (combined == NULL) {
      return -1;
    }
    
    if (gsl_blas_dgemm(CblasNoTrans,
		       CblasNoTrans,
		       1.0,
		       c,
		       r,
		       0.0,
		       combined) < 0) {
      return -1;
    }
    
    for (i = 0; i < combined->size1; i ++) {
      for (j = 0; j < combined->size2; j ++) {
	
	gsl_matrix_set(m, i, j, gsl_matrix_get(combined, i, j));
	
      }
    }
    
    gsl_matrix_free(combined);
  }
  
  gsl_matrix_free(r);
  gsl_matrix_free(c);

  return 0;
}  

int
generic_matrix_create_inverse2d_step(int j_max,
				     int nscaling,
				     const double *scaling,
				     const int *soffset,
				     int nwavelet,
				     const double *wavelet,
				     const int *woffset,
				     gsl_matrix **_m,
				     generic_matrix_edge_func_t edge)
{
  int rowstride;
  int size;
  
  gsl_matrix *m;
  
  rowstride = 1 << j_max;
  size = rowstride * rowstride;

  m = gsl_matrix_alloc(size, size);
  if (m == NULL) {
    ERROR("failed to allocate matrix");
    return -1;
  }

  gsl_matrix_set_all(m, 0.0);

  if (generic_matrix_fill_inverse2d_step(j_max,
					 nscaling,
					 scaling,
					 soffset,
					 nwavelet,
					 wavelet,
					 woffset,
					 m,
					 edge) < 0) {
    ERROR("failed to fill matrix");
    return -1;
  }

  *_m = m;
						      
  return 0;
}  

int
generic_matrix_create_inverse2d_substep(int j_max,
					int j0,
					int nscaling,
					const double *scaling,
					const int *soffset,
					int nwavelet,
					const double *wavelet,
					const int *woffset,
					gsl_matrix **_m,
					generic_matrix_edge_func_t edge)
{
  int rowstride;
  int size;
  int rowstride0;
  int size0;
  int i;
  
  gsl_matrix *m;
  
  rowstride = 1 << j_max;
  size = rowstride * rowstride;
  rowstride0 = 1 << j0;
  size0 = rowstride0 * rowstride0;

  m = gsl_matrix_alloc(size, size);
  if (m == NULL) {
    ERROR("failed to allocate matrix");
    return -1;
  }

  gsl_matrix_set_all(m, 0.0);

  if (generic_matrix_fill_inverse2d_step(j0,
					 nscaling,
					 scaling,
					 soffset,
					 nwavelet,
					 wavelet,
					 woffset,
					 m,
					 edge) < 0) {
    ERROR("failed to fill matrix (%d %d)", j_max, j0);
    return -1;
  }

  for (i = size0; i < size; i ++) {
    gsl_matrix_set(m, i, i, 1.0);
  }

  *_m = m;
						      
  return 0;
}  


int
generic_matrix_create_inverse2d(int j_max,
				int nscaling,
				const double *scaling,
				const int *soffset,
				int nwavelet,
				const double *wavelet,
				const int *woffset,
				gsl_matrix **_m,
				generic_matrix_edge_func_t edge)
{
  gsl_matrix *m;
  gsl_matrix *a;
  gsl_matrix *b;

  int j;

  if (generic_matrix_create_inverse2d_substep(j_max,
					      1,
					      nscaling,
					      scaling,
					      soffset,
					      nwavelet,
					      wavelet,
					      woffset,
					      &m,
					      edge) < 0) {
    return -1;
  }

  for (j = 2; j <= j_max; j ++) {
    if (generic_matrix_create_inverse2d_substep(j_max,
						j,
						nscaling,
						scaling,
						soffset,
						nwavelet,
						wavelet,
						woffset,
						&a,
						edge) < 0) {
      return -1;
    }

    b = m;

    m = gsl_matrix_alloc(b->size1, b->size2);
    
    if (gsl_blas_dgemm(CblasNoTrans,
		       CblasNoTrans,
		       1.0,
		       a,
		       b,
		       0.0,
		       m) < 0) {
      return -1;
    }
    
    gsl_matrix_free(a);
    gsl_matrix_free(b);
  }

  *_m = m;
  return 0;
}
