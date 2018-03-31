
#include <stdio.h>
#include <stdlib.h>
#include <check.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include "generic_matrix.h"
#include "boundary.h"

START_TEST (test_wavelet_boundary_reflect)
{
  ck_assert(wavelet_boundary_reflect(0, 8) == 0);
  ck_assert(wavelet_boundary_reflect(-1, 8) == 1);
  ck_assert(wavelet_boundary_reflect(-2, 8) == 2);

  ck_assert(wavelet_boundary_reflect(7, 8) == 7);
  ck_assert(wavelet_boundary_reflect(8, 8) == 6);
  ck_assert(wavelet_boundary_reflect(9, 8) == 5);
  
  ck_assert(wavelet_boundary_reflect(8, 2) == 0);
  ck_assert(wavelet_boundary_reflect(9, 2) == 1);
}
END_TEST

START_TEST (test_wavelet_boundary_periodic)
{
  ck_assert(wavelet_boundary_periodic(0, 8) == 0);
  ck_assert(wavelet_boundary_periodic(-1, 8) == 7);
  ck_assert(wavelet_boundary_periodic(-2, 8) == 6);

  ck_assert(wavelet_boundary_periodic(7, 8) == 7);
  ck_assert(wavelet_boundary_periodic(8, 8) == 0);
  ck_assert(wavelet_boundary_periodic(9, 8) == 1);
  
  ck_assert(wavelet_boundary_periodic(8, 2) == 0);
  ck_assert(wavelet_boundary_periodic(9, 2) == 1);
}
END_TEST

START_TEST (test_generic_matrix_interlace_index)
{
  ck_assert(generic_matrix_interlace_index(0, 4) == 0);
  ck_assert(generic_matrix_interlace_index(1, 4) == 2);
  ck_assert(generic_matrix_interlace_index(2, 4) == 1);
  ck_assert(generic_matrix_interlace_index(3, 4) == 3);
}
END_TEST

START_TEST (test_generic_matrix_fill_forward2d_col_A)
{
  static const int JMAX = 2;
  
  static const double SCALING[3] = {0.1, 1.0, 0.2};
  static const int SOFFSET[3] = {-1, 0, 1};

  gsl_matrix *m;
  int rowstride;
  int size;
  int i;

  rowstride = 1 << JMAX;
  size = rowstride * rowstride;
  m = gsl_matrix_alloc(size, size);
  ck_assert(m != NULL);

  gsl_matrix_set_all(m, 0.0);
  
  ck_assert(generic_matrix_fill_forward2d_col_A(JMAX,
						0,
						3,
						SCALING,
						SOFFSET,
						m,
						0,
						0,
						wavelet_boundary_reflect) >= 0);

  for (i = 0; i < rowstride; i ++) {
    ck_assert(gsl_matrix_get(m, i, i) == SCALING[1]);
    ck_assert(gsl_matrix_get(m, i, rowstride + i) == (SCALING[0] + SCALING[2]));
  }

  gsl_matrix_set_all(m, 0.0);
  
  ck_assert(generic_matrix_fill_forward2d_col_A(JMAX,
						1,
						3,
						SCALING,
						SOFFSET,
						m,
						0,
						0,
						wavelet_boundary_reflect) >= 0);

  for (i = 0; i < rowstride; i ++) {
    ck_assert(gsl_matrix_get(m, i, rowstride + i) == SCALING[0]);
    ck_assert(gsl_matrix_get(m, i, rowstride*2 + i) == SCALING[1]);
    ck_assert(gsl_matrix_get(m, i, rowstride*3 + i) == SCALING[2]);
  }

  gsl_matrix_free(m);
}
END_TEST

START_TEST (test_generic_matrix_fill_forward2d_col_B)
{
  static const int JMAX = 2;
  
  static const double WAVELET[3] = {-0.1, -1.0, -0.1};
  static const int WOFFSET[3] = {-1, 0, 1};

  gsl_matrix *m;
  int rowstride;
  int size;
  int i;

  rowstride = 1 << JMAX;
  size = rowstride * rowstride;
  m = gsl_matrix_alloc(size, size);
  ck_assert(m != NULL);
  
  ck_assert(generic_matrix_fill_forward2d_col_B(JMAX,
						0,
						3,
						WAVELET,
						WOFFSET,
						m,
						0,
						0,
						wavelet_boundary_reflect) >= 0);

  for (i = 0; i < rowstride; i ++) {
    ck_assert(gsl_matrix_get(m, i, i) == WAVELET[0]);
    ck_assert(gsl_matrix_get(m, i, rowstride + i) == WAVELET[1]);
    ck_assert(gsl_matrix_get(m, i, 2*rowstride + i) == WAVELET[2]);
  }

  gsl_matrix_set_all(m, 0.0);
  
  ck_assert(generic_matrix_fill_forward2d_col_B(JMAX,
						rowstride/2 - 1,
						3,
						WAVELET,
						WOFFSET,
						m,
						0,
						0,
						wavelet_boundary_reflect) >= 0);

  for (i = 0; i < rowstride; i ++) {
    ck_assert(gsl_matrix_get(m, i, rowstride*rowstride - rowstride + i) == WAVELET[1]);
    ck_assert(gsl_matrix_get(m, i, rowstride*rowstride - 2*rowstride + i) == (WAVELET[0] + WAVELET[2]));
  }

  gsl_matrix_free(m);
}
END_TEST

START_TEST (test_generic_matrix_fill_inverse2d_col_A)
{
  static const int JMAX = 2;
  
  static const double SCALING[3] = {0.1, 1.0, 0.2};
  static const int SOFFSET[3] = {-1, 0, 1};

  gsl_matrix *m;
  int rowstride;
  int size;
  int i;

  rowstride = 1 << JMAX;
  size = rowstride * rowstride;
  m = gsl_matrix_alloc(size, size);
  ck_assert(m != NULL);

  gsl_matrix_set_all(m, 0.0);
  
  ck_assert(generic_matrix_fill_inverse2d_col_A(JMAX,
						0,
						3,
						SCALING,
						SOFFSET,
						m,
						0,
						0,
						wavelet_boundary_reflect) >= 0);

  for (i = 0; i < rowstride; i ++) {
    ck_assert(gsl_matrix_get(m, i, i) == SCALING[1]);
    ck_assert(gsl_matrix_get(m, i, rowstride*rowstride/2 + i) == (SCALING[0] + SCALING[2]));
  }

  gsl_matrix_set_all(m, 0.0);
  
  ck_assert(generic_matrix_fill_inverse2d_col_A(JMAX,
						1,
						3,
						SCALING,
						SOFFSET,
						m,
						0,
						0,
						wavelet_boundary_reflect) >= 0);

  for (i = 0; i < rowstride; i ++) {
    ck_assert(gsl_matrix_get(m, i, rowstride + i) == SCALING[1]);
    ck_assert(gsl_matrix_get(m, i, rowstride*rowstride/2 + i) == SCALING[0]);
    ck_assert(gsl_matrix_get(m, i, rowstride*rowstride/2 + rowstride + i) == SCALING[2]);
  }

  gsl_matrix_free(m);
}
END_TEST

START_TEST (test_generic_matrix_fill_inverse2d_col_B)
{
  static const int JMAX = 2;
  
  static const double WAVELET[3] = {-0.1, -1.0, -0.1};
  static const int WOFFSET[3] = {-1, 0, 1};

  gsl_matrix *m;
  int rowstride;
  int size;
  int i;
  
  rowstride = 1 << JMAX;
  size = rowstride * rowstride;
  m = gsl_matrix_alloc(size, size);
  ck_assert(m != NULL);
  
  ck_assert(generic_matrix_fill_inverse2d_col_B(JMAX,
						0,
						3,
						WAVELET,
						WOFFSET,
						m,
						0,
						0,
						wavelet_boundary_reflect) >= 0);

  for (i = 0; i < rowstride; i ++) {
    ck_assert(gsl_matrix_get(m, i, i) == WAVELET[0]);
    ck_assert(gsl_matrix_get(m, i, rowstride*rowstride/2 + i) == WAVELET[1]);
    ck_assert(gsl_matrix_get(m, i, rowstride + i) == WAVELET[2]);
  }

  gsl_matrix_set_all(m, 0.0);
  
  ck_assert(generic_matrix_fill_inverse2d_col_B(JMAX,
						rowstride/2 - 1,
						3,
						WAVELET,
						WOFFSET,
						m,
						0,
						0,
						wavelet_boundary_reflect) >= 0);
  
  for (i = 0; i < rowstride; i ++) {
    
    ck_assert(gsl_matrix_get(m, i, rowstride*rowstride - rowstride + i) == WAVELET[1]);
    ck_assert(gsl_matrix_get(m, i, rowstride*rowstride/2 - rowstride + i) == (WAVELET[0] + WAVELET[2]));
  }

  gsl_matrix_free(m);
}
END_TEST

START_TEST (test_generic_matrix_fill_forward2d_row_A)
{
  static const int JMAX = 2;
  
  static const double SCALING[3] = {0.1, 1.0, 0.2};
  static const int SOFFSET[3] = {-1, 0, 1};

  gsl_matrix *m;
  int rowstride;
  int size;
  int i;

  rowstride = 1 << JMAX;
  size = rowstride * rowstride;
  m = gsl_matrix_alloc(size, size);
  ck_assert(m != NULL);

  gsl_matrix_set_all(m, 0.0);
  
  ck_assert(generic_matrix_fill_forward2d_row_A(JMAX,
						3,
						SCALING,
						SOFFSET,
						m,
						0,
						0,
						wavelet_boundary_reflect) >= 0);

  for (i = 0; i < rowstride/2; i ++) {
    ck_assert(gsl_matrix_get(m, i, 2*i) == SCALING[1]);
    if (i > 0) {
      ck_assert(gsl_matrix_get(m, i, 2*i - 1) == SCALING[0]);
      ck_assert(gsl_matrix_get(m, i, 2*i + 1) == SCALING[2]);
    } else {
      ck_assert(gsl_matrix_get(m, i, 1) == (SCALING[0] + SCALING[2]));
    }
  }

  gsl_matrix_free(m);
}
END_TEST

START_TEST (test_generic_matrix_fill_forward2d_row_B)
{
  static const int JMAX = 2;
  
  static const double WAVELET[3] = {-0.1, -1.0, -0.1};
  static const int WOFFSET[3] = {-1, 0, 1};

  gsl_matrix *m;
  int rowstride;
  int size;
  int i;

  rowstride = 1 << JMAX;
  size = rowstride * rowstride;
  m = gsl_matrix_alloc(size, size);
  ck_assert(m != NULL);
  
  ck_assert(generic_matrix_fill_forward2d_row_B(JMAX,
						3,
						WAVELET,
						WOFFSET,
						m,
						0,
						0,
						wavelet_boundary_reflect) >= 0);

  for (i = 0; i < rowstride/2; i ++) {
    ck_assert(gsl_matrix_get(m, i, 2*i + 1) == WAVELET[1]);
    if (i < (rowstride/2 - 1)) {
      ck_assert(gsl_matrix_get(m, i, 2*i) == WAVELET[0]);
      ck_assert(gsl_matrix_get(m, i, 2*i + 2) == WAVELET[2]);
    } else {
      ck_assert(gsl_matrix_get(m, i, rowstride - 2) == (WAVELET[0] + WAVELET[2]));
    }
  }

  gsl_matrix_free(m);
}
END_TEST

START_TEST (test_generic_matrix_fill_inverse2d_row_A)
{
  static const int JMAX = 2;
  
  static const double SCALING[3] = {0.1, 1.0, 0.2};
  static const int SOFFSET[3] = {-1, 0, 1};

  gsl_matrix *m;
  int rowstride;
  int size;
  int i;

  rowstride = 1 << JMAX;
  size = rowstride * rowstride;
  m = gsl_matrix_alloc(size, size);
  ck_assert(m != NULL);

  gsl_matrix_set_all(m, 0.0);
  
  ck_assert(generic_matrix_fill_inverse2d_row_A(JMAX,
						3,
						SCALING,
						SOFFSET,
						m,
						0,
						0,
						wavelet_boundary_reflect) >= 0);

  for (i = 0; i < rowstride/2; i ++) {
    ck_assert(gsl_matrix_get(m, 2*i, i) == SCALING[1]);
    if (i == 0) {
      ck_assert(gsl_matrix_get(m, 0, size/4) == (SCALING[0] + SCALING[2]));
    } else {
      ck_assert(gsl_matrix_get(m, 2*i, size/4 + i - 1) == SCALING[0]);
      ck_assert(gsl_matrix_get(m, 2*i, size/4 + i) == SCALING[2]);
    }
  }

  gsl_matrix_free(m);
}
END_TEST

START_TEST (test_generic_matrix_fill_inverse2d_row_B)
{
  static const int JMAX = 2;
  
  static const double WAVELET[3] = {-0.1, -1.0, -0.1};
  static const int WOFFSET[3] = {-1, 0, 1};

  gsl_matrix *m;
  int rowstride;
  int size;
  int i;
  
  rowstride = 1 << JMAX;
  size = rowstride * rowstride;
  m = gsl_matrix_alloc(size, size);
  ck_assert(m != NULL);
  
  ck_assert(generic_matrix_fill_inverse2d_row_B(JMAX,
						3,
						WAVELET,
						WOFFSET,
						m,
						0,
						0,
						wavelet_boundary_reflect) >= 0);

  for (i = 0; i < rowstride/2; i ++) {
    ck_assert(gsl_matrix_get(m, 2*i + 1, size/4 + i) == WAVELET[1]);
    if (i == (rowstride/2 - 1)) { 
      ck_assert(gsl_matrix_get(m, rowstride - 1, rowstride/2 - 1) == (WAVELET[0] + WAVELET[2])); 
    } else {
       ck_assert(gsl_matrix_get(m, 2*i + 1, i) == WAVELET[0]); 
       ck_assert(gsl_matrix_get(m, 2*i + 1, i + 1) == WAVELET[2]); 
    } 
  }

  gsl_matrix_free(m);
}
END_TEST

START_TEST (test_generic_matrix_create_forward2d_col_step_1)
{
  /*
   * Tests a single component wavelet to check location of h0,g0 coefficients in matrix
   */
  static const int JMAX = 2;
  
  static const double SCALING[1] = {1.0};
  static const int SOFFSET[1] = {0};

  static const double WAVELET[1] = {-1.0};
  static const int WOFFSET[1] = {0};

  gsl_matrix *m;
  int width;
  int size;
  int i;
  int j;
  int k;
  int o;

  ck_assert(generic_matrix_create_forward2d_col_step(JMAX,
						     1,
						     SCALING,
						     SOFFSET,
						     1,
						     WAVELET,
						     WOFFSET,
						     &m,
						     wavelet_boundary_reflect) == 0);

  width = 1 << JMAX;
  size = width * width;
  ck_assert(m->size1 == size);
  ck_assert(m->size2 == size);

  /* printf("Forward Col 1\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   for (j = 0; j < size; j ++) { */

  /*     printf("%4.1f ", gsl_matrix_get(m, i, j)); */

  /*   } */
  /*   printf("\n"); */
  /* } */
  /* printf("\n"); */

  /*
   * Check Scaling Rows
   */
  for (i = 0; i < width/2; i ++) {
    o = i*2*width;
    for (j = 0; j < width; j ++) {
      for (k = 0; k < size; k ++) {
	if (k == (o + j)) {
	  ck_assert(gsl_matrix_get(m, width*i + j, k) == 1.0);
	} else {
	  ck_assert(gsl_matrix_get(m, width*i + j, k) == 0.0);
	}
      }
    }
  }

  /*
   * Check Wavelet Rows
   */
  for (i = 0; i < width/2; i ++) {
    o = (2*i + 1) * width;
    for (j = 0; j < width; j ++) {
      for (k = 0; k < size; k ++) {
	if (k == (o + j)) {
	  ck_assert(gsl_matrix_get(m, width/2*width + width*i + j, k) == -1.0);
	} else {
	  ck_assert(gsl_matrix_get(m, width/2*width + width*i + j, k) == 0.0);
	}
      }
    }
  }

  gsl_matrix_free(m);
}
END_TEST

START_TEST (test_generic_matrix_create_inverse2d_col_step_1)
{
  /*
   * Tests a single component wavelet to check location of h0,g0 coefficients in matrix
   */
  static const int JMAX = 2;
  
  static const double SCALING[1] = {1.0};
  static const int SOFFSET[1] = {0};

  static const double WAVELET[1] = {-1.0};
  static const int WOFFSET[1] = {0};

  gsl_matrix *m;
  int width;
  int size;
  int i;
  int j;
  int k;
  int o;

  ck_assert(generic_matrix_create_inverse2d_col_step(JMAX,
						     1,
						     SCALING,
						     SOFFSET,
						     1,
						     WAVELET,
						     WOFFSET,
						     &m,
						     wavelet_boundary_reflect) == 0);

  width = 1 << JMAX;
  size = width * width;
  ck_assert(m->size1 == size);
  ck_assert(m->size2 == size);

  /* printf("Inverse Col 1\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   for (j = 0; j < size; j ++) { */

  /*     printf("%4.1f ", gsl_matrix_get(m, i, j)); */

  /*   } */
  /*   printf("\n"); */
  /* } */
  /* printf("\n"); */

  /*
   * Check Scaling Rows
   */
  for (i = 0; i < width/2; i ++) {
    o = i*width;
    for (j = 0; j < width; j ++) {
      for (k = 0; k < size; k ++) {
	if (k == (o + j)) {
	  ck_assert(gsl_matrix_get(m, 2*width*i + j, k) == 1.0);
	} else {
	  ck_assert(gsl_matrix_get(m, 2*width*i + j, k) == 0.0);
	}
      }
    }
  }

  /*
   * Check Wavelet Rows
   */
  for (i = 0; i < width/2; i ++) {
    o = width/2*width + i*width;
    for (j = 0; j < width; j ++) {
      for (k = 0; k < size; k ++) {
	if (k == (o + j)) {
	  ck_assert(gsl_matrix_get(m, (2*i + 1)*width + j, k) == -1.0);
	} else {
	  ck_assert(gsl_matrix_get(m, (2*i + 1)*width + j, k) == 0.0);
	}
      }
    }
  }

  gsl_matrix_free(m);
}
END_TEST

START_TEST (test_generic_matrix_create_forward2d_col_step_2)
{
  /*
   * Tests a symmetric wavelet transform with more than one coefficient 
   * to verify wrapping
   */
  static const int JMAX = 2;
  
  static const double SCALING[3] = {0.1, 1.0, 0.2};
  static const int SOFFSET[3] = {-1, 0, 1};

  static const double WAVELET[3] = {-0.1, -1.0, -0.2};
  static const int WOFFSET[3] = {-1, 0, 1};

  gsl_matrix *m;
  int width;
  int size;
  int i;
  int j;
  int k;
  int o;

  ck_assert(generic_matrix_create_forward2d_col_step(JMAX,
						     3,
						     SCALING,
						     SOFFSET,
						     3,
						     WAVELET,
						     WOFFSET,
						     &m,
						     wavelet_boundary_reflect) == 0);

  width = 1 << JMAX;
  size = width * width;
  ck_assert(m->size1 == size);
  ck_assert(m->size2 == size);

  /* printf("Forward Col 2\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   for (j = 0; j < size; j ++) { */

  /*     printf("%4.1f ", gsl_matrix_get(m, i, j)); */

  /*   } */
  /*   printf("\n"); */
  /* } */

  /*
   * Check Scaling Rows
   */
  for (i = 0; i < width/2; i ++) {
    o = i*2*width;
    for (j = 0; j < width; j ++) {
      for (k = 0; k < size; k ++) {
	if (k == (o + j)) {
	  ck_assert(gsl_matrix_get(m, width*i + j, k) == SCALING[1]);
	} else if (o < width && (k - width) == (o + j)) {
	  ck_assert(gsl_matrix_get(m, width*i + j, k) == (SCALING[0] + SCALING[2]));
	} else if ((k - width) == (o + j)) {
	  ck_assert(gsl_matrix_get(m, width*i + j, k) == SCALING[2]);
	} else if ((k + width) == (o + j)) {
	  ck_assert(gsl_matrix_get(m, width*i + j, k) == SCALING[0]);
	} else {
	  ck_assert(gsl_matrix_get(m, width*i + j, k) == 0.0);
	}
      }
    }
  }

  /*
   * Check Wavelet Rows
   */
  for (i = 0; i < width/2; i ++) {
    o = (2*i + 1) * width;
    for (j = 0; j < width; j ++) {
      for (k = 0; k < size; k ++) {
	if (k == (o + j)) {
	  ck_assert(gsl_matrix_get(m, width/2*width + width*i + j, k) == WAVELET[1]);
	} else if (o >= (size - width) && (k + width) == (o + j)) {
	  ck_assert(gsl_matrix_get(m, width/2*width + width*i + j, k) == (WAVELET[0] + WAVELET[2]));
	} else if ((k - width) == (o + j)) {
	  ck_assert(gsl_matrix_get(m, width/2*width + width*i + j, k) == WAVELET[2]);
	} else if ((k + width) == (o + j)) {
	  ck_assert(gsl_matrix_get(m, width/2*width + width*i + j, k) == WAVELET[0]);
	} else {
	  ck_assert(gsl_matrix_get(m, width/2*width + width*i + j, k) == 0.0);
	}
      }
    }
  }

  gsl_matrix_free(m);
}
END_TEST

START_TEST (test_generic_matrix_create_inverse2d_col_step_2)
{
  /*
   * Tests a symmetric wavelet transform with more than one coefficient 
   * to verify wrapping
   */
  static const int JMAX = 3;
  
  static const double SCALING[3] = {0.1, 1.0, 0.2};
  static const int SOFFSET[3] = {-1, 0, 1};

  static const double WAVELET[3] = {-0.1, -1.0, -0.2};
  static const int WOFFSET[3] = {-1, 0, 1};

  gsl_matrix *m;
  int width;
  int size;
  int i;
  int j;
  int k;
  int o;

  ck_assert(generic_matrix_create_inverse2d_col_step(JMAX,
						     3,
						     SCALING,
						     SOFFSET,
						     3,
						     WAVELET,
						     WOFFSET,
						     &m,
						     wavelet_boundary_reflect) == 0);

  width = 1 << JMAX;
  size = width * width;
  ck_assert(m->size1 == size);
  ck_assert(m->size2 == size);

  /* printf("Inverse Col 2\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   for (j = 0; j < size; j ++) { */

  /*     printf("%4.1f ", gsl_matrix_get(m, i, j)); */

  /*   } */
  /*   printf("\n"); */
  /* } */
  /* printf("\n"); */
  
  /*
   * Check Scaling Row 1
   */
  i = 0;
  o = width/2*width;
  for (j = 0; j < width; j ++) {
    for (k = 0; k < size; k ++) {
      if (k == (i + j)) {
	ck_assert(gsl_matrix_get(m, j, k) == SCALING[1]);
	
      } else if (k == (o + j)) {
	ck_assert(gsl_matrix_get(m, j, k) == (SCALING[0] + SCALING[2]));
      } else {
	ck_assert(gsl_matrix_get(m, j, k) == 0.0);
      }
    }
  }

  /*
   * Check Scaling Row 2
   */
  i = width;
  o = width/2*width;
  for (j = 0; j < width; j ++) {
    for (k = 0; k < size; k ++) {
      if (k == (i + j)) {
	ck_assert(gsl_matrix_get(m, j + 2*width, k) == SCALING[1]);
	
      } else if (k == (o + j)) {
	ck_assert(gsl_matrix_get(m, j + 2*width, k) == SCALING[0]);
      } else if (k == (o + j + width)) {
 	ck_assert(gsl_matrix_get(m, j + 2*width, k) == SCALING[2]);
      } else {
	ck_assert(gsl_matrix_get(m, j + 2*width, k) == 0.0);
      }
    }
  }

  /*
   * Check Wavelet Row 1
   */
  i = width/2*width;
  o = 0;
  for (j = 0; j < width; j ++) {
    for (k = 0; k < size; k ++) {
      if (k == (i + j)) {
	ck_assert(gsl_matrix_get(m, j + width, k) == WAVELET[1]);
	
      } else if (k == (o + j)) {
	ck_assert(gsl_matrix_get(m, j + width, k) == WAVELET[0]);
      } else if (k == (o + j + width)) {
	ck_assert(gsl_matrix_get(m, j + width, k) == WAVELET[2]);
      } else {
	ck_assert(gsl_matrix_get(m, j + width, k) == 0.0);
      }
    }
  }

  /*
   * Check Wavelet Last Row
   */
  i = width/2*width + width;
  o = width;
  for (j = 0; j < width; j ++) {
    for (k = 0; k < size; k ++) {
      if (k == (size - width + j)) {
	ck_assert(gsl_matrix_get(m, size - width + j, k) == WAVELET[1]);
      } else if (k == (size/2 - width + j)) {
	ck_assert(gsl_matrix_get(m, size - width + j, k) == (WAVELET[0] + WAVELET[2]));
      } else {
	ck_assert(gsl_matrix_get(m, size - width + j, k) == 0.0);
      }
    }
  }
  
  gsl_matrix_free(m);
}
END_TEST

START_TEST (test_generic_matrix_create_forward2d_row_step_1)
{
  /*
   * Tests a single component wavelet to check location of h0,g0 coefficients in matrix
   */
  static const int JMAX = 2;
  
  static const double SCALING[1] = {1.0};
  static const int SOFFSET[1] = {0};

  static const double WAVELET[1] = {-1.0};
  static const int WOFFSET[1] = {0};

  gsl_matrix *m;
  int width;
  int size;
  int i;
  int j;
  int k;
  int o;

  ck_assert(generic_matrix_create_forward2d_row_step(JMAX,
						     1,
						     SCALING,
						     SOFFSET,
						     1,
						     WAVELET,
						     WOFFSET,
						     &m,
						     wavelet_boundary_reflect) == 0);

  width = 1 << JMAX;
  size = width * width;
  ck_assert(m->size1 == size);
  ck_assert(m->size2 == size);

  /* printf("Forward Row 1\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   for (j = 0; j < size; j ++) { */

  /*     printf("%4.1f ", gsl_matrix_get(m, i, j)); */

  /*   } */
  /*   printf("\n"); */
  /* } */
  /* printf("\n"); */
  
  /*
   * Check Scaling Rows 1
   */
  for (j = 0; j < width; j ++) {
    o = j*2;
    for (k = 0; k < size; k ++) {
      if (k == o) {
	ck_assert(gsl_matrix_get(m, j, k) == SCALING[0]);
      } else {
	ck_assert(gsl_matrix_get(m, j, k) == 0.0);
      }
    }
  }

  /*
   * Check Wavelet Rows 1
   */
  for (i = 0; i < width; i ++) {
    o = 2*i + 1;
    for (k = 0; k < size; k ++) {
      if (k == o) {
	ck_assert(gsl_matrix_get(m, i + size/4, k) == WAVELET[0]);
      } else {
	ck_assert(gsl_matrix_get(m, i + size/4, k) == 0.0);
      }
    }
  }

  /*
   * Check Scaling Rows 2
   */
  for (i = 0; i < width; i ++) {
    o = width*width/2 + 2*i;
    for (k = 0; k < size; k ++) {
      if (k == o) {
	ck_assert(gsl_matrix_get(m, i + width*width/2, k) == SCALING[0]);
      } else {
	ck_assert(gsl_matrix_get(m, i + width*width/2, k) == 0.0);
      }
    }
  }

  /*
   * Check Wavelet Rows 1
   */
  for (i = 0; i < width; i ++) {
    o = width*width/2 + 2*i + 1;
    for (k = 0; k < size; k ++) {
      if (k == o) {
	ck_assert(gsl_matrix_get(m, i + 3*size/4, k) == WAVELET[0]);
      } else {
	ck_assert(gsl_matrix_get(m, i + 3*size/4, k) == 0.0);
      }
    }
  }

  gsl_matrix_free(m);
}
END_TEST

START_TEST (test_generic_matrix_create_inverse2d_row_step_1)
{
  /*
   * Tests a single component wavelet to check location of h0,g0 coefficients in matrix
   */
  static const int JMAX = 4;
  
  static const double SCALING[1] = {1.0};
  static const int SOFFSET[1] = {0};

  static const double WAVELET[1] = {-1.0};
  static const int WOFFSET[1] = {0};

  gsl_matrix *m;
  int width;
  int size;
  int i;
  int j;
  int k;
  int o;

  ck_assert(generic_matrix_create_inverse2d_row_step(JMAX,
						     1,
						     SCALING,
						     SOFFSET,
						     1,
						     WAVELET,
						     WOFFSET,
						     &m,
						     wavelet_boundary_reflect) == 0);

  width = 1 << JMAX;
  size = width * width;
  ck_assert(m->size1 == size);
  ck_assert(m->size2 == size);

  /* printf("Forward Row 1\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   for (j = 0; j < size; j ++) { */

  /*     printf("%4.1f ", gsl_matrix_get(m, i, j)); */

  /*   } */
  /*   printf("\n"); */
  /* } */
  /* printf("\n"); */
  
  /*
   * Check Scaling Rows 1
   */
  for (j = 0; j < width; j ++) {
    o = j;
    for (k = 0; k < size; k ++) {
      if (k == o) {
	ck_assert(gsl_matrix_get(m, 2*j, k) == SCALING[0]);
      } else {
	ck_assert(gsl_matrix_get(m, 2*j, k) == 0.0);
      }
    }
  }

  /*
   * Check Wavelet Rows 1
   */
  for (i = 0; i < width; i ++) {
    o = size/4 + i;
    for (k = 0; k < size; k ++) {
      if (k == o) {
	ck_assert(gsl_matrix_get(m, 2*i + 1, k) == WAVELET[0]);
      } else {
	ck_assert(gsl_matrix_get(m, 2*i + 1, k) == 0.0);
      }
    }
  }

  /*
   * Check Scaling Rows 2
   */
  for (i = 0; i < width; i ++) {
    o = size/2 + i;
    for (k = 0; k < size; k ++) {
      if (k == o) {
	ck_assert(gsl_matrix_get(m, size/2 + 2*i, k) == SCALING[0]);
      } else {
	ck_assert(gsl_matrix_get(m, size/2 + 2*i, k) == 0.0);
      }
    }
  }

  /*
   * Check Wavelet Rows 1
   */
  for (i = 0; i < width; i ++) {
    o = 3*size/4 + i;
    for (k = 0; k < size; k ++) {
      if (k == o) {
	ck_assert(gsl_matrix_get(m, size/2 + 2*i + 1, k) == WAVELET[0]);
      } else {
	ck_assert(gsl_matrix_get(m, size/2 + 2*i + 1, k) == 0.0);
      }
    }
  }

  gsl_matrix_free(m);
}
END_TEST

START_TEST (test_generic_matrix_create_forward2d_row_step_2)
{
  /*
   * Tests a single component wavelet to check location of h0,g0 coefficients in matrix
   */
  static const int JMAX = 2;
  
  static const double SCALING[3] = {0.1, 1.0, 0.2};
  static const int SOFFSET[3] = {-1, 0, 1};

  static const double WAVELET[3] = {-0.1, -1.0, -0.2};
  static const int WOFFSET[3] = {-1, 0, 1};

  gsl_matrix *m;
  int width;
  int size;
  int i;
  int j;
  int k;
  int o;

  ck_assert(generic_matrix_create_forward2d_row_step(JMAX,
						     3,
						     SCALING,
						     SOFFSET,
						     3,
						     WAVELET,
						     WOFFSET,
						     &m,
						     wavelet_boundary_reflect) == 0);

  width = 1 << JMAX;
  size = width * width;
  ck_assert(m->size1 == size);
  ck_assert(m->size2 == size);

  /* printf("Forward Row 2\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   for (j = 0; j < size; j ++) { */

  /*     printf("%4.1f ", gsl_matrix_get(m, i, j)); */

  /*   } */
  /*   printf("\n"); */
  /* } */
  /* printf("\n"); */

  /*
   * Check Scaling Rows 1
   */
  for (j = 0; j < width/2; j ++) {
    o = j*2;
    for (k = 0; k < size; k ++) {
      if (k == o) {
	ck_assert(gsl_matrix_get(m, j, k) == SCALING[1]);
      } else {
	if (j == 0) {
	  if (k == 1) {
	    ck_assert(gsl_matrix_get(m, j, k) == (SCALING[0] + SCALING[2]));
	  } else {
	    ck_assert(gsl_matrix_get(m, j, k) == 0.0);
	  }
	} else if (k == (o - 1)) {
	  ck_assert(gsl_matrix_get(m, j, k) == SCALING[0]);
	} else if (k == (o + 1)) {
	  ck_assert(gsl_matrix_get(m, j, k) == SCALING[2]);
	} else {
	  ck_assert(gsl_matrix_get(m, j, k) == 0.0);
	}
      }
    }
  }

  /*
   * Check Wavelet Rows 1
   */
  for (i = 0; i < width/2; i ++) {
    o = 2*i + 1;
    for (k = 0; k < size; k ++) {
      if (k == o) {
	ck_assert(gsl_matrix_get(m, i + size/4, k) == WAVELET[1]);
      } else {
	if (i == (width/2 - 1)) {
	  if (k == (o - 1)) {
	    ck_assert(gsl_matrix_get(m, i + size/4, k) == (WAVELET[0] + WAVELET[2]));
	  } else {
	    ck_assert(gsl_matrix_get(m, i + size/4, k) == 0.0);
	  }
	} else if (k == (o - 1)) {
	  ck_assert(gsl_matrix_get(m, i + size/4, k) == WAVELET[0]);
	} else if (k == (o + 1)) {
	  ck_assert(gsl_matrix_get(m, i + size/4, k) == WAVELET[2]);
	} else {
	  ck_assert(gsl_matrix_get(m, i + size/4, k) == 0.0);
	}
      }
    }
  }


  gsl_matrix_free(m);
}
END_TEST

START_TEST (test_generic_matrix_create_inverse2d_row_step_2)
{
  /*
   * Tests a single component wavelet to check location of h0,g0 coefficients in matrix
   */
  static const int JMAX = 2;
  
  static const double SCALING[3] = {0.1, 1.0, 0.1};
  static const int SOFFSET[3] = {-1, 0, 1};

  static const double WAVELET[3] = {-0.1, -1.0, -0.1};
  static const int WOFFSET[3] = {-1, 0, 1};

  gsl_matrix *m;
  int width;
  int size;
  int i;
  int j;
  int k;
  int o;

  ck_assert(generic_matrix_create_inverse2d_row_step(JMAX,
						     3,
						     SCALING,
						     SOFFSET,
						     3,
						     WAVELET,
						     WOFFSET,
						     &m,
						     wavelet_boundary_reflect) == 0);

  width = 1 << JMAX;
  size = width * width;
  ck_assert(m->size1 == size);
  ck_assert(m->size2 == size);

  /* printf("Row Inverse 2\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   for (j = 0; j < size; j ++) { */

  /*     printf("%4.1f ", gsl_matrix_get(m, i, j)); */

  /*   } */
  /*   printf("\n"); */
  /* } */
  /* printf("\n"); */
  
  /*
   * Check Scaling Rows 1
   */
  for (j = 0; j < width/2; j ++) {
    o = j;
    for (k = 0; k < size; k ++) {
      if (k == o) {
	ck_assert(gsl_matrix_get(m, 2*j, k) == SCALING[1]);
      } else if (k == o + size/4) {
	if (j == 0) {
	  ck_assert(gsl_matrix_get(m, 2*j, k) == (SCALING[0] + SCALING[2]));
	} else {
	  ck_assert(gsl_matrix_get(m, 2*j, k) == SCALING[2]);
	}
      } else if (k == o + size/4 - 1) {
	if (j != 0) {
	  ck_assert(gsl_matrix_get(m, 2*j, k) == SCALING[0]);
	} else {
	  ck_assert(gsl_matrix_get(m, 2*j, k) == 0.0);
	}
      } else {
	ck_assert(gsl_matrix_get(m, 2*j, k) == 0.0);
      }
    }
  }

  /*
   * Check Wavelet Rows 1
   */
  for (i = 0; i < width/2; i ++) {
    o = size/4 + i;
    for (k = 0; k < size; k ++) {
      if (k == o) {
	ck_assert(gsl_matrix_get(m, 2*i + 1, k) == WAVELET[1]);
      } else if (k == (o - size/4)) {
	if (i == (width/2 - 1)) {
	  ck_assert(gsl_matrix_get(m, 2*i + 1, k) == (WAVELET[0] + WAVELET[2]));
	} else {
	  ck_assert(gsl_matrix_get(m, 2*i + 1, k) == WAVELET[0]);
	}
      } else if (k == (o - size/4 + 1)) {
	if (i != (width/2 - 1)) {
	  ck_assert(gsl_matrix_get(m, 2*i + 1, k) == WAVELET[2]);
	} else {
	  ck_assert(gsl_matrix_get(m, 2*i + 1, k) == 0.0);
	}
      } else {
	ck_assert(gsl_matrix_get(m, 2*i + 1, k) == 0.0);
      }
    }
  }

  gsl_matrix_free(m);
}
END_TEST


Suite *
cdf97_suite (void)
{
  Suite *s = suite_create ("Generic Matrix");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  
  tcase_add_test (tc_core, test_wavelet_boundary_reflect);
  tcase_add_test (tc_core, test_wavelet_boundary_periodic);
  tcase_add_test (tc_core, test_generic_matrix_interlace_index);

  tcase_add_test (tc_core, test_generic_matrix_fill_forward2d_col_A);
  tcase_add_test (tc_core, test_generic_matrix_fill_forward2d_col_B);
  
  tcase_add_test (tc_core, test_generic_matrix_fill_inverse2d_col_A);
  tcase_add_test (tc_core, test_generic_matrix_fill_inverse2d_col_B);
	  
  tcase_add_test (tc_core, test_generic_matrix_fill_forward2d_row_A);
  tcase_add_test (tc_core, test_generic_matrix_fill_forward2d_row_B);
  
  tcase_add_test (tc_core, test_generic_matrix_fill_inverse2d_row_A);
  tcase_add_test (tc_core, test_generic_matrix_fill_inverse2d_row_B);

  tcase_add_test (tc_core, test_generic_matrix_create_forward2d_col_step_1);
  tcase_add_test (tc_core, test_generic_matrix_create_inverse2d_col_step_1);

  tcase_add_test (tc_core, test_generic_matrix_create_forward2d_col_step_2);
  tcase_add_test (tc_core, test_generic_matrix_create_inverse2d_col_step_2);

  tcase_add_test (tc_core, test_generic_matrix_create_forward2d_row_step_1);
  tcase_add_test (tc_core, test_generic_matrix_create_inverse2d_row_step_1);

  tcase_add_test (tc_core, test_generic_matrix_create_forward2d_row_step_2);
  tcase_add_test (tc_core, test_generic_matrix_create_inverse2d_row_step_2);

  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = cdf97_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
