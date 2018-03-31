
#include <stdio.h>
#include <stdlib.h>
#include <check.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include "cdf97_matrix.h"
#include "cdf97_lift.h"

START_TEST (test_cdf97_matrix_col_step_versus_lift)
{
  static const int J_MAX = 3;
  static const double INITIAL_VALUE = 3.1415;
  
  gsl_matrix *c;

  gsl_vector *n;
  gsl_vector *v;

  int width;
  int size;
  int i;
  int j;

  double *data;
  double *workspace;

  double t;
  
  width = 1 << J_MAX;
  size = width * width;

  v = gsl_vector_alloc(size);
  ck_assert(v != NULL);

  n = gsl_vector_alloc(size);
  ck_assert(n != NULL);

  /*
   * Create an individual row steps
   */
  ck_assert(cdf97_matrix_forward2d_create_col_step(J_MAX, J_MAX, &c) >= 0);
  ck_assert(c->size1 == size);
  ck_assert(c->size2 == size);

  /* printf("Forward 2D Col Step\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   for (j = 0; j < size; j ++) { */
  /*     fprintf(stderr, "%6.4f ", gsl_matrix_get(c, i, j)); */
  /*   } */
  /*   fprintf(stderr, "\n"); */
  /* } */
  
  data = malloc(sizeof(double) * size);
  ck_assert(data != NULL);

  workspace = malloc(sizeof(double) * width);
  ck_assert(workspace != NULL);


  /*
   * First Test a Constant value
   */
  for (i = 0; i < size; i ++) {
    data[i] = INITIAL_VALUE;
  }
  gsl_vector_set_all(v, INITIAL_VALUE);

  /*
   * Col matrix transform
   */
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, c, v, 0.0, n) >= 0);

  /*
   * Equivalent lifting steps
   */
  for (i = 0; i < width; i ++) {
    ck_assert(cdf97_lift_forward1d_cdf97_step(data + i,
					      width,
					      width,
					      workspace) >= 0);
  }

  for (i = 0; i < size; i ++) {
    ck_assert(fabs(gsl_vector_get(n, i) - data[i]) < 1.0e-6);
    /* printf("%3d %10.6f %10.6f\n", i, gsl_vector_get(n, i), data[i]); */
  }

  /*
   * Next test a sinusoidal function
   */
  for (j = 0; j < width; j ++) {
    for (i = 0; i < width; i ++) {
      t = sin((double)i/64.0 * M_PI * 2.0) * cos((double)j/8.0 * M_PI * 2.0);
      gsl_vector_set(v, j * width + i, t);
      data[j*width + i] = t;
    }
  }

  /*
   * Col matrix transform
   */
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, c, v, 0.0, n) >= 0);

  /*
   * Equivalent lifting steps
   */
  for (i = 0; i < width; i ++) {
    ck_assert(cdf97_lift_forward1d_cdf97_step(data + i,
					      width,
					      width,
					      workspace) >= 0);
  }

  for (i = 0; i < size; i ++) {
    /* printf("%3d %10.6f %10.6f\n", i, gsl_vector_get(n, i), data[i]); */
    ck_assert(fabs(gsl_vector_get(n, i) - data[i]) < 1.0e-6);
  }
  
  free(data);
  free(workspace);

  gsl_vector_free(v);
  gsl_vector_free(n);

  gsl_matrix_free(c);
}
END_TEST

START_TEST (test_cdf97_matrix_inverse_col_step_versus_lift)
{
  static const int J_MAX = 3;
  static const double INITIAL_VALUE = 3.1415;
  
  gsl_matrix *ci;

  gsl_vector *n;
  gsl_vector *v;

  int width;
  int size;
  int i;
  int j;

  double *data;
  double *workspace;

  double t;
  
  width = 1 << J_MAX;
  size = width * width;

  v = gsl_vector_alloc(size);
  ck_assert(v != NULL);

  n = gsl_vector_alloc(size);
  ck_assert(n != NULL);

  /*
   * Create an individual row steps
   */
  ck_assert(cdf97_matrix_inverse2d_create_col_step(J_MAX, J_MAX, &ci) >= 0);
  ck_assert(ci->size1 == size);
  ck_assert(ci->size2 == size);

  /* printf("Forward 2D Col Step\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   for (j = 0; j < size; j ++) { */
  /*     fprintf(stderr, "%6.4f ", gsl_matrix_get(ci, i, j)); */
  /*   } */
  /*   fprintf(stderr, "\n"); */
  /* } */
  
  data = malloc(sizeof(double) * size);
  ck_assert(data != NULL);

  workspace = malloc(sizeof(double) * width);
  ck_assert(workspace != NULL);


  /*
   * First Test a Constant value
   */
  for (j = 0; j < width/2; j ++) {
    for (i = 0; i < width; i ++) {
    
      data[j*width + i] = INITIAL_VALUE;
      data[(j + width/2)*width + i] = 0;
      
      gsl_vector_set(v, j*width + i, INITIAL_VALUE);
      gsl_vector_set(v, (j + width/2)*width + i, 0.0);
    }
  }

  /*
   * Col matrix transform
   */
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, ci, v, 0.0, n) >= 0);

  /*
   * Equivalent lifting steps
   */
  for (i = 0; i < width; i ++) {
    ck_assert(cdf97_lift_inverse1d_cdf97_step(data + i,
					      width,
					      width,
					      workspace) >= 0);
  }

  /* printf("Inverse Col Results: constant\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   printf("%3d %10.6f %10.6f\n", i, gsl_vector_get(n, i), data[i]); */
  /* } */
  
  for (i = 0; i < size; i ++) {
    ck_assert(fabs(gsl_vector_get(n, i) - data[i]) < 1.0e-6);
  }

  /*
   * Next test a sinusoidal function
   */
  for (j = 0; j < width; j ++) {
    for (i = 0; i < width; i ++) {
      t = sin((double)i/64.0 * M_PI * 2.0) * cos((double)j/8.0 * M_PI * 2.0);
      gsl_vector_set(v, j * width + i, t);
      data[j*width + i] = t;
    }
  }

  /*
   * Col matrix transform
   */
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, ci, v, 0.0, n) >= 0);

  /*
   * Equivalent lifting steps
   */
  for (i = 0; i < width; i ++) {
    ck_assert(cdf97_lift_inverse1d_cdf97_step(data + i,
					      width,
					      width,
					      workspace) >= 0);
  }

  /* printf("Inverse Col Results: sinusoidal\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   printf("%3d %10.6f %10.6f\n", i, gsl_vector_get(n, i), data[i]); */
  /* } */

  for (i = 0; i < size; i ++) {
    ck_assert(fabs(gsl_vector_get(n, i) - data[i]) < 1.0e-6);
  }
  
  free(data);
  free(workspace);

  gsl_vector_free(v);
  gsl_vector_free(n);

  gsl_matrix_free(ci);
}
END_TEST

START_TEST (test_cdf97_matrix_row_step_versus_lift)
{
  static const int J_MAX = 3;
  static const double INITIAL_VALUE = 3.1415;
  
  gsl_matrix *r;

  gsl_vector *n;
  gsl_vector *v;

  int width;
  int size;
  int i;
  int j;
  int k;
  
  double *data;
  double *workspace;

  double t;

  int vi;
  int li;
  
  width = 1 << J_MAX;
  size = width * width;

  v = gsl_vector_alloc(size);
  ck_assert(v != NULL);

  n = gsl_vector_alloc(size);
  ck_assert(n != NULL);

  /*
   * Create an individual row steps
   */
  ck_assert(cdf97_matrix_forward2d_create_row_step(J_MAX, J_MAX, &r) >= 0);
  ck_assert(r->size1 == size);
  ck_assert(r->size2 == size);

  /* printf("Forward 2D Row Step\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   for (j = 0; j < size; j ++) { */
  /*     fprintf(stderr, "%6.4f ", gsl_matrix_get(r, i, j)); */
  /*   } */
  /*   fprintf(stderr, "\n"); */
  /* } */
  
  data = malloc(sizeof(double) * size);
  ck_assert(data != NULL);

  workspace = malloc(sizeof(double) * width);
  ck_assert(workspace != NULL);

  /*
   * First Test a Constant value
   */
  for (i = 0; i < size; i ++) {
    data[i] = INITIAL_VALUE;
  }
  gsl_vector_set_all(v, INITIAL_VALUE);

  /*
   * Col matrix transform
   */
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, r, v, 0.0, n) >= 0);

  /*
   * Equivalent lifting steps
   */
  for (i = 0; i < width; i ++) {
    ck_assert(cdf97_lift_forward1d_cdf97_step(data + i*width,
					      width,
					      1,
					      workspace) >= 0);
  }

  for (k = 0; k < 4; k ++) {
    for (j = 0; j < width/2; j ++) {
      for (i = 0; i < width/2; i ++) {
	vi = k*width*width/4 + j*width/2 + i;
	li = k/2 * width + k%2 * width/2 + j*width + i;

	/* printf("%d %d %d: %10.6f %10.6f\n", k, j, i, gsl_vector_get(n, vi), data[li]); */
	ck_assert(fabs(gsl_vector_get(n, vi) - data[li]) < 1.0e-6);
	
      }
    }
  }
      

  /*
   * Next test a sinusoidal function
   */
  for (j = 0; j < width; j ++) {
    for (i = 0; i < width; i ++) {
      t = sin((double)i/64.0 * M_PI * 2.0) * cos((double)j/8.0 * M_PI * 2.0);
      gsl_vector_set(v, j * width + i, t);
      data[j*width + i] = t;
    }
  }

  /*
   * Row matrix transform
   */
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, r, v, 0.0, n) >= 0);

  /*
   * Equivalent lifting steps
   */
  for (i = 0; i < width; i ++) {
    ck_assert(cdf97_lift_forward1d_cdf97_step(data + i*width,
					      width,
					      1,
					      workspace) >= 0);
  }

  for (k = 0; k < 4; k ++) {
    for (j = 0; j < width/2; j ++) {
      for (i = 0; i < width/2; i ++) {
	vi = k*width*width/4 + j*width/2 + i;
	li = k/2*width*width/2 + k%2 * width/2 + j*width + i;

	/* printf("%d %d %d %3d %3d: %10.6f %10.6f\n", k, j, i, vi, li, gsl_vector_get(n, vi), data[li]);  */
	ck_assert(fabs(gsl_vector_get(n, vi) - data[li]) < 1.0e-6);
	
      }
    }
  }
  
  free(data);
  free(workspace);

  gsl_vector_free(v);
  gsl_vector_free(n);

  gsl_matrix_free(r);
}
END_TEST

START_TEST (test_cdf97_matrix_inverse_row_step_versus_lift)
{
  static const int J_MAX = 3;
  static const double INITIAL_VALUE = 3.1415;
  
  gsl_matrix *ri;

  gsl_vector *n;
  gsl_vector *v;

  int width;
  int size;
  int i;
  int j;
  int k;
  
  double *data;
  double *workspace;

  double t;

  int vi;
  int li;
  
  width = 1 << J_MAX;
  size = width * width;

  v = gsl_vector_alloc(size);
  ck_assert(v != NULL);

  n = gsl_vector_alloc(size);
  ck_assert(n != NULL);

  /*
   * Create an individual row steps
   */
  ck_assert(cdf97_matrix_inverse2d_create_row_step(J_MAX, J_MAX, &ri) >= 0);
  ck_assert(ri->size1 == size);
  ck_assert(ri->size2 == size);

  /* printf("Forward 2D Row Step\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   for (j = 0; j < size; j ++) { */
  /*     fprintf(stderr, "%6.4f ", gsl_matrix_get(r, i, j)); */
  /*   } */
  /*   fprintf(stderr, "\n"); */
  /* } */
  
  data = malloc(sizeof(double) * size);
  ck_assert(data != NULL);

  workspace = malloc(sizeof(double) * width);
  ck_assert(workspace != NULL);

  /*
   * First Test a Constant value
   */
  for (k = 0; k < 4; k ++) {
    for (j = 0; j < width/2; j ++) {
      for (i = 0; i < width/2; i ++) {
	vi = k*width*width/4 + j*width/2 + i;
	li = k/2 * width * width/2 + k%2 * width/2 + j*width + i;

	if (k % 2 == 0) {
	  data[li] = INITIAL_VALUE;
	  gsl_vector_set(v, vi, INITIAL_VALUE);
	} else {
	  data[li] = 0.0;
	  gsl_vector_set(v, vi, 0.0);
	}
      }
    }
  }

  /*
   * Row matrix transform
   */
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, ri, v, 0.0, n) >= 0);

  /*
   * Equivalent lifting steps
   */
  for (i = 0; i < width; i ++) {
    ck_assert(cdf97_lift_inverse1d_cdf97_step(data + i*width,
					      width,
					      1,
					      workspace) >= 0);
  }

  for (i = 0; i < size; i ++) {
    /* printf("%3d: %10.6f %10.6f\n", i, gsl_vector_get(n, vi), data[li]); */
    ck_assert(fabs(gsl_vector_get(n, i) - data[i]) < 1.0e-6);
	
  }
      
  /*
   * Next test a sinusoidal function
   */
  for (k = 0; k < 4; k ++) {
    for (j = 0; j < width/2; j ++) {
      for (i = 0; i < width/2; i ++) {
	vi = k*width*width/4 + j*width/2 + i;
	li = k/2 * width * width/2 + k%2 * width/2 + j*width + i;

	t = sin((double)i/64.0 * M_PI * 2.0) * cos((double)j/8.0 * M_PI * 2.0);
	
	if (k % 2 == 0) {
	  data[li] = t;
	  gsl_vector_set(v, vi, t);
	} else {
	  data[li] = 0.0;
	  gsl_vector_set(v, vi, 0.0);
	}
      }
    }
  }

  /*
   * Row matrix transform
   */
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, ri, v, 0.0, n) >= 0);

  /*
   * Equivalent lifting steps
   */
  for (i = 0; i < width; i ++) {
    ck_assert(cdf97_lift_inverse1d_cdf97_step(data + i*width,
					      width,
					      1,
					      workspace) >= 0);
  }

  for (i = 0; i < size; i ++) {
    /* printf("%3d: %10.6f %10.6f\n", i, gsl_vector_get(n, vi), data[li]); */
    ck_assert(fabs(gsl_vector_get(n, i) - data[i]) < 1.0e-6);
  }
  
  free(data);
  free(workspace);

  gsl_vector_free(v);
  gsl_vector_free(n);

  gsl_matrix_free(ri);
}
END_TEST

START_TEST (test_cdf97_matrix_step_constant)
{
  static const int J_MAX = 3;
  static const double INITIAL_VALUE = 3.1415;
  
  gsl_matrix *r;
  gsl_matrix *ri;
  
  gsl_matrix *c;
  gsl_matrix *ci;

  gsl_matrix *m;
  gsl_matrix *mi;
  
  gsl_vector *n;
  gsl_vector *v;
  int width;
  int size;
  int i;
  int j;
  
  width = 1 << J_MAX;
  size = width * width;

  v = gsl_vector_alloc(size);
  ck_assert(v != NULL);

  n = gsl_vector_alloc(size);
  ck_assert(n != NULL);

  m = gsl_matrix_alloc(size, size);
  ck_assert(m != NULL);
  
  /*
   * Create an individual row steps
   */
  ck_assert(cdf97_matrix_forward2d_create_row_step(J_MAX, J_MAX, &r) >= 0);
  ck_assert(r->size1 == size);
  ck_assert(r->size2 == size);

  /* printf("Forward 2D Row Step\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   for (j = 0; j < size; j ++) { */
  /*     fprintf(stderr, "%6.4f ", gsl_matrix_get(r, i, j)); */
  /*   } */
  /*   fprintf(stderr, "\n"); */
  /* } */
  
  ck_assert(cdf97_matrix_inverse2d_create_row_step(J_MAX, J_MAX, &ri) >= 0);
  ck_assert(ri->size1 == size);
  ck_assert(ri->size2 == size);

  /* printf("Inverse 2D Row Step\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   for (j = 0; j < size; j ++) { */
  /*     fprintf(stderr, "%6.4f ", gsl_matrix_get(ri, i, j)); */
  /*   } */
  /*   fprintf(stderr, "\n"); */
  /* } */

  /*
   * Test product of transform matrices (Not identity matrix as expected?)
   */
  ck_assert(gsl_blas_dgemm(CblasNoTrans,
			   CblasNoTrans,
			   1.0,
			   r,
			   ri,
			   0.0,
			   m) >= 0);
  /* printf("Product 2D Row Step\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   for (j = 0; j < size; j ++) { */
  /*     fprintf(stderr, "%6.4f ", gsl_matrix_get(m, i, j)); */
  /*   } */
  /*   fprintf(stderr, "\n"); */
  /* } */
  
  
  /*
   * Test doing a row transformation forward and reverse
   */
  gsl_vector_set_all(v, INITIAL_VALUE);
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, r, v, 0.0, n) >= 0);

  /* printf("Forward Row Transform Result\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   printf("%3d %10.6f\n", i, gsl_vector_get(n, i)); */
  /* } */
  
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, ri, n, 0.0, v) >= 0);

  /* printf("Inverse Row Transform Result\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   printf("%3d %10.6f\n", i, gsl_vector_get(v, i)); */
  /* } */
  for (i = 0; i < size; i ++) {
    ck_assert(fabs(gsl_vector_get(v, i) - INITIAL_VALUE) < 1.0e-6); 
  }

  /*
   * Create an individual col steps
   */
  ck_assert(cdf97_matrix_forward2d_create_col_step(J_MAX, J_MAX, &c) >= 0);
  ck_assert(c->size1 == size);
  ck_assert(c->size2 == size);

  
  ck_assert(cdf97_matrix_inverse2d_create_col_step(J_MAX, J_MAX, &ci) >= 0);
  ck_assert(ci->size1 == size);
  ck_assert(ci->size2 == size);

  /*
   * Test doing a row transformation forward and reverse
   */
  gsl_vector_set_all(v, INITIAL_VALUE);
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, c, v, 0.0, n) >= 0);

  /* printf("Forward Col Transform Result\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   printf("%d %f\n", i, gsl_vector_get(n, i)); */
  /* } */
  
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, ci, n, 0.0, v) >= 0);

  /* printf("Inverse Col Transform Result\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   printf("%d %f\n", i, gsl_vector_get(v, i)); */
  /* } */
  
  for (i = 0; i < size; i ++) {
    ck_assert(fabs(gsl_vector_get(v, i) - INITIAL_VALUE) < 1.0e-6);
  }

  /*
   * Create combined steps
   */
  gsl_matrix_free(m);
  ck_assert(cdf97_matrix_forward2d_create_step(J_MAX, J_MAX, &m) >= 0);
  ck_assert(cdf97_matrix_inverse2d_create_step(J_MAX, J_MAX, &mi) >= 0);

  /*
   * Test doing a transform forward and reverse
   */
  gsl_vector_set_all(v, INITIAL_VALUE);
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, m, v, 0.0, n) >= 0);

  /* printf("Forward Transform Result\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   printf("%d %f\n", i, gsl_vector_get(n, i)); */
  /* } */
  
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, mi, n, 0.0, v) >= 0);

  /* printf("Inverse Transform Result\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   printf("%d %f\n", i, gsl_vector_get(v, i)); */
  /* } */

  for (i = 0; i < size; i ++) {
    ck_assert(fabs(gsl_vector_get(v, i) - INITIAL_VALUE) < 1.0e-6);
  }

  gsl_matrix_free(r);
  gsl_matrix_free(ri);
  gsl_matrix_free(c);
  gsl_matrix_free(ci);
  gsl_matrix_free(m);
  gsl_matrix_free(mi);

  gsl_vector_free(n);
  gsl_vector_free(v);
}
END_TEST

START_TEST(test_cdf97_matrix_step_sinusoid)
{
  static const int J_MAX = 2;

  double t;
  
  int i;
  int j;

  gsl_vector *v;
  gsl_vector *c;
  gsl_vector *n;
  gsl_matrix *W;
  gsl_matrix *Wi;

  gsl_matrix *r;
  gsl_matrix *ri;

  gsl_matrix *C;
  gsl_matrix *Ci;

  int width;
  int size;

  width = 1 << J_MAX;
  size = width * width;

  ck_assert(cdf97_matrix_forward2d_create_step(J_MAX, J_MAX, &W) >= 0);
  ck_assert(cdf97_matrix_inverse2d_create_step(J_MAX, J_MAX, &Wi) >= 0);

  ck_assert(cdf97_matrix_forward2d_create_row_step(J_MAX, J_MAX, &r) >= 0);
  ck_assert(cdf97_matrix_inverse2d_create_row_step(J_MAX, J_MAX, &ri) >= 0);

  ck_assert(cdf97_matrix_forward2d_create_col_step(J_MAX, J_MAX, &C) >= 0);
  ck_assert(cdf97_matrix_inverse2d_create_col_step(J_MAX, J_MAX, &Ci) >= 0);

  v = gsl_vector_alloc(size);
  ck_assert(v != NULL);

  c = gsl_vector_alloc(size);
  ck_assert(c != NULL);

  n = gsl_vector_alloc(size);
  ck_assert(n != NULL);

  for (j = 0; j < width; j ++) {
    
    for (i = 0; i < width; i ++) {
      t = sin((double)i/64.0 * M_PI * 2.0) * cos((double)j/8.0 * M_PI * 2.0);
      gsl_vector_set(v, j * width + i, t);
    }

  }

  /*
   * Forward row transform on sinusoid signal
   */
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, r, v, 0.0, c) >= 0);

  /* printf("Forward Row Result Sinusoid\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   printf("%3d %10.6f\n", i, gsl_vector_get(c, i)); */
  /* } */

  /*
   * Inverse row transform on sinusoid signal
   */
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, ri, c, 0.0, n) >= 0);

  /* printf("Inverse Row Result Sinusoid\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   printf("%3d %10.6f %10.6f\n", i, gsl_vector_get(n, i), gsl_vector_get(v, i)); */
  /* } */
  
  for (i = 0; i < size; i ++) {
    /* ck_assert(fabs(gsl_vector_get(n, i) - gsl_vector_get(v, i)) < 1.0e-6); */
  }

  /*
   * Forward col transform on sinusoid signal
   */
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, C, v, 0.0, c) >= 0);

  /* printf("Forward Col Result Sinusoid\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   printf("%3d %10.6f\n", i, gsl_vector_get(c, i)); */
  /* } */

  /*
   * Inverse col transform on sinusoid signal
   */
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, Ci, c, 0.0, n) >= 0);

  /* printf("Inverse Col Result Sinusoid\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   printf("%3d %10.6f %10.6f\n", i, gsl_vector_get(n, i), gsl_vector_get(v, i)); */
  /* } */
  
  for (i = 0; i < size; i ++) {
    ck_assert(fabs(gsl_vector_get(n, i) - gsl_vector_get(v, i)) < 1.0e-6);
  }


  /*
   * Forward transform on sinusoid signal
   */
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, W, v, 0.0, c) >= 0);

  /* printf("Forward Result Sinusoid\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   printf("%3d %10.6f\n", i, gsl_vector_get(c, i)); */
  /* } */

  /*
   * Inverse transform on sinusoid signal
   */
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, Wi, c, 0.0, n) >= 0);

  /* printf("Inverse Result Sinusoid\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   printf("%3d %10.6f\n", i, gsl_vector_get(n, i)); */
  /* } */
  
  for (i = 0; i < size; i ++) {
    ck_assert(fabs(gsl_vector_get(n, i) - gsl_vector_get(v, i)) < 1.0e-6);
  }

      
  gsl_vector_free(n);
  gsl_vector_free(v);
  gsl_vector_free(c);
  gsl_matrix_free(W);
}
END_TEST

START_TEST (test_cdf97_matrix_constant)
{
  static const int J_MAX = 3;
  static const double INITIAL_VALUE = 3.1415;
  
  gsl_matrix *m;
  gsl_matrix *mi;
  
  gsl_vector *n;
  gsl_vector *v;
  int width;
  int size;
  int i;
  
  width = 1 << J_MAX;
  size = width * width;
  
  ck_assert(cdf97_matrix_forward2d_create(J_MAX, &m) >= 0);
  ck_assert(cdf97_matrix_inverse2d_create(J_MAX, &mi) >= 0);

  ck_assert(m->size1 == size);
  ck_assert(m->size2 == size);

  v = gsl_vector_alloc(size);
  ck_assert(v != NULL);

  gsl_vector_set_all(v, INITIAL_VALUE);

  n = gsl_vector_alloc(size);
  ck_assert(n != NULL);

  /*
   * Forward transform on constant signal
   */
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, m, v, 0.0, n) >= 0);

  ck_assert(fabs(gsl_vector_get(n, 0) - INITIAL_VALUE) < 1.0e-6);
  for (i = 1; i < size; i ++) {
    ck_assert(fabs(gsl_vector_get(n, i)) < 1.0e-6);
    /* printf("%d %g\n", i, gsl_vector_get(n, i)); */
  }

  /*
   * Inverse transform to recover constant signal
   */
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, mi, n, 0.0, v) >= 0);
  for (i = 0; i < size; i ++) {
    ck_assert(fabs(gsl_vector_get(v, i) - INITIAL_VALUE) < 1.0e-6);
    /* printf("%d %g\n", i, gsl_vector_get(v, i)); */
  }

  gsl_matrix_free(m);
  gsl_matrix_free(mi);
  gsl_vector_free(v);
  gsl_vector_free(n);
  
}
END_TEST

START_TEST(test_cdf97_matrix_sinusoid)
{
  static const int J_MAX = 3;

  double t;
  
  int i;
  int j;

  gsl_vector *v;
  gsl_vector *c;
  gsl_vector *n;
  gsl_matrix *W;
  gsl_matrix *Wi;
  
  int width;
  int size;

  width = 1 << J_MAX;
  size = width * width;

  ck_assert(cdf97_matrix_forward2d_create(J_MAX, &W) >= 0);
  ck_assert(cdf97_matrix_inverse2d_create(J_MAX, &Wi) >= 0);

  
  v = gsl_vector_alloc(size);
  ck_assert(v != NULL);

  c = gsl_vector_alloc(size);
  ck_assert(c != NULL);

  n = gsl_vector_alloc(size);
  ck_assert(n != NULL);

  for (j = 0; j < width; j ++) {
    
    for (i = 0; i < width; i ++) {
      t = sin((double)i/64.0 * M_PI * 2.0) * cos((double)j/8.0 * M_PI * 2.0);
      gsl_vector_set(v, j * width + i, t);
    }

  }

  /*
   * Forward transform on sinusoid signal
   */
  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, W, v, 0.0, c) >= 0);

  /* printf("Forward Transform Result:\n"); */
  /* printf("0 %g\n", gsl_vector_get(n, 0)); */
  /* for (i = 1; i < size; i ++) { */
  /*   printf("%d %g\n", i, gsl_vector_get(c, i)); */
  /* } */

  /*
   * Inverse transform to recover constant signal
   */

  ck_assert(gsl_blas_dgemv(CblasNoTrans, 1.0, Wi, c, 0.0, n) >= 0);

  /* printf("Inverse Transform Result:\n"); */
  /* for (i = 0; i < size; i ++) { */
  /*   printf("%3d %10.6f %10.6f\n", i, gsl_vector_get(n, i), gsl_vector_get(v, i)); */
  /* } */
  
  for (i = 0; i < size; i ++) {
    ck_assert(fabs(gsl_vector_get(v, i) - gsl_vector_get(n, i)) < 1.0e-6);
  }

  gsl_vector_free(v);
  gsl_vector_free(c);
  gsl_vector_free(n);
  gsl_matrix_free(W);
  gsl_matrix_free(Wi);
}
END_TEST

Suite *
cdf97_suite (void)
{
  Suite *s = suite_create ("CDF97 Matrix");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  
  tcase_add_test(tc_core, test_cdf97_matrix_col_step_versus_lift);
  tcase_add_test(tc_core, test_cdf97_matrix_inverse_col_step_versus_lift);
  tcase_add_test(tc_core, test_cdf97_matrix_row_step_versus_lift);
  tcase_add_test(tc_core, test_cdf97_matrix_inverse_row_step_versus_lift);
  
  tcase_add_test (tc_core, test_cdf97_matrix_step_constant);
  tcase_add_test (tc_core, test_cdf97_matrix_step_sinusoid);

  tcase_add_test (tc_core, test_cdf97_matrix_constant);
  tcase_add_test (tc_core, test_cdf97_matrix_sinusoid);
  
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
