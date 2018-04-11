
#include <stdio.h>
#include <stdlib.h>
#include <check.h>

#include <math.h>

#include <gsl/gsl_rng.h>

#include "manifold.h"
#include "octahedron.h"
#include "icosahedron.h"
#include "face_subdivision.h"

START_TEST(test_octahedron_single_step)
{
  static const int DEGREE = 2;
  gsl_rng *rng;
  
  manifold_t *o;
  double *input_coeff;
  double *recon_coeff;
  int exact;
  int sstart;
  int send;

  int i;
  
  rng = gsl_rng_alloc(gsl_rng_taus);
  ck_assert(rng != NULL);
  
  o = octahedron_create(DEGREE);

  ck_assert(o != NULL);
  ck_assert(manifold_valid(o));

  input_coeff = malloc(sizeof(double) * o->ntotaltriangles);
  ck_assert(input_coeff != NULL);

  recon_coeff = malloc(sizeof(double) * o->ntotaltriangles);
  ck_assert(recon_coeff != NULL);

  /*
   * We only really reconstruct the max level so only test these coefficients
   */
  sstart = 0;
  for (i = 0; i < DEGREE; i ++) {
    sstart += o->ntriangles[i];
  }
  send = sstart + o->ntriangles[DEGREE];
  
  for (i = 0; i < o->ntotaltriangles; i ++) {
    input_coeff[i] = gsl_rng_uniform(rng);
    recon_coeff[i] = input_coeff[i];
  }

  ck_assert(face_subdivision_forward_step(o, recon_coeff, DEGREE) >= 0);

  exact = -1;
  for (i = 0; i < o->ntotaltriangles; i ++) {
    /* Check for numerical errors */
    ck_assert(isnormal(recon_coeff[i]) || recon_coeff[i] == 0.0);
    
    if (recon_coeff[i] != input_coeff[i]) {
      exact = 0;
    }
  }
  /* Make sure some of the values have changed */
  ck_assert(!exact);

  ck_assert(face_subdivision_inverse_step(o, recon_coeff, DEGREE) >= 0);
  for (i = sstart; i < send; i ++) {
    /* Check for numerical errors */
    ck_assert(isnormal(recon_coeff[i]) || recon_coeff[i] == 0.0);
    ck_assert(fabs(recon_coeff[i] - input_coeff[i]) < 1.0e-9);
  }
  

  free(input_coeff);
  free(recon_coeff);
  gsl_rng_free(rng);
  manifold_destroy(o);
}
END_TEST

START_TEST(test_octahedron_full)
{
  static const int DEGREE = 7;
  
  gsl_rng *rng;
  
  manifold_t *o;
  double *input_coeff;
  double *recon_coeff;
  int exact;
  int sstart;
  int send;

  int i;
  
  rng = gsl_rng_alloc(gsl_rng_taus);
  ck_assert(rng != NULL);
  
  o = octahedron_create(DEGREE);

  ck_assert(o != NULL);
  ck_assert(manifold_valid(o));

  input_coeff = malloc(sizeof(double) * o->ntotaltriangles);
  ck_assert(input_coeff != NULL);

  recon_coeff = malloc(sizeof(double) * o->ntotaltriangles);
  ck_assert(recon_coeff != NULL);

  printf("Octahedron: %d vertices\n", o->ntotaltriangles);
  
  sstart = 0;
  for (i = 0; i < DEGREE; i ++) {
    sstart += o->ntriangles[i];
  }
  send = sstart + o->ntriangles[DEGREE];

  for (i = 0; i < o->ntotaltriangles; i ++) {
    input_coeff[i] = gsl_rng_uniform(rng);
    recon_coeff[i] = input_coeff[i];
  }

  ck_assert(face_subdivision_forward(o, recon_coeff) >= 0);

  exact = -1;
  for (i = 0; i < o->ntotaltriangles; i ++) {
    /* Check for numerical errors */
    ck_assert(isnormal(recon_coeff[i]) || recon_coeff[i] == 0.0);
    
    if (recon_coeff[i] != input_coeff[i]) {
      exact = 0;
    }
  }
  /* Make sure some of the values have changed */
  ck_assert(!exact);

  ck_assert(face_subdivision_inverse(o, recon_coeff) >= 0);
  for (i = sstart; i < send; i ++) {
    /* Check for numerical errors */
    ck_assert(isnormal(recon_coeff[i]) || recon_coeff[i] == 0.0);

    ck_assert(fabs(recon_coeff[i] - input_coeff[i]) < 1.0e-9);
  }
  

  free(input_coeff);
  free(recon_coeff);
  gsl_rng_free(rng);
  manifold_destroy(o);
}
END_TEST

START_TEST(test_octahedron_forward_constant)
{
  static const int DEGREE = 7;
  static const double CONST_VALUE = 3.1415;
  gsl_rng *rng;
  
  manifold_t *o;
  double *input_coeff;
  double *recon_coeff;
  int sstart;
  int send;

  int i;
  
  rng = gsl_rng_alloc(gsl_rng_taus);
  ck_assert(rng != NULL);
  
  o = octahedron_create(DEGREE);

  ck_assert(o != NULL);
  ck_assert(manifold_valid(o));

  input_coeff = malloc(sizeof(double) * o->ntotaltriangles);
  ck_assert(input_coeff != NULL);

  recon_coeff = malloc(sizeof(double) * o->ntotaltriangles);
  ck_assert(recon_coeff != NULL);

  printf("Octahedron: %d faces\n", o->ntotaltriangles);
  
  sstart = 0;
  for (i = 0; i < DEGREE; i ++) {
    sstart += o->ntriangles[i];
  }
  send = sstart + o->ntriangles[DEGREE];

  for (i = 0; i < o->ntotaltriangles; i ++) {
    input_coeff[i] = CONST_VALUE;
    recon_coeff[i] = input_coeff[i];
  }

  ck_assert(face_subdivision_forward(o, recon_coeff) >= 0);

  /*
   * For constant surface, first 8 faces  
   * should be the constant value and the rest should be 0.
   */
  for (i = 0; i < 8; i ++) {
    ck_assert(isnormal(recon_coeff[i]) || recon_coeff[i] == 0.0);
    ck_assert(fabs(recon_coeff[i] - CONST_VALUE) < 1.0e-9);
  }
  
  for (i = 8; i < o->ntotaltriangles; i ++) {
    /* Check for numerical errors */
    ck_assert(isnormal(recon_coeff[i]) || recon_coeff[i] == 0.0);
    ck_assert(fabs(recon_coeff[i]) < 1.0e-9);
  }

  ck_assert(face_subdivision_inverse(o, recon_coeff) >= 0);
  for (i = sstart; i < send; i ++) {
    /* Check for numerical errors */
    ck_assert(isnormal(recon_coeff[i]) || recon_coeff[i] == 0.0);

    ck_assert(fabs(recon_coeff[i] - input_coeff[i]) < 1.0e-9);
  }
  

  free(input_coeff);
  free(recon_coeff);
  gsl_rng_free(rng);
  manifold_destroy(o);
}
END_TEST

#if 0
START_TEST(test_octahedron_shell_full)
{
  static const int DEGREE = 5;
  
  gsl_rng *rng;
  
  manifold_t *o;
  double *input_coeff;
  double *recon_coeff;
  int exact;
  double *workspace;
  int radial_size;
  int total_size;

  int i;
  
  rng = gsl_rng_alloc(gsl_rng_taus);
  ck_assert(rng != NULL);
  
  o = octahedron_create(DEGREE);

  ck_assert(o != NULL);
  ck_assert(manifold_valid(o));

  radial_size = 1 << DEGREE;
  workspace = malloc(sizeof(double) * radial_size);
  ck_assert(workspace != NULL);

  total_size = o->nvertices * radial_size;
  
  input_coeff = malloc(sizeof(double) * total_size);
  ck_assert(input_coeff != NULL);

  recon_coeff = malloc(sizeof(double) * total_size);
  ck_assert(recon_coeff != NULL);

  printf("Octahedron Shell Degree %d: %d vertices\n", DEGREE, total_size);
  
  for (i = 0; i < total_size; i ++) {
    input_coeff[i] = gsl_rng_uniform(rng);
    recon_coeff[i] = input_coeff[i];
  }

  ck_assert(face_subdivision_shell_forward(o, recon_coeff, total_size, workspace) >= 0);

  exact = -1;
  for (i = 0; i < total_size; i ++) {
    /* Check for numerical errors */
    ck_assert(isnormal(recon_coeff[i]) || recon_coeff[i] == 0.0);
    
    if (recon_coeff[i] != input_coeff[i]) {
      exact = 0;
    }
  }
  /* Make sure some of the values have changed */
  ck_assert(!exact);

  ck_assert(face_subdivision_shell_inverse(o, recon_coeff, total_size, workspace) >= 0);
  for (i = 0; i < total_size; i ++) {
    /* Check for numerical errors */
    ck_assert(isnormal(recon_coeff[i]) || recon_coeff[i] == 0.0);

    ck_assert(fabs(recon_coeff[i] - input_coeff[i]) < 1.0e-9);
  }
  
  free(workspace);
  free(input_coeff);
  free(recon_coeff);
  gsl_rng_free(rng);
  manifold_destroy(o);
}
END_TEST
#endif /* 0 */

START_TEST(test_icosahedron_single_step)
{
  static const int DEGREE = 2;
  gsl_rng *rng;
  
  manifold_t *o;
  double *input_coeff;
  double *recon_coeff;
  int exact;
  int sstart;
  int send;

  int i;
  
  rng = gsl_rng_alloc(gsl_rng_taus);
  ck_assert(rng != NULL);
  
  o = icosahedron_create(DEGREE);

  ck_assert(o != NULL);
  ck_assert(manifold_valid(o));

  input_coeff = malloc(sizeof(double) * o->ntotaltriangles);
  ck_assert(input_coeff != NULL);

  recon_coeff = malloc(sizeof(double) * o->ntotaltriangles);
  ck_assert(recon_coeff != NULL);

  for (i = 0; i < o->ntotaltriangles; i ++) {
    input_coeff[i] = gsl_rng_uniform(rng);
    recon_coeff[i] = input_coeff[i];
  }

  sstart = 0;
  for (i = 0; i < DEGREE; i ++) {
    sstart += o->ntriangles[i];
  }
  send = sstart + o->ntriangles[DEGREE];

  ck_assert(face_subdivision_forward_step(o, recon_coeff, DEGREE) >= 0);

  exact = -1;
  for (i = 0; i < o->ntotaltriangles; i ++) {
    /* Check for numerical errors */
    ck_assert(isnormal(recon_coeff[i]) || recon_coeff[i] == 0.0);
    
    if (recon_coeff[i] != input_coeff[i]) {
      exact = 0;
    }
  }
  /* Make sure some of the values have changed */
  ck_assert(!exact);

  ck_assert(face_subdivision_inverse_step(o, recon_coeff, DEGREE) >= 0);
  for (i = sstart; i < send; i ++) {
    /* Check for numerical errors */
    ck_assert(isnormal(recon_coeff[i]) || recon_coeff[i] == 0.0);

    ck_assert(fabs(recon_coeff[i] - input_coeff[i]) < 1.0e-9);
  }
  

  free(input_coeff);
  free(recon_coeff);
  gsl_rng_free(rng);
  manifold_destroy(o);
}
END_TEST

START_TEST(test_icosahedron_full)
{
  static const int DEGREE = 6;
  
  gsl_rng *rng;
  
  manifold_t *o;
  double *input_coeff;
  double *recon_coeff;
  int exact;

  int sstart;
  int send;
  
  int i;
  
  rng = gsl_rng_alloc(gsl_rng_taus);
  ck_assert(rng != NULL);
  
  o = icosahedron_create(DEGREE);

  ck_assert(o != NULL);
  ck_assert(manifold_valid(o));

  input_coeff = malloc(sizeof(double) * o->ntotaltriangles);
  ck_assert(input_coeff != NULL);

  recon_coeff = malloc(sizeof(double) * o->ntotaltriangles);
  ck_assert(recon_coeff != NULL);

  printf("Icosahedron: %d vertices\n", o->ntotaltriangles);
  
  for (i = 0; i < o->ntotaltriangles; i ++) {
    input_coeff[i] = gsl_rng_uniform(rng);
    recon_coeff[i] = input_coeff[i];
  }

  ck_assert(face_subdivision_forward(o, recon_coeff) >= 0);

  sstart = 0;
  for (i = 0; i < DEGREE; i ++) {
    sstart += o->ntriangles[i];
  }
  send = sstart + o->ntriangles[DEGREE];

  exact = -1;
  for (i = 0; i < o->ntotaltriangles; i ++) {
    /* Check for numerical errors */
    ck_assert(isnormal(recon_coeff[i]) || recon_coeff[i] == 0.0);
    
    if (recon_coeff[i] != input_coeff[i]) {
      exact = 0;
    }
  }
  /* Make sure some of the values have changed */
  ck_assert(!exact);

  ck_assert(face_subdivision_inverse(o, recon_coeff) >= 0);
  for (i = sstart; i < send; i ++) {
    /* Check for numerical errors */
    ck_assert(isnormal(recon_coeff[i]) || recon_coeff[i] == 0.0);

    ck_assert(fabs(recon_coeff[i] - input_coeff[i]) < 1.0e-9);
  }
  

  free(input_coeff);
  free(recon_coeff);
  gsl_rng_free(rng);
  manifold_destroy(o);
}
END_TEST

#if 0
START_TEST(test_icosahedron_shell_full)
{
  static const int DEGREE = 5;
  
  gsl_rng *rng;
  
  manifold_t *o;
  double *input_coeff;
  double *recon_coeff;
  int exact;
  double *workspace;
  int radial_size;
  int total_size;

  int i;
  
  rng = gsl_rng_alloc(gsl_rng_taus);
  ck_assert(rng != NULL);
  
  o = icosahedron_create(DEGREE);

  ck_assert(o != NULL);
  ck_assert(manifold_valid(o));

  radial_size = 1 << DEGREE;
  workspace = malloc(sizeof(double) * radial_size);
  ck_assert(workspace != NULL);

  total_size = o->ntotaltriangles * radial_size;
  
  input_coeff = malloc(sizeof(double) * total_size);
  ck_assert(input_coeff != NULL);

  recon_coeff = malloc(sizeof(double) * total_size);
  ck_assert(recon_coeff != NULL);

  printf("Icosahedron Shell Degree %d: %d vertices\n", DEGREE, total_size);
  
  for (i = 0; i < total_size; i ++) {
    input_coeff[i] = gsl_rng_uniform(rng);
    recon_coeff[i] = input_coeff[i];
  }

  ck_assert(face_subdivision_shell_forward(o, recon_coeff, total_size, workspace) >= 0);

  exact = -1;
  for (i = 0; i < total_size; i ++) {
    /* Check for numerical errors */
    ck_assert(isnormal(recon_coeff[i]) || recon_coeff[i] == 0.0);
    
    if (recon_coeff[i] != input_coeff[i]) {
      exact = 0;
    }
  }
  /* Make sure some of the values have changed */
  ck_assert(!exact);

  ck_assert(face_subdivision_shell_inverse(o, recon_coeff, total_size, workspace) >= 0);
  for (i = 0; i < total_size; i ++) {
    /* Check for numerical errors */
    ck_assert(isnormal(recon_coeff[i]) || recon_coeff[i] == 0.0);

    ck_assert(fabs(recon_coeff[i] - input_coeff[i]) < 1.0e-9);
  }
  
  free(workspace);
  free(input_coeff);
  free(recon_coeff);
  gsl_rng_free(rng);
  manifold_destroy(o);
}
END_TEST
#endif


Suite *
face_subdivision_suite (void)
{
  Suite *s = suite_create ("Face Subdivision");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");

  tcase_add_test (tc_core, test_octahedron_single_step);
  tcase_add_test (tc_core, test_octahedron_full);
  tcase_add_test (tc_core, test_octahedron_forward_constant);

  tcase_add_test (tc_core, test_icosahedron_single_step);
  tcase_add_test (tc_core, test_icosahedron_full);

  /* tcase_add_test (tc_core, test_octahedron_shell_full); */
  /* tcase_add_test (tc_core, test_icosahedron_shell_full); */
  
  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = face_subdivision_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
