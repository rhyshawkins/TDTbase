

#include <stdio.h>
#include <stdlib.h>
#include <check.h>
#include <math.h>

#include "subdivisiontree2d.h"

START_TEST (test_subdivisiontree2d_toplevel_0overlap)
{
  #define IM_DEGREE 3
  #define IM_W 8
  #define IM_SIZE 64

  subdivisiontree2d_t *s;
  int i;
  int j;

  double img[IM_W*IM_W];

  double x;
  double y;
  double dx;

  double xp;
  double yp;
  double dxp;

  double v;
  double expected;


  s = subdivisiontree2d_create(IM_DEGREE, 0.0, 0.0, SUBDIVISION_BASIS_PYRAMID);
  ck_assert(s != NULL);

  ck_assert(subdivisiontree2d_get_size(s) == (IM_SIZE));

  /*
   * Initialization
   */
  ck_assert(subdivisiontree2d_initialize(s, 1.0) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 1);
  
  /*
   * Map root node to check top level pyramid
   */
  ck_assert(subdivisiontree2d_map_to_array(s, img, IM_SIZE) >= 0);

  dx = 2.0/((double)IM_W);
  for (j = 0; j < (IM_W); j ++) {
    for (i = 0; i < (IM_W); i ++) {
      
      x = -1.0 + dx/2 + (double)i*dx;
      y = -1.0 + dx/2 + (double)j*dx;
      
      ck_assert(img[j*IM_W + i] == basis_pyramid(x, y, 1.0, 0.0));
    }
  }

  /*
   * Add at index 1 and commit
   */
  ck_assert(subdivisiontree2d_propose_birth(s, 1, 1, 1.0) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 2);
  ck_assert(subdivisiontree2d_commit(s) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 2);
  ck_assert(subdivisiontree2d_coeff_count_scan(s) == 2);

  /*
   * Check Top-Left quadrant is a pyramid
   */
  ck_assert(subdivisiontree2d_map_to_array(s, img, IM_SIZE) >= 0);
  
  dx = 2.0/((double)IM_W/2);
  dxp = 2.0/((double)IM_W);
  
  for (j = 0; j < (IM_W/2); j ++) {
    for (i = 0; i < (IM_W/2); i ++) {

      x = -1.0 + dx/2 + (double)i*dx;
      y = -1.0 + dx/2 + (double)j*dx;

      xp = -1.0 + dxp/2 + (double)i*dxp;
      yp = -1.0 + dxp/2 + (double)j*dxp;

      v = img[j*IM_W + i];
      expected = basis_pyramid(xp, yp, 1.0, 0.0) + 
	basis_pyramid(x, y, 1.0, 0.0);
	/* ((1.0 - fabs(xp))*(1.0 - fabs(yp)) +  */
	/*  (1.0 - fabs(x))*(1.0 - fabs(y))); */
      
      ck_assert(v == expected);

      xp = -1.0 + dxp/2 + (double)(i + IM_W/2)*dxp;
      yp = -1.0 + dxp/2 + (double)j*dxp;
      expected = basis_pyramid(xp, yp, 1.0, 0.0);
      ck_assert(img[j*IM_W + i + IM_W/2] == expected);
      
      xp = -1.0 + dxp/2 + (double)i*dxp;
      yp = -1.0 + dxp/2 + (double)(j + IM_W/2)*dxp;
      expected = basis_pyramid(xp, yp, 1.0, 0.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i ] == expected);

      xp = -1.0 + dxp/2 + (double)(i + IM_W/2)*dxp;
      yp = -1.0 + dxp/2 + (double)(j + IM_W/2)*dxp;
      expected = basis_pyramid(xp, yp, 1.0, 0.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i + IM_W/2] == expected);
    }
  }
}
END_TEST

START_TEST (test_subdivisiontree2d_toplevel_halfoverlap)
{
  #define IM_DEGREE 3
  #define IM_W 8
  #define IM_SIZE 64

  subdivisiontree2d_t *s;
  int i;
  int j;

  double img[IM_W*IM_W];

  double x;
  double y;
  double dx;
  
  double xp;
  double yp;
  double dxp;


  s = subdivisiontree2d_create(IM_DEGREE, 0.0, 0.5, SUBDIVISION_BASIS_PYRAMID);
  ck_assert(s != NULL);

  ck_assert(subdivisiontree2d_get_size(s) == (IM_SIZE));

  /*
   * Initialization
   */
  ck_assert(subdivisiontree2d_initialize(s, 1.0) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 1);

  /*
   * Map root node and check base pyramid
   */
  ck_assert(subdivisiontree2d_map_to_array(s, img, IM_SIZE) >= 0);

  dx = 2.0/((double)IM_W);
  for (j = 0; j < (IM_W); j ++) {
    for (i = 0; i < (IM_W); i ++) {
      
      x = -1.0 + dx/2 + (double)i*dx;
      y = -1.0 + dx/2 + (double)j*dx;
      
      ck_assert(img[j*IM_W + i] == basis_pyramid(x, y, 1.0, 0.5));
    }
  }

  /*
   * Add at index 1 and commit
   */
  ck_assert(subdivisiontree2d_propose_birth(s, 1, 1, 1.0) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 2);
  ck_assert(subdivisiontree2d_commit(s) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 2);
  ck_assert(subdivisiontree2d_coeff_count_scan(s) == 2);

  /*
   * Check Top-Left quadrant is a pyramid
   */
  ck_assert(subdivisiontree2d_map_to_array(s, img, IM_SIZE) >= 0);
  
  dx = 2.0/((double)IM_W/2);
  dxp = 2.0/((double)IM_W);
  
  for (j = 0; j < (IM_W/2); j ++) {
    for (i = 0; i < (IM_W/2); i ++) {

      x = -1.0 + dx/2 + (double)i*dx;
      y = -1.0 + dx/2 + (double)j*dx;

      xp = -1.0 + dxp/2 + (double)i*dxp;
      yp = -1.0 + dxp/2 + (double)j*dxp;

      if (fabs(x) <= 2.0 && fabs(y) <= 2.0) {
        ck_assert(img[j*IM_W + i] == (basis_pyramid(x, y, 1.0, 0.5) + 
				      basis_pyramid(xp, yp, 1.0, 0.5)));
      } else {
        ck_assert(img[j*IM_W + i] == basis_pyramid(xp, yp, 1.0, 0.5));
      }
    }
  }
}
END_TEST

START_TEST (test_subdivisiontree2d_toplevel_1overlap)
{
  #define IM_DEGREE 3
  #define IM_W 8
  #define IM_SIZE 64

  subdivisiontree2d_t *s;
  int i;
  int j;

  double img[IM_W*IM_W];

  double x;
  double y;
  double dx;

  double xp;
  double yp;
  double dxp;

  s = subdivisiontree2d_create(IM_DEGREE, 0.0, 1.0, SUBDIVISION_BASIS_PYRAMID);
  ck_assert(s != NULL);

  ck_assert(subdivisiontree2d_get_size(s) == (IM_SIZE));

  /*
   * Initialization
   */
  ck_assert(subdivisiontree2d_initialize(s, 1.0) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 1);

  /*
   * Map root node to get constant image
   */
  ck_assert(subdivisiontree2d_map_to_array(s, img, IM_SIZE) >= 0);

  dx = 2.0/((double)IM_W);
  for (j = 0; j < (IM_W); j ++) {
    for (i = 0; i < (IM_W); i ++) {
      
      x = -1.0 + dx/2 + (double)i*dx;
      y = -1.0 + dx/2 + (double)j*dx;
      
      ck_assert(img[j*IM_W + i] == basis_pyramid(x, y, 1.0, 1.0));
    }
  }

  /*
   * Add at index 1 and commit
   */
  ck_assert(subdivisiontree2d_propose_birth(s, 1, 1, 1.0) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 2);
  ck_assert(subdivisiontree2d_commit(s) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 2);
  ck_assert(subdivisiontree2d_coeff_count_scan(s) == 2);

  /*
   * Check Top-Left quadrant is a pyramid
   */
  ck_assert(subdivisiontree2d_map_to_array(s, img, IM_SIZE) >= 0);
  
  dx = 2.0/((double)IM_W/2);
  dxp = 2.0/((double)IM_W);
  
  for (j = 0; j < (IM_W/2); j ++) {
    for (i = 0; i < (IM_W/2); i ++) {

      x = -1.0 + dx/2 + (double)i*dx;
      y = -1.0 + dx/2 + (double)j*dx;

      xp = -1.0 + dxp/2 + (double)i*dxp;
      yp = -1.0 + dxp/2 + (double)j*dxp;

      if (fabs(x) <= 3.0 && fabs(y) <= 3.0) {
        ck_assert(img[j*IM_W + i] == (basis_pyramid(x, y, 1.0, 1.0) + 
				      basis_pyramid(xp, yp, 1.0, 1.0)));
      } else {
        ck_assert(img[j*IM_W + i] == basis_pyramid(xp, yp, 1.0, 1.0));
      }
    }
  }
}
END_TEST

Suite *
subdivisiontree2d_suite (void)
{
  Suite *s = suite_create ("Subdivision Tree 2D: Pyramid Basis");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_subdivisiontree2d_toplevel_0overlap);
  tcase_add_test (tc_core, test_subdivisiontree2d_toplevel_halfoverlap);
  tcase_add_test (tc_core, test_subdivisiontree2d_toplevel_1overlap);
  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = subdivisiontree2d_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
