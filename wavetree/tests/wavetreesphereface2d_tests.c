
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <check.h>

#include "manifold.h"
#include "icosahedron.h"
#include "wavetreesphereface2d.h"

START_TEST (test_wavetreesphereface2d_create)
{
  static const int DEGREE = 1;
  wavetreesphereface2d_t *s;
  manifold_t *manifold;

  ck_assert(DEGREE > 0);
  
  manifold = icosahedron_create(DEGREE - 1);
  ck_assert(manifold != NULL);
  
  s = wavetreesphereface2d_create(manifold, 
				  0.0);
  ck_assert(s != NULL);

  wavetreesphereface2d_destroy(s);
  manifold_destroy(manifold);
}
END_TEST

START_TEST (test_wavetreesphereface2d_depth_base)
{
  static const int DEGREE = 5;
  wavetreesphereface2d_t *s;
  manifold_t *manifold;

  ck_assert(DEGREE > 0);
  
  manifold = icosahedron_create(DEGREE - 1);
  ck_assert(manifold != NULL);
  
  s = wavetreesphereface2d_create(manifold, 
				  0.0);
  ck_assert(s != NULL);

  ck_assert(manifold->ntrianglesatdepth(0) == 20);

  ck_assert(wavetreesphereface2d_depth_base(s, 0) == 0);
  ck_assert(wavetreesphereface2d_depth_base(s, 1) == 1);
  ck_assert(wavetreesphereface2d_depth_base(s, 2) == 21);
  ck_assert(wavetreesphereface2d_depth_base(s, 3) == 101);
  ck_assert(wavetreesphereface2d_depth_base(s, 4) == 421);
  ck_assert(wavetreesphereface2d_depth_base(s, 5) == 1701);
  ck_assert(wavetreesphereface2d_depth_base(s, 6) == 6821);
  

  wavetreesphereface2d_destroy(s);
  manifold_destroy(manifold);
}
END_TEST

START_TEST (test_wavetreesphereface2d_child_indices)
{
  static const int DEGREE = 5;
  wavetreesphereface2d_t *s;
  manifold_t *manifold;
  int indices[20];
  int i;
  int n;

  ck_assert(DEGREE > 0);
  
  manifold = icosahedron_create(DEGREE - 1);
  ck_assert(manifold != NULL);
  
  s = wavetreesphereface2d_create(manifold, 
				  0.0);
  ck_assert(s != NULL);

  ck_assert(wavetreesphereface2d_child_indices(s, 0, 0, indices, &n, 20) >= 0);
  ck_assert(n == 20);

  for (i = 0; i < n; i ++) {

    ck_assert(wavetreesphereface2d_parent_index(s, indices[i]) == 0);

  }

  ck_assert(wavetreesphereface2d_child_indices(s, 1, 1, indices, &n, 20) >= 0);
  ck_assert(n == 3);

  for (i = 0; i < n; i ++) {

    ck_assert(wavetreesphereface2d_parent_index(s, indices[i]) == 1);

  }
  
  ck_assert(wavetreesphereface2d_child_indices(s, 21, 2, indices, &n, 20) >= 0);
  ck_assert(n == 3);
  
  for (i = 0; i < n; i ++) {

    ck_assert(wavetreesphereface2d_parent_index(s, indices[i]) == 21);

  }

  wavetreesphereface2d_destroy(s);
  manifold_destroy(manifold);
}
END_TEST

Suite *
wavetreesphereface2d_suite (void)
{
  Suite *s = suite_create ("Wave Tree Face Wavelet 2D");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");

  tcase_add_test (tc_core, test_wavetreesphereface2d_create);

  tcase_add_test (tc_core, test_wavetreesphereface2d_depth_base);

  tcase_add_test (tc_core, test_wavetreesphereface2d_child_indices);

  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = wavetreesphereface2d_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
