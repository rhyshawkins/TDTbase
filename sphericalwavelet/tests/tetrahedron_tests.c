
#include <stdio.h>
#include <stdlib.h>
#include <check.h>

#include "tetrahedron.h"

START_TEST(test_tetrahedron_counting)
{
  ck_assert(tetrahedron_nvertices(0) == 4);
  ck_assert(tetrahedron_nedges(0) == 6);
  ck_assert(tetrahedron_ntriangles(0) == 4);

  ck_assert(tetrahedron_nvertices(1) == 10);
  ck_assert(tetrahedron_nedges(1) == 24);
  ck_assert(tetrahedron_ntriangles(1) == 16);

  ck_assert(tetrahedron_nvertices(2) == 34);
  ck_assert(tetrahedron_nedges(2) == 96);
  ck_assert(tetrahedron_ntriangles(2) == 64);

  ck_assert(tetrahedron_nvertices(3) == 130);
  ck_assert(tetrahedron_nedges(3) == 384);
  ck_assert(tetrahedron_ntriangles(3) == 256);
}
END_TEST

START_TEST(test_tetrahedron_create)
{
  manifold_t *o;

  o = tetrahedron_create(2);

  ck_assert(o != NULL);

  manifold_destroy(o);
}
END_TEST


Suite *
tetrahedron_suite (void)
{
  Suite *s = suite_create ("Tetrahedron");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");

  tcase_add_test (tc_core, test_tetrahedron_counting);
  tcase_add_test (tc_core, test_tetrahedron_create);

  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = tetrahedron_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
