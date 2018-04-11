
#include <stdio.h>
#include <stdlib.h>
#include <check.h>

#include "manifold.h"
#include "octahedron.h"

START_TEST(test_octahedron_counting)
{
  ck_assert(octahedron_nvertices(0) == 6);
  ck_assert(octahedron_nedges(0) == 12);
  ck_assert(octahedron_ntriangles(0) == 8);

  ck_assert(octahedron_nvertices(1) == 18);
  ck_assert(octahedron_nedges(1) == 48);
  ck_assert(octahedron_ntriangles(1) == 32);

  ck_assert(octahedron_nvertices(2) == 66);
  ck_assert(octahedron_nedges(2) == 192);
  ck_assert(octahedron_ntriangles(2) == 128);

  ck_assert(octahedron_nvertices(3) == 258);
  ck_assert(octahedron_nedges(3) == 768);
  ck_assert(octahedron_ntriangles(3) == 512);
}
END_TEST

START_TEST(test_octahedron_create)
{
  manifold_t *o;

  o = octahedron_create(2);

  ck_assert(o != NULL);

  ck_assert(manifold_valid(o));

  manifold_destroy(o);
}
END_TEST

Suite *
octahedron_suite (void)
{
  Suite *s = suite_create ("Octahedron");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  /* tcase_add_test (tc_core, test_octahedron_create); */
  /* tcase_add_test (tc_core, test_octahedron_subdivide1); */
  /* tcase_add_test (tc_core, test_octahedron_subdivide2); */

  tcase_add_test(tc_core, test_octahedron_counting);
  tcase_add_test(tc_core, test_octahedron_create);

  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = octahedron_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
