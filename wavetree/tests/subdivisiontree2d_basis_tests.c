

#include <stdio.h>
#include <stdlib.h>
#include <check.h>
#include <math.h>

#include "subdivisiontree2d.h"

START_TEST (test_subdivisiontree2d_zero)
{
  ck_assert(basis_constant(0.0, 0.0, 1.0, 0.0) == 1.0);
}
END_TEST

START_TEST (test_subdivisiontree2d_pyramid)
{
  /*
   * With overlap 0
   */

  /* Peak should be 1 */
  ck_assert(basis_pyramid(0.0, 0.0, 1.0, 0.0) == 1.0);

  /* Left should be 0 */
  ck_assert(basis_pyramid(-1.0, 0.0, 1.0, 0.0) == 0.0);

  /* Right should be 0 */
  ck_assert(basis_pyramid(1.0, 0.0, 1.0, 0.0) == 0.0);

  /* Mid left should be half */
  ck_assert(basis_pyramid(-0.5, 0.0, 1.0, 0.0) == 0.5);

  /* Mid Right should be half */
  ck_assert(basis_pyramid(0.5, 0.0, 1.0, 0.0) == 0.5);


  /*
   * With overlap 0.5
   */

  /* Peak should be 1 */
  ck_assert(basis_pyramid(0.0, 0.0, 1.0, 0.5) == 1.0);

  /* Left should be 0 */
  ck_assert(basis_pyramid(-2.0, 0.0, 1.0, 0.5) == 0.0);

  /* Right should be 0 */
  ck_assert(basis_pyramid(2.0, 0.0, 1.0, 0.5) == 0.0);

  /* Mid left should be half */
  ck_assert(basis_pyramid(-1.0, 0.0, 1.0, 0.5) == 0.5);

  /* Mid Right should be half */
  ck_assert(basis_pyramid(1.0, 0.0, 1.0, 0.5) == 0.5);

  /*
   * With overlap 1.0
   */

  /* Peak should be 1 */
  ck_assert(basis_pyramid(0.0, 0.0, 1.0, 1.0) == 1.0);

  /* Left should be 0 */
  ck_assert(basis_pyramid(-3.0, 0.0, 1.0, 1.0) == 0.0);

  /* Right should be 0 */
  ck_assert(basis_pyramid(3.0, 0.0, 1.0, 1.0) == 0.0);

  /* Mid left should be half */
  ck_assert(basis_pyramid(-1.5, 0.0, 1.0, 1.0) == 0.5);

  /* Mid Right should be half */
  ck_assert(basis_pyramid(1.5, 0.0, 1.0, 1.0) == 0.5);

}
END_TEST

START_TEST (test_subdivisiontree2d_lanczos)
{
}
END_TEST

Suite *
subdivisiontree2d_suite (void)
{
  Suite *s = suite_create ("Subdivision Tree 2D: Pyramid Basis");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_subdivisiontree2d_zero);
  tcase_add_test (tc_core, test_subdivisiontree2d_pyramid);
  tcase_add_test (tc_core, test_subdivisiontree2d_lanczos);
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
