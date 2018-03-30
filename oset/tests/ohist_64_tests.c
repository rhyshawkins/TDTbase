
#include <stdio.h>
#include <stdlib.h>
#include <check.h>

#include "ohist_64.h"

START_TEST (test_ohist_64_create)
{
  ohist_64_t *set;

  set = ohist_64_create();
  ck_assert_ptr_ne(set, NULL);

  ohist_64_destroy(set);
}
END_TEST

START_TEST(test_ohist_64_insert)
{
  uint64_t sequence[] = {123, /* x2 */
			 54,  /* x3 */
			 232, /* x2 */
			 97,  /* x1 */
			 103, /* x1 */
			 54,
			 232,
			 54,
			 123};
  int i;
  int k;

  uint64_t idx;
  int dep;

  ohist_64_t *set;

  set = ohist_64_create();
  ck_assert_ptr_ne(set, NULL);

  k = 1;

  for (i = 0; i < sizeof(sequence)/sizeof(uint64_t); i ++) {

    ck_assert(ohist_64_insert(set, sequence[i], k, 1) >= 0);

  }

  ck_assert(ohist_64_nelements(set, k) == 5);

  ck_assert(ohist_64_nth_element(set, k, 0, &idx, &dep) == 0);
  ck_assert(idx == 54);
  ck_assert(dep == 3);

  ck_assert(ohist_64_nth_element(set, k, 1, &idx, &dep) == 0);
  ck_assert(idx == 97);
  ck_assert(dep == 1);

  ck_assert(ohist_64_nth_element(set, k, 2, &idx, &dep) == 0);
  ck_assert(idx == 103);
  ck_assert(dep == 1);

  ck_assert(ohist_64_nth_element(set, k, 3, &idx, &dep) == 0);
  ck_assert(idx == 123);
  ck_assert(dep == 2);

  ck_assert(ohist_64_nth_element(set, k, 4, &idx, &dep) == 0);
  ck_assert(idx == 232);
  ck_assert(dep == 2);

  ck_assert(ohist_64_nth_element(set, k, 5, &idx, &dep) == -1);

  ohist_64_destroy(set);
}
END_TEST

Suite *
ohist_64_suite (void)
{
  Suite *s = suite_create ("OHist 64");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_ohist_64_create);
  tcase_add_test (tc_core, test_ohist_64_insert);
  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = ohist_64_suite ();
  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
