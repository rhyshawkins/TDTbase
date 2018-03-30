
#include <stdio.h>
#include <stdlib.h>
#include <check.h>

#include "oset_int.h"

START_TEST (test_oset_int_create)
{
  oset_int_t *set;

  set = oset_int_create();
  ck_assert_ptr_ne(set, NULL);

  oset_int_destroy(set);
}
END_TEST

START_TEST(test_oset_int_insert)
{
  int indices[] = {9, 4, 7, 8, 1, 3, 6};
  int depths[] = {3, 1, 2, 2, 1, 1, 2};
  int i;

  int idx;
  int dep;

  oset_int_t *set;

  set = oset_int_create();
  ck_assert_ptr_ne(set, NULL);

  for (i = 0; i < sizeof(indices)/sizeof(int); i ++) {

    ck_assert(oset_int_insert(set, indices[i], depths[i]) > 0);

  }

  ck_assert(oset_int_count(set) == 7);

  ck_assert(oset_int_nth_element(set, 0, &idx, &dep) == 1);
  ck_assert(idx == 1);
  ck_assert(dep == 1);

  ck_assert(oset_int_nth_element(set, 1, &idx, &dep) == 3);
  ck_assert(idx == 3);
  ck_assert(dep == 1);

  ck_assert(oset_int_nth_element(set, 2, &idx, &dep) == 4);
  ck_assert(idx == 4);
  ck_assert(dep == 1);

  ck_assert(oset_int_nth_element(set, 3, &idx, &dep) == 6);
  ck_assert(idx == 6);
  ck_assert(dep == 2);

  ck_assert(oset_int_nth_element(set, 4, &idx, &dep) == 7);
  ck_assert(idx == 7);
  ck_assert(dep == 2);

  ck_assert(oset_int_nth_element(set, 5, &idx, &dep) == 8);
  ck_assert(idx == 8);
  ck_assert(dep == 2);

  ck_assert(oset_int_nth_element(set, 6, &idx, &dep) == 9);
  ck_assert(idx == 9);
  ck_assert(dep == 3);

  for (i = 0; i < sizeof(indices)/sizeof(int); i ++) {
    ck_assert(oset_int_is_element(set, indices[i]));
  }

  oset_int_destroy(set);
}
END_TEST


START_TEST(test_oset_int_remove)
{
  int indices[] = {9, 4, 7, 8, 1, 3, 6};
  int depths[] = {3, 1, 2, 2, 1, 1, 2};
  int i;

  oset_int_t *set;

  set = oset_int_create();
  ck_assert_ptr_ne(set, NULL);

  ck_assert(oset_int_remove(set, 8) == 0);

  ck_assert(oset_int_count(set) == 0);

  for (i = 0; i < sizeof(indices)/sizeof(int); i ++) {

    ck_assert(oset_int_insert(set, indices[i], depths[i]) > 0);

  }

  ck_assert(oset_int_remove(set, 5) == 0);

  ck_assert(oset_int_count(set) == 7);

  for (i = 0; i < sizeof(indices)/sizeof(int); i ++) {
    ck_assert(oset_int_is_element(set, indices[i]));
  }

  /*
   * Remove middle
   */
  ck_assert(oset_int_remove(set, 8) == 1);

  ck_assert(oset_int_count(set) == 6);
  
  ck_assert(oset_int_is_element(set, 8) == 0);

  for (i = 0; i < sizeof(indices)/sizeof(int); i ++) {
    if (indices[i] != 8) {
      ck_assert(oset_int_is_element(set, indices[i]));
    }
  }

  /* 2nd remove should return 0, ie nothing removed */
  ck_assert(oset_int_remove(set, 8) == 0); 

  ck_assert(oset_int_insert(set, 8, 2) > 0);


  /*
   * Remove first
   */
  ck_assert(oset_int_remove(set, 1) == 1);

  ck_assert(oset_int_count(set) == 6);
  
  ck_assert(oset_int_is_element(set, 1) == 0);

  for (i = 0; i < sizeof(indices)/sizeof(int); i ++) {
    if (indices[i] != 1) {
      ck_assert(oset_int_is_element(set, indices[i]));
    }
  }

  /* 2nd remove should return 0, ie nothing removed */
  ck_assert(oset_int_remove(set, 1) == 0); 

  ck_assert(oset_int_insert(set, 1, 1) > 0);

  /*
   * Remove last
   */
  ck_assert(oset_int_remove(set, 9) == 1);

  ck_assert(oset_int_count(set) == 6);
  
  ck_assert(oset_int_is_element(set, 9) == 0);

  for (i = 0; i < sizeof(indices)/sizeof(int); i ++) {
    if (indices[i] != 9) {
      ck_assert(oset_int_is_element(set, indices[i]));
    }
  }

  /* 2nd remove should return 0, ie nothing removed */
  ck_assert(oset_int_remove(set, 9) == 0); 

  oset_int_destroy(set);
}
END_TEST

START_TEST(test_oset_int_depth_restriction)
{
  int indices[] = {9, 4, 7, 8, 1, 3, 6};
  int depths[] = {3, 1, 2, 2, 1, 1, 2};
  int i;

  oset_int_t *set;

  set = oset_int_create();
  ck_assert_ptr_ne(set, NULL);

  ck_assert(oset_int_remove(set, 8) == 0);

  ck_assert(oset_int_count(set) == 0);

  for (i = 0; i < sizeof(indices)/sizeof(int); i ++) {

    ck_assert(oset_int_insert(set, indices[i], depths[i]) > 0);

  }

  ck_assert(oset_int_remove(set, 5) == 0);

  ck_assert(oset_int_count(set) == 7);

  for (i = 0; i < sizeof(indices)/sizeof(int); i ++) {
    ck_assert(oset_int_is_element(set, indices[i]));
  }

  ck_assert(oset_int_weighted_sum(set, 0.0, 3) == 7.0);
  ck_assert(oset_int_weighted_sum(set, 0.0, 2) == 6.0);
  ck_assert(oset_int_weighted_sum(set, 0.0, 1) == 3.0);

  ck_assert(oset_int_inverse_weighted_sum(set, 0.0, 3) == 7.0);
  ck_assert(oset_int_inverse_weighted_sum(set, 0.0, 2) == 6.0);
  ck_assert(oset_int_inverse_weighted_sum(set, 0.0, 1) == 3.0);
}
END_TEST

Suite *
oset_int_suite (void)
{
  Suite *s = suite_create ("Oset Int");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_oset_int_create);
  tcase_add_test (tc_core, test_oset_int_insert);
  tcase_add_test (tc_core, test_oset_int_remove);
  tcase_add_test (tc_core, test_oset_int_depth_restriction);
  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = oset_int_suite ();
  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
