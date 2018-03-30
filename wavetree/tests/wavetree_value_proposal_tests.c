
#include <stdio.h>
#include <stdlib.h>
#include <check.h>
#include <math.h>

#include "wavetree_value_proposal.h"

START_TEST(test_global_gaussian)
{
  wavetree_prior_t *prior;
  wavetree_value_t *value;

  int i;
  double coeff;
  double prob;
  int error_count;

  prior = wavetree_prior_create_globally_uniform(-1.0, 1.0, 1234);
  ck_assert(prior != NULL);

  value = wavetree_value_create_global_gaussian(1.0,
						1234);
  ck_assert(value != NULL);

  ck_assert(value->init(value->user) == 0);
  
  for (i = 0; i < 100; i ++) {

    coeff = 0.0;
    ck_assert(value->perturb(prior,
			     value->user,
			     1, 2, 3,
			     8,
			     8,
			     0.0,
			     &coeff,
			     &prob) == 0);

  }

  ck_assert(value->errors(value->user, &error_count) == 0);

  /*
   * Error count should be greater than 0 as std dev is same as prior width.
   */
  ck_assert(error_count > 0); 

  wavetree_value_destroy(value);

  wavetree_prior_destroy(prior);
}
END_TEST

START_TEST(test_depth_gaussian)
{
  wavetree_prior_t *prior;
  wavetree_value_t *value;
  
  #define NDEPTHS 8
  double std[NDEPTHS] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  int i;
  double coeff;
  double prob;
  int error_count;

  prior = wavetree_prior_create_globally_uniform(-1.0, 1.0, 1234);
  ck_assert(prior != NULL);

  value = wavetree_value_create_depth_gaussian(NDEPTHS,
					       std,
					       1234);
  ck_assert(value != NULL);

  ck_assert(value->init(value->user) == 0);
  
  for (i = 0; i < 100; i ++) {

    coeff = 0.0;
    ck_assert(value->perturb(prior,
			     value->user,
			     1, 2, 3,
			     i % NDEPTHS,
			     8,
			     0.0,
			     &coeff,
			     &prob) == 0);

  }

  ck_assert(value->errors(value->user, &error_count) == 0);

  /*
   * Error count should be greater than 0 as std dev is same as prior width.
   */
  ck_assert(error_count > 0); 

  wavetree_value_destroy(value);

  wavetree_prior_destroy(prior);
}
END_TEST

Suite *
prior_suite (void)
{
  Suite *s = suite_create ("Value Proposal Tests");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_global_gaussian);
  tcase_add_test (tc_core, test_depth_gaussian);
  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = prior_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
