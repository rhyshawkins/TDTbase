
#include <stdio.h>
#include <stdlib.h>
#include <check.h>
#include <math.h>

#include "wavetree_prior.h"

START_TEST(test_globally_uniform)
{
  wavetree_prior_t *prior;
  int i;
  double coeff;
  double mean;
  

  prior = wavetree_prior_create_globally_uniform(-1.0, 1.0,
						 1234);
  ck_assert(prior != NULL);

  /*
   * Check generation works
   */
  mean = 0.0;
  for (i = 0; i < 1000; i ++) {
    ck_assert(prior->sample(prior->user,
			    1, 2, 3,
			    5, 
			    6,
			    0.0,
			    &coeff) >= 0);

    ck_assert(coeff >= -1.0);
    ck_assert(coeff <= 1.0);

    mean = mean + coeff;
    
    ck_assert(prior->prob(prior->user,
			  1, 2, 3,
			  5, 
			  6,
			  0.0,
			  coeff) == 0.5);

    ck_assert(prior->valid(prior->user,
			   1, 2, 3,
			   5, 
			   6,
			   0.0,
			   coeff));
  }

  /*
   * Ensure that the mean is near 0
   */
  mean /= (double)1000;
  ck_assert(fabs(mean) < 0.1);

  /*
   * Check values outside of prior return invalid an 0 prob
   */
  ck_assert(prior->prob(prior->user,
			1, 2, 3,
			5, 
			6,
			0.0,
			2.0) == 0.0);
  ck_assert(prior->valid(prior->user,
			 1, 2, 3,
			 5, 
			 6,
			 0.0,
			 2.0) == 0);

  ck_assert(prior->prob(prior->user,
			1, 2, 3,
			5, 
			6,
			0.0,
			-2.0) == 0.0);
  ck_assert(prior->valid(prior->user,
			 1, 2, 3,
			 5, 
			 6,
			 0.0,
			 -2.0) == 0);

  wavetree_prior_free_globally_uniform(prior);
}
END_TEST

Suite *
prior_suite (void)
{
  Suite *s = suite_create ("Prior Tests");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_globally_uniform);
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
