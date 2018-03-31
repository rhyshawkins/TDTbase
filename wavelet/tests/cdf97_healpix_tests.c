
#include <stdio.h>
#include <stdlib.h>
#include <check.h>
#include <math.h>

#include "cdf97_healpix.h"

static inline double within1pc(double a, double b)
{
  double delta;

  delta = fabs(a - b);
  if (delta < 1.0e-6) {
    return -1;
  }

  if (fabs(a) > fabs(b)) {
    if ((delta/fabs(a)) < 0.01) {
      return -1;
    } else {
      fprintf(stderr, "within1p: %f %f\n", a, b);
      return 0;
    }
  } else {
    if ((delta/fabs(b)) < 0.01) {
      return -1;
    } else {
      fprintf(stderr, "within1p: %f %f\n", a, b);
      return 0;
    }
  }
}

START_TEST (test_cdf97_healpix_traverse)
{
  const int WIDTH = 8;
  const int ROW = 2;

  cdf97_healpix_t *c;

  int next_tile;
  int next_type;
  int next_index;
  int next_dir;

  c = cdf97_healpix_create(WIDTH);
  ck_assert(c != NULL);

  /*
   * From tile 0 tests
   */
  ck_assert(cdf97_healpix_traverse_row(c, 
				       WIDTH,
				       0,
				       ROW,
				       -1, 
				       &next_tile,
				       &next_type,
				       &next_index,
				       &next_dir) == 0);
  ck_assert(next_tile == 4);
  ck_assert(next_type == 0);
  ck_assert(next_index == ROW);
  ck_assert(next_dir == -1);
	    
  ck_assert(cdf97_healpix_traverse_row(c, 
				       WIDTH,
				       0,
				       ROW,
				       1, 
				       &next_tile,
				       &next_type,
				       &next_index,
				       &next_dir) == 0);
  ck_assert(next_tile == 1);
  ck_assert(next_type == 1);
  ck_assert(next_index == WIDTH - 1 - ROW);
  ck_assert(next_dir == 1);

  ck_assert(cdf97_healpix_traverse_col(c, 
				       WIDTH,
				       0,
				       ROW,
				       -1, 
				       &next_tile,
				       &next_type,
				       &next_index,
				       &next_dir) == 0);
  ck_assert(next_tile == 3);
  ck_assert(next_type == 0);
  ck_assert(next_index == WIDTH - 1 - ROW);
  ck_assert(next_dir == -1);
	    
  ck_assert(cdf97_healpix_traverse_col(c, 
				       WIDTH,
				       0,
				       ROW,
				       1, 
				       &next_tile,
				       &next_type,
				       &next_index,
				       &next_dir) == 0);
  ck_assert(next_tile == 5);
  ck_assert(next_type == 1);
  ck_assert(next_index == ROW);
  ck_assert(next_dir == 1);
  
  cdf97_healpix_destroy(c);
}
END_TEST

START_TEST (test_cdf97_healpix_wrap)
{
}
END_TEST

 
Suite *
cdf97_suite (void)
{
  Suite *s = suite_create ("CDF97");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_cdf97_healpix_traverse);
  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = cdf97_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
