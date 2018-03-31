
#include <stdio.h>
#include <stdlib.h>
#include <check.h>
#include <math.h>

#include "daubechies.h"

#define WIDTH 256

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

START_TEST (test_daub4_constant)
{
  double data[WIDTH * WIDTH];
  double work[WIDTH];
  int i;

  for (i = 0; i < (WIDTH * WIDTH); i ++) {
    data[i] = 3.14;
  }

  ck_assert(daubechies2d_forward_d4(data, 
				    WIDTH,
				    WIDTH,
				    WIDTH, 
				    work) >= 0);

  for (i = 1; i < (WIDTH * WIDTH); i ++) {
    ck_assert(isnormal(data[i]));
  }

  ck_assert(within1pc(data[0], 3.14));

  for (i = 1; i < (WIDTH * WIDTH); i ++) {
    ck_assert(within1pc(data[i], 0.0));
  }
  
  ck_assert(daubechies2d_inverse_d4(data, 
				    WIDTH,
				    WIDTH,
				    WIDTH, 
				    work) >= 0);

  for (i = 0; i < (WIDTH * WIDTH); i ++) {
    ck_assert(within1pc(data[i], 3.14));
  }
}
END_TEST

START_TEST (test_daub4_sinusoid)
{
  double data[WIDTH * WIDTH];
  double recon[WIDTH * WIDTH];
  double work[WIDTH];
  int i;
  int j;

  for (j = 0; j < WIDTH; j ++) {
    
    for (i = 0; i < WIDTH; i ++) {
      data[j * WIDTH + i] = sin((double)i/64.0 * M_PI * 2.0) * cos((double)j/8.0 * M_PI * 2.0);
      recon[j * WIDTH + i] = data[j * WIDTH + i];
    }

  }

  for (i = 0; i < (WIDTH * WIDTH); i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  ck_assert(daubechies2d_forward_d4(recon, 
				    WIDTH,
				    WIDTH,
				    WIDTH, 
				    work) >= 0);

  for (i = 0; i < (WIDTH * WIDTH); i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  ck_assert(daubechies2d_inverse_d4(recon, 
				    WIDTH,
				    WIDTH,
				    WIDTH, 
				    work) >= 0);

  for (i = 0; i < (WIDTH * WIDTH); i ++) {
    ck_assert(isnormal(recon[i]) || isfinite(recon[i]));
  }

  for (i = 0; i < (WIDTH * WIDTH); i ++) {
    ck_assert(within1pc(data[i], recon[i]));
  }
}
END_TEST
 
Suite *
daubechies_suite (void)
{
  Suite *s = suite_create ("Daubechies");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_daub4_1d_simple_step);
  
  tcase_add_test (tc_core, test_daub4_1d_constant_step);
  tcase_add_test (tc_core, test_daub4_1d_sinusoid_step);

  tcase_add_test (tc_core, test_daub4_constant);
  tcase_add_test (tc_core, test_daub4_sinusoid);

  tcase_add_test (tc_core, test_daub4_2d_constant);
  tcase_add_test (tc_core, test_daub4_2d_sinusoid);
  tcase_add_test (tc_core, test_daub4_2d_constant_nonsquare1);
  tcase_add_test (tc_core, test_daub4_2d_constant_nonsquare2);
  tcase_add_test (tc_core, test_daub4_2d_sinusoid_nonsquare1);
  tcase_add_test (tc_core, test_daub4_2d_sinusoid_nonsquare2);
  tcase_add_test (tc_core, test_daub4_3d_sinusoid);
  tcase_add_test (tc_core, test_daub4_3d_sinusoid_nonsquare1);
  tcase_add_test (tc_core, test_daub4_3d_sinusoid_nonsquare2);


  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = daubechies_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
