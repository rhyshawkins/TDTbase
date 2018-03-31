
#include <stdio.h>
#include <stdlib.h>
#include <check.h>
#include <math.h>

#include "cdf97_lift.h"
#include "cdf97_lift_impulse.h"

START_TEST (test_cdf97_impulse_size)
{
  ck_assert(cdf97_lift_impulse_1dsize(8, 8) == 9);
  ck_assert(cdf97_lift_impulse_1dsize(8, 7) == 23);
  ck_assert(cdf97_lift_impulse_1dsize(8, 6) == 51);
  ck_assert(cdf97_lift_impulse_1dsize(8, 5) == 107);
  ck_assert(cdf97_lift_impulse_1dsize(8, 4) == 219);
}
END_TEST

START_TEST (test_cdf97_impulse_1d_top)
{
  /*
   * Twice as large as needed to check for errors
   */
  #define SIZE_1D_TOP 18
  #define SIZE_1D_WAVE 64

  double coeff[SIZE_1D_TOP];
  double wave[SIZE_1D_WAVE];
  double workspace[SIZE_1D_WAVE];

  int i;

  int coeff_centre;
  int wave_centre;
  

  ck_assert(cdf97_lift_impulse_1d(8,
				  8,
				  1.0,
				  coeff,
				  SIZE_1D_TOP,
				  SIZE_1D_TOP/2) >= 0);

  memset(wave, 0, sizeof(double) * SIZE_1D_WAVE);
  wave[(SIZE_1D_WAVE * 3)/4] = 1.0;

  ck_assert(cdf97_lift_inverse1d_cdf97(wave, SIZE_1D_WAVE, 1, workspace) >= 0);
  
  coeff_centre = SIZE_1D_TOP/2;
  wave_centre = SIZE_1D_WAVE/2 + 1;

  ck_assert(coeff[coeff_centre] == wave[wave_centre]);

  /* This assumes that the size of the coefficient array is smaller than the wavelet */

  /* Backwards */
  for (i = 1; (coeff_centre - i) >= 0; i ++) {
    /* printf(" %d %f %d %f\n", coeff_centre - i, coeff[coeff_centre - i], wave_centre - i, wave[wave_centre - i]);  */
    ck_assert(coeff[coeff_centre - i] == wave[wave_centre - i]);
  }

  /* Forwards */
  for (i = 1; (coeff_centre + i) < SIZE_1D_TOP; i ++) {
    /* printf(" %d %f %d %f\n", coeff_centre + i, coeff[coeff_centre + i], wave_centre + i, wave[wave_centre + i]);  */
    ck_assert(coeff[coeff_centre + i] == wave[wave_centre + i]);
  }
}
END_TEST


Suite *
cdf97_suite (void)
{
  Suite *s = suite_create ("CDF97");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");

  tcase_add_test (tc_core, test_cdf97_impulse_size);
  tcase_add_test (tc_core, test_cdf97_impulse_1d_top);

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
