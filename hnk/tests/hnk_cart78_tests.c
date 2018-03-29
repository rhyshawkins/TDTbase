//
//    HNK Library : A library for computing combinations of arrangements of
//    general trees for the Trans-dimensional Tree algorithm. See
//
//      R Hawkins and M Sambridge, "Geophysical imaging using trans-dimensional trees",
//      Geophysical Journal International, 2015, 203:2, 972 - 1000,
//      https://doi.org/10.1093/gji/ggv326
//    
//    Copyright (C) 2014 - 2018 Rhys Hawkins
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//
#include <stdio.h>
#include <stdlib.h>
#include <check.h>
#include <math.h>

#include "hnk.h"
#include "hnk_standard.h"
#include "hnk_cartesian.h"

START_TEST (test_hnk_cart78_create)
{
  hnk_t *cart78;

  cart78 = hnk_cartesian_78_create(4, 100);
  ck_assert_ptr_ne(cart78, NULL);
  hnk_destroy(cart78);

}
END_TEST


START_TEST (test_hnk_cart78_h3)
{
  hnk_t *cart78;

  int kmaxs[] = {1, 8, 64, 512, -1};

  int hnks[] = { 1, 
		 1,
		 7,
		 77,
		 1015,
		 11179,
		 115563,
		 1155707,
		 11191895,
		 105454216,
		 969258381 };
  int i;
  mpz_t a;

  cart78 = hnk_cartesian_78_create(3, 
				   100);
  ck_assert_ptr_ne(cart78, NULL);
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {

    ck_assert(hnk_get_maxk_at_h(cart78, i) == kmaxs[i]);

  }

  for (i = 0; i < sizeof(hnks)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(cart78, 3, i, a) == 0);

    //gmp_printf("hnk: %d %d\n", mpz_get_ui(a), hnks[i]);
    ck_assert(mpz_get_ui(a) == hnks[i]);

  }

  mpz_clear(a);
  hnk_destroy(cart78);
}
END_TEST

START_TEST (test_hnk_cart78_h4)
{
  hnk_t *cart78;

  int kmaxs[] = {1, 8, 64, 512, 4096, -1};

  int hnks[] = { 1, 
		 1,
		 7,
		 77,
		 1015,
		 14763,
		 199787,
		 2585723,
		 32673495,
		 406123144 };

  int i;
  mpz_t a;

  cart78 = hnk_cartesian_78_create(4, 
				   100);
  ck_assert_ptr_ne(cart78, NULL);
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {

    //printf("kmax: %d %d\n", hnk_get_maxk_at_h(cart78, i), kmaxs[i]);
    ck_assert(hnk_get_maxk_at_h(cart78, i) == kmaxs[i]);

  }

  for (i = 0; i < sizeof(hnks)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(cart78, 4, i, a) == 0);

    //gmp_printf("hnk: %d %d\n", mpz_get_ui(a), hnks[i]);
    ck_assert(mpz_get_ui(a) == hnks[i]);

  }

  mpz_clear(a);
  hnk_destroy(cart78);
}
END_TEST

START_TEST (test_hnk_cart78_h5)
{
  hnk_t *cart78;

  int kmaxs[] = {1, 8, 64, 512, 4096, 32768, -1};

  int hnks[] = {1, 
		1,
		7,
		77,
		1015,
		14763,
		228459,
		3460219,
		51037911,
		739040904 };

  int i;

  mpz_t a;

  cart78 = hnk_cartesian_78_create(5, 
				   100);
  ck_assert_ptr_ne(cart78, NULL);
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {

    //printf("kmax: %d %d\n", hnk_get_maxk_at_h(cart78, i), kmaxs[i]);
    ck_assert(hnk_get_maxk_at_h(cart78, i) == kmaxs[i]);

  }

  for (i = 0; i < sizeof(hnks)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(cart78, 5, i, a) == 0);

    //gmp_printf("hnk: %d %d\n", mpz_get_ui(a), hnks[i]);
    ck_assert(mpz_get_ui(a) == hnks[i]);

  }

  mpz_clear(a);
  hnk_destroy(cart78);
}
END_TEST
 
START_TEST (test_hnk_cart78_h6)
{
  hnk_t *cart78;

  int kmaxs[] = {1, 8, 64, 512, 4096, 32768, 262144, -1};

  int hnks[] = {1, 1, 3, 15, 91, 612, 4389, 29818, 194571};

  int i;

  mpz_t a;

  cart78 = hnk_cartesian_78_create(6, 
				   100);
  ck_assert_ptr_ne(cart78, NULL);
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {

    printf("kmax: %d %d\n", hnk_get_maxk_at_h(cart78, i), kmaxs[i]);
    ck_assert(hnk_get_maxk_at_h(cart78, i) == kmaxs[i]);

  }

  for (i = 0; i < sizeof(hnks)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(cart78, 6, i, a) == 0);

    gmp_printf("hnk: %d %d\n", mpz_get_ui(a), hnks[i]);
    //ck_assert(mpz_get_ui(a) == hnks[i]);

  }

  mpz_clear(a);
  hnk_destroy(cart78);
}
END_TEST

Suite *
hnk_cart78_suite (void)
{
  Suite *s = suite_create ("HNK Cart78 Int");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_hnk_cart78_create);
  tcase_add_test (tc_core, test_hnk_cart78_h3);
  tcase_add_test (tc_core, test_hnk_cart78_h4);
  tcase_add_test (tc_core, test_hnk_cart78_h5);
  tcase_add_test (tc_core, test_hnk_cart78_h6);
  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = hnk_cart78_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
