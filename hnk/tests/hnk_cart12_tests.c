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

/*
 * Counts for this test come from:
 *
 *   http://oeis.org/A137560
 * 
 * And the raw table:
 *
 *   http://oeis.org/A137560/b137560.txt
 */
START_TEST (test_hnk_cart12_create)
{
  hnk_t *cart12;

  cart12 = hnk_cartesian_12_create(4, 100);
  ck_assert_ptr_ne(cart12, NULL);
  hnk_destroy(cart12);

}
END_TEST

START_TEST (test_hnk_cart12_h3)
{
  hnk_t *cart12;

  int kmaxs[] = {1, 2, 4, 8, -1};

  int hnks[] = {1, 
		1,
		1,
		2,
		5,
		6,
		6,
		4,
		1,
		0};

  int i;

  mpz_t a;

  cart12 = hnk_cartesian_12_create(3, 
				   100);
  ck_assert_ptr_ne(cart12, NULL);
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {

    ck_assert(hnk_get_maxk_at_h(cart12, i) == kmaxs[i]);

  }
  
  for (i = 0; i < sizeof(hnks)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(cart12, 3, i, a) == 0);

    //gmp_printf("hnk: %d: %Zd %d\n", i, a, hnks[i]);
    ck_assert(mpz_get_ui(a) == hnks[i]);

  }

  mpz_clear(a);
  hnk_destroy(cart12);
}
END_TEST

START_TEST (test_hnk_cart12_h4)
{
  hnk_t *cart12;

  int kmaxs[] = {1, 2, 4, 8, 16, -1};

  int hnks[] = {1, 
		1,
		1,
		2,
		5,
		14,
		26,
		44,
		69,
		94,
		114,
		116,
		94,
		60,
		28,
		8,
		1,
		0};

  int i;

  mpz_t a;

  cart12 = hnk_cartesian_12_create(4, 
				   100);
  ck_assert_ptr_ne(cart12, NULL);
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {

    //printf("kmax: %d %d\n", hnk_get_maxk_at_h(cart12, i), kmaxs[i]);
    ck_assert(hnk_get_maxk_at_h(cart12, i) == kmaxs[i]);

  }

  for (i = 0; i < sizeof(hnks)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(cart12, 4, i, a) == 0);
    //gmp_printf("hnk: %d: %Zd %d\n", i, a, hnks[i]);
    ck_assert(mpz_get_ui(a) == hnks[i]);

  }

  mpz_clear(a);
  hnk_destroy(cart12);
}
END_TEST

START_TEST (test_hnk_cart12_h5)
{
  hnk_t *cart12;

  int kmaxs[] = {1, 2, 4, 8, 16, 32, -1};

  int hnks[] = { 1,
		 1,
		 1,
		 2,
		 5,
		 14,
		 42,
		 100,
		 221,
		 470,
		 958,
		 1860,
		 3434,
		 6036,
		 10068,
		 15864,
		 23461,
		 32398,
		 41658,
		 49700,
		 54746,
		 55308,
		 50788,
		 41944,
		 30782,
		 19788,
		 10948,
		 5096,
		 1932,
		 568,
		 120,
		 16,
		 1,
		 0 };

  int i;
  mpz_t a;

  //printf("hnk12 5\n");
  cart12 = hnk_cartesian_12_create(5, 
				   100);
  ck_assert_ptr_ne(cart12, NULL);
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {

    //printf("kmax: %d %d\n", hnk_get_maxk_at_h(cart12, i), kmaxs[i]);
    ck_assert(hnk_get_maxk_at_h(cart12, i) == kmaxs[i]);

  }

  for (i = 0; i < sizeof(hnks)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(cart12, 5, i, a) == 0);

    //gmp_printf("hnk: %d %Zd %d\n", i, a, hnks[i]);
    ck_assert(mpz_get_ui(a) == hnks[i]);

  }

  mpz_clear(a);
  hnk_destroy(cart12);
}
END_TEST
 
START_TEST (test_hnk_cart12_h6)
{
  hnk_t *cart12;

  int kmaxs[] = {1, 2, 4, 8, 16, 32, 64, -1};

  int hnks[] = { 1,
		 1,
		 1,
		 2,
		 5,
		 14,
		 42,
		 132,
		 365,
		 950,
		 2398,
		 5916,
		 14290,
		 33708,
		 77684,
		 175048,
		 385741,
		 831014,
		 1749654,
		 3598964,
		 7228014,
		 14162220,
		 27049196,
		 50323496,
		 91143114,
		 160617860,
		 275276716,
		 458591432,
		 742179284,
		 1166067016,
		 1777171560 };

  int i;
  mpz_t a;

  cart12 = hnk_cartesian_12_create(6, 
				   100);
  ck_assert_ptr_ne(cart12, NULL);
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {

    //printf("kmax: %d %d\n", hnk_get_maxk_at_h(cart12, i), kmaxs[i]);
    ck_assert(hnk_get_maxk_at_h(cart12, i) == kmaxs[i]);

  }

  for (i = 0; i < sizeof(hnks)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(cart12, 6, i, a) == 0);

    //gmp_printf("hnk: %d %Zd %d\n", i, a, hnks[i]);
    ck_assert(mpz_get_ui(a) == hnks[i]);

  }

  mpz_clear(a);
  hnk_destroy(cart12);
}
END_TEST

Suite *
hnk_cart12_suite (void)
{
  Suite *s = suite_create ("HNK Cart12 Int");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_hnk_cart12_create);
  tcase_add_test (tc_core, test_hnk_cart12_h3);
  tcase_add_test (tc_core, test_hnk_cart12_h4);
  tcase_add_test (tc_core, test_hnk_cart12_h5);
  tcase_add_test (tc_core, test_hnk_cart12_h6);
  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = hnk_cart12_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
