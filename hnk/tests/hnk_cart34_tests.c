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

START_TEST (test_hnk_cart34_create)
{
  hnk_t *cart34;

  cart34 = hnk_cartesian_34_create(4, 100);
  ck_assert_ptr_ne(cart34, NULL);
  hnk_destroy(cart34);

}
END_TEST

START_TEST (test_hnk_cart34_h3)
{
  hnk_t *cart34;

  int kmaxs[] = {1, 4, 16, 64, -1};

  int hnks[] = {1, 
		1,
		3,
		15,
		91,
		420,
		1797,
		7354,
		28635,
		107052,
		384076,
		1320567,
		4357437,
		13821268,
		42199062,
		124126764,
		351761603,
		959756628};

  double ratios[] = 
    {1.0, 3.0, 5.0, 91.0/15.0, 420.0/91.0};

  int i;

  mpz_t a;
  double r;

  cart34 = hnk_cartesian_34_create(3, 
				   100);
  ck_assert_ptr_ne(cart34, NULL);
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {

    ck_assert(hnk_get_maxk_at_h(cart34, i) == kmaxs[i]);

  }
  
  for (i = 0; i < sizeof(hnks)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(cart34, 3, i, a) == 0);

    //gmp_printf("hnk: %d: %Zd %d\n", i, a, hnks[i]);
    ck_assert(mpz_get_ui(a) == hnks[i]);

  }

  for (i = 0; i < sizeof(ratios)/sizeof(double); i ++) {

    ck_assert(hnk_get_kplus1_ratio(cart34, 
				   3, 
				   i, 
				   &r) == 0);

    ck_assert(fabs(r - ratios[i]) < 1.0e-9);
  }
  
  mpz_clear(a);
  hnk_destroy(cart34);
}
END_TEST

START_TEST (test_hnk_cart34_h4)
{
  hnk_t *cart34;

  int kmaxs[] = {1, 4, 16, 64, 256, -1};

  int hnks[] = {1, 
		1,
		3,
		15,
		91,
		612,
		3621,
		20218,
		109707,
		582156,
		3028780,
		15474231,
		77752605,
		384661772,
		1875096126};

  int i;

  mpz_t a;

  cart34 = hnk_cartesian_34_create(4, 
				   100);
  ck_assert_ptr_ne(cart34, NULL);
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {

    //printf("kmax: %d %d\n", hnk_get_maxk_at_h(cart34, i), kmaxs[i]);
    ck_assert(hnk_get_maxk_at_h(cart34, i) == kmaxs[i]);

  }

  for (i = 0; i < sizeof(hnks)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(cart34, 4, i, a) == 0);
    //gmp_printf("hnk: %d: %d %d\n", i, mpz_get_ui(a), hnks[i]);
    ck_assert(mpz_get_ui(a) == hnks[i]);

  }

  mpz_clear(a);
  hnk_destroy(cart34);
}
END_TEST

START_TEST (test_hnk_cart34_h5)
{
  hnk_t *cart34;

  int kmaxs[] = {1, 4, 16, 64, 256, 1024, -1};

  int hnks[] = { 1, 
		 1,
		 3,
		 15,
		 91,
		 612,
		 4389,
		 29818,
		 194571,
		 1240140,
		 7783276,
		 48276855,
		 296371677,
		 1802463884 };

  int i;
  mpz_t a;

  //printf("hnk34 5\n");
  cart34 = hnk_cartesian_34_create(5, 
				   100);
  ck_assert_ptr_ne(cart34, NULL);
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {

    //printf("kmax: %d %d\n", hnk_get_maxk_at_h(cart34, i), kmaxs[i]);
    ck_assert(hnk_get_maxk_at_h(cart34, i) == kmaxs[i]);

  }

  for (i = 0; i < sizeof(hnks)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(cart34, 5, i, a) == 0);

    //gmp_printf("hnk: %d %d %d\n", i, mpz_get_ui(a), hnks[i]);
    ck_assert(mpz_get_ui(a) == hnks[i]);

  }

  mpz_clear(a);
  hnk_destroy(cart34);
}
END_TEST
 
START_TEST (test_hnk_cart34_h6)
{
  hnk_t *cart34;

  int kmaxs[] = {1, 4, 16, 64, 256, 1024, 4096, -1};

  int hnks[] = {1, 
		1,
		3,
		15,
		91,
		612,
		4389,
		32890,
		242187,
		1740876,
		12289132,
		85656951,
		591596253};

  int i;
  mpz_t a;

  cart34 = hnk_cartesian_34_create(6, 
				   100);
  ck_assert_ptr_ne(cart34, NULL);
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {

    //printf("kmax: %d %d\n", hnk_get_maxk_at_h(cart34, i), kmaxs[i]);
    ck_assert(hnk_get_maxk_at_h(cart34, i) == kmaxs[i]);

  }

  for (i = 0; i < sizeof(hnks)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(cart34, 6, i, a) == 0);

    //gmp_printf("hnk: %d %d\n", mpz_get_ui(a), hnks[i]);
    ck_assert(mpz_get_ui(a) == hnks[i]);

  }

  mpz_clear(a);
  hnk_destroy(cart34);
}
END_TEST

Suite *
hnk_cart34_suite (void)
{
  Suite *s = suite_create ("HNK Cart34 Int");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_hnk_cart34_create);
  tcase_add_test (tc_core, test_hnk_cart34_h3);
  tcase_add_test (tc_core, test_hnk_cart34_h4);
  tcase_add_test (tc_core, test_hnk_cart34_h5);
  tcase_add_test (tc_core, test_hnk_cart34_h6);
  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = hnk_cart34_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
