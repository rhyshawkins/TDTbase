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

#include "hnk.h"
#include "hnk_s2.h"

START_TEST (test_hnk_s2_create)
{
  hnk_t *t;

  /*
   * Icosahedron 
   */
  t = hnk_s2_icosahedron_create(4,
				100);
  ck_assert(t != NULL);

  hnk_destroy(t);
}
END_TEST

START_TEST (test_hnk_s2_saveload)
{
  hnk_t *t;
  hnk_t *u;
  mpz_t a;
  mpz_t b;
  
  t = hnk_s2_icosahedron_create(4,
				100);
  ck_assert(t != NULL);
  
  mpz_init(a);
  
  ck_assert(hnk_is_hnk_memoized(t, 4, 20) == 0);
  ck_assert(hnk_get_hnk(t, 4, 20, a) == 0);
  ck_assert(hnk_is_hnk_memoized(t, 4, 20) == 1);
  
  ck_assert(hnk_save(t, "hnk_s2.data") == 0);

  u = hnk_s2_icosahedron_create(4,
				100);
  ck_assert(u != NULL);

  ck_assert(hnk_is_hnk_memoized(u, 4, 20) == 0);
  ck_assert(hnk_restore(u, "hnk_s2.data") == 0);
  ck_assert(hnk_is_hnk_memoized(u, 4, 20) == 1);

  mpz_init(b);

  ck_assert(hnk_get_hnk(u, 4, 20, b) == 0);
  ck_assert(mpz_cmp(a, b) == 0);

  
  ck_assert(hnk_get_hnk(t, 4, 21, a) == 0);
  ck_assert(hnk_get_hnk(u, 4, 21, b) == 0);
  ck_assert(mpz_cmp(a, b) == 0);
  

  hnk_destroy(t);
  hnk_destroy(u);
  mpz_clear(a);
  mpz_clear(b);
}
END_TEST

START_TEST (test_hnk_s2_icosahedron_h1)
{
  hnk_t *s2;

  int kmaxs[] = {1, 21, -1};
  
  /* Verified by hand (complete) */
  int hnks1[] = {1, 1, 20, 190, 1140, 4845, 15504, 38760, 77520, 125970, 167960, 184756,
		 167960, 125970, 77520, 38760, 15504, 4845, 1140, 190, 20, 1, 0};

  int i;

  mpz_t a;

  s2 = hnk_s2_icosahedron_create(1,   
				 100);
  ck_assert_ptr_ne(s2, NULL);
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {
    printf("maxk: %d %d %d\n", i, hnk_get_maxk_at_h(s2, i), kmaxs[i]);
    /* ck_assert(hnk_get_maxk_at_h(s2, i) == kmaxs[i]); */
  }

  for (i = 0; i < sizeof(hnks1)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(s2, 1, i, a) == 0);
    gmp_printf("hnk1: %d %Zd %d\n", i, a, hnks1[i]);
    /* ck_assert(mpz_get_ui(a) == hnks1[i]); */

  }

  mpz_clear(a);
  hnk_destroy(s2);
}
END_TEST

Suite *
hnk_s2_suite (void)
{
  Suite *s = suite_create ("HNK S2");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");

  tcase_add_test (tc_core, test_hnk_s2_create);
  tcase_add_test (tc_core, test_hnk_s2_saveload);

  tcase_add_test (tc_core, test_hnk_s2_icosahedron_h1);

  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = hnk_s2_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
