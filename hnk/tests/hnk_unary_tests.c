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
#include "hnk_standard.h"

START_TEST (test_hnk_utree_create)
{
  hnk_t *utree;

  utree = hnk_create(4, 
		     100,
		     2,
		     unary_tree_maxk_at_h,
		     unary_tree_hnk,
		     NULL);
  ck_assert_ptr_ne(utree, NULL);

  hnk_destroy(utree);
}
END_TEST

START_TEST (test_hnk_utree_h3)
{
  hnk_t *utree;

  int kmaxs[] = {1, 2, 3, 4, -1};

  int i;

  mpz_t a;

  utree = hnk_create(3,    /* Max Height */
		     100,  /* Max k */
		     1,    /* Max no. of splits */
		     unary_tree_maxk_at_h,
		     unary_tree_hnk,
		     NULL);/* Sub tree */
  ck_assert_ptr_ne(utree, NULL);
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {
    /* printf("maxk: %d %d %d\n", i, hnk_get_maxk_at_h(utree, i), kmaxs[i]); */
    ck_assert(hnk_get_maxk_at_h(utree, i) == kmaxs[i]);
  }

  ck_assert(hnk_get_hnk(utree, 1, 2, a) == 0);
  /* gmp_printf("hnk 1, 2 = %Zd\n", a); */

  for (i = 0; i <= 3; i ++) {

    ck_assert(hnk_get_hnk(utree, 2, i, a) == 0);

    /* gmp_printf("hnk: %d %Zd\n", i, a); */
    ck_assert(mpz_get_ui(a) == 1);
  }

  mpz_clear(a);
  hnk_destroy(utree);
}
END_TEST
 
Suite *
hnk_utree_suite (void)
{
  Suite *s = suite_create ("HNK Utree Int");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_hnk_utree_create);
  tcase_add_test (tc_core, test_hnk_utree_h3);
  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = hnk_utree_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
