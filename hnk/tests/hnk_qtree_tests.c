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

START_TEST (test_hnk_qtree_create)
{
  hnk_t *qtree;

  qtree = hnk_create(4, 
		     100,
		     2,
		     binary_tree_maxk_at_h,
		     binary_tree_hnk,
		     NULL);
  ck_assert_ptr_ne(qtree, NULL);

  hnk_destroy(qtree);
}
END_TEST

START_TEST (test_hnk_qtree_h3)
{
  hnk_t *qtree;

  int kmaxs[] = {1, 5, 21, 85, -1};
  int hnks[] = {1, 1, 4, 22, 76, 233, 620, 1420, 2876, 5156, 8112, 11182, 13420, 
		13750, 11704, 8056, 4372, 1820, 560, 120, 16, 1, 0};
  double ratios[] = 
    {1.0, 4.0, 11.0/2.0, 38.0/11.0, 233.0/76.0,
     620.0/233.0, 1420.0/620.0, 2876.0/1420.0, 5156.0/2876.0};

  int i;

  mpz_t a;
  double r;

  qtree = hnk_create(3, 
		     100,
		     4,
		     quaternary_tree_maxk_at_h,
		     quaternary_tree_hnk,
		     NULL);
  ck_assert_ptr_ne(qtree, NULL);
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {

    ck_assert(hnk_get_maxk_at_h(qtree, i) == kmaxs[i]);

  }

  for (i = 0; i < sizeof(hnks)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(qtree, 2, i, a) == 0);
    ck_assert(mpz_get_ui(a) == hnks[i]);

  }

  for (i = 0; i < sizeof(ratios)/sizeof(double); i ++) {

    ck_assert(hnk_get_kplus1_ratio(qtree, 
				   2, 
				   i, 
				   &r) == 0);

    ck_assert(fabs(r - ratios[i]) < 1.0e-9);
  }
  
  mpz_clear(a);
  hnk_destroy(qtree);
}
END_TEST
 
Suite *
hnk_qtree_suite (void)
{
  Suite *s = suite_create ("HNK Qtree Int");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_hnk_qtree_create);
  tcase_add_test (tc_core, test_hnk_qtree_h3);
  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = hnk_qtree_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
