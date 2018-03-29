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

START_TEST (test_hnk_octree_create)
{
  hnk_t *octree;

  octree = hnk_create(4, 
		      100,
		      8,
		      octary_tree_maxk_at_h,
		      octary_tree_hnk,
		     NULL);
  ck_assert_ptr_ne(octree, NULL);

  hnk_destroy(octree);
}
END_TEST

START_TEST (test_hnk_octree_h3)
{
  hnk_t *octree;

  int kmaxs[] = {1, 9, 73, 585, -1};

  int hnk_i[] = {0, 1, 2, 3, 72, 73, 74};
  int hnk_c[] = {1, 1, 8, 92, 64, 1, 0};

  double ratios[] = 
    {1.0, 8.0, 92.0/8.0};

  int i;

  mpz_t a;
  double r;

  octree = hnk_create(3, 
		      100,
		      8,
		      octary_tree_maxk_at_h,
		      octary_tree_hnk,
		     NULL);
  ck_assert_ptr_ne(octree, NULL);
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {

    ck_assert(hnk_get_maxk_at_h(octree, i) == kmaxs[i]);

  }

  for (i = 0; i < sizeof(hnk_i)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(octree, 2, hnk_i[i], a) == 0);

    ck_assert(mpz_get_ui(a) == hnk_c[i]);

  }

  for (i = 0; i < sizeof(ratios)/sizeof(double); i ++) {

    ck_assert(hnk_get_kplus1_ratio(octree, 
				   2, 
				   i, 
				   &r) == 0);

    ck_assert(fabs(r - ratios[i]) < 1.0e-9);
  }
  
  mpz_clear(a);
  hnk_destroy(octree);
}
END_TEST

START_TEST (test_hnk_aggregate_octree_h3)
{
  hnk_t *octree;
  hnk_t *agg;
  hnk_t *subtree;

  int kmaxs[] = {1, 9, 73, 585, -1};

  int hnk_i[] = {0, 1, 2, 3, 72, 73, 74};
  int hnk_c[] = {1, 1, 8, 92, 64, 1, 0};

  int i;

  mpz_t a;

  octree = hnk_create(3, 
		      100,
		      8,
		      octary_tree_maxk_at_h,
		      octary_tree_hnk,
		     NULL);
  ck_assert_ptr_ne(octree, NULL);

  subtree = hnk_create_octary_tree(2,
				   100);
  ck_assert(subtree != NULL);

  agg = hnk_create_aggregate(3, 100, 8, 8,
			     subtree, subtree, subtree, subtree,
			     subtree, subtree, subtree, subtree);
  ck_assert(agg != NULL);
  
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {

    ck_assert(hnk_get_maxk_at_h(octree, i) == kmaxs[i]);
    ck_assert(hnk_get_maxk_at_h(agg, i) == kmaxs[i]);
    
  }

  for (i = 0; i < sizeof(hnk_i)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(octree, 2, hnk_i[i], a) == 0);

    ck_assert(mpz_get_ui(a) == hnk_c[i]);

    ck_assert(hnk_get_hnk(agg, 2, hnk_i[i], a) == 0);

    ck_assert(mpz_get_ui(a) == hnk_c[i]);
  }

  mpz_clear(a);
  hnk_destroy(octree);
  hnk_destroy(agg);
}
END_TEST

Suite *
hnk_octree_suite (void)
{
  Suite *s = suite_create ("HNK Octree Int");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_hnk_octree_create);
  tcase_add_test (tc_core, test_hnk_octree_h3);
  tcase_add_test (tc_core, test_hnk_aggregate_octree_h3);
  
  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = hnk_octree_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
