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

START_TEST (test_hnk_btree_create)
{
  hnk_t *btree;

  btree = hnk_create(4, 
		     100,
		     2,
		     binary_tree_maxk_at_h,
		     binary_tree_hnk,
		     NULL);
  ck_assert_ptr_ne(btree, NULL);

  hnk_destroy(btree);
}
END_TEST

START_TEST (test_hnk_btree_saveload)
{
  hnk_t *t;
  mpz_t a;
  mpz_t b;

  t = hnk_create(4,
		 100,
		 2,
		 binary_tree_maxk_at_h,
		 binary_tree_hnk,
		 NULL);

  ck_assert(t != NULL);

  mpz_init(a);
  mpz_init(b);

  ck_assert(hnk_is_hnk_memoized(t, 4, 20) == 0);
  ck_assert(hnk_highest_memoized_k(t, 4) == 0);
  ck_assert(hnk_get_hnk(t, 4, 20, a) == 0);
  ck_assert(hnk_is_hnk_memoized(t, 4, 20) == 1);
  ck_assert(hnk_highest_memoized_k(t, 4) == 20);

  ck_assert(hnk_save(t, "hnk_btree.data") == 0);

  hnk_destroy(t);
  
  t = hnk_create(4,
		 100,
		 2,
		 binary_tree_maxk_at_h,
		 binary_tree_hnk,
		 NULL);

  printf("Loading\n");
  ck_assert(hnk_highest_memoized_k(t, 4) == 0);
  ck_assert(hnk_is_hnk_memoized(t, 4, 20) == 0);
  ck_assert(hnk_restore(t, "hnk_btree.data") == 0);
  ck_assert(hnk_is_hnk_memoized(t, 4, 20) == 1);
  ck_assert(hnk_highest_memoized_k(t, 4) == 20);
  ck_assert(hnk_get_hnk(t, 4, 20, b) == 0);

  gmp_printf("hnk 4, 20 = %Zd, %Zd\n", a, b);
  ck_assert(mpz_cmp(a, b) == 0);

  
  hnk_destroy(t);
  mpz_clear(a);
  mpz_clear(b);
}
END_TEST


START_TEST (test_hnk_btree_h3)
{
  hnk_t *btree;

  int kmaxs[] = {1, 3, 7, 15, -1};
  int hnks2[] = {1, 1, 2, 5,  6,  6,  4,  1,  0};
  int hnks3[] = {1, 1, 2, 5, 14, 26, 44, 69, 94, 114, 116, 94, 60, 28, 8, 1, 0};

  double ratios2[] = {1.0, 2.0, 2.5, 6.0/5.0, 1.0, 2.0/3.0, 0.25, 0.0};

  int i;

  mpz_t a;
  double r;

  btree = hnk_create(3,    /* Max Height */
		     100,  /* Max k */
		     2,    /* Max no. of splits */
		     binary_tree_maxk_at_h,
		     binary_tree_hnk,
		     NULL);/* Sub tree */
  ck_assert_ptr_ne(btree, NULL);
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {
    /* printf("maxk: %d %d %d\n", i, hnk_get_maxk_at_h(btree, i), kmaxs[i]); */
    ck_assert(hnk_get_maxk_at_h(btree, i) == kmaxs[i]);
  }

  ck_assert(hnk_get_hnk(btree, 1, 2, a) == 0);
  gmp_printf("hnk 1, 2 = %Zd\n", a);

  for (i = 0; i < sizeof(hnks2)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(btree, 2, i, a) == 0);
    ck_assert(mpz_get_ui(a) == hnks2[i]);

  }

  for (i = 0; i < sizeof(ratios2)/sizeof(double); i ++) {

    ck_assert(hnk_get_kplus1_ratio(btree, 
				   2, 
				   i, 
				   &r) == 0);
    ck_assert(r == ratios2[i]);
  }
  
  for (i = 0; i < sizeof(hnks3)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(btree, 3, i, a) == 0);
    ck_assert(mpz_get_ui(a) == hnks3[i]);

  }

  mpz_clear(a);
  hnk_destroy(btree);
}
END_TEST
 
Suite *
hnk_btree_suite (void)
{
  Suite *s = suite_create ("HNK Btree Int");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_hnk_btree_create);
  tcase_add_test (tc_core, test_hnk_btree_saveload);
  tcase_add_test (tc_core, test_hnk_btree_h3);
  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = hnk_btree_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
