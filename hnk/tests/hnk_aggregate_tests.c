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

START_TEST (test_hnk_aggregate_indexing)
{
  int i;
  int splits[] =   {2, 3, 4, 5, 6, 7, 8,  9, 10, 11, 12, 13, 14, 15, 16, 17};
  int expected[] = {2, 4, 4, 8, 8, 8, 8, 16, 16, 16, 16, 16, 16, 16, 16, 32};

  for (i = 0; i < sizeof(splits)/sizeof(int); i ++) {
    /* printf("%d: %d (%d)\n", i, hnk_naggregatesplits(splits[i]), expected[i]); */
    ck_assert(hnk_naggregatesplits(splits[i]) == expected[i]);
  }

}
END_TEST

START_TEST (test_hnk_aggregate_create)
{
  hnk_t *aggregate;
  hnk_t *binary;

  int i;

  binary = hnk_create(5, 100, 2, binary_tree_maxk_at_h, binary_tree_hnk, NULL);
  ck_assert(binary != NULL);

  aggregate = hnk_create_aggregate(6, 
				   100,
				   2,
				   2,
				   binary,
				   binary);
  ck_assert_ptr_ne(aggregate, NULL);

  for (i = 0; i < 5; i ++) {
    /* printf("  %d: %d %d\n", i,  */
    /* 	   hnk_get_maxk_at_h(aggregate, i), */
    /* 	   hnk_get_maxk_at_h(binary, i)); */
    ck_assert(hnk_get_maxk_at_h(aggregate, i) ==
	      hnk_get_maxk_at_h(binary, i));
  }

  hnk_destroy(aggregate);
}
END_TEST

START_TEST (test_hnk_aggregate_binary)
{
  hnk_t *aggregate;
  hnk_t *binary;

  int kmaxs[] = {1, 3, 7, 15, 31, 63, 127, -1};
  int hnks2[] = {1, 1, 2, 5,  6,  6,  4,  1,  0};
  int hnks3[] = {1, 1, 2, 5, 14, 26, 44, 69, 94, 114, 116, 94, 60, 28, 8, 1, 0};

  int i;

  mpz_t a;

  binary = hnk_create(5, 100, 2, binary_tree_maxk_at_h, binary_tree_hnk, NULL);
  ck_assert(binary != NULL);

  aggregate = hnk_create_aggregate(6, 
				   100,
				   2,
				   2,
				   binary,
				   binary);
  ck_assert_ptr_ne(aggregate, NULL);

  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {
    /* printf("maxk: %d %d %d\n", i, hnk_get_maxk_at_h(aggregate, i), kmaxs[i]);  */
    ck_assert(hnk_get_maxk_at_h(aggregate, i) == kmaxs[i]);
  }

  ck_assert(hnk_get_hnk(aggregate, 1, 2, a) == 0);
  /* gmp_printf("hnk 1, 2 = %Zd\n", a); */

  for (i = 0; i < sizeof(hnks2)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(aggregate, 2, i, a) == 0);
    /* gmp_printf("h2: %d: %Zd %d\n", i, a, hnks2[i]); */
    ck_assert(mpz_get_ui(a) == hnks2[i]);

  }

  for (i = 0; i < sizeof(hnks3)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(aggregate, 3, i, a) == 0);
    /* gmp_printf("h3: %d: %Zd %d\n", i, a, hnks3[i]); */
    ck_assert(mpz_get_ui(a) == hnks3[i]);

  }

  mpz_clear(a);
  hnk_destroy(aggregate);
}
END_TEST
 
START_TEST (test_hnk_aggregate_ternary)
{
  hnk_t *aggregate;
  hnk_t *ternary;

  int kmaxs[] = {1, 4, 13, 40, 121, 364, 1093, -1};
  int hnks1[] = {1, 1, 3, 3, 1, 0};
  int hnks2[] = {1, 1, 3, 12, 28, 57, 96, 129, 144};
  int hnks3[] = {1, 1, 3, 12, 55, 192, 618, 1893, 5436, 14772};

  int i;

  mpz_t a;
  mpz_t b;

  ternary = hnk_create(5, 100, 3, ternary_tree_maxk_at_h, ternary_tree_hnk, NULL);
  ck_assert(ternary != NULL);

  aggregate = hnk_create_aggregate(6, 
				   100,
				   3,
				   3,
				   ternary,
				   ternary,
				   ternary);
  ck_assert_ptr_ne(aggregate, NULL);

  mpz_init(a);
  mpz_init(b);

  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {
    /* printf("maxk: %d %d %d\n", i, hnk_get_maxk_at_h(aggregate, i), kmaxs[i]); */
    ck_assert(hnk_get_maxk_at_h(aggregate, i) == kmaxs[i]);
  }

  for (i = 0; i < sizeof(hnks1)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(aggregate, 1, i, a) == 0);
    ck_assert(hnk_get_hnk(ternary, 1, i, b) == 0);
    /* gmp_printf("h1: %d: %Zd (%Zd) %d\n", i, a, b, hnks1[i]);  */
    ck_assert(mpz_get_ui(a) == hnks1[i]);

  }
  
  for (i = 0; i < sizeof(hnks2)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(aggregate, 2, i, a) == 0);
    ck_assert(hnk_get_hnk(ternary, 2, i, b) == 0);
    /* gmp_printf("h2: %d: %Zd (%Zd) %d\n", i, a, b, hnks2[i]); */
    ck_assert(mpz_get_ui(a) == hnks2[i]);

  }

  for (i = 0; i < sizeof(hnks3)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(aggregate, 3, i, a) == 0);
    ck_assert(hnk_get_hnk(ternary, 3, i, b) == 0);
    /* gmp_printf("h3: %d: %Zd (%Zd) %d\n", i, a, b, hnks3[i]);  */
    ck_assert(mpz_get_ui(a) == hnks3[i]);

  }

  mpz_clear(a);
  mpz_clear(b);
  hnk_destroy(aggregate);
}
END_TEST

Suite *
hnk_aggregate_suite (void)
{
  Suite *s = suite_create ("HNK Aggregate Int");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_hnk_aggregate_indexing);
  tcase_add_test (tc_core, test_hnk_aggregate_create);
  tcase_add_test (tc_core, test_hnk_aggregate_binary);
  tcase_add_test (tc_core, test_hnk_aggregate_ternary);
  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = hnk_aggregate_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
