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
#include "hnk_butterfly.h"

START_TEST (test_hnk_butterfly_create)
{
  hnk_t *butterfly;

  /*
   * Tetrahedron
   */
  butterfly = hnk_butterfly_tetrahedron_create(4, 
					       100);
  ck_assert(butterfly != NULL);

  hnk_destroy(butterfly);

  /* 
   * Octahedron
   */
  butterfly = hnk_butterfly_octahedron_create(4,
					      100);
  ck_assert(butterfly != NULL);

  hnk_destroy(butterfly);

  /*
   * Icosahedron
   */
  butterfly = hnk_butterfly_icosahedron_create(4,
					       100);
  ck_assert(butterfly != NULL);

  hnk_destroy(butterfly);
}
END_TEST

START_TEST (test_hnk_butterfly_tetrahedron_h3)
{
  hnk_t *butterfly;

  int kmaxs[] = {1, 5, 11, 35, -1};

  /* Verified by hand (complete) */
  int hnks1[] = {1, 1, 4, 6, 4, 1, 0};

  /* Verified by hand (k = 4) */
  int hnks2[] = {1, 1, 4, 12, 28, 48, 66, 72, 56, 28, 8, 1, 0};

  /* Verified by hand (k = 3) */
  int hnks3[] = {1, 1, 4, 12, 52, 204, 678, 2082, 6010, 16414};

  int i;

  mpz_t a;

  butterfly = hnk_butterfly_tetrahedron_create(3,   
					       100);
  ck_assert_ptr_ne(butterfly, NULL);
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {
    /* printf("maxk: %d %d %d\n", i, hnk_get_maxk_at_h(butterfly, i), kmaxs[i]);  */
    ck_assert(hnk_get_maxk_at_h(butterfly, i) == kmaxs[i]);
  }

  for (i = 0; i < sizeof(hnks1)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(butterfly, 1, i, a) == 0);
    /* gmp_printf("hnk1: %d %Zd %d\n", i, a, hnks1[i]); */
    ck_assert(mpz_get_ui(a) == hnks1[i]);

  }

  for (i = 0; i < sizeof(hnks2)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(butterfly, 2, i, a) == 0);
    /* gmp_printf("hnk2: %d %Zd %d\n", i, a, hnks2[i]); */
    ck_assert(mpz_get_ui(a) == hnks2[i]);

  }

  for (i = 0; i < sizeof(hnks3)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(butterfly, 3, i, a) == 0);
    /* gmp_printf("hnk3: %d %Zd %d\n", i, a, hnks3[i]); */
    ck_assert(mpz_get_ui(a) == hnks3[i]);

  }

  mpz_clear(a);
  hnk_destroy(butterfly);
}
END_TEST
 
START_TEST (test_hnk_butterfly_octahedron_h3)
{
  hnk_t *butterfly;

  int kmaxs[] = {1, 7, 19, 67, -1};

  /* Verified by hand (complete)*/
  int hnks1[] = {1, 1, 6, 15, 20, 15, 6, 1, 0};

  /* Verified by hand (k = 3) */
  int hnks2[] = {1, 1, 6, 27, 92, 253, 590, 1175, 2020};

  /* Verified by hand (k = 3) */
  int hnks3[] = {1, 1, 6, 27, 140, 661, 2774, 10979, 41432};

  int i;

  mpz_t a;

  butterfly = hnk_butterfly_octahedron_create(3,   
					      100);
  ck_assert_ptr_ne(butterfly, NULL);
  mpz_init(a);

  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {
    /* printf("maxk: %d %d %d\n", i, hnk_get_maxk_at_h(butterfly, i), kmaxs[i]);  */
    ck_assert(hnk_get_maxk_at_h(butterfly, i) == kmaxs[i]);
  }

  for (i = 0; i < sizeof(hnks1)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(butterfly, 1, i, a) == 0);
    /* gmp_printf("hnk1: %d %Zd %d\n", i, a, hnks1[i]); */
    ck_assert(mpz_get_ui(a) == hnks1[i]);

  }

  for (i = 0; i < sizeof(hnks2)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(butterfly, 2, i, a) == 0);
    /* gmp_printf("hnk2: %d %Zd %d\n", i, a, hnks2[i]); */
    ck_assert(mpz_get_ui(a) == hnks2[i]);

  }

  for (i = 0; i < sizeof(hnks3)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(butterfly, 3, i, a) == 0);
    /* gmp_printf("hnk3: %d %Zd %d\n", i, a, hnks3[i]); */
    ck_assert(mpz_get_ui(a) == hnks3[i]);

  }

  mpz_clear(a);
  hnk_destroy(butterfly);
}
END_TEST

START_TEST (test_hnk_butterfly_icosahedron_h3)
{
  hnk_t *butterfly;

  int kmaxs[] = {1, 13, 43, 163, -1};

  /* Verified by hand (complete)*/
  int hnks1[] = {1, 1, 6, 15, 20, 15, 6, 1, 0};

  /* Verified by hand (k = 3) */
  int hnks2[] = {1, 1, 6, 27, 92, 253, 590, 1175, 2020};

  /* Verified by hand (k = 3) */
  int hnks3[] = {1, 1, 6, 27, 140, 661, 2774, 10979, 41432};

  int i;

  mpz_t a;

  butterfly = hnk_butterfly_icosahedron_create(3,   
					       100);
  ck_assert_ptr_ne(butterfly, NULL);
  mpz_init(a);

  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {
    printf("maxk: %d %d %d\n", i, hnk_get_maxk_at_h(butterfly, i), kmaxs[i]);
    ck_assert(hnk_get_maxk_at_h(butterfly, i) == kmaxs[i]);
  }

  for (i = 0; i < sizeof(hnks1)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(butterfly, 1, i, a) == 0);
    gmp_printf("tthnk1: %d %Zd %d\n", i, a, hnks1[i]); 
    ck_assert(mpz_get_ui(a) == hnks1[i]);

  }

  for (i = 0; i < sizeof(hnks2)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(butterfly, 2, i, a) == 0);
    /* gmp_printf("hnk2: %d %Zd %d\n", i, a, hnks2[i]); */
    ck_assert(mpz_get_ui(a) == hnks2[i]);

  }

  for (i = 0; i < sizeof(hnks3)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(butterfly, 3, i, a) == 0);
    /* gmp_printf("hnk3: %d %Zd %d\n", i, a, hnks3[i]); */
    ck_assert(mpz_get_ui(a) == hnks3[i]);

  }

  mpz_clear(a);
  hnk_destroy(butterfly);
}
END_TEST


START_TEST (test_hnk_butterfly_tetrahedron_shell_h3)
{
  hnk_t *butterfly;

  int kmaxs[] = {1, 5, 21, 137, -1};

  /* Verified by hand (complete) */
  int hnks1[] = {1, 1, 4, 6, 4, 1, 0};

  /* Verified by hand (k = 3) */
  int hnks2[] = {1, 1, 4, 22, 94, 323, 956, 2521, 5782, 11195};

  /* Verified by hand (k = 3) */
  int hnks3[] = {1, 1, 4, 22, 210, 1737, 13074, 92257, 615504};

  int i;

  mpz_t a;

  butterfly = hnk_butterfly_tetrahedron_shell_create(3,   
						     100,
						     2);
  ck_assert_ptr_ne(butterfly, NULL);
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {
    /* printf("maxk: %d %d %d\n", i, hnk_get_maxk_at_h(butterfly, i), kmaxs[i]); */
    ck_assert(hnk_get_maxk_at_h(butterfly, i) == kmaxs[i]);
  }

  for (i = 0; i < sizeof(hnks1)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(butterfly, 1, i, a) == 0);
    /* gmp_printf("hnk1: %d %Zd %d\n", i, a, hnks1[i]); */
    ck_assert(mpz_get_ui(a) == hnks1[i]);

  }

  for (i = 0; i < sizeof(hnks2)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(butterfly, 2, i, a) == 0);
    /* gmp_printf("hnk2: %d %Zd %d\n", i, a, hnks2[i]); */
    ck_assert(mpz_get_ui(a) == hnks2[i]);

  }

  for (i = 0; i < sizeof(hnks3)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(butterfly, 3, i, a) == 0);
    /* gmp_printf("hnk3: %d %Zd %d\n", i, a, hnks3[i]); */
    ck_assert(mpz_get_ui(a) == hnks3[i]);

  }

  mpz_clear(a);
  hnk_destroy(butterfly);
}
END_TEST

START_TEST (test_hnk_butterfly_octahedron_shell_h3)
{
  hnk_t *butterfly;

  int kmaxs[] = {1, 7, 37, 265, -1};

  /* Verified by hand (complete) */
  int hnks1[] = {1, 1, 6, 15, 20, 15, 6, 1, 0};

  /* Verified by hand (k = 3) */
  int hnks2[] = {1, 1, 6, 45, 254, 1226, 5322, 20863, 74122, 240924};

  /* Verified by hand (k = 3) */
  int hnks3[] = {1, 1, 6, 45, 482, 4496, 38528, 315841, 2495934};

  int i;

  mpz_t a;

  butterfly = hnk_butterfly_octahedron_shell_create(3,   
						     100,
						     2);
  ck_assert_ptr_ne(butterfly, NULL);
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {
    /* printf("maxk: %d %d %d\n", i, hnk_get_maxk_at_h(butterfly, i), kmaxs[i]); */
    ck_assert(hnk_get_maxk_at_h(butterfly, i) == kmaxs[i]);
  }

  for (i = 0; i < sizeof(hnks1)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(butterfly, 1, i, a) == 0);
    /* gmp_printf("hnk1: %d %Zd %d\n", i, a, hnks1[i]); */
    ck_assert(mpz_get_ui(a) == hnks1[i]);

  }

  for (i = 0; i < sizeof(hnks2)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(butterfly, 2, i, a) == 0);
    /* gmp_printf("hnk2: %d %Zd %d\n", i, a, hnks2[i]); */
    ck_assert(mpz_get_ui(a) == hnks2[i]);

  }

  for (i = 0; i < sizeof(hnks3)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(butterfly, 3, i, a) == 0);
    /* gmp_printf("hnk3: %d %Zd %d\n", i, a, hnks3[i]); */
    ck_assert(mpz_get_ui(a) == hnks3[i]);

  }

  mpz_clear(a);
  hnk_destroy(butterfly);
}
END_TEST

START_TEST (test_hnk_butterfly_icosahedron_shell_h3)
{
  hnk_t *butterfly;

  int kmaxs[] = {1, 7, 37, 265, -1};

  /* Verified by hand (complete) */
  int hnks1[] = {1, 1, 6, 15, 20, 15, 6, 1, 0};

  /* Verified by hand (k = 3) */
  int hnks2[] = {1, 1, 6, 45, 254, 1226, 5322, 20863, 74122, 240924};

  /* Verified by hand (k = 3) */
  int hnks3[] = {1, 1, 6, 45, 482, 4496, 38528, 315841, 2495934};

  int i;

  mpz_t a;

  butterfly = hnk_butterfly_icosahedron_shell_create(4,   
						     100,
						     3);
  ck_assert_ptr_ne(butterfly, NULL);
  mpz_init(a);
  
  for (i = 0; i < sizeof(kmaxs)/sizeof(int); i ++) {
    printf("maxk: %d %d %d\n", i, hnk_get_maxk_at_h(butterfly, i), kmaxs[i]);
    /* ck_assert(hnk_get_maxk_at_h(butterfly, i) == kmaxs[i]); */
  }

  for (i = 0; i < sizeof(hnks1)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(butterfly, 1, i, a) == 0);
    gmp_printf("hnk1: %d %Zd %d\n", i, a, hnks1[i]);
    /* ck_assert(mpz_get_ui(a) == hnks1[i]); */

  }

  for (i = 0; i < sizeof(hnks2)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(butterfly, 2, i, a) == 0);
    gmp_printf("hnk2: %d %Zd %d\n", i, a, hnks2[i]);
    /* ck_assert(mpz_get_ui(a) == hnks2[i]); */

  }

  for (i = 0; i < sizeof(hnks3)/sizeof(int); i ++) {

    ck_assert(hnk_get_hnk(butterfly, 3, i, a) == 0);
    /* gmp_printf("hnk3: %d %Zd %d\n", i, a, hnks3[i]); */
    ck_assert(mpz_get_ui(a) == hnks3[i]);

  }

  mpz_clear(a);
  hnk_destroy(butterfly);
}
END_TEST

Suite *
hnk_butterfly_suite (void)
{
  Suite *s = suite_create ("HNK Butterfly Int");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_hnk_butterfly_create);
  tcase_add_test (tc_core, test_hnk_butterfly_tetrahedron_h3);
  tcase_add_test (tc_core, test_hnk_butterfly_octahedron_h3);
  tcase_add_test (tc_core, test_hnk_butterfly_icosahedron_h3);

  tcase_add_test (tc_core, test_hnk_butterfly_tetrahedron_shell_h3);
  tcase_add_test (tc_core, test_hnk_butterfly_octahedron_shell_h3);
  tcase_add_test (tc_core, test_hnk_butterfly_icosahedron_shell_h3);

  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = hnk_butterfly_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
