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
#include "hnk_cartesian_nonsquare.h"

START_TEST (test_hnk_cartesian_nonsquare_2D_aspect_2_create)
{
  hnk_t *cartesian_nonsquare;

  cartesian_nonsquare = hnk_cartesian_nonsquare_2D_create(4,
							  3,
							  100);
  ck_assert_ptr_ne(cartesian_nonsquare, NULL);
  hnk_destroy(cartesian_nonsquare);

  cartesian_nonsquare = hnk_cartesian_nonsquare_2D_create(3,
							  4,
							  100);
  ck_assert_ptr_ne(cartesian_nonsquare, NULL);
  hnk_destroy(cartesian_nonsquare);
}
END_TEST

START_TEST (test_hnk_cartesian_nonsquare_2D_aspect_4_create)
{
  hnk_t *cartesian_nonsquare;

  cartesian_nonsquare = hnk_cartesian_nonsquare_2D_create(5,
							  3,
							  100);
  ck_assert_ptr_ne(cartesian_nonsquare, NULL);
  hnk_destroy(cartesian_nonsquare);

  cartesian_nonsquare = hnk_cartesian_nonsquare_2D_create(3,
							  5,
							  100);
  ck_assert_ptr_ne(cartesian_nonsquare, NULL);
  hnk_destroy(cartesian_nonsquare);
}
END_TEST

START_TEST (test_hnk_cartesian_nonsquare_2D_aspect_2_maxk)
{
  hnk_t *cartesian_nonsquare;

  cartesian_nonsquare = hnk_cartesian_nonsquare_2D_create(4,
							  3,
							  100);
  ck_assert_ptr_ne(cartesian_nonsquare, NULL);

  /* 4, 3 represents a 16*8 grid, maxk will follow same pattern as square grid until
   * the restriction (at h = 4) is reached  */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 0) == 1);   /*  1 x 1 */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 1) == 4);   /*  2 x 2 */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 2) == 16);  /*  4 x 4 */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 3) == 64);  /*  8 x 8 */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 4) == 128); /* 16 x 8 */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 5) == -1);

  hnk_destroy(cartesian_nonsquare);
}
END_TEST
  
START_TEST (test_hnk_cartesian_nonsquare_2D_aspect_2_hnk)
{
  hnk_t *cartesian_nonsquare;
  mpz_t a;
  int k;
  
  mpz_init(a);
  
  /* cartesian_nonsquare = hnk_cartesian_nonsquare_2D_create(4, */
  /* 							  4, */
  /* 							  100); */
  cartesian_nonsquare = hnk_cartesian_34_create(3,
						100);
  ck_assert_ptr_ne(cartesian_nonsquare, NULL);

  /* 4, 3 represents a 16*8 grid, maxk will follow same pattern as square grid until
   * the restriction (at h = 4) is reached  */
  for (k = 1; k < 64; k ++) {
    ck_assert(hnk_get_hnk(cartesian_nonsquare, 3, k, a) == 0);
    /* gmp_printf("%3d %Zd\n", k, a); */
  }

  mpz_clear(a);
  hnk_destroy(cartesian_nonsquare);
}
END_TEST

START_TEST (test_hnk_cartesian_nonsquare_2D_aspect_4_maxk)
{
  hnk_t *cartesian_nonsquare;

  cartesian_nonsquare = hnk_cartesian_nonsquare_2D_create(5,
							  3,
							  100);
  ck_assert_ptr_ne(cartesian_nonsquare, NULL);

  /* 5, 3 represents a 32*8 grid, maxk will follow same pattern as square grid until
   * the restriction (at h = 3) is reached  */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 0) == 1);   /*  1 x 1 */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 1) == 4);   /*  2 x 2 */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 2) == 16);  /*  4 x 4 */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 3) == 64);  /*  8 x 8 */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 4) == 128); /* 16 x 8 */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 5) == 256); /* 32 x 8 */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 6) == -1);

  hnk_destroy(cartesian_nonsquare);
}
END_TEST

START_TEST (test_hnk_cartesian_nonsquare_2D_sub_aspect_2_maxk)
{
  hnk_t *cartesian_nonsquare;

  cartesian_nonsquare = hnk_cartesian_nonsquare_2D_create_sub(4,
							      3,
							      100);
  ck_assert_ptr_ne(cartesian_nonsquare, NULL);

  /* 4, 3 represents a 16*8 grid, maxk will follow same pattern as square grid until
   * the restriction (at h = 4) is reached  */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 0) == 1);   /*  1 x 1 */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 1) == 3);   /*  2 x 1 + 1 */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 2) == 9);   /*  4 x 2 + 1 */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 3) == 33);  /*  8 x 4 + 1 */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 4) == 129); /* 16 x 8 + 1 */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 5) == -1);

  hnk_destroy(cartesian_nonsquare);
}
END_TEST
  
START_TEST (test_hnk_cartesian_nonsquare_2D_sub_aspect_4_maxk)
{
  hnk_t *cartesian_nonsquare;

  cartesian_nonsquare = hnk_cartesian_nonsquare_2D_create_sub(5,
							      3,
							      100);
  ck_assert_ptr_ne(cartesian_nonsquare, NULL);

  /* 5, 3 represents a 32*8 grid, maxk will follow same pattern as square grid until
   * the restriction (at h = 3) is reached  */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 0) == 1);   /*  1 x 1 */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 1) == 5);   /*  4 x 1 + 1 */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 2) == 17);  /*  8 x 2 + 1 */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 3) == 65);  /* 16 x 4 + 1 */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 4) == 257); /* 32 x 8 + 1 */
  ck_assert(hnk_get_maxk_at_h(cartesian_nonsquare, 5) == -1);

  hnk_destroy(cartesian_nonsquare);
}
END_TEST


/*
 * 3D Non-square
 */

START_TEST (test_hnk_cartesian_nonsquare_3D_aspect_2_create)
{
  hnk_t *h;

  h = hnk_cartesian_nonsquare_3D_create(4, 4, 3, 100);
  ck_assert_ptr_ne(h, NULL);

  hnk_destroy(h);

  h = hnk_cartesian_nonsquare_3D_create(4, 3, 3, 100);
  ck_assert_ptr_ne(h, NULL);

  hnk_destroy(h);

  h = hnk_cartesian_nonsquare_3D_create(4, 3, 2, 100);
  ck_assert_ptr_ne(h, NULL);

  hnk_destroy(h);
  
}
END_TEST

START_TEST (test_hnk_cartesian_nonsquare_3D_aspect_2_maxk)
{
  hnk_t *h;

  h = hnk_cartesian_nonsquare_3D_create(4,
					4,
					3,
					100);
  ck_assert_ptr_ne(h, NULL);

  /* 5, 3 represents a 16*16*8 grid, maxk will follow same pattern as square grid until
   * the restriction (at h = 4) is reached  */
  ck_assert(hnk_get_maxk_at_h(h, 0) == 1);    /*  1x 1x 1 */
  ck_assert(hnk_get_maxk_at_h(h, 1) == 8);    /*  2x 2x 2 */
  ck_assert(hnk_get_maxk_at_h(h, 2) == 64);   /*  4x 4x 4 */
  ck_assert(hnk_get_maxk_at_h(h, 3) == 512);  /*  8x 8x 8 */
  ck_assert(hnk_get_maxk_at_h(h, 4) == 2048); /* 16x16x 8 */
  ck_assert(hnk_get_maxk_at_h(h, 5) == -1);

  hnk_destroy(h);

  h = hnk_cartesian_nonsquare_3D_create(4,
					3,
					3,
					100);
  ck_assert_ptr_ne(h, NULL);

  /* 5, 3 represents a 16*8*8 grid, maxk will follow same pattern as square grid until
   * the restriction (at h = 4) is reached  */
  ck_assert(hnk_get_maxk_at_h(h, 0) == 1);    /*  1x 1x 1 */
  ck_assert(hnk_get_maxk_at_h(h, 1) == 8);    /*  2x 2x 2 */
  ck_assert(hnk_get_maxk_at_h(h, 2) == 64);   /*  4x 4x 4 */
  ck_assert(hnk_get_maxk_at_h(h, 3) == 512);  /*  8x 8x 8 */
  ck_assert(hnk_get_maxk_at_h(h, 4) == 1024); /* 16x 8x 8 */
  ck_assert(hnk_get_maxk_at_h(h, 5) == -1);

  hnk_destroy(h);
}
END_TEST

START_TEST (test_hnk_cartesian_nonsquare_3D_aspect_4_maxk)
{
  hnk_t *h;

  /*
   * This represents a 32 x 16 x 8 grid
   */
  h = hnk_cartesian_nonsquare_3D_create(5,
					4,
					3,
					100);
  ck_assert(h != NULL);

  ck_assert(hnk_get_maxk_at_h(h, 0) == 1);    /*  1 x 1 x 1 */
  ck_assert(hnk_get_maxk_at_h(h, 1) == 8);    /*  2 x 2 x 2 */
  ck_assert(hnk_get_maxk_at_h(h, 2) == 64);   /*  4 x 4 x 4 */
  ck_assert(hnk_get_maxk_at_h(h, 3) == 512);  /*  8 x 8 x 8 */
  ck_assert(hnk_get_maxk_at_h(h, 4) == 2048); /* 16 x16 x 8 */
  ck_assert(hnk_get_maxk_at_h(h, 5) == 4096); /* 32 x16 x 8 */
  ck_assert(hnk_get_maxk_at_h(h, 6) == -1);

  hnk_destroy(h);

  
  /*
   * This represents a 32 x 8 x 8 grid
   */
  h = hnk_cartesian_nonsquare_3D_create(5,
					3,
					3,
					100);
  ck_assert(h != NULL);

  ck_assert(hnk_get_maxk_at_h(h, 0) == 1);    /*  1 x 1 x 1 */
  ck_assert(hnk_get_maxk_at_h(h, 1) == 8);    /*  2 x 2 x 2 */
  ck_assert(hnk_get_maxk_at_h(h, 2) == 64);   /*  4 x 4 x 4 */
  ck_assert(hnk_get_maxk_at_h(h, 3) == 512);  /*  8 x 8 x 8 */
  ck_assert(hnk_get_maxk_at_h(h, 4) == 1024); /* 16 x 8 x 8 */
  ck_assert(hnk_get_maxk_at_h(h, 5) == 2048); /* 32 x 8 x 8 */
  ck_assert(hnk_get_maxk_at_h(h, 6) == -1);

  hnk_destroy(h);
  
  /*
   * This represents a 32 x 32 x 8 grid
   */
  h = hnk_cartesian_nonsquare_3D_create(5,
					5,
					3,
					100);
  ck_assert(h != NULL);

  ck_assert(hnk_get_maxk_at_h(h, 0) == 1);    /*  1 x 1 x 1 */
  ck_assert(hnk_get_maxk_at_h(h, 1) == 8);    /*  2 x 2 x 2 */
  ck_assert(hnk_get_maxk_at_h(h, 2) == 64);   /*  4 x 4 x 4 */
  ck_assert(hnk_get_maxk_at_h(h, 3) == 512);  /*  8 x 8 x 8 */
  ck_assert(hnk_get_maxk_at_h(h, 4) == 2048); /* 16 x16 x 8 */
  ck_assert(hnk_get_maxk_at_h(h, 5) == 8192); /* 32 x32 x 8 */
  ck_assert(hnk_get_maxk_at_h(h, 6) == -1);

  hnk_destroy(h);

  /* 
   * This represents a 256 x 128 x 32 grid
   */
  h = hnk_cartesian_nonsquare_3D_create(8, 7, 5, 100);
  ck_assert(h != NULL);

  ck_assert(hnk_get_maxk_at_h(h, 0) == 1);       /*   1 x   1 x  1 */
  ck_assert(hnk_get_maxk_at_h(h, 1) == 8);       /*   2 x   2 x  2 */
  ck_assert(hnk_get_maxk_at_h(h, 2) == 64);      /*   4 x   4 x  4 */
  ck_assert(hnk_get_maxk_at_h(h, 3) == 512);     /*   8 x   8 x  8 */
  ck_assert(hnk_get_maxk_at_h(h, 4) == 4096);    /*  16 x  16 x 16 */
  ck_assert(hnk_get_maxk_at_h(h, 5) == 32768);   /*  32 x  32 x 32 */
  ck_assert(hnk_get_maxk_at_h(h, 6) == 131072);  /*  64 x  64 x 32 */
  ck_assert(hnk_get_maxk_at_h(h, 7) == 524288);  /* 128 x 128 x 32 */
  ck_assert(hnk_get_maxk_at_h(h, 8) == 1048576); /* 256 x 128 x 32 */
  ck_assert(hnk_get_maxk_at_h(h, 9) == -1);    

  hnk_destroy(h);

}
END_TEST

START_TEST (test_hnk_cartesian_nonsquare_3D_sub_aspect_4_maxk)
{
  hnk_t *h;
  mpz_t a;
  int i;

  mpz_init(a);
  
  /*
   * This represents a 32 x 16 x 8 grid with a starting subdivision of 4 x 2 x 1
   */
  h = hnk_cartesian_nonsquare_3D_create_sub(5,
					    4,
					    3,
					    100);
  ck_assert(h != NULL);

  ck_assert(hnk_get_maxk_at_h(h, 0) == 1);    /*   1 x 1 x 1 */
  ck_assert(hnk_get_maxk_at_h(h, 1) == 9);    /*   4 x 2 x 1 + 1*/
  ck_assert(hnk_get_maxk_at_h(h, 2) == 65);   /*   8 x 4 x 2 + 1*/
  ck_assert(hnk_get_maxk_at_h(h, 3) == 513);  /*  16 x 8 x 4 + 1*/
  ck_assert(hnk_get_maxk_at_h(h, 4) == 4097); /*  32 x16 x 8 + 1*/
  ck_assert(hnk_get_maxk_at_h(h, 5) == -1);

  for (i = 0; i < 10; i ++) {
    ck_assert(hnk_get_hnk(h, 4, i, a) == 0);
  }

  hnk_destroy(h);

  /*
   * This represents a 32 x 8 x 8 grid with a starting subdivision of 4 x 1 x 1
   */
  h = hnk_cartesian_nonsquare_3D_create_sub(5,
					    3,
					    3,
					    100);
  ck_assert(h != NULL);

  ck_assert(hnk_get_maxk_at_h(h, 0) == 1);    /*  1 x 1 x 1 */
  ck_assert(hnk_get_maxk_at_h(h, 1) == 5);    /*  4 x 1 x 1 + 1 */
  ck_assert(hnk_get_maxk_at_h(h, 2) == 33);   /*  8 x 2 x 2 + 1 */
  ck_assert(hnk_get_maxk_at_h(h, 3) == 257);  /* 16 x 4 x 4 + 1 */
  ck_assert(hnk_get_maxk_at_h(h, 4) == 2049); /* 32 x 8 x 8 + 1 */
  ck_assert(hnk_get_maxk_at_h(h, 5) == -1);

  hnk_destroy(h);

  /*
   * This represents a 32 x 32 x 8 grid with a starting subdivision of 4 x 4 x 1
   */
  h = hnk_cartesian_nonsquare_3D_create_sub(5,
					    5,
					    3,
					    100);
  ck_assert(h != NULL);

  ck_assert(hnk_get_maxk_at_h(h, 0) == 1);    /*  1 x  1 x 1 */
  ck_assert(hnk_get_maxk_at_h(h, 1) == 17);   /*  4 x  4 x 1 + 1 */
  ck_assert(hnk_get_maxk_at_h(h, 2) == 129);  /*  8 x  8 x 2 + 1 */
  ck_assert(hnk_get_maxk_at_h(h, 3) == 1025); /* 16 x 16 x 4 + 1 */
  ck_assert(hnk_get_maxk_at_h(h, 4) == 8193); /* 32 x 32 x 8 + 1 */
  ck_assert(hnk_get_maxk_at_h(h, 5) == -1);

  hnk_destroy(h);

  /* 
   * This represents a 256 x 128 x 32 grid with a starting subdivision of 8 x 4 x 1
   */
  h = hnk_cartesian_nonsquare_3D_create_sub(8, 7, 5, 100);
  ck_assert(h != NULL);

  ck_assert(hnk_get_maxk_at_h(h, 0) == 1);       /*   1 x   1 x  1 */
  ck_assert(hnk_get_maxk_at_h(h, 1) == 33);      /*   8 x   4 x  1 + 1 */
  ck_assert(hnk_get_maxk_at_h(h, 2) == 257);     /*  16 x   8 x  2 + 1 */
  ck_assert(hnk_get_maxk_at_h(h, 3) == 2049);    /*  32 x  16 x  4 + 1 */
  ck_assert(hnk_get_maxk_at_h(h, 4) == 16385);   /*  64 x  32 x  8 + 1 */
  ck_assert(hnk_get_maxk_at_h(h, 5) == 131073);  /* 128 x  64 x 16 + 1 */
  ck_assert(hnk_get_maxk_at_h(h, 6) == 1048577); /* 258 x 128 x 32 + 1 */
  ck_assert(hnk_get_maxk_at_h(h, 7) == -1);    

  hnk_destroy(h);

  mpz_clear(a);
}
END_TEST

START_TEST (test_hnk_cartesian_nonsquare_3D_spectral)
{
  hnk_t *h;
  int i;
  double ratio;
  mpz_t a;
  
  h = hnk_cartesian_nonsquare_3D_create_spectral(4,
						 4,
						 1000);
  ck_assert(h != NULL);

  ck_assert(hnk_get_maxk_at_h(h, 0) == 1);   
  ck_assert(hnk_get_maxk_at_h(h, 1) == 5);   
  ck_assert(hnk_get_maxk_at_h(h, 2) == 37);  
  ck_assert(hnk_get_maxk_at_h(h, 3) == 293); 
  ck_assert(hnk_get_maxk_at_h(h, 4) == 2341);
  ck_assert(hnk_get_maxk_at_h(h, 5) == -1);    

  mpz_init(a);

  for (i = 0; i < 50; i ++) {
    ck_assert(hnk_get_hnk(h, 2, i, a) >= 0);
    gmp_printf("%3d %Zd\n", i, a);
  }
    
  
  for (i = 0; i < 100; i ++) {
    ck_assert(hnk_get_kplus1_ratio(h, 4, i, &ratio) >= 0);
    
    /* printf("%2d %f\n", i, ratio); */
  }
  
  hnk_destroy(h);
  mpz_clear(a);
}
END_TEST

Suite *
hnk_cartesian_nonsquare_suite (void)
{
  Suite *s = suite_create ("HNK Cartesian Nonsquare");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_hnk_cartesian_nonsquare_2D_aspect_2_create);
  tcase_add_test (tc_core, test_hnk_cartesian_nonsquare_2D_aspect_4_create);
  tcase_add_test (tc_core, test_hnk_cartesian_nonsquare_2D_aspect_2_maxk);
  tcase_add_test (tc_core, test_hnk_cartesian_nonsquare_2D_aspect_2_hnk);
  tcase_add_test (tc_core, test_hnk_cartesian_nonsquare_2D_aspect_4_maxk);

  tcase_add_test (tc_core, test_hnk_cartesian_nonsquare_2D_sub_aspect_2_maxk);
  tcase_add_test (tc_core, test_hnk_cartesian_nonsquare_2D_sub_aspect_4_maxk);

  tcase_add_test (tc_core, test_hnk_cartesian_nonsquare_3D_aspect_2_create);
  tcase_add_test (tc_core, test_hnk_cartesian_nonsquare_3D_aspect_2_maxk);
  tcase_add_test (tc_core, test_hnk_cartesian_nonsquare_3D_aspect_4_maxk);

  tcase_add_test (tc_core, test_hnk_cartesian_nonsquare_3D_sub_aspect_4_maxk);

  tcase_add_test (tc_core, test_hnk_cartesian_nonsquare_3D_spectral);
  
  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = hnk_cartesian_nonsquare_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
