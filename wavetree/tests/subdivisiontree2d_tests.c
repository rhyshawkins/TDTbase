
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <check.h>

#include "subdivisiontree2d.h"

START_TEST (test_subdivisiontree2d_ncoefficients)
{
  ck_assert(subdivisiontree2d_total_coefficients(0) == 1);
  ck_assert(subdivisiontree2d_total_coefficients(1) == 5);
  ck_assert(subdivisiontree2d_total_coefficients(2) == 21);
  ck_assert(subdivisiontree2d_total_coefficients(3) == 85);
  ck_assert(subdivisiontree2d_total_coefficients(4) == 341);
  ck_assert(subdivisiontree2d_total_coefficients(5) == 1365);
  ck_assert(subdivisiontree2d_total_coefficients(6) == 5461);
}
END_TEST

START_TEST (test_subdivisiontree2d_2dindices)
{
  subdivisiontree2d_t *s;
  int i;
  int j;
  int ii;
  int ij;
  int rs;

  int index;

  /*
   * Create a 128*128 subdivision grid
   */
  s = subdivisiontree2d_create(7, 1.0, 0.0, SUBDIVISION_BASIS_CONSTANT);
  ck_assert(s != NULL);

  /*
   * Depth 1
   */
  ck_assert(subdivisiontree2d_2dindices(s, 0, &ii, &ij, &rs) == 0);
  ck_assert(ii == 0);
  ck_assert(ij == 0);
  ck_assert(rs == 1);
  ck_assert(subdivisiontree2d_from_2dindices(s, 0, 0, 1) == 0);

  /*
   * Depth 2
   */
  ck_assert(subdivisiontree2d_2dindices(s, 1, &ii, &ij, &rs) == 0);
  ck_assert(ii == 0);
  ck_assert(ij == 0);
  ck_assert(rs == 2);
  ck_assert(subdivisiontree2d_from_2dindices(s, 0, 0, 2) == 1);

  ck_assert(subdivisiontree2d_2dindices(s, 2, &ii, &ij, &rs) == 0);
  ck_assert(ii == 1);
  ck_assert(ij == 0);
  ck_assert(rs == 2);
  ck_assert(subdivisiontree2d_from_2dindices(s, 1, 0, 2) == 2);

  ck_assert(subdivisiontree2d_2dindices(s, 3, &ii, &ij, &rs) == 0);
  ck_assert(ii == 0);
  ck_assert(ij == 1);
  ck_assert(rs == 2);
  ck_assert(subdivisiontree2d_from_2dindices(s, 0, 1, 2) == 3);

  ck_assert(subdivisiontree2d_2dindices(s, 4, &ii, &ij, &rs) == 0);
  ck_assert(ii == 1);
  ck_assert(ij == 1);
  ck_assert(rs == 2);
  ck_assert(subdivisiontree2d_from_2dindices(s, 1, 1, 2) == 4);


  /*
   * Depth 3
   */
  for (j = 0; j < 4; j ++) {
    for (i = 0; i < 4; i ++) {

      index = j*4 + i + 5;
      ck_assert(subdivisiontree2d_2dindices(s, index, &ii, &ij, &rs) == 0);
      ck_assert(ii == i);
      ck_assert(ij == j);
      ck_assert(rs == 4);

      ck_assert(subdivisiontree2d_from_2dindices(s, i, j, 4) == index);
    }
  }

  /*
   * Depth 4
   */
  for (j = 0; j < 8; j ++) {
    for (i = 0; i < 8; i ++) {

      index = j*8 + i + 21;

      ck_assert(subdivisiontree2d_2dindices(s, index, &ii, &ij, &rs) == 0);
      ck_assert(ii == i);
      ck_assert(ij == j);
      ck_assert(rs == 8);

      ck_assert(subdivisiontree2d_from_2dindices(s, i, j, 8) == index);
    }
  }

  subdivisiontree2d_destroy(s);
}
END_TEST

START_TEST (test_subdivisiontree2d_childindices)
{
  subdivisiontree2d_t *s;

  /*
   * Create a 128*128 subdivision grid
   */
  s = subdivisiontree2d_create(7, 1.0, 0.0, SUBDIVISION_BASIS_CONSTANT);
  ck_assert(s != NULL);

  ck_assert(subdivisiontree2d_get_width(s) == 128);
  
#define PARENT_CHILDCHECK(s, p, tl, tr, bl, br) \
  ck_assert(subdivisiontree2d_TL(s, p) == tl); \
  ck_assert(subdivisiontree2d_parent_index(s, tl) == p); \
  ck_assert(subdivisiontree2d_TR(s, p) == tr); \
  ck_assert(subdivisiontree2d_parent_index(s, tr) == p); \
  ck_assert(subdivisiontree2d_BL(s, p) == bl); \
  ck_assert(subdivisiontree2d_parent_index(s, bl) == p); \
  ck_assert(subdivisiontree2d_BR(s, p) == br); \
  ck_assert(subdivisiontree2d_parent_index(s, br) == p);


  /*
   * Depth 1
   */
  PARENT_CHILDCHECK(s, 0, 1, 2, 3, 4);


  /*
   * Depth 2
   */
  PARENT_CHILDCHECK(s, 1, 5, 6, 9, 10);
  PARENT_CHILDCHECK(s, 2, 7, 8, 11, 12);
  PARENT_CHILDCHECK(s, 3, 13, 14, 17, 18);
  PARENT_CHILDCHECK(s, 4, 15, 16, 19, 20);
  
  /*
   * Depth 3
   */
  PARENT_CHILDCHECK(s, 5, 21, 22, 29, 30);
  PARENT_CHILDCHECK(s, 6, 23, 24, 31, 32);
  PARENT_CHILDCHECK(s, 7, 25, 26, 33, 34);
  PARENT_CHILDCHECK(s, 8, 27, 28, 35, 36);

  PARENT_CHILDCHECK(s, 9, 37, 38, 45, 46);
  PARENT_CHILDCHECK(s, 10, 39, 40, 47, 48);
  PARENT_CHILDCHECK(s, 11, 41, 42, 49, 50);
  PARENT_CHILDCHECK(s, 12, 43, 44, 51, 52);

  PARENT_CHILDCHECK(s, 13, 53, 54, 61, 62);
  PARENT_CHILDCHECK(s, 14, 55, 56, 63, 64);
  PARENT_CHILDCHECK(s, 15, 57, 58, 65, 66);
  PARENT_CHILDCHECK(s, 16, 59, 60, 67, 68);
  
  PARENT_CHILDCHECK(s, 17, 69, 70, 77, 78);
  PARENT_CHILDCHECK(s, 18, 71, 72, 79, 80);
  PARENT_CHILDCHECK(s, 19, 73, 74, 81, 82);
  PARENT_CHILDCHECK(s, 20, 75, 76, 83, 84);

  subdivisiontree2d_destroy(s);
}
END_TEST

START_TEST (test_subdivisiontree2d_depth)
{
  subdivisiontree2d_t *s;
  int i;

  /*
   * Create a 128*128 subdivision grid
   */
  s = subdivisiontree2d_create(7, 1.0, 0.0, SUBDIVISION_BASIS_CONSTANT);
  ck_assert(s != NULL);

  ck_assert(subdivisiontree2d_get_width(s) == 128);
  
  ck_assert(subdivisiontree2d_depthofindex(s, 0) == 0);

  /* Indices 1..4 should be depth 1 */
  for (i = 1; i <= 4; i ++) {
    ck_assert(subdivisiontree2d_depthofindex(s, i) == 1);
  }

  /* Indices 5..20 should be depth 2 */
  for (i = 5; i <= 20; i ++) {
    ck_assert(subdivisiontree2d_depthofindex(s, i) == 2);
  }

  /* Indices 21..84 should be depth 3 */
  for (i = 21; i <= 84; i ++) {
    ck_assert(subdivisiontree2d_depthofindex(s, i) == 3);
  }

  /* Indices 85..340 should be depth 4 */
  for (i = 85; i <= 340; i ++) {
    ck_assert(subdivisiontree2d_depthofindex(s, i) == 4);
  }

  subdivisiontree2d_destroy(s);
}
END_TEST

START_TEST(test_subdivisiontree2d_birth)
{
  subdivisiontree2d_t *s;

  /*
   * Create a 128*128 subdivision grid
   */
  s = subdivisiontree2d_create(7, 1.0, 0.0, SUBDIVISION_BASIS_CONSTANT);
  ck_assert(s != NULL);

  ck_assert(subdivisiontree2d_get_width(s) == 128);

  /*
   * Initialization
   */
  ck_assert(subdivisiontree2d_initialize(s, 0.0) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 1);

  /*
   * Add at index 1 and commit
   */
  ck_assert(subdivisiontree2d_propose_birth(s, 1, 1, 1.0) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 2);
  ck_assert(subdivisiontree2d_commit(s) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 2);
  ck_assert(subdivisiontree2d_coeff_count_scan(s) == 2);

  /*
   * Add at index 2 and undo
   */
  ck_assert(subdivisiontree2d_propose_birth(s, 2, 1, 2.0) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 3);
  ck_assert(subdivisiontree2d_undo(s) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 2);
  ck_assert(subdivisiontree2d_coeff_count_scan(s) == 2);

  /*
   * Add at index 3 and commit
   */
  ck_assert(subdivisiontree2d_propose_birth(s, 3, 1, 3.0) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 3);
  ck_assert(subdivisiontree2d_commit(s) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 3);
  ck_assert(subdivisiontree2d_coeff_count_scan(s) == 3);

  subdivisiontree2d_destroy(s);
}
END_TEST

START_TEST(test_subdivisiontree2d_image_mapping)
{
  #define IM_DEGREE 3
  #define IM_W 8
  #define IM_SIZE 64
  subdivisiontree2d_t *s;
  int i;
  int j;

  double img[IM_W*IM_W];
  
  /*
   * Create a 8x8 subdivision grid
   */
  s = subdivisiontree2d_create(IM_DEGREE, 0.0, 0.0, SUBDIVISION_BASIS_CONSTANT);
  ck_assert(s != NULL);

  ck_assert(subdivisiontree2d_get_width(s) == (IM_W));
  ck_assert(subdivisiontree2d_get_size(s) == (IM_SIZE));

  /*
   * Initialization
   */
  ck_assert(subdivisiontree2d_initialize(s, 1.0) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 1);

  /*
   * Map root node to get constant image
   */
  ck_assert(subdivisiontree2d_map_to_array(s, img, IM_SIZE) >= 0);

  for (i = 0; i < (IM_SIZE); i ++) {
    ck_assert(img[i] == 1.0);
  }

  /*
   * Add at index 1 and commit
   */
  ck_assert(subdivisiontree2d_propose_birth(s, 1, 1, 0.5) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 2);
  ck_assert(subdivisiontree2d_commit(s) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 2);
  ck_assert(subdivisiontree2d_coeff_count_scan(s) == 2);
  ck_assert(subdivisiontree2d_valid(s));

  ck_assert(subdivisiontree2d_map_to_array(s, img, IM_SIZE) >= 0);

  /*
   * Check TL quadrant is correct
   */
  for (j = 0; j < (IM_W/2); j ++) {
    for (i = 0; i < (IM_W/2); i ++) {

      ck_assert(img[j*IM_W + i] == 1.5);

      ck_assert(img[j*IM_W + i + IM_W/2] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i ] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i + IM_W/2] == 1.0);
    }
  }

  /*
   * Add at index 5 and commit
   */
  ck_assert(subdivisiontree2d_propose_birth(s, 5, 2, 0.25) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 3);
  ck_assert(subdivisiontree2d_commit(s) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 3);
  ck_assert(subdivisiontree2d_coeff_count_scan(s) == 3);
  ck_assert(subdivisiontree2d_valid(s));

  ck_assert(subdivisiontree2d_map_to_array(s, img, IM_SIZE) >= 0);

  /*
   * Check TL quadrant is correct
   */
  for (j = 0; j < (IM_W/2); j ++) {
    for (i = 0; i < (IM_W/2); i ++) {

      if (j < 2 && i < 2) {
	ck_assert(img[j*IM_W + i] == 1.75);
      } else {
	ck_assert(img[j*IM_W + i] == 1.5);
      }

      ck_assert(img[j*IM_W + i + IM_W/2] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i ] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i + IM_W/2] == 1.0);
    }
  }

  /*
   * Add at index 21 and commit
   */
  ck_assert(subdivisiontree2d_coeff_count(s) == 3);
  ck_assert(subdivisiontree2d_propose_birth(s, 21, 3, 0.125) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 4);
  ck_assert(subdivisiontree2d_commit(s) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 4);
  ck_assert(subdivisiontree2d_coeff_count_scan(s) == 4);
  ck_assert(subdivisiontree2d_valid(s));

  ck_assert(subdivisiontree2d_map_to_array(s, img, IM_SIZE) >= 0);

  /*
   * Check TL quadrant is correct
   */
  for (j = 0; j < (IM_W/2); j ++) {
    for (i = 0; i < (IM_W/2); i ++) {

      if (j < 2 && i < 2) {
	if (j == 0 && i == 0) {
	  ck_assert(img[0] == 1.875);
	} else {
	  ck_assert(img[j*IM_W + i] == 1.75);
	}
      } else {
	ck_assert(img[j*IM_W + i] == 1.5);
      }

      ck_assert(img[j*IM_W + i + IM_W/2] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i ] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i + IM_W/2] == 1.0);
    }
  }

  subdivisiontree2d_destroy(s);
}
END_TEST

START_TEST(test_subdivisiontree2d_incremental_image_mapping)
{
  #define IM_DEGREE 3
  #define IM_W 8
  #define IM_SIZE 64
  subdivisiontree2d_t *s;
  int i;
  int j;

  double img[IM_W*IM_W];
  double old_value;
  
  /*
   * Create a 8x8 subdivision grid
   */
  s = subdivisiontree2d_create(IM_DEGREE, 0.0, 0.0, SUBDIVISION_BASIS_CONSTANT);
  ck_assert(s != NULL);

  ck_assert(subdivisiontree2d_get_width(s) == (IM_W));
  ck_assert(subdivisiontree2d_get_size(s) == (IM_SIZE));

  /*
   * Initialization
   */
  ck_assert(subdivisiontree2d_initialize(s, 1.0) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 1);

  /*
   * Map root node to get constant image
   */
  ck_assert(subdivisiontree2d_map_to_array(s, img, IM_SIZE) >= 0);

  for (i = 0; i < (IM_SIZE); i ++) {
    ck_assert(img[i] == 1.0);
  }

  /*
   * Add at index 1 and commit
   */
  ck_assert(subdivisiontree2d_propose_birth(s, 1, 1, 0.5) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 2);
  ck_assert(subdivisiontree2d_propose_to_array(s, img, IM_SIZE) >= 0);
  ck_assert(subdivisiontree2d_commit(s) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 2);
  ck_assert(subdivisiontree2d_coeff_count_scan(s) == 2);
  ck_assert(subdivisiontree2d_valid(s));

  /*
   * Check TL quadrant is correct
   */
  for (j = 0; j < (IM_W/2); j ++) {
    for (i = 0; i < (IM_W/2); i ++) {

      ck_assert(img[j*IM_W + i] == 1.5);

      ck_assert(img[j*IM_W + i + IM_W/2] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i ] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i + IM_W/2] == 1.0);
    }
  }

  /*
   * Add at index 5 and commit
   */
  ck_assert(subdivisiontree2d_propose_birth(s, 5, 2, 0.25) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 3);
  ck_assert(subdivisiontree2d_propose_to_array(s, img, IM_SIZE) >= 0);
  ck_assert(subdivisiontree2d_commit(s) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 3);
  ck_assert(subdivisiontree2d_coeff_count_scan(s) == 3);
  ck_assert(subdivisiontree2d_valid(s));

  /*
   * Check TL quadrant is correct
   */
  for (j = 0; j < (IM_W/2); j ++) {
    for (i = 0; i < (IM_W/2); i ++) {

      if (j < 2 && i < 2) {
	ck_assert(img[j*IM_W + i] == 1.75);
      } else {
	ck_assert(img[j*IM_W + i] == 1.5);
      }

      ck_assert(img[j*IM_W + i + IM_W/2] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i ] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i + IM_W/2] == 1.0);
    }
  }

  /*
   * Add at index 21 and commit
   */
  ck_assert(subdivisiontree2d_coeff_count(s) == 3);
  ck_assert(subdivisiontree2d_propose_birth(s, 21, 3, 0.125) >= 0);
  ck_assert(subdivisiontree2d_propose_to_array(s, img, IM_SIZE) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 4);
  ck_assert(subdivisiontree2d_commit(s) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 4);
  ck_assert(subdivisiontree2d_coeff_count_scan(s) == 4);
  ck_assert(subdivisiontree2d_valid(s));

  /*
   * Check TL quadrant is correct
   */
  for (j = 0; j < (IM_W/2); j ++) {
    for (i = 0; i < (IM_W/2); i ++) {

      if (j < 2 && i < 2) {
	if (j == 0 && i == 0) {
	  ck_assert(img[0] == 1.875);
	} else {
	  ck_assert(img[j*IM_W + i] == 1.75);
	}
      } else {
	ck_assert(img[j*IM_W + i] == 1.5);
      }

      ck_assert(img[j*IM_W + i + IM_W/2] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i ] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i + IM_W/2] == 1.0);
    }
  }

  /*
   * Remove index 21 and undo
   */
  ck_assert(subdivisiontree2d_propose_death(s, 21, 3, &old_value) >= 0);
  ck_assert(old_value == 0.125);
  ck_assert(subdivisiontree2d_propose_to_array(s, img, IM_SIZE) >= 0);

  for (j = 0; j < (IM_W/2); j ++) {
    for (i = 0; i < (IM_W/2); i ++) {

      if (j < 2 && i < 2) {
	ck_assert(img[j*IM_W + i] == 1.75);
      } else {
	ck_assert(img[j*IM_W + i] == 1.5);
      }

      ck_assert(img[j*IM_W + i + IM_W/2] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i ] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i + IM_W/2] == 1.0);
    }
  }

  ck_assert(subdivisiontree2d_revert_to_array(s, img, IM_SIZE) >= 0);

  for (j = 0; j < (IM_W/2); j ++) {
    for (i = 0; i < (IM_W/2); i ++) {

      if (j < 2 && i < 2) {
	if (j == 0 && i == 0) {
	  ck_assert(img[0] == 1.875);
	} else {
	  ck_assert(img[j*IM_W + i] == 1.75);
	}
      } else {
	ck_assert(img[j*IM_W + i] == 1.5);
      }

      ck_assert(img[j*IM_W + i + IM_W/2] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i ] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i + IM_W/2] == 1.0);
    }
  }

  ck_assert(subdivisiontree2d_undo(s) >= 0);

  /*
   * Change value of index 5 
   */
  ck_assert(subdivisiontree2d_propose_value(s, 5, 2, 0.35) >= 0);
  ck_assert(subdivisiontree2d_propose_to_array(s, img, IM_SIZE) >= 0);

  /*
   * Check TL quadrant is correct
   */
  for (j = 0; j < (IM_W/2); j ++) {
    for (i = 0; i < (IM_W/2); i ++) {

      if (j < 2 && i < 2) {
	if (j == 0 && i == 0) {
	  ck_assert(img[0] == 1.975);
	} else {
	  ck_assert(img[j*IM_W + i] == 1.85);
	}
      } else {
	ck_assert(img[j*IM_W + i] == 1.5);
      }

      ck_assert(img[j*IM_W + i + IM_W/2] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i ] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i + IM_W/2] == 1.0);
    }
  }

  ck_assert(subdivisiontree2d_revert_to_array(s, img, IM_SIZE) >= 0);

  /*
   * Check TL quadrant is correct
   */
  for (j = 0; j < (IM_W/2); j ++) {
    for (i = 0; i < (IM_W/2); i ++) {

      if (j < 2 && i < 2) {
	if (j == 0 && i == 0) {
	  ck_assert(img[0] == 1.875);
	} else {
	  ck_assert(img[j*IM_W + i] == 1.75);
	}
      } else {
	ck_assert(img[j*IM_W + i] == 1.5);
      }

      ck_assert(img[j*IM_W + i + IM_W/2] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i ] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i + IM_W/2] == 1.0);
    }
  }

  ck_assert(subdivisiontree2d_undo(s) >= 0);

  subdivisiontree2d_destroy(s);
}
END_TEST

START_TEST(test_subdivisiontree2d_incremental_image_mapping_with_overlap)
{
  #define IM_DEGREE 3
  #define IM_W 8
  #define IM_SIZE 64
  subdivisiontree2d_t *s;
  int i;
  int j;

  double img[IM_W*IM_W];
  double img_true[IM_W*IM_W];
  double old_value;
  
  /*
   * Create a 8x8 subdivision grid
   */
  s = subdivisiontree2d_create(IM_DEGREE, 0.0, 0.5, SUBDIVISION_BASIS_PYRAMID);
  ck_assert(s != NULL);

  ck_assert(subdivisiontree2d_get_width(s) == (IM_W));
  ck_assert(subdivisiontree2d_get_size(s) == (IM_SIZE));

  /*
   * Initialization
   */
  ck_assert(subdivisiontree2d_initialize(s, 1.0) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 1);

  /*
   * Map root node to get constant image
   */
  ck_assert(subdivisiontree2d_map_to_array(s, img, IM_SIZE) >= 0);

  /*
   * Add at index 1 and commit
   */
  ck_assert(subdivisiontree2d_propose_birth(s, 1, 1, 0.5) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 2);
  ck_assert(subdivisiontree2d_propose_to_array(s, img, IM_SIZE) >= 0);
  ck_assert(subdivisiontree2d_commit(s) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 2);
  ck_assert(subdivisiontree2d_coeff_count_scan(s) == 2);
  ck_assert(subdivisiontree2d_valid(s));

  ck_assert(subdivisiontree2d_map_to_array(s, img_true, IM_SIZE) >= 0);

  /*
   * Check images match to within a small error
   */
  for (i = 0; i < IM_SIZE; i ++) {
    ck_assert(fabs(img[i] - img_true[i]) < 1.0e-14);
  }

  /*
   * Add at index 5 and commit
   */
  ck_assert(subdivisiontree2d_propose_birth(s, 5, 2, 0.25) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 3);
  ck_assert(subdivisiontree2d_propose_to_array(s, img, IM_SIZE) >= 0);
  ck_assert(subdivisiontree2d_commit(s) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 3);
  ck_assert(subdivisiontree2d_coeff_count_scan(s) == 3);
  ck_assert(subdivisiontree2d_valid(s));

  ck_assert(subdivisiontree2d_map_to_array(s, img_true, IM_SIZE) >= 0);

  /*
   * Check images match to within a small error
   */
  for (i = 0; i < IM_SIZE; i ++) {
    ck_assert(fabs(img[i] - img_true[i]) < 1.0e-14);
  }
  
  return;

  /*
   * Add at index 21 and commit
   */
  ck_assert(subdivisiontree2d_coeff_count(s) == 3);
  ck_assert(subdivisiontree2d_propose_birth(s, 21, 3, 0.125) >= 0);
  ck_assert(subdivisiontree2d_propose_to_array(s, img, IM_SIZE) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 4);
  ck_assert(subdivisiontree2d_commit(s) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 4);
  ck_assert(subdivisiontree2d_coeff_count_scan(s) == 4);
  ck_assert(subdivisiontree2d_valid(s));

  /*
   * Check TL quadrant is correct
   */
  for (j = 0; j < (IM_W/2); j ++) {
    for (i = 0; i < (IM_W/2); i ++) {

      if (j < 2 && i < 2) {
	if (j == 0 && i == 0) {
	  ck_assert(img[0] == 1.875);
	} else {
	  ck_assert(img[j*IM_W + i] == 1.75);
	}
      } else {
	ck_assert(img[j*IM_W + i] == 1.5);
      }

      ck_assert(img[j*IM_W + i + IM_W/2] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i ] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i + IM_W/2] == 1.0);
    }
  }

  /*
   * Remove index 21 and undo
   */
  ck_assert(subdivisiontree2d_propose_death(s, 21, 3, &old_value) >= 0);
  ck_assert(old_value == 0.125);
  ck_assert(subdivisiontree2d_propose_to_array(s, img, IM_SIZE) >= 0);

  for (j = 0; j < (IM_W/2); j ++) {
    for (i = 0; i < (IM_W/2); i ++) {

      if (j < 2 && i < 2) {
	ck_assert(img[j*IM_W + i] == 1.75);
      } else {
	ck_assert(img[j*IM_W + i] == 1.5);
      }

      ck_assert(img[j*IM_W + i + IM_W/2] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i ] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i + IM_W/2] == 1.0);
    }
  }

  ck_assert(subdivisiontree2d_revert_to_array(s, img, IM_SIZE) >= 0);

  for (j = 0; j < (IM_W/2); j ++) {
    for (i = 0; i < (IM_W/2); i ++) {

      if (j < 2 && i < 2) {
	if (j == 0 && i == 0) {
	  ck_assert(img[0] == 1.875);
	} else {
	  ck_assert(img[j*IM_W + i] == 1.75);
	}
      } else {
	ck_assert(img[j*IM_W + i] == 1.5);
      }

      ck_assert(img[j*IM_W + i + IM_W/2] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i ] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i + IM_W/2] == 1.0);
    }
  }

  ck_assert(subdivisiontree2d_undo(s) >= 0);

  /*
   * Change value of index 5 
   */
  ck_assert(subdivisiontree2d_propose_value(s, 5, 2, 0.35) >= 0);
  ck_assert(subdivisiontree2d_propose_to_array(s, img, IM_SIZE) >= 0);

  /*
   * Check TL quadrant is correct
   */
  for (j = 0; j < (IM_W/2); j ++) {
    for (i = 0; i < (IM_W/2); i ++) {

      if (j < 2 && i < 2) {
	if (j == 0 && i == 0) {
	  ck_assert(img[0] == 1.975);
	} else {
	  ck_assert(img[j*IM_W + i] == 1.85);
	}
      } else {
	ck_assert(img[j*IM_W + i] == 1.5);
      }

      ck_assert(img[j*IM_W + i + IM_W/2] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i ] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i + IM_W/2] == 1.0);
    }
  }

  ck_assert(subdivisiontree2d_revert_to_array(s, img, IM_SIZE) >= 0);

  /*
   * Check TL quadrant is correct
   */
  for (j = 0; j < (IM_W/2); j ++) {
    for (i = 0; i < (IM_W/2); i ++) {

      if (j < 2 && i < 2) {
	if (j == 0 && i == 0) {
	  ck_assert(img[0] == 1.875);
	} else {
	  ck_assert(img[j*IM_W + i] == 1.75);
	}
      } else {
	ck_assert(img[j*IM_W + i] == 1.5);
      }

      ck_assert(img[j*IM_W + i + IM_W/2] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i ] == 1.0);
      ck_assert(img[(j + IM_W/2)*IM_W + i + IM_W/2] == 1.0);
    }
  }

  ck_assert(subdivisiontree2d_undo(s) >= 0);

  subdivisiontree2d_destroy(s);
}
END_TEST

START_TEST(test_subdivisiontree2d_saveload)
{
  subdivisiontree2d_t *s;
  double c;

  /*
   * Create a 128*128 subdivision grid
   */
  s = subdivisiontree2d_create(7, 0.0, 0.0, SUBDIVISION_BASIS_CONSTANT);
  ck_assert(s != NULL);

  /*
   * Initialization
   */
  ck_assert(subdivisiontree2d_initialize(s, 0.0) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 1);

  /*
   * Add at index 1 and commit
   */
  ck_assert(subdivisiontree2d_propose_birth(s, 1, 1, 1.0) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 2);
  ck_assert(subdivisiontree2d_commit(s) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 2);
  ck_assert(subdivisiontree2d_coeff_count_scan(s) == 2);

  /*
   * Add at index 2 and undo
   */
  ck_assert(subdivisiontree2d_propose_birth(s, 2, 2, 2.0) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 3);
  ck_assert(subdivisiontree2d_commit(s) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 3);
  ck_assert(subdivisiontree2d_coeff_count_scan(s) == 3);

  /*
   * Add at index 3 and commit
   */
  ck_assert(subdivisiontree2d_propose_birth(s, 3, 2, 3.0) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 4);
  ck_assert(subdivisiontree2d_commit(s) >= 0);
  ck_assert(subdivisiontree2d_coeff_count(s) == 4);
  ck_assert(subdivisiontree2d_coeff_count_scan(s) == 4);

  /*
   * Save & Destroy
   */
  ck_assert(subdivisiontree2d_save(s, "subdivisiontree2d_tests_saveload.txt") >= 0);
  subdivisiontree2d_destroy(s);

  /* 
   * Recreate 
   */
  s = subdivisiontree2d_create(7, 0.0, 0.0, SUBDIVISION_BASIS_CONSTANT);
  ck_assert(s != NULL);

  /*
   * Load
   */
  ck_assert(subdivisiontree2d_load(s, "subdivisiontree2d_tests_saveload.txt") >= 0);

  ck_assert(subdivisiontree2d_coeff_count(s) == 4);
  ck_assert(subdivisiontree2d_coeff_count_scan(s) == 4);

  ck_assert(subdivisiontree2d_get_coeff(s, 1, &c) >= 0);
  ck_assert(c == 1.0);

  ck_assert(subdivisiontree2d_get_coeff(s, 2, &c) >= 0);
  ck_assert(c == 2.0);

  ck_assert(subdivisiontree2d_get_coeff(s, 3, &c) >= 0);
  ck_assert(c == 3.0);

  subdivisiontree2d_destroy(s);

  /*
   * Recreate a higher level 
   */
  s = subdivisiontree2d_create(8, 0.0, 0.0, SUBDIVISION_BASIS_CONSTANT);
  ck_assert(s != NULL);

  /*
   * Load
   */
  ck_assert(subdivisiontree2d_load_promote(s, "subdivisiontree2d_tests_saveload.txt") >= 0);

  ck_assert(subdivisiontree2d_coeff_count(s) == 4);
  ck_assert(subdivisiontree2d_coeff_count_scan(s) == 4);

  ck_assert(subdivisiontree2d_get_coeff(s, 1, &c) >= 0);
  ck_assert(c == 1.0);

  ck_assert(subdivisiontree2d_get_coeff(s, 2, &c) >= 0);
  ck_assert(c == 2.0);

  ck_assert(subdivisiontree2d_get_coeff(s, 3, &c) >= 0);
  ck_assert(c == 3.0);

  subdivisiontree2d_destroy(s);
}
END_TEST

Suite *
subdivisiontree2d_suite (void)
{
  Suite *s = suite_create ("Subdivision Tree 2D");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_subdivisiontree2d_ncoefficients);
  tcase_add_test (tc_core, test_subdivisiontree2d_2dindices);
  tcase_add_test (tc_core, test_subdivisiontree2d_childindices);
  tcase_add_test (tc_core, test_subdivisiontree2d_depth);
  tcase_add_test (tc_core, test_subdivisiontree2d_birth);
  tcase_add_test (tc_core, test_subdivisiontree2d_image_mapping);
  tcase_add_test (tc_core, test_subdivisiontree2d_incremental_image_mapping);
  tcase_add_test (tc_core, test_subdivisiontree2d_incremental_image_mapping_with_overlap);
  tcase_add_test (tc_core, test_subdivisiontree2d_saveload);
  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = subdivisiontree2d_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
