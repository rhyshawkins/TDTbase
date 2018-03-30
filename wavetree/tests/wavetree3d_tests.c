
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <check.h>

#include "wavetree3d.h"

START_TEST (test_wavetree3d_ncoefficients)
{
  wavetree3d_t *s;
  int i;
  
  for (i = 0; i <= 8; i ++) {
    s = wavetree3d_create(i, i, i, 0.0);
    ck_assert(s != NULL);
    ck_assert(wavetree3d_get_size(s) == 1 << (3*i));
    wavetree3d_destroy(s);
  }
}
END_TEST

START_TEST (test_wavetree3d_3dindices)
{
  wavetree3d_t *s;
  int i;
  int j;
  int k;
  int ii;
  int ij;
  int ik;

  int index;
  int width;

  /*
   * Create a 32x32x32 wave grid
   */
  s = wavetree3d_create(5, 5, 5, 0.0);
  ck_assert(s != NULL);

  width = wavetree3d_get_width(s);
  ck_assert(width == 32);

  for (k = 0; k < width; k ++) {
    for (j = 0; j < width; j ++) {
      for (i = 0; i < width; i ++) {
	
	index = (k*width + j)*width + i;
	ck_assert(wavetree3d_3dindices(s, index, &ii, &ij, &ik) == 0);
	ck_assert(ii == i);
	ck_assert(ij == j);
	ck_assert(ik == k);
	
	ck_assert(wavetree3d_from_3dindices(s, i, j, k) == index);
      }
    }
  }

  wavetree3d_destroy(s);
}
END_TEST

START_TEST (test_wavetree3d_childindices)
{
  wavetree3d_t *s;
  int width;
  int rowstride;
  int depthstride;

  /*
   * Create a 32x32x32 wave grid
   */
  s = wavetree3d_create(5, 5, 5, 0.0);
  ck_assert(s != NULL);

  width = wavetree3d_get_width(s);
  ck_assert(width == 32);

  rowstride = width;
  depthstride = width*width;
  
#define PARENT_CHILDCHECK(s, p, utl, utr, ubl, ubr, ltl, ltr, lbl, lbr)	\
  ck_assert(wavetree3d_UTL(s, p) == (utl));				\
  if (utl > 0) {							\
    ck_assert(wavetree3d_parent_index(s, (utl)) == p);			\
  }									\
  ck_assert(wavetree3d_UTR(s, p) == (utr));				\
  ck_assert(wavetree3d_parent_index(s, (utr)) == p);			\
  ck_assert(wavetree3d_UBL(s, p) == (ubl));				\
  ck_assert(wavetree3d_parent_index(s, (ubl)) == p);			\
  ck_assert(wavetree3d_UBR(s, p) == (ubr));				\
  ck_assert(wavetree3d_parent_index(s, (ubr)) == p);			\
  ck_assert(wavetree3d_LTL(s, p) == (ltl));				\
  ck_assert(wavetree3d_parent_index(s, (ltl)) == p);			\
  ck_assert(wavetree3d_LTR(s, p) == (ltr));				\
  ck_assert(wavetree3d_parent_index(s, (ltr)) == p);			\
  ck_assert(wavetree3d_LBL(s, p) == (lbl));				\
  ck_assert(wavetree3d_parent_index(s, (lbl)) == p);			\
  ck_assert(wavetree3d_LBR(s, p) == (lbr));				\
  ck_assert(wavetree3d_parent_index(s, (lbr)) == p);			\


  /*
   * Depth 1
   */
  PARENT_CHILDCHECK(s, 0, 
		    -1, 1, rowstride, rowstride + 1, 
		    depthstride, depthstride + 1, depthstride + rowstride, depthstride + rowstride + 1);


  /*
   * Depth 2 (check each orthogonal direction)
   */
  /* 1, 0, 0 */
  PARENT_CHILDCHECK(s, 1,
		    2, 3, rowstride + 2, rowstride + 3,
		    depthstride + 2, depthstride + 3, depthstride + rowstride + 2, depthstride + rowstride + 3);

  /* 0, 1, 0 */
  PARENT_CHILDCHECK(s, rowstride,
		    2*rowstride, 2*rowstride + 1, 3*rowstride, 3*rowstride + 1,
		    depthstride + 2*rowstride, depthstride + 2*rowstride + 1, 
		    depthstride + 3*rowstride, depthstride + 3*rowstride + 1);

  /* 0, 0, 1 */
  PARENT_CHILDCHECK(s, depthstride, 
		    2*depthstride, 2*depthstride + 1, 2*depthstride + rowstride, 2*depthstride + rowstride + 1,
		    3*depthstride, 3*depthstride + 1, 3*depthstride + rowstride, 3*depthstride + rowstride + 1);

  wavetree3d_destroy(s);
}
END_TEST

START_TEST (test_wavetree3d_depth)
{
  wavetree3d_t *s;
  int i;
  int j;
  int width;
  
  /*
   * Create a 32x32x32 wave grid
   */
  s = wavetree3d_create(5, 5, 5, 0.0);
  ck_assert(s != NULL);

  width = wavetree3d_get_width(s);
  ck_assert(width == 32);
  
  ck_assert(wavetree3d_depthofindex(s, 0) == 0);

  /* Indices 1..4 should be depth 1 */
  for (i = 1, j = 1; i < width; i *= 2, j ++) {
    /*
     * Along width
     */
    ck_assert(wavetree3d_depthofindex(s, i) == j);

    /*
     * Along height
     */
    ck_assert(wavetree3d_depthofindex(s, i*width) == j);

    /*
     * Along depth
     */
    ck_assert(wavetree3d_depthofindex(s, i*width*width) == j);

    /* 
     * Diagonal
     */
    ck_assert(wavetree3d_depthofindex(s, (i*width + i)*width + i) == j);
  }

  wavetree3d_destroy(s);
}
END_TEST

START_TEST(test_wavetree3d_birth)
{
  wavetree3d_t *s;
  char buffer[256];
  uint64_t binary;

  /*
   * Create a 32x32x32 wave grid
   */
  s = wavetree3d_create(5, 5, 5, 0.0);
  ck_assert(s != NULL);

  /*
   * Initialization
   */
  ck_assert(wavetree3d_initialize(s, 0.0) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 1);

  /*
   * Add at index 1 and commit
   */
  ck_assert(wavetree3d_propose_birth(s, 1, 1, 1.0) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 2);
  ck_assert(wavetree3d_commit(s) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 2);

  /*
   * Add at index 2 and undo
   */
  ck_assert(wavetree3d_propose_birth(s, 2, 2, 2.0) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 3);
  ck_assert(wavetree3d_undo(s) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 2);

  /*
   * Add at index 3 and commit
   */
  ck_assert(wavetree3d_propose_birth(s, 3, 2, 3.0) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 3);
  ck_assert(wavetree3d_commit(s) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 3);

  /*
   * Check Dyck word is correct (we use a pseudo Dyck with "." == "()" to
   * aid readability.
   *
   * 0 = (.......)
   * 0, 1 = ((........)......)
   * 0, 1, 2 = (((........).......)......)
   * 0, 1, 3 = ((.(........)......)......)
   */
  ck_assert(wavetree3d_generate_dyck_word(s, buffer, 256) == 0);
  ck_assert(strcmp(buffer, "((.(........)......)......)") == 0);
    
  /* Also check binary version.
   * 0, 1, 3 in Reverse binary:
   * 0001 0010 1010 1010 1010 1101 0101 0101 0110 1010 1010 1011
   *    8    4    5    5    5    b    a    a    6    5    5    d
   */
  ck_assert(wavetree3d_generate_dyck_binary(s, &binary) == 0);
  ck_assert(binary == 0xd556aab55548L);

  wavetree3d_destroy(s);
}
END_TEST

START_TEST(test_wavetree3d_death)
{
  wavetree3d_t *s;
  double old_value;
  /*
   * Create a 32x32x32 wave grid
   */
  s = wavetree3d_create(5, 5, 5, 0.0);
  ck_assert(s != NULL);

  /*
   * Initialization
   */
  ck_assert(wavetree3d_initialize(s, 0.0) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 1);

  /*
   * Add at index 1 and commit
   */
  ck_assert(wavetree3d_propose_birth(s, 1, 1, 1.0) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 2);
  ck_assert(wavetree3d_commit(s) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 2);

  /*
   * Add at index 2 and commit
   */
  ck_assert(wavetree3d_propose_birth(s, 2, 2, 2.0) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 3);
  ck_assert(wavetree3d_commit(s) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 3);

  /*
   * Add at index 3 and commit
   */
  ck_assert(wavetree3d_propose_birth(s, 3, 2, 3.0) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 4);
  ck_assert(wavetree3d_commit(s) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 4);

  /*
   * Remove index 3 and undo 
   */
  ck_assert(wavetree3d_propose_death(s, 3, 2, &old_value) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 3);
  ck_assert(old_value == 3.0);
  ck_assert(wavetree3d_undo(s) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 4);

  /*
   * Remove index 2 and commit
   */
  ck_assert(wavetree3d_propose_death(s, 2, 2, &old_value) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 3);
  ck_assert(old_value == 2.0);
  ck_assert(wavetree3d_commit(s) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 3);

  wavetree3d_destroy(s);
}
END_TEST

START_TEST(test_wavetree3d_value)
{
  wavetree3d_t *s;
  double cvalue;

  /*
   * Create a 32x32x32 wave grid
   */
  s = wavetree3d_create(5, 5, 5, 0.0);
  ck_assert(s != NULL);

  /*
   * Initialization
   */
  ck_assert(wavetree3d_initialize(s, 0.0) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 1);

  /*
   * Add at index 1 and commit
   */
  ck_assert(wavetree3d_propose_birth(s, 1, 1, 1.0) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 2);
  ck_assert(wavetree3d_commit(s) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 2);

  /*
   * Add at index 2 and commit
   */
  ck_assert(wavetree3d_propose_birth(s, 2, 2, 2.0) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 3);
  ck_assert(wavetree3d_commit(s) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 3);

  /*
   * Add at index 3 and commit
   */
  ck_assert(wavetree3d_propose_birth(s, 3, 2, 3.0) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 4);
  ck_assert(wavetree3d_commit(s) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 4);

  /*
   * Change index 3 and undo 
   */
  ck_assert(wavetree3d_propose_value(s, 3, 2, 3.5) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 4);
  ck_assert(wavetree3d_get_coeff(s, 3, 2, &cvalue) >= 0);
  ck_assert(cvalue == 3.5);
  ck_assert(wavetree3d_undo(s) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 4);
  ck_assert(wavetree3d_get_coeff(s, 3, 2, &cvalue) >= 0);
  ck_assert(cvalue == 3.0);

  /*
   * Remove index 2 and commit
   */
  ck_assert(wavetree3d_propose_value(s, 2, 2, 2.5) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 4);
  ck_assert(wavetree3d_get_coeff(s, 2, 2, &cvalue) >= 0);
  ck_assert(cvalue == 2.5);
  ck_assert(wavetree3d_commit(s) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 4);
  ck_assert(wavetree3d_get_coeff(s, 2, 2, &cvalue) >= 0);
  ck_assert(cvalue == 2.5);

  wavetree3d_destroy(s);
}
END_TEST

START_TEST(test_wavetree3d_image_mapping)
{
  #define IM_DEGREE 4
  #define IM_W 16
  #define IM_SIZE 4096
  wavetree3d_t *s;
  int i;
  int width;

  double img[IM_W*IM_W*IM_W];
  
  /*
   * Create a 8x8 wave grid
   */
  s = wavetree3d_create(IM_DEGREE, IM_DEGREE, IM_DEGREE, 0.0);
  ck_assert(s != NULL);

  width = wavetree3d_get_width(s);
  ck_assert(width == 16);

  ck_assert(wavetree3d_get_size(s) == (IM_SIZE));

  /*
   * Initialization
   */
  ck_assert(wavetree3d_initialize(s, 1.0) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 1);

  /*
   * Map root node to get constant image
   */
  memset(img, 0, sizeof(img));
  ck_assert(wavetree3d_map_to_array(s, img, IM_SIZE) >= 0);

  ck_assert(img[0] == 1.0);
  for (i = 1; i < (IM_SIZE); i ++) {
    ck_assert(img[i] == 0.0);
  }

  /*
   * Add at index 1 and commit
   */
  ck_assert(wavetree3d_propose_birth(s, 1, 1, 0.5) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 2);
  ck_assert(wavetree3d_commit(s) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 2);
  ck_assert(wavetree3d_valid(s));

  memset(img, 0, sizeof(img));
  ck_assert(wavetree3d_map_to_array(s, img, IM_SIZE) >= 0);

  /*
   * Check map is correct
   */
  ck_assert(img[0] == 1.0);
  ck_assert(img[1] == 0.5);
  for (i = 2; i < (IM_SIZE); i ++) {
    ck_assert(img[i] == 0.0);
  }

  /*
   * Add at index width and commit
   */
  ck_assert(wavetree3d_propose_birth(s, width, 1, 0.25) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 3);
  ck_assert(wavetree3d_commit(s) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 3);
  ck_assert(wavetree3d_valid(s));

  memset(img, 0, sizeof(img));
  ck_assert(wavetree3d_map_to_array(s, img, IM_SIZE) >= 0);

  /*
   * Check TL quadrant is correct
   */
  ck_assert(img[0] == 1.0);
  ck_assert(img[1] == 0.5);
  ck_assert(img[IM_W] == 0.25);
  for (i = 2; i < (IM_W); i ++) {
    ck_assert(img[i] == 0.0);
  }
  for (i = IM_W + 1; i < (IM_SIZE); i ++) {
    ck_assert(img[i] == 0.0);
  }
  
  wavetree3d_destroy(s);
}
END_TEST

START_TEST(test_wavetree3d_saveload)
{
  wavetree3d_t *s;
  double c;

  /*
   * Create a 32x32x32 wave grid
   */
  s = wavetree3d_create(5, 5, 5, 0.0);
  ck_assert(s != NULL);

  /*
   * Initialization
   */
  ck_assert(wavetree3d_initialize(s, 0.0) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 1);

  /*
   * Add at index 1 and commit
   */
  ck_assert(wavetree3d_propose_birth(s, 1, 1, 1.0) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 2);
  ck_assert(wavetree3d_commit(s) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 2);

  /*
   * Add at index 2 and commit
   */
  ck_assert(wavetree3d_propose_birth(s, 2, 2, 2.0) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 3);
  ck_assert(wavetree3d_commit(s) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 3);

  /*
   * Add at index 3 and commit
   */
  ck_assert(wavetree3d_propose_birth(s, 3, 2, 3.0) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 4);
  ck_assert(wavetree3d_commit(s) >= 0);
  ck_assert(wavetree3d_coeff_count(s) == 4);

  /*
   * Save & Destroy
   */
  ck_assert(wavetree3d_save(s, "wavetree3d_tests_saveload.txt") >= 0);
  wavetree3d_destroy(s);

  /* 
   * Recreate 
   */
  s = wavetree3d_create(5, 5, 5, 0.0);
  ck_assert(s != NULL);

  /*
   * Load
   */
  ck_assert(wavetree3d_load(s, "wavetree3d_tests_saveload.txt") >= 0);

  ck_assert(wavetree3d_coeff_count(s) == 4);

  ck_assert(wavetree3d_get_coeff(s, 1, 1, &c) >= 0);
  ck_assert(c == 1.0);

  ck_assert(wavetree3d_get_coeff(s, 2, 2, &c) >= 0);
  ck_assert(c == 2.0);

  ck_assert(wavetree3d_get_coeff(s, 3, 2, &c) >= 0);
  ck_assert(c == 3.0);

  wavetree3d_destroy(s);

#if 0
  /*
   * Recreate a higher level 
   */
  s = wavetree3d_create(5, 5, 5, 0.0);
  ck_assert(s != NULL);

  /*
   * Load
   */
  ck_assert(wavetree3d_load_promote(s, "wavetree3d_tests_saveload.txt") >= 0);

  ck_assert(wavetree3d_coeff_count(s) == 4);

  ck_assert(wavetree3d_get_coeff(s, 1, &c) >= 0);
  ck_assert(c == 1.0);

  ck_assert(wavetree3d_get_coeff(s, 2, &c) >= 0);
  ck_assert(c == 2.0);

  ck_assert(wavetree3d_get_coeff(s, 3, &c) >= 0);
  ck_assert(c == 3.0);

  wavetree3d_destroy(s);
#endif 
}
END_TEST

START_TEST(test_wavetree3d_nonsquare)
{
  wavetree3d_t *s;
  int rowstride;
  int slicestride;

  int ii, ij, ik;

  s = wavetree3d_create(5, 4, 3, 0.0);
  ck_assert(s != NULL);

  ck_assert(wavetree3d_get_width(s) == 32);
  ck_assert(wavetree3d_get_height(s) == 16);
  ck_assert(wavetree3d_get_depth(s) == 8);
  ck_assert(wavetree3d_get_size(s) == (32*16*8));

  rowstride = 32;
  slicestride = 32*16;
    
  /*
   * Check height terminates correctly
   */
  ck_assert(wavetree3d_UBL(s, 0) == rowstride);

  ck_assert(wavetree3d_3dindices(s, rowstride, &ii, &ij, &ik) >= 0);
  ck_assert(ii == 0);
  ck_assert(ij == 1);
  ck_assert(ik == 0);
  
  ck_assert(wavetree3d_UTL(s, rowstride) == (rowstride * 2));
  ck_assert(wavetree3d_UBL(s, rowstride) == (rowstride * 3));
  
  ck_assert(wavetree3d_UTL(s, 2*rowstride) == (rowstride * 4));
  ck_assert(wavetree3d_UBL(s, 2*rowstride) == (rowstride * 5));

  ck_assert(wavetree3d_UTL(s, 3*rowstride) == (rowstride * 6));
  ck_assert(wavetree3d_UBL(s, 3*rowstride) == (rowstride * 7));

  ck_assert(wavetree3d_UTL(s, 4*rowstride) == (rowstride * 8));
  ck_assert(wavetree3d_UBL(s, 4*rowstride) == (rowstride * 9));

  ck_assert(wavetree3d_UTL(s, 5*rowstride) == (rowstride * 10));
  ck_assert(wavetree3d_UBL(s, 5*rowstride) == (rowstride * 11));

  ck_assert(wavetree3d_UTL(s, 6*rowstride) == (rowstride * 12));
  ck_assert(wavetree3d_UBL(s, 6*rowstride) == (rowstride * 13));

  ck_assert(wavetree3d_UTL(s, 7*rowstride) == (rowstride * 14));
  ck_assert(wavetree3d_UBL(s, 7*rowstride) == (rowstride * 15));

  ck_assert(wavetree3d_UTL(s,  8*rowstride) == (-1));
  ck_assert(wavetree3d_UTL(s,  9*rowstride) == (-1));
  ck_assert(wavetree3d_UTL(s, 10*rowstride) == (-1));
  ck_assert(wavetree3d_UTL(s, 11*rowstride) == (-1));
  ck_assert(wavetree3d_UTL(s, 12*rowstride) == (-1));
  ck_assert(wavetree3d_UTL(s, 13*rowstride) == (-1));
  ck_assert(wavetree3d_UTL(s, 14*rowstride) == (-1));
  ck_assert(wavetree3d_UTL(s, 15*rowstride) == (-1));
  
  /*
   * Check depth terminates correctly
   */
  ck_assert(wavetree3d_LTL(s, 0) == slicestride);

  ck_assert(wavetree3d_3dindices(s, slicestride, &ii, &ij, &ik) >= 0);
  ck_assert(ii == 0);
  ck_assert(ij == 0);
  ck_assert(ik == 1);

  ck_assert(wavetree3d_UTL(s, slicestride) == (slicestride * 2));
  ck_assert(wavetree3d_LTL(s, slicestride) == (slicestride * 3));

  ck_assert(wavetree3d_UTL(s, 2*slicestride) == (slicestride * 4));
  ck_assert(wavetree3d_LTL(s, 2*slicestride) == (slicestride * 5));
  
  ck_assert(wavetree3d_UTL(s, 3*slicestride) == (slicestride * 6));
  ck_assert(wavetree3d_LTL(s, 3*slicestride) == (slicestride * 7));

  ck_assert(wavetree3d_UTL(s, 4*slicestride) == (-1));
  ck_assert(wavetree3d_UTL(s, 5*slicestride) == (-1));
  ck_assert(wavetree3d_UTL(s, 6*slicestride) == (-1));
  ck_assert(wavetree3d_UTL(s, 7*slicestride) == (-1));
  
  ck_assert(wavetree3d_LTL(s, 4*slicestride) == (-1));
  ck_assert(wavetree3d_LTL(s, 5*slicestride) == (-1));
  ck_assert(wavetree3d_LTL(s, 6*slicestride) == (-1));
  ck_assert(wavetree3d_LTL(s, 7*slicestride) == (-1));

  wavetree3d_destroy(s);
}
END_TEST

void recurse_indices(wavetree3d_t *s, int index, int *indices)
{
  if (index < 0) {
    return;
  }

  ck_assert(indices[index] == 0);
  indices[index] = 1;

  recurse_indices(s, wavetree3d_UTL(s, index), indices);
  recurse_indices(s, wavetree3d_UTR(s, index), indices);
  recurse_indices(s, wavetree3d_UBL(s, index), indices);
  recurse_indices(s, wavetree3d_UBR(s, index), indices);

  recurse_indices(s, wavetree3d_LTL(s, index), indices);
  recurse_indices(s, wavetree3d_LTR(s, index), indices);
  recurse_indices(s, wavetree3d_LBL(s, index), indices);
  recurse_indices(s, wavetree3d_LBR(s, index), indices);
}
START_TEST(test_wavetree3d_nonsquare_coverage)
{
  wavetree3d_t *s;
  int *indices;
  int size;
  int i;
  
  s = wavetree3d_create(8, 7, 5, 0.0);
  ck_assert(s != NULL);

  ck_assert(wavetree3d_get_width(s) == 256);
  ck_assert(wavetree3d_get_height(s) == 128);
  ck_assert(wavetree3d_get_depth(s) == 32);

  size = 256 * 128 * 32;
  ck_assert(wavetree3d_get_size(s) == size);

  indices = malloc(sizeof(int) * size);
  ck_assert(indices != NULL);

  memset(indices, 0, sizeof(int) * size);

  recurse_indices(s, 0, indices);

  for (i = 0; i < size; i ++) {
    ck_assert(indices[i] == 1);
  }

  free(indices);
  wavetree3d_destroy(s);
}
END_TEST
     


Suite *
wavetree3d_suite (void)
{
  Suite *s = suite_create ("Wave Tree 3D");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_wavetree3d_ncoefficients);
  tcase_add_test (tc_core, test_wavetree3d_3dindices);
  tcase_add_test (tc_core, test_wavetree3d_childindices); 
  tcase_add_test (tc_core, test_wavetree3d_depth);
  tcase_add_test (tc_core, test_wavetree3d_birth);
  tcase_add_test (tc_core, test_wavetree3d_value);
  tcase_add_test (tc_core, test_wavetree3d_death);
  tcase_add_test (tc_core, test_wavetree3d_image_mapping);
  tcase_add_test (tc_core, test_wavetree3d_saveload);

  tcase_add_test (tc_core, test_wavetree3d_nonsquare);
  tcase_add_test (tc_core, test_wavetree3d_nonsquare_coverage);

  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = wavetree3d_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
