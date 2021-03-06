
#include <stdio.h>
#include <stdlib.h>
#include <check.h>

#include "wavetree2d_sub.h"

START_TEST (test_wavetree2d_sub_ncoefficients)
{
  wavetree2d_sub_t *s;
  int i;
  
  for (i = 0; i <= 8; i ++) {
    s = wavetree2d_sub_create(i, i, 0.0);
    ck_assert(s != NULL);
    ck_assert(wavetree2d_sub_get_size(s) == 1 << (2*i));
    wavetree2d_sub_destroy(s);
  }
}
END_TEST

START_TEST (test_wavetree2d_sub_2dindices)
{
  wavetree2d_sub_t *s;
  int i;
  int j;
  int ii;
  int ij;

  int index;
  int width;
  int height;
  
  /*
   * Create a 128*128 wave grid
   */
  s = wavetree2d_sub_create(7, 7, 0.0);
  ck_assert(s != NULL);

  width = wavetree2d_sub_get_width(s);
  ck_assert(width == 128);

  height = wavetree2d_sub_get_height(s);
  ck_assert(height == 128);

  for (j = 0; j < height; j ++) {
    for (i = 0; i < width; i ++) {

      index = j*width + i;
      ck_assert(wavetree2d_sub_2dindices(s, index, &ii, &ij) == 0);
      ck_assert(ii == i);
      ck_assert(ij == j);

      ck_assert(wavetree2d_sub_from_2dindices(s, i, j) == index);
    }
  }

  wavetree2d_sub_destroy(s);
}
END_TEST

START_TEST (test_wavetree2d_sub_2dindices_rectangle)
{
  wavetree2d_sub_t *s;
  int i;
  int j;
  int ii;
  int ij;

  int index;
  int width;
  int height;
  
  /*
   * Create a 128*128 wave grid
   */
  s = wavetree2d_sub_create(7, 5, 0.0);
  ck_assert(s != NULL);

  width = wavetree2d_sub_get_width(s);
  ck_assert(width == 128);

  height = wavetree2d_sub_get_height(s);
  ck_assert(height == 32);

  for (j = 0; j < height; j ++) {
    for (i = 0; i < width; i ++) {

      index = j*width + i + 1;
      ck_assert(wavetree2d_sub_2dindices(s, index, &ii, &ij) == 0);

      ck_assert(ii == i);
      ck_assert(ij == j);

      ck_assert(wavetree2d_sub_from_2dindices(s, i, j) == index);
    }
  }

  wavetree2d_sub_destroy(s);
}
END_TEST

START_TEST (test_wavetree2d_sub_childindices)
{
  wavetree2d_sub_t *s;
  int width;

  /*
   * Create a 128*128 wave grid
   */
  s = wavetree2d_sub_create(7, 7, 0.0);
  ck_assert(s != NULL);

  width = wavetree2d_sub_get_width(s);
  ck_assert(width == 128);
  
#define PARENT_CHILDCHECK(s, p, tl, tr, bl, br) \
  ck_assert(wavetree2d_sub_TL(s, p) == tl); \
  ck_assert(wavetree2d_sub_parent_index(s, tl) == p); \
  ck_assert(wavetree2d_sub_TR(s, p) == tr); \
  ck_assert(wavetree2d_sub_parent_index(s, tr) == p); \
  ck_assert(wavetree2d_sub_BL(s, p) == bl); \
  ck_assert(wavetree2d_sub_parent_index(s, bl) == p); \
  ck_assert(wavetree2d_sub_BR(s, p) == br); \
  ck_assert(wavetree2d_sub_parent_index(s, br) == p);


  /*
   * Depth 1
   */
  ck_assert(wavetree2d_sub_TL(s, 0) == -1);
  ck_assert(wavetree2d_sub_TR(s, 0) == 1);
  ck_assert(wavetree2d_sub_parent_index(s, 1) == 0);
  ck_assert(wavetree2d_sub_BL(s, 0) == width);
  ck_assert(wavetree2d_sub_parent_index(s, width) == 0);
  ck_assert(wavetree2d_sub_BR(s, 0) == width + 1);
  ck_assert(wavetree2d_sub_parent_index(s, width + 1) == 0);

  wavetree2d_sub_destroy(s);
}
END_TEST

START_TEST (test_wavetree2d_sub_depth)
{
  wavetree2d_sub_t *s;
  int i;
  int j;
  int width;
  

  /*
   * Create a 128*128 wave grid
   */
  s = wavetree2d_sub_create(7, 7, 0.0);
  ck_assert(s != NULL);

  width = wavetree2d_sub_get_width(s);
  ck_assert(width == 128);
  
  ck_assert(wavetree2d_sub_depthofindex(s, 0) == 0);

  /* Indices 1..4 should be depth 1 */
  for (i = 1, j = 1; i < width; i *= 2, j ++) {
    /*
     * Horizontal
     */
    ck_assert(wavetree2d_sub_depthofindex(s, i) == j);

    /* 
     * Diagonal
     */
    ck_assert(wavetree2d_sub_depthofindex(s, i*width + i) == j);
    
    /*
     * Vertical
     */
    ck_assert(wavetree2d_sub_depthofindex(s, i*width) == j);
  }

  wavetree2d_sub_destroy(s);
}
END_TEST

START_TEST(test_wavetree2d_sub_birth)
{
  wavetree2d_sub_t *s;
  char buffer[256];
  uint64_t binary;

  /*
   * Create a 128*128 wave grid
   */
  s = wavetree2d_sub_create(7, 7, 0.0);
  ck_assert(s != NULL);

  /*
   * Initialization
   */
  ck_assert(wavetree2d_sub_initialize(s, 0.0) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 1);

  /*
   * Add at index 1 and commit
   */
  ck_assert(wavetree2d_sub_propose_birth(s, 1, 1, 1.0) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 2);
  ck_assert(wavetree2d_sub_commit(s) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 2);

  /*
   * Add at index 2 and undo
   */
  ck_assert(wavetree2d_sub_propose_birth(s, 2, 2, 2.0) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 3);
  ck_assert(wavetree2d_sub_undo(s) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 2);

  /*
   * Add at index 3 and commit
   */
  ck_assert(wavetree2d_sub_propose_birth(s, 3, 2, 3.0) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 3);
  ck_assert(wavetree2d_sub_commit(s) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 3);

  /*
   * Check Dyck word is correct (we use a pseudo Dyck with "." == "()" to
   * aid readability.
   *
   * 0 = (...)
   * 0, 1 = ((....)..)
   * 0, 1, 2 = (((....)...)..)
   * 0, 1, 3 = ((.(....)..)..)
   */
  ck_assert(wavetree2d_sub_generate_dyck_word(s, buffer, 256) == 0);
  ck_assert(strcmp(buffer, "((.(....)..)..)") == 0);
    
  /* Also check binary version.
   * 0, 1, 3 in Reverse binary:
   * 0001 0010 1010 1101 0110 1011 = 845B6D = d6b548
   */
  ck_assert(wavetree2d_sub_generate_dyck_binary(s, &binary) == 0);
  ck_assert(binary == 0xd6b548);

  wavetree2d_sub_destroy(s);
}
END_TEST

START_TEST(test_wavetree2d_sub_move)
{
  wavetree2d_sub_t *s;

  int indices[10];
  int nindices;

  /*
   * Create a 128*128 wave grid
   */
  s = wavetree2d_sub_create(7, 7, 0.0);
  ck_assert(s != NULL);

  /*
   * Initialization
   */
  ck_assert(wavetree2d_sub_initialize(s, 0.0) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 1);

  /*
   * Add at index 1 and commit
   */
  ck_assert(wavetree2d_sub_propose_birth(s, 1, 1, 1.0) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 2);
  ck_assert(wavetree2d_sub_commit(s) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 2);

  /*
   * Add at index 128 and commit
   */
  ck_assert(wavetree2d_sub_propose_birth(s, 128, 1, 2.0) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 3);
  ck_assert(wavetree2d_sub_commit(s) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 3);

  /*
   * Check potential move destinations
   */
  ck_assert(wavetree2d_sub_move_available_siblings(s, 1, 1, indices, &nindices) == 0);
  /* for (i = 0; i < nindices; i ++) { */
  /*   printf("%d %d\n", i, indices[i]); */
  /* } */
  ck_assert(nindices == 1);
  ck_assert(indices[0] == 129);

  ck_assert(wavetree2d_sub_move_available_siblings(s, 1, 128, indices, &nindices) == 0);
  ck_assert(nindices == 1);
  ck_assert(indices[0] == 129);

  /*
   * Propose move of 128 to 129
   */
  ck_assert(wavetree2d_sub_propose_move(s, 128, 1, 129, 3.0) == 0);

  nindices = 10;
  ck_assert(wavetree2d_sub_get_indices(s, indices, &nindices) == 0);
  ck_assert(nindices == 3);

  ck_assert(indices[0] == 0);
  ck_assert(indices[1] == 1);
  ck_assert(indices[2] == 129);

  /* for (i = 0; i < nindices; i ++) { */
  /*   printf("%d %d\n", i, indices[i]); */
  /* } */

  wavetree2d_sub_destroy(s);
}
END_TEST

START_TEST(test_wavetree2d_sub_image_mapping)
{
  static const int IM_DEGREE = 4;
  static const int IM_W = 16;
  static const int IM_SIZE = 256;
  wavetree2d_sub_t *s;
  int i;
  int width;

  double *img;

  img = malloc(sizeof(double) * IM_SIZE);
  ck_assert(img);
  
  /*
   * Create a 8x8 wave grid
   */
  s = wavetree2d_sub_create(IM_DEGREE, IM_DEGREE, 0.0);
  ck_assert(s != NULL);

  width = wavetree2d_sub_get_width(s);
  ck_assert(width == 16);

  ck_assert(wavetree2d_sub_get_size(s) == (IM_SIZE));

  /*
   * Initialization
   */
  ck_assert(wavetree2d_sub_initialize(s, 1.0) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 1);

  /*
   * Map root node to get constant image
   */
  for (i = 0; i < (IM_SIZE); i ++) {
    img[i] = 0.0;
  }
				    
  ck_assert(wavetree2d_sub_map_to_array(s, img, IM_SIZE) >= 0);

  ck_assert(img[0] == 1.0);
  for (i = 1; i < (IM_SIZE); i ++) {
    ck_assert(img[i] == 0.0);
  }

  /*
   * Add at index 1 and commit
   */
  ck_assert(wavetree2d_sub_propose_birth(s, 1, 1, 0.5) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 2);
  ck_assert(wavetree2d_sub_commit(s) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 2);
  ck_assert(wavetree2d_sub_valid(s));

  ck_assert(wavetree2d_sub_map_to_array(s, img, IM_SIZE) >= 0);

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
  ck_assert(wavetree2d_sub_propose_birth(s, width, 1, 0.25) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 3);
  ck_assert(wavetree2d_sub_commit(s) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 3);
  ck_assert(wavetree2d_sub_valid(s));

  ck_assert(wavetree2d_sub_map_to_array(s, img, IM_SIZE) >= 0);

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
  
  wavetree2d_sub_destroy(s);

  free(img);
}
END_TEST

START_TEST(test_wavetree2d_sub_image_mapping_rectangular)
{
  static const int IM_DEGREEX = 4;
  static const int IM_DEGREEY = 5;
  static const int IM_W = 16;
  static const int IM_H = 32;
  static const int IM_SIZE = 512;
  
  wavetree2d_sub_t *s;
  int i;
  int width;
  int height;

  double *img;

  img = malloc(sizeof(double) * IM_SIZE);
  ck_assert(img != NULL);
  
  /*
   * Create a 8x8 wave grid
   */
  s = wavetree2d_sub_create(IM_DEGREEX, IM_DEGREEY, 0.0);
  ck_assert(s != NULL);

  width = wavetree2d_sub_get_width(s);
  ck_assert(width == IM_W);

  height = wavetree2d_sub_get_height(s);
  ck_assert(height == IM_H);

  ck_assert(wavetree2d_sub_get_size(s) == (IM_SIZE));

  /*
   * Initialization
   */
  ck_assert(wavetree2d_sub_initialize(s, 1.0) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 1);

  /*
   * Map root node to get constant image
   */
  for (i = 0; i < (IM_SIZE); i ++) {
    img[i] = 0.0;
  }
  
  ck_assert(wavetree2d_sub_map_to_array(s, img, IM_SIZE) >= 0);

  for (i = 0; i < (IM_SIZE); i ++) {
    if (i == 0 || i == IM_W) {
      ck_assert(img[i] == 1.0);
    } else {
      ck_assert(img[i] == 0.0);
    }
  }

  /*
   * Add at index 1 and commit
   */
  ck_assert(wavetree2d_sub_propose_birth(s, 1, 1, 0.5) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 2);
  ck_assert(wavetree2d_sub_commit(s) >= 0);

  /* wavetree2d_sub_dump_sets(s); */

  ck_assert(wavetree2d_sub_coeff_count(s) == 2);
  ck_assert(wavetree2d_sub_valid(s));


  ck_assert(wavetree2d_sub_map_to_array(s, img, IM_SIZE) >= 0);

  /*
   * Check map is correct
   */
  for (i = 0; i < (IM_SIZE); i ++) {
    if (i == 0) {
      ck_assert(img[i] == 1.5);
    } else if (i == IM_W) {
      ck_assert(img[i] == 1.0);
    } else {
      ck_assert(img[i] == 0.0);
    }
  }

  /*
   * Add at index width and commit
   */
  ck_assert(wavetree2d_sub_propose_birth(s, IM_W + 1, 1, 0.25) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 3);
  ck_assert(wavetree2d_sub_commit(s) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 3);
  ck_assert(wavetree2d_sub_valid(s));

  ck_assert(wavetree2d_sub_map_to_array(s, img, IM_SIZE) >= 0);

  /*
   * Check TL quadrant is correct
   */
  for (i = 0; i < (IM_SIZE); i ++) {
    if (i == 0) {
      ck_assert(img[i] == 1.5);
    } else if (i == IM_W) {
      ck_assert(img[i] == 1.25);
    } else {
      ck_assert(img[i] == 0.0);
    }
  }
  
  wavetree2d_sub_destroy(s);

  free(img);
}
END_TEST

START_TEST(test_wavetree2d_sub_saveload)
{
  wavetree2d_sub_t *s;
  double c;

  /*
   * Create a 128*128 wave grid
   */
  s = wavetree2d_sub_create(7, 7, 0.0);
  ck_assert(s != NULL);

  /*
   * Initialization
   */
  ck_assert(wavetree2d_sub_initialize(s, 0.0) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 1);

  /*
   * Add at index 1 and commit
   */
  ck_assert(wavetree2d_sub_propose_birth(s, 1, 1, 1.0) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 2);
  ck_assert(wavetree2d_sub_commit(s) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 2);

  /*
   * Add at index 2 and undo
   */
  ck_assert(wavetree2d_sub_propose_birth(s, 2, 2, 2.0) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 3);
  ck_assert(wavetree2d_sub_commit(s) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 3);

  /*
   * Add at index 3 and commit
   */
  ck_assert(wavetree2d_sub_propose_birth(s, 3, 2, 3.0) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 4);
  ck_assert(wavetree2d_sub_commit(s) >= 0);
  ck_assert(wavetree2d_sub_coeff_count(s) == 4);

  /*
   * Save & Destroy
   */
  ck_assert(wavetree2d_sub_save(s, "wavetree2d_sub_tests_saveload.txt") >= 0);
  wavetree2d_sub_destroy(s);

  /* 
   * Recreate 
   */
  s = wavetree2d_sub_create(7, 7, 0.0);
  ck_assert(s != NULL);

  /*
   * Load
   */
  ck_assert(wavetree2d_sub_load(s, "wavetree2d_sub_tests_saveload.txt") >= 0);

  ck_assert(wavetree2d_sub_coeff_count(s) == 4);

  ck_assert(wavetree2d_sub_get_coeff(s, 1, &c) >= 0);
  ck_assert(c == 1.0);

  ck_assert(wavetree2d_sub_get_coeff(s, 2, &c) >= 0);
  ck_assert(c == 2.0);

  ck_assert(wavetree2d_sub_get_coeff(s, 3, &c) >= 0);
  ck_assert(c == 3.0);

  wavetree2d_sub_destroy(s);

  /*
   * Recreate a higher level 
   */
  s = wavetree2d_sub_create(8, 8, 0.0);
  ck_assert(s != NULL);

  /*
   * Load
   */
  ck_assert(wavetree2d_sub_load_promote(s, "wavetree2d_sub_tests_saveload.txt") >= 0);

  ck_assert(wavetree2d_sub_coeff_count(s) == 4);

  ck_assert(wavetree2d_sub_get_coeff(s, 1, &c) >= 0);
  ck_assert(c == 1.0);

  ck_assert(wavetree2d_sub_get_coeff(s, 2, &c) >= 0);
  ck_assert(c == 2.0);

  ck_assert(wavetree2d_sub_get_coeff(s, 3, &c) >= 0);
  ck_assert(c == 3.0);

  wavetree2d_sub_destroy(s);
}
END_TEST

START_TEST(test_wavetree2d_sub_dyck_duplicate)
{
  wavetree2d_sub_t *s1;
  wavetree2d_sub_t *s2;

  static const int K = 5;
  int indices1[] = {0, 65, 130, 131, 263};
  int indices2[] = {0, 1, 66, 67, 135};
  //0x5aa4d5ff

  uint64_t binary1;
  uint64_t binary2;

  char dyck1[80];
  char dyck2[80];

  int rindices1[10];
  int rk1;
  int rindices2[10];
  int rk2;
  int i;
  int d;

  s1 = wavetree2d_sub_create(6, 6, 0.0);
  ck_assert(s1 != NULL);

  ck_assert(wavetree2d_sub_initialize(s1, 0.0) >= 0);

  for (i = 1; i < K; i ++) {
    ck_assert(wavetree2d_sub_propose_birth(s1, indices1[i], wavetree2d_sub_depthofindex(s1, indices1[i]), 0.0) >= 0);
    ck_assert(wavetree2d_sub_commit(s1) >= 0);
  }
  ck_assert(wavetree2d_sub_coeff_count(s1) == K);

  s2 = wavetree2d_sub_create(6, 6, 0.0);
  ck_assert(s2 != NULL);

  ck_assert(wavetree2d_sub_initialize(s2, 0.0) >= 0);

  for (i = 1; i < K; i ++) {
    ck_assert(wavetree2d_sub_propose_birth(s2, indices2[i], wavetree2d_sub_depthofindex(s2, indices2[i]), 0.0) >= 0);
    ck_assert(wavetree2d_sub_commit(s2) >= 0);
  }
  ck_assert(wavetree2d_sub_coeff_count(s2) == K);

  ck_assert(wavetree2d_sub_get_indices(s1, rindices1, &rk1) >= 0);
  ck_assert(rk1 == K);
  for (i = 0; i < K; i ++) {
    ck_assert(rindices1[i] == indices1[i]);
  }

  ck_assert(wavetree2d_sub_get_indices(s2, rindices2, &rk2) >= 0);
  ck_assert(rk2 == K);
  for (i = 0; i < K; i ++) {
    ck_assert(rindices2[i] == indices2[i]);
  }

  d = 0;
  for (i = 0; i < K; i ++) {
    d += abs(rindices1[i] - rindices2[i]);
  }
  ck_assert(d > 0);

  /* for(i = 0; i < K; i ++) { */
  /*   printf("%d %d\n", rindices1[i], rindices2[i]); */
  /* } */

  ck_assert(wavetree2d_sub_generate_dyck_binary(s1, &binary1) >= 0);
  ck_assert(wavetree2d_sub_generate_dyck_binary(s2, &binary2) >= 0);

  ck_assert(wavetree2d_sub_generate_dyck_word(s1, dyck1, 80) >= 0);
  ck_assert(wavetree2d_sub_generate_dyck_word(s2, dyck2, 80) >= 0);

  /* printf("%lx %lx\n", binary1, binary2); */
  /* printf("%s\n%s\n", dyck1, dyck2); */

  ck_assert(binary1 == 0xeb5aa4d514);
  ck_assert(binary2 == 0xd75aa4d528);

  ck_assert(binary1 != binary2);
}
END_TEST

int recurse_indices(wavetree2d_sub_t *s, int index, int depth, int *indices)
{
  #define MAXCHILDREN 32
  int children[MAXCHILDREN];
  int nchildren;
  int i;
  int j;
  
  if (index < 0) {
    return 0;
  }

  if (indices[index] != 0) {
    printf("index already set %d\n", index);
    return -1;
  }
  
  ck_assert(indices[index] == 0);
  indices[index] = 1;

  ck_assert(wavetree2d_sub_child_indices(s, index, depth, children, &nchildren, MAXCHILDREN) >= 0);

  for (i = 0; i < nchildren; i ++) {

    ck_assert(wavetree2d_sub_parent_index(s, children[i]) == index);
    
    if (recurse_indices(s, children[i], depth + 1, indices) < 0) {
      printf("  error from index %d -> %d\n", index, children[i]);
      return -1;
    }
  }

  return 0;
}

START_TEST(test_wavetree2d_sub_nonsquare_coverage)
{
  wavetree2d_sub_t *s;
  int *indices;
  int size;
  int i;
  int nindices;
  
  s = wavetree2d_sub_create(8, 7, 0.0);
  ck_assert(s != NULL);

  ck_assert(wavetree2d_sub_get_width(s) == 256);
  ck_assert(wavetree2d_sub_get_height(s) == 128);


  ck_assert(wavetree2d_sub_max_child_count(s) == 4);

  size = 256 * 128;
  ck_assert(wavetree2d_sub_get_size(s) == size);
  size ++;

  indices = malloc(sizeof(int) * size);
  ck_assert(indices != NULL);

  ck_assert(wavetree2d_sub_child_indices(s, 0, 0, indices, &nindices, 8) >= 0);
  ck_assert(nindices == 2);
  ck_assert(indices[0] == 1);
  ck_assert(indices[1] == 2);
  
  ck_assert(wavetree2d_sub_child_indices(s, 1, 1, indices, &nindices, 8) >= 0);
  ck_assert(nindices == 3);

  ck_assert(indices[0] == 3);
  ck_assert(indices[1] == 257);
  ck_assert(indices[2] == 259);
  
  ck_assert(wavetree2d_sub_child_indices(s, 2, 1, indices, &nindices, 8) >= 0);
  ck_assert(nindices == 3);

  ck_assert(indices[0] == 4);
  ck_assert(indices[1] == 258);
  ck_assert(indices[2] == 260);
  
  memset(indices, 0, sizeof(int) * size);

  ck_assert(recurse_indices(s, 0, 0, indices) == 0);

  for (i = 0; i < size; i ++) {
    if (indices[i] != 1) {
      printf("%d = %d\n", i, indices[i]);
    }
    ck_assert(indices[i] == 1);
  }

  free(indices);
  wavetree2d_sub_destroy(s);
}
END_TEST

START_TEST(test_wavetree2d_sub_map_from_array)
{
  wavetree2d_sub_t *s;
  double input_image[32*32];
  double output_image[32*32];
  int i;
  int j;
  
  s = wavetree2d_sub_create(5, 5, 0.0);
  ck_assert(s != NULL);

  ck_assert(wavetree2d_sub_get_width(s) == 32);
  ck_assert(wavetree2d_sub_get_height(s) == 32);

  for (j = 0; j < 32; j ++) {
    for (i = 0; i < 32; i ++) {
      input_image[j*32 + i] = 3.14;
    }
  }

  ck_assert(wavetree2d_sub_create_from_array_with_threshold(s, input_image, 32*32, 0.0) == 0);

  ck_assert(wavetree2d_sub_coeff_count(s) == (32*32));

  ck_assert(wavetree2d_sub_map_to_array(s, output_image, 32*32) == 0);

  for (j = 0; j < 32; j ++) {
    for (i = 0; i < 32; i ++) {
      ck_assert(input_image[j*32 + i] == output_image[j*32 + i]);
    }
  }

  /*
   * Now check thresholding
   */
  for (j = 0; j < 4; j ++) {
    for (i = 0; i < 4; i ++) {
      input_image[j*32 + i] = 6.28;
    }
  }
  
  ck_assert(wavetree2d_sub_create_from_array_with_threshold(s, input_image, 32*32, 5.0) == 0);

  ck_assert(wavetree2d_sub_coeff_count(s) == (4*4));

  for (j = 0; j < 32; j ++) {
    for (i = 0; i < 32; i ++) {
      output_image[j*32 + i] = 0.0;
    }
  }
  
  ck_assert(wavetree2d_sub_map_to_array(s, output_image, 32*32) == 0);

  for (j = 0; j < 32; j ++) {
    for (i = 0; i < 32; i ++) {

      if (j < 4 && i < 4) {
	ck_assert(output_image[j*32 + i] == input_image[j*32 + i]);
      } else {
	ck_assert(output_image[j*32 + i] == 0.0);
      }
      
    }
  }
  
}
END_TEST

START_TEST(test_wavetree2d_sub_map_from_array_nonsquare)
{
  wavetree2d_sub_t *s;
  double input_image[32*16];
  double output_image[32*16];
  int i;
  int j;
  
  s = wavetree2d_sub_create(5, 4, 0.0);
  ck_assert(s != NULL);

  ck_assert(wavetree2d_sub_get_width(s) == 32);
  ck_assert(wavetree2d_sub_get_height(s) == 16);

  for (j = 0; j < 16; j ++) {
    for (i = 0; i < 32; i ++) {
      input_image[j*32 + i] = 3.14;
    }
  }

  ck_assert(wavetree2d_sub_create_from_array_with_threshold(s, input_image, 32*16, 0.0) == 0);

  ck_assert(wavetree2d_sub_coeff_count(s) == (32*16 + 1));

  ck_assert(wavetree2d_sub_map_to_array(s, output_image, 32*16) == 0);

  for (j = 0; j < 16; j ++) {
    for (i = 0; i < 32; i ++) {
      ck_assert(input_image[j*32 + i] == output_image[j*32 + i]);
    }
  }

  /*
   * Now check thresholding
   */
  for (j = 0; j < 4; j ++) {
    for (i = 0; i < 8; i ++) {
      input_image[j*32 + i] = 6.28;
    }
  }
  
  ck_assert(wavetree2d_sub_create_from_array_with_threshold(s, input_image, 32*16, 5.0) == 0);

  ck_assert(wavetree2d_sub_coeff_count(s) == (8*4 + 1));

  for (j = 0; j < 16; j ++) {
    for (i = 0; i < 32; i ++) {
      output_image[j*32 + i] = 0.0;
    }
  }
  
  ck_assert(wavetree2d_sub_map_to_array(s, output_image, 32*16) == 0);

  for (j = 0; j < 16; j ++) {
    for (i = 0; i < 32; i ++) {

      if (j < 4 && i < 8) {
	ck_assert(output_image[j*32 + i] == input_image[j*32 + i]);
      } else {
	ck_assert(output_image[j*32 + i] == 0.0);
      }
      
    }
  }
  
}
END_TEST

START_TEST(test_wavetree2d_sub_map_impulse)
{
  
  wavetree2d_sub_t *s;
  double image[16*8];
  int i;
  int j;
  int idx;
  int size;

  size = 16 * 8;
  
  s = wavetree2d_sub_create(4, 3, 0.0);
  ck_assert(s != NULL);

  ck_assert(wavetree2d_sub_get_width(s) == 16);
  ck_assert(wavetree2d_sub_get_height(s) == 8);

  //
  // DC
  //
  memset(image, 0, size * sizeof(double));
  ck_assert(wavetree2d_sub_map_impulse_to_array(s, 0, image, size) == 0);

  for (j = 0; j < 8; j ++) {
    for (i = 0; i < 16; i ++) {

      idx = j * 16 + i;

      if (idx == 0 || idx == 1) {
	ck_assert(image[idx] == 1.0);
      } else {
	ck_assert(image[idx] == 0.0);
      }
    }
  }

  memset(image, 0, size * sizeof(double));
  ck_assert(wavetree2d_sub_map_impulse_to_array(s, 15, image, size) == 0);

  for (j = 0; j < 8; j ++) {
    for (i = 0; i < 16; i ++) {

      idx = j * 16 + i;

      if (idx == 14) {
	ck_assert(image[idx] == 1.0);
      } else {
	ck_assert(image[idx] == 0.0);
      }
    }
  }

  memset(image, 0, size * sizeof(double));
  ck_assert(wavetree2d_sub_map_impulse_to_array(s, 128, image, size) == 0);

  for (j = 0; j < 8; j ++) {
    for (i = 0; i < 16; i ++) {

      idx = j * 16 + i;

      if (idx == 127) {
	ck_assert(image[idx] == 1.0);
      } else {
	ck_assert(image[idx] == 0.0);
      }
    }
  }

  wavetree2d_sub_destroy(s);
}
END_TEST

Suite *
wavetree2d_sub_suite (void)
{
  Suite *s = suite_create ("Wave Tree 2D");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_wavetree2d_sub_ncoefficients);
  tcase_add_test (tc_core, test_wavetree2d_sub_2dindices);
  tcase_add_test (tc_core, test_wavetree2d_sub_2dindices_rectangle);
  tcase_add_test (tc_core, test_wavetree2d_sub_childindices);
  tcase_add_test (tc_core, test_wavetree2d_sub_depth);
  tcase_add_test (tc_core, test_wavetree2d_sub_birth);
  tcase_add_test (tc_core, test_wavetree2d_sub_image_mapping);
  tcase_add_test (tc_core, test_wavetree2d_sub_image_mapping_rectangular);
  tcase_add_test (tc_core, test_wavetree2d_sub_saveload);

  tcase_add_test (tc_core, test_wavetree2d_sub_dyck_duplicate);

  tcase_add_test (tc_core, test_wavetree2d_sub_move);

  tcase_add_test (tc_core, test_wavetree2d_sub_nonsquare_coverage);

  tcase_add_test (tc_core, test_wavetree2d_sub_map_from_array);
  tcase_add_test (tc_core, test_wavetree2d_sub_map_from_array_nonsquare);

  tcase_add_test (tc_core, test_wavetree2d_sub_map_impulse);
  
  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = wavetree2d_sub_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
