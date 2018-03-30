
#include <stdio.h>
#include <stdlib.h>

#include "subdivisiontree2d.h"

int save_image(const char *filename,
	       int width,
	       const double *img);

int main(int argc, char *argv[])
{
  subdivisiontree2d_t *s;

  #define DEGREE 8
  #define IM_W 256
  #define IM_SIZE 65536

  double img[IM_SIZE];
  int width;
  int size = IM_SIZE;

  int i;

  s = subdivisiontree2d_create(DEGREE, 0.0, 0.0, SUBDIVISION_BASIS_LANCZOS);
  if (s == NULL) {
    fprintf(stderr, "error: failed to create tree\n");
    return -1;
  }

  width = subdivisiontree2d_get_width(s);
  if (width != IM_W) {
    fprintf(stderr, "error: width mismatch\n");
    return -1;
  }

  if (subdivisiontree2d_initialize(s, 1.0) < 0) {
    fprintf(stderr, "error: failed to initialize\n");
    return -1;
  }

  if (subdivisiontree2d_map_to_array(s, img, size) < 0) {
    fprintf(stderr, "error: failed to map depth 0 to array\n");
    return -1;
  }

  if (save_image("lanczos_image_0.txt", width, img) < 0) {
    fprintf(stderr, "error: failed to save depth 0 image\n");
    return -1;
  }

  /*
   * Depth 1
   */
  if (subdivisiontree2d_propose_value(s, 0, 0, 0.0) < 0 ||
      subdivisiontree2d_commit(s) < 0) {
    fprintf(stderr, "error: failed to set component\n");
    return -1;
  }

  for (i = subdivisiontree2d_total_coefficients(0); 
       i < subdivisiontree2d_total_coefficients(1); i ++) {
    if (subdivisiontree2d_propose_birth(s, i, 1, 1.0) < 0 ||
	subdivisiontree2d_commit(s) < 0) {
      fprintf(stderr, "error: failed to add component\n");
      return -1;
    }
  }
  
  if (subdivisiontree2d_map_to_array(s, img, size) < 0) {
    fprintf(stderr, "error: failed to map depth 0 to array\n");
    return -1;
  }

  if (save_image("lanczos_image_1.txt", width, img) < 0) {
    fprintf(stderr, "error: failed to save depth 1 image\n");
    return -1;
  }

  /*
   * Depth 2
   */
  for (i = 1; i < 5; i ++) {
    if (subdivisiontree2d_propose_value(s, i, 1, 0.0) < 0 ||
	subdivisiontree2d_commit(s) < 0) {
      fprintf(stderr, "error: failed to set component\n");
      return -1;
    }
  }

  for (i = subdivisiontree2d_total_coefficients(1); i < subdivisiontree2d_total_coefficients(2); i ++) {
    if (subdivisiontree2d_propose_birth(s, i, 2, 1.0) < 0 ||
	subdivisiontree2d_commit(s) < 0) {
      fprintf(stderr, "error: failed to add component\n");
      return -1;
    }
  }
  
  if (subdivisiontree2d_map_to_array(s, img, size) < 0) {
    fprintf(stderr, "error: failed to map depth 0 to array\n");
    return -1;
  }

  if (save_image("lanczos_image_2.txt", width, img) < 0) {
    fprintf(stderr, "error: failed to save depth 1 image\n");
    return -1;
  }

  subdivisiontree2d_destroy(s);
  return 0;
}

int save_image(const char *filename,
	       int width,
	       const double *img)
{
  FILE *fp;
  int i;
  int j;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    fprintf(stderr, "save_image: failed to create file\n");
    return -1;
  }
  
  for (j = 0; j < width; j ++) {
    for (i = 0; i < width; i ++) {
      fprintf(fp, "%f ", img[j*width + i]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
  return 0;
}

