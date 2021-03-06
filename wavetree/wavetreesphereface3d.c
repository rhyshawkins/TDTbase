//
//    Wavetree Library : A library for performed trans-dimensional tree inversion,
//    See
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


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "wavetreesphereface3d.h"

#include "multiset_int.h"
#include "multiset_int_double.h"
#include "oset_int.h"

#include "slog.h"

typedef enum {
  UNDO_NONE = 0,
  UNDO_VALUE,
  UNDO_BIRTH,
  UNDO_DEATH
} wavetreesphereface3d_undo_t;

/*
 * Derived from icosahedron
 */
#define MAX_CHILD_INDICES 12

struct _wavetreesphereface3d {

  manifold_t *manifold; /* Manifold subdivision structure: icosahedron etc */
  int base;             /* Base no. of triangles: 20 for icosahedron */
  
  int degree;           /* Degree of subdivision, 0 = no subdivision */
  int radial_degree;    /* Radial degree, no. slices = 2^radial_degree */
  int rowstride;        /* Total no. of triangles in each slice */
  int depth;            

  int toffset[16];      /* Array offset of triangles for each depth */
  
  int coeff_size;
  int image_size;

  double alpha;

  multiset_int_double_t *S_v;
  multiset_int_t *S_b;
  multiset_int_t *S_d;

  wavetreesphereface3d_undo_t undo;
  int u_i;
  int u_d;
  double u_v;

  int child_indices[MAX_CHILD_INDICES];
};

static int add_node(wavetreesphereface3d_t *t, int i, int d, double coeff);
static int remove_node(wavetreesphereface3d_t *t, int i, int d);

wavetreesphereface3d_t *
wavetreesphereface3d_create(manifold_t *manifold, int max_depth, double alpha)
{
  wavetreesphereface3d_t *r;
  int i;

  /*
   * This tree is a little different. Since the manifold starts with some initial
   * number of vertices > 1, we have encode the root as the mean of these 
   * initial values and the vertices as differences from that mean value to 
   * construct a single rooted tree. This means that our manifold degree must
   * be one less than the degree of the wavetreesphereface3d object. For sanity,
   * we also don't support a 0 degree (which would just be a mean of the entire
   * globe so not much use. See the map_to_array function for how the image
   * is created from the mean + difference values.
   */
  if (max_depth < 1) {
    ERROR("max depth must be 1 or greater");
    return NULL;
  }
  
  if (manifold->degree != (max_depth - 1)) {
    ERROR("inconsistent degree of manifold and wavetreesphereface3d (%d != %d)",
	  manifold->degree, max_depth);
    return NULL;
  }
  
  r = malloc(sizeof(wavetreesphereface3d_t));
  if (r == NULL) {
    return NULL;
  }

  r->manifold = manifold;
  r->base = manifold->ntrianglesatdepth(0);
  r->rowstride = manifold->ntotaltriangles;
  
  r->degree = max_depth;

  /* Radial degree fixed for now */
  r->radial_degree = max_depth - 1;
  r->depth = 1 << r->radial_degree;

  r->toffset[0] = 0;
  for (i = 1; i <= r->degree; i ++) {
    r->toffset[i] = r->toffset[i - 1] + manifold->ntrianglesatdepth(i - 1);
  }
  /* for (i = 0; i < r->degree; i ++) { */
  /*   printf("toffset: %d %d\n", i, r->toffset[i]); */
  /* } */
  

  r->image_size = r->rowstride * r->depth;
  r->coeff_size = r->image_size + 1;

  r->alpha = alpha;
  

  r->S_v = multiset_int_double_create();
  if (r->S_v == NULL) {
    ERROR("failed to allocate Set v");
    return NULL;
  }

  r->S_b = multiset_int_create();
  if (r->S_b == NULL) {
    ERROR("failed to allocate Set b");
    return NULL;
  }

  r->S_d = multiset_int_create();
  if (r->S_d == NULL) {
    ERROR("failed to allocate Set d");
    return NULL;
  }

  /* Initialse undo information */
  r->undo = UNDO_NONE;
  r->u_i = 0;
  r->u_d = 0;
  r->u_v = 0.0;

  return r;
}

void wavetreesphereface3d_destroy(wavetreesphereface3d_t *t)
{
  if (t != NULL) {
    multiset_int_double_destroy(t->S_v);
    multiset_int_destroy(t->S_d);
    multiset_int_destroy(t->S_b);
    free(t);
  }
}

int
wavetreesphereface3d_save(const wavetreesphereface3d_t *t,
			  const char *filename)
{
  FILE *fp;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    ERROR("failed to save wavetree");
    return -1;
  }

  fprintf(fp, "%d\n", t->degree);
  fprintf(fp, "%d %d %d\n", t->rowstride, t->depth, t->image_size);
  fprintf(fp, "%.10g\n", t->alpha);

  if (multiset_int_double_write(t->S_v, fp) < 0) {
    ERROR("failed to save S_v");
    return -1;
  }

  if (multiset_int_write(t->S_b, fp) < 0) {
    ERROR("failed to write S_b");
    return -1;
  }

  if (multiset_int_write(t->S_d, fp) < 0) {
    ERROR("failed to write S_d");
    return -1;
  }

  fclose(fp);

  return 0;
}

int 
wavetreesphereface3d_load(wavetreesphereface3d_t *t,
			  const char *filename)
{
  FILE *fp;

  int degree;
  int rowstride;
  int depth;
  int size;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    ERROR("failed to open file");
    return -1;
  }

  if (fscanf(fp, "%d\n", &degree) != 1) {
    ERROR("failed to read degree");
    return -1;
  }

  if (degree != t->degree) {
    ERROR("degree mismatch");
    return -1;
  }

  if (fscanf(fp, "%d %d %d\n", &rowstride, &depth, &size) != 3) {
    ERROR("failed to read header");
    return -1;
  }

  if (rowstride != t->rowstride ||
      depth != t->depth ||
      size != t->image_size) {
    ERROR("size mismatch");
    return -1;
  }

  if (fscanf(fp, "%lf\n", &(t->alpha)) != 1) {
    ERROR("failed to read alpha");
    return -1;
  }

  /* Clear everything first */

  multiset_int_double_clear(t->S_v);
  multiset_int_clear(t->S_b);
  multiset_int_clear(t->S_d);

  if (multiset_int_double_read(t->S_v, fp) < 0) {
    ERROR("failed to read S_v");
    return -1;
  }

  if (multiset_int_read(t->S_b, fp) < 0) {
    ERROR("failed to read S_b");
    return -1;
  }

  if (multiset_int_read(t->S_d, fp) < 0) {
    ERROR("failed to read S_d");
    return -1;
  }

  return 0;
}

int
wavetreesphereface3d_get_rowstride(wavetreesphereface3d_t *t)
{
  return t->rowstride;
}

int
wavetreesphereface3d_get_width(wavetreesphereface3d_t *t)
{
  ERROR("not implemented");
  return -1;
}

int
wavetreesphereface3d_get_height(wavetreesphereface3d_t *t)
{
  ERROR("not implemented");
  return -1;
}

int
wavetreesphereface3d_get_depth(wavetreesphereface3d_t *t)
{
  if (t != NULL) {
    return t->depth;
  } 
  
  return -1;
}

int
wavetreesphereface3d_get_image_size(wavetreesphereface3d_t *t)
{
  if (t != NULL) {
    return t->image_size;
  }

  return -1;
}

int
wavetreesphereface3d_get_coeff_size(wavetreesphereface3d_t *t)
{
  if (t != NULL) {
    return t->coeff_size;
  }

  return -1;
}

int
wavetreesphereface3d_initialize(wavetreesphereface3d_t *t,
				double dc)
{
  int d;
  int pc;
  int j;

  /*
   * Insert the dc value
   */
  d = wavetreesphereface3d_depthofindex(t, 0);
  multiset_int_double_insert(t->S_v,
			     0,
			     d,
			     dc);
  
  /*
   * Add the potential birthing nodes
   */
  pc = wavetreesphereface3d_get_child_indices(t, 0, t->child_indices, MAX_CHILD_INDICES);
  if (pc < 0) {
    ERROR("failed to get child indices\n");
    return -1;
  }

  for (j = 0; j < pc; j ++) {
    multiset_int_insert(t->S_b,
			t->child_indices[j],
			d + 1);
  }
  
  if (multiset_int_double_total_count(t->S_v) != 1) {
    ERROR("failed to initialize S_v\n");
    return -1;
  }

  if (multiset_int_total_count(t->S_b) != pc) {
    ERROR("failed to initialize S_b: %d\n",
	  multiset_int_total_count(t->S_b));
    return -1;
  }

  return 0;
}

double wavetreesphereface3d_dc(wavetreesphereface3d_t *t)
{
  double dc;

  if (multiset_int_double_get(t->S_v, 0, 0, &dc) < 0) {
    return -1.0;
  }

  return dc;
}

int wavetreesphereface3d_prunable_leaves(const wavetreesphereface3d_t *t)
{
  return multiset_int_total_count(t->S_d);
}


int wavetreesphereface3d_attachable_branches(const wavetreesphereface3d_t *t)
{
  return multiset_int_total_count(t->S_b);
}

int wavetreesphereface3d_coeff_count(const wavetreesphereface3d_t *t)
{
  return multiset_int_double_total_count(t->S_v);
}

int wavetreesphereface3d_map_to_array(wavetreesphereface3d_t *wt, 
				      double *model,
				      int n)
{
  int d;
  int c;

  int i;
  int index;
  double value;
  double root;
  
  /*
   * Get the root element
   */
  if (multiset_int_double_nth_element(wt->S_v, 0, 0, &index, &root) < 0) {
    ERROR("failed to get root element");
    return -1;
  }

  for (i = 0; i < wt->manifold->nverticesatdepth(0); i ++) {
    model[i] = root;
  }

  /* 
   * First level (encoded as difference from root)
   */
  d = 1;
  c = multiset_int_double_depth_count(wt->S_v, d);
  if (c > 0) {
  
    for (i = 0; i < c; i ++) {
      
      if (multiset_int_double_nth_element(wt->S_v, d, i, &index, &value) < 0) {
	ERROR("failed to get nth element");
	return -1;
      }

      model[index - 1] -= value;
    }
  }

  for (d = 2; d < wt->degree; d ++) {
    
    c = multiset_int_double_depth_count(wt->S_v, d);
    if (c > 0) {

      for (i = 0; i < c; i ++) {

	if (multiset_int_double_nth_element(wt->S_v, d, i, &index, &value) < 0) {
	  ERROR("failed to get nth element");
	  return -1;
	}

	model[index - 1] = value;
      }
    }
  }

  return 0;
}

int
wavetreesphereface3d_propose_value(wavetreesphereface3d_t *t,
			       int i,
			       int d,
			       double value)
{
  double old_value;

  if (t == NULL || 
      i < 0 ||
      i >= t->coeff_size) {
    ERROR("invalid parameters");
    return -1;
  }

  if (multiset_int_double_get(t->S_v, i, d, &old_value) < 0) {
    ERROR("failed to get old value (index %d, depth %d)", i, d);
    return -1;
  }

  /*
   * Store the old value in undo information
   */
  t->undo = UNDO_VALUE;
  t->u_i = i;
  t->u_d = d;
  t->u_v = old_value;

  if (multiset_int_double_set(t->S_v, i, d, value) < 0) {
    ERROR("failed to set new value");
    return -1;
  }

  return 0;
}

int
wavetreesphereface3d_propose_birth(wavetreesphereface3d_t *t,
			       int i,
			       int d,
			       double value)
{
  if (t == NULL || 
      i < 0 ||
      i >= t->coeff_size) {
    ERROR("invalid parameters");
    return -1;
  }

  /*
   * Store the undo information
   */
  t->undo = UNDO_BIRTH;
  t->u_i = i;
  t->u_d = d;

  return add_node(t, i, d, value);
}

int 
wavetreesphereface3d_propose_death(wavetreesphereface3d_t *t,
			       int i,
			       int d,
			       double *old_value)
{
  if (t == NULL || 
      i < 0 ||
      i >= t->coeff_size) {
    ERROR("invalid parameters");
    return -1;
  }

  if (multiset_int_double_get(t->S_v, i, d, old_value) < 0) {
    ERROR("failed to get old value (index %d, depth %d)", i, d);
    multiset_int_double_dump(t->S_v);
    return -1;
  }

  t->undo = UNDO_DEATH;
  t->u_i = i;
  t->u_d = d;
  t->u_v = *old_value;

  return remove_node(t, i, d);
}

int 
wavetreesphereface3d_undo(wavetreesphereface3d_t *t)
{
  switch (t->undo) {
  case UNDO_NONE:
    ERROR("nothing to undo");
    return -1;

  case UNDO_VALUE:
    if (multiset_int_double_set(t->S_v, t->u_i, t->u_d, t->u_v) < 0) {
      ERROR("failed to undo value");
      return -1;
    }
    break;

  case UNDO_BIRTH:
    if (remove_node(t, t->u_i, t->u_d) < 0) {
      ERROR("failed to undo birth");
      return -1;
    }
    break;

  case UNDO_DEATH:
    if (add_node(t, t->u_i, t->u_d, t->u_v) < 0) {
      ERROR("failed to undo death");
      return -1;
    }
    break;

  default:
    ERROR("invalid undo type");
    return -1;
  }

  t->undo = UNDO_NONE;
  t->u_i = 0;
  t->u_d = 0;
  t->u_v = 0.0;

  return 0;
}

int
wavetreesphereface3d_commit(wavetreesphereface3d_t *t)
{
  if (t->undo == UNDO_NONE) {
    ERROR("nothing to commit");
    return -1;
  }

  t->undo = UNDO_NONE;
  t->u_i = 0;
  t->u_d = 0;
  t->u_v = 0.0;
  
  return 0;
}

int wavetreesphereface3d_valid(wavetreesphereface3d_t *t)
{
  int error_count;
  int d;
  int n;
  int i;

  int index;
  double value;

  int pi;

  error_count = 0;

  /*
   * Check S_v, that every element has its parent in the list, index 0 exists
   */
  if (!multiset_int_double_is_element(t->S_v, 0, 0)) {
    error_count ++;
    ERROR("missing 0,0 element");
  }

  for (d = 1; d <= t->degree; d ++) {
    n = multiset_int_double_depth_count(t->S_v, d);
    for (i = 0; i < n; i ++) {

      if (multiset_int_double_nth_element(t->S_v, d, i, &index, &value) < 0) {
	ERROR("failed to get nth element of S_v");
	return 0;
      }

      if (wavetreesphereface3d_depthofindex(t, index) != d) {
	error_count ++;
	ERROR("index %d in S_v at incorrect depth %d", index, d);
      }

      pi = wavetreesphereface3d_parent_index(t, index);
      if (!multiset_int_double_is_element(t->S_v, pi, d - 1)) {
	error_count ++;
	ERROR("parent %d not set for index %d", pi, index);
      }
    }
  }

  /*
   * Check S_b, elements in S_b must not be in S_v or S_d
   */
  for (d = 1; d <= t->degree; d ++) {
    n = multiset_int_depth_count(t->S_b, d);
    for (i = 0; i < n; i ++) {

      if (multiset_int_nth_element(t->S_b, d, i, &index) < 0) {
	ERROR("failed to get nth element of S_b");
	return 0;
      }

      if (wavetreesphereface3d_depthofindex(t, index) != d) {
	error_count ++;
	ERROR("index %d in S_b at incorrect depth %d (%d,%d)", index, d, i, n);
	return 0;
      }

      if (multiset_int_is_element(t->S_d, index, d)) {
	error_count ++;
	ERROR("index %d in S_b and S_d", index);
      }

      if (multiset_int_double_is_element(t->S_v, index, d)) {
	error_count ++;
	ERROR("index %d in S_b and S_d", index);
      }
    }
  }

  /*
   * Check S_d, elements in S_d must be in S_v and not in S_b
   */
  for (d = 1; d <= t->degree; d ++) {
    n = multiset_int_depth_count(t->S_d, d);
    for (i = 0; i < n; i ++) {

      if (multiset_int_nth_element(t->S_d, d, i, &index) < 0) {
	ERROR("failed to get nth element of S_b");
	return 0;
      }

      if (wavetreesphereface3d_depthofindex(t, index) != d) {
	error_count ++;
	ERROR("index %d in S_d at incorrect depth %d", index, d);
      }
      
      if (multiset_int_is_element(t->S_b, index, d)) {
	error_count ++;
	ERROR("index %d in S_d and S_b", index);
      }

      if (!multiset_int_double_is_element(t->S_v, index, d)) {
	error_count ++;
	ERROR("index %d in S_d but not in S_v", index);
      }
    }
  }

  return (error_count == 0);
}

void wavetreesphereface3d_print_setinfo(wavetreesphereface3d_t *t)
{
  printf("Nb %d Nd %d Nv %d\n", 
         multiset_int_total_count(t->S_b), 
         multiset_int_total_count(t->S_d), 
         multiset_int_double_total_count(t->S_v));
}

void wavetreesphereface3d_dump_sets(wavetreesphereface3d_t *t)
{
  printf("Sb:\n  ");
  multiset_int_dump(t->S_b);
  printf("Sd:\n  ");
  multiset_int_dump(t->S_d);
  printf("Sv:\n  ");
  multiset_int_double_dump(t->S_v);
}

void wavetreesphereface3d_dump_coeffs(wavetreesphereface3d_t *t)
{
  multiset_int_double_dump(t->S_v);
}

int wavetreesphereface3d_depth(wavetreesphereface3d_t *t) 
{
  int d;
  
  for (d = 1; d <= t->degree; d ++) {
    if (multiset_int_double_depth_count(t->S_v, d) == 0) {
      return d - 1;
    }
  }

  return -1;
}

int wavetreesphereface3d_maxdepth(wavetreesphereface3d_t *t)
{
  return t->degree;
}

int wavetreesphereface3d_update_histogram(const wavetreesphereface3d_t *t, coefficient_histogram_t *hist)
{
  int d;
  int n;
  int i;
  int index;
  double value;

  for (d = 0; d < t->degree; d ++) {
    n = multiset_int_double_depth_count(t->S_v, d);
    for (i = 0; i < n; i ++) {

      if (multiset_int_double_nth_element(t->S_v, d, i, &index, &value) < 0) {
	ERROR("failed to get nth element");
	return -1;
      }

      coefficient_histogram_sample(hist, index, value);
    }
  }

  return 0;
}

int wavetreesphereface3d_depth_filter(void *wt, int subset, int index)
{
  wavetreesphereface3d_t *t = (wavetreesphereface3d_t *)wt;

  return wavetreesphereface3d_depthofindex(t, index) == subset;
}

/*
 * Functions for setting up a birth proposal
 */
int 
wavetreesphereface3d_choose_birth_depth(const wavetreesphereface3d_t *t, double u, int maxdepth, int *depth, double *prob)
{
  int ndepths;

  if (t == NULL) {
    return -1;
  }

  if (multiset_int_choose_depth(t->S_b, u, maxdepth, depth, &ndepths) < 0) {
    return -1;
  }

  *prob = 1.0/((double)ndepths);
  return 0;
}

int wavetreesphereface3d_reverse_birth_depth(const wavetreesphereface3d_t *t, int depth, int maxdepth, double *prob)
{
  if (t == NULL) {
    return -1;
  }

  *prob = 1.0/((double)multiset_int_nonempty_count(t->S_d, maxdepth));
  return 0;
}

int wavetreesphereface3d_choose_birth(const wavetreesphereface3d_t *t, int depth, double u, int *coeff, double *prob)
{
  int nindices;

  if (t == NULL) {
    return -1;
  }

  if (multiset_int_choose_index(t->S_b, depth, u, coeff, &nindices) < 0) {
    return -1;
  }

  if (*coeff <= 0) {
    ERROR("-ve coeff in S_b");
    multiset_int_dump(t->S_b);
    return -1;
  }
  
  *prob = 1.0/(double)(nindices);
  return 0;
}

int wavetreesphereface3d_choose_birth_global(const wavetreesphereface3d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob)
{
  int nindices;

  if (t == NULL) {
    return -1;
  }

  if (t->alpha == 0.0) {
    if (multiset_int_choose_index_globally(t->S_b, u, maxdepth, coeff, depth, &nindices) < 0) {
      return -1;
    }
    *prob = 1.0/(double)(nindices);
  } else {
    if (multiset_int_choose_index_weighted(t->S_b, u, maxdepth, t->alpha, depth, coeff, prob) < 0) {
      return -1;
    }
  }

  if (*coeff <= 0) {
    ERROR("-ve coeff in S_b");
    return -1;
  }

  return 0;
}

int wavetreesphereface3d_reverse_birth(const wavetreesphereface3d_t *t, int depth, int coeff, double *prob)
{
  if (t == NULL) {
    return -1;
  }
  
  *prob = 1.0/(double)(multiset_int_depth_count(t->S_d, depth));
  return 0;
}

int wavetreesphereface3d_reverse_birth_global(const wavetreesphereface3d_t *t, int maxdepth, int depth, int coeff, double *prob)
{
  if (t == NULL) {
    return -1;
  }

  if (t->alpha == 0.0) {
    *prob = 1.0/(double)(multiset_int_restricted_total_count(t->S_d, maxdepth));
  } else {
    if (multiset_int_reverse_choose_index_weighted(t->S_d, maxdepth, t->alpha, coeff, depth, prob) < 0) {
      return -1;
    }
  }

  return 0;
}

/*
 * Functions for setting up a death proposal
 */
int wavetreesphereface3d_choose_death_depth(const wavetreesphereface3d_t *t, double u, int maxdepth, int *depth, double *prob)
{
  int ndepths;

  if (t == NULL) {
    return -1;
  }

  if (multiset_int_choose_depth(t->S_d, u, maxdepth, depth, &ndepths) < 0) {
    return -1;
  }

  *prob = 1.0/((double)ndepths);
  return 0;
}

int wavetreesphereface3d_reverse_death_depth(const wavetreesphereface3d_t *t, int depth, int maxdepth, double *prob)
{
  if (t == NULL) {
    return -1;
  }

  *prob = 1.0/((double)multiset_int_nonempty_count(t->S_b, maxdepth));
  return 0;
}

int wavetreesphereface3d_choose_death(const wavetreesphereface3d_t *t, int depth, double u, int *coeff, double *prob)
{
  int nindices;

  if (t == NULL) {
    return -1;
  }

  if (multiset_int_choose_index(t->S_d, depth, u, coeff, &nindices) < 0) {
    return -1;
  }

  *prob = 1.0/(double)(nindices);
  return 0;
}


int wavetreesphereface3d_choose_death_global(const wavetreesphereface3d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob)
{
  int nindices;

  if (t == NULL) {
    return -1;
  }

  if (t->alpha == 0.0) {
    if (multiset_int_choose_index_globally(t->S_d, u, maxdepth, coeff, depth, &nindices) < 0) {
      return -1;
    }
    *prob = 1.0/(double)(nindices);
  } else {
    if (multiset_int_choose_index_weighted(t->S_d, u, maxdepth, t->alpha, depth, coeff, prob) < 0) {
      return -1;
    }
  }

  return 0;
}

int wavetreesphereface3d_reverse_death(const wavetreesphereface3d_t *t, int depth, int coeff, double *prob)
{
  if (t == NULL) {
    return -1;
  }
  
  *prob = 1.0/(double)(multiset_int_depth_count(t->S_b, depth));
  return 0;
}

int wavetreesphereface3d_reverse_death_global(const wavetreesphereface3d_t *t, int maxdepth, int depth, int coeff, double *prob)
{
  if (t == NULL) {
    return -1;
  }

  if (t->alpha == 0.0) {
    *prob = 1.0/(double)(multiset_int_restricted_total_count(t->S_b, maxdepth));
  } else {
    if (multiset_int_reverse_choose_index_weighted(t->S_b, maxdepth, t->alpha, coeff, depth, prob) < 0) {
      return -1;
    }
  }

  return 0;
}

/*
 * Functions for setting up a value proposal
 */
int wavetreesphereface3d_choose_value_depth(const wavetreesphereface3d_t *t, double u, int maxdepth, int *depth, double *prob)
{
  int ndepths;

  if (t == NULL) {
    return -1;
  }

  if (multiset_int_double_choose_depth(t->S_v, u, maxdepth, depth, &ndepths) < 0) {
    return -1;
  }

  *prob = 1.0/((double)ndepths);
  return 0;
}

int wavetreesphereface3d_choose_value(const wavetreesphereface3d_t *t, int depth, double u, int *coeff, double *prob)
{
  int nindices;

  if (t == NULL) {
    return -1;
  }

  if (multiset_int_double_choose_index(t->S_v, depth, u, coeff, &nindices) < 0) {
    return -1;
  }

  *prob = 1.0/(double)(nindices);
  return 0;
}


int wavetreesphereface3d_choose_value_global(const wavetreesphereface3d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob)
{
  int nindices;

  if (t == NULL) {
    return -1;
  }

  if (multiset_int_double_choose_index_globally(t->S_v, u, maxdepth, coeff, depth, &nindices) < 0) {
    return -1;
  }

  *prob = 1.0/(double)(nindices);
  return 0;
}

int wavetreesphereface3d_perturb(wavetreesphereface3d_t *t,
				 wavetreesphereface3d_perturb_func_t f,
				 void *user,
				 double *prior_ratio)
{
  int d;
  int n;
  int i;
  int index;
  double value;
  double old_value;

  int it;
  int ir;
  int id;

  int pi;
  double pvalue;

  for (d = 0; d < t->degree; d ++) {
    n = multiset_int_double_depth_count(t->S_v, d);
    for (i = 0; i < n; i ++) {

      if (multiset_int_double_nth_element(t->S_v, d, i, &index, &value) < 0) {
	ERROR("failed to get nth element");
	return -1;
      }

      if (wavetreesphereface3d_index_to_sphereindices(t, i, &it, &ir, &id) < 0) {
	ERROR("failed to get sphere3d indices");
	return -1;
      }

      /* Get parent coeff */
      if (index == 0) {
	pvalue = 0.0;
      } else {
	pi = wavetreesphereface3d_parent_index(t, index);
	if (wavetreesphereface3d_get_coeff(t, pi, d - 1, &pvalue) < 0) {
	  ERROR("failed to get parent coeff");
	  return -1;
	}
      }
       
      old_value = value;
      if (f(user, it, 0, ir, d, t->degree, pvalue, &value, prior_ratio) < 0) {
	return -1;
      }

      /* Update the value if required */
      if (old_value != value) {
	if (multiset_int_double_set(t->S_v, index, d, value) < 0) {
	  ERROR("failed to set value");
	  return -1;
	}
      }

    }
  }

  return 0;
}

int wavetreesphereface3d_parent_index(wavetreesphereface3d_t *t, int c)
{
  int it;
  int ir;
  int itd;
  int ird;
  int id;
  int pt;
  int pr;
   
  if (c <= 0) {
    /* 0th has no parent */
    return -1;
  }

  if (wavetreesphereface3d_index_to_separatesphereindices(t, c, &it, &ir, &itd, &ird, &id) < 0) {
    ERROR("failed to get 2d indices");
    return -1;
  }

  if (it < t->base) {
    if (id == 1) {
      /* Special early termination for first level nodes */
      return 0;
    } else {
      pt = it;
    }
  } else {
    /* printf("parent %d it %d id %d offset %d offset -1 %d\n", c, it, id, t->toffset[id - 1], t->toffset[id - 2]); */
    pt = (it - t->toffset[itd])/4 + t->toffset[itd - 1];
  }

  pr = ir/2;
    
  return wavetreesphereface3d_sphereindices_to_index(t, pt, pr, id - 1);  
}

int
wavetreesphereface3d_index_to_separatesphereindices(wavetreesphereface3d_t *wt,
						    int i,
						    int *_t,
						    int *_r,
						    int *_td,
						    int *_rd,
						    int *_d)
{
  int t, r;
  int td;
  int rd;
  int d;
  
  if (i < 0) {
    return -1;
  }
  
  if (i == 0) {
    *_t = 0;
    *_r = 0;
    *_d = 0;

    return 0;
  }

  /* Remove the root of the tree */
  i = i - 1;
  
  t = i % wt->rowstride;
  r = (i - t)/wt->rowstride;

  *_t = t;
  *_r = r;
  
  rd = 0;
  while (r > 0) {
    rd ++;
    r /= 2;
  }

  td = 0;
  while (wt->toffset[td + 1] <= t) {
    td ++;
  }

  *_rd = rd;
  *_td = td;

  d = rd;
  if (td > d) {
    d = td;
  }

  /* We need to add one for the root of the tree (previously we used d for the
   * depth in the manifold subdivision */
  *_d = d + 1;
  
  return 0;
}
						    
int 
wavetreesphereface3d_index_to_sphereindices(wavetreesphereface3d_t *wt,
					    int i,
					    int *_t,
					    int *_r,
					    int *_d)
{
  int td, rd;
  
  if (wavetreesphereface3d_index_to_separatesphereindices(wt, i, _t, _r, &td, &rd, _d) < 0) {
    return -1;
  }
  
  return 0;
}

int 
wavetreesphereface3d_sphereindices_to_index(wavetreesphereface3d_t *wt,
					    int t,
					    int r,
					    int d)
{
  if (d == 0) {
    return 0;
  } else {
    return t + r * wt->rowstride + 1;
  }
}

int wavetreesphereface3d_depthofindex(wavetreesphereface3d_t *wt,
				      int i)
{
  int t, r, d;

  if (wavetreesphereface3d_index_to_sphereindices(wt, i, &t, &r, &d) < 0) {
    return -1;
  }

  return d;
}

int wavetreesphereface3d_get_coeff(wavetreesphereface3d_t *t,
				   int i,
				   int d,
				   double *coeff)
{
  if (multiset_int_double_get(t->S_v, i, d, coeff) < 0) {
    return -1;
  }

  return 0;
}

int
wavetreesphereface3d_get_child_indices(wavetreesphereface3d_t *t,
				       int index,
				       int *child_indices,
				       int maxindices)
{
  int i;
  int it, ir, itd, ird, id;
  int n;
  
  if (index == 0) {
    if (t->base > maxindices) {
      ERROR("indices overflow %d > %d",
	    t->base, maxindices);
      return -1;
    }
    
    for (i = 0; i < t->base; i ++) {
      child_indices[i] = i + 1;
    }
    return t->base;
  } else {

    if (wavetreesphereface3d_index_to_separatesphereindices(t, index, &it, &ir, &itd, &ird, &id) < 0) {
      ERROR("failed to get sphere indices");
      return -1;
    }

    n = 0;
    if (id < t->degree) {
      
      /* 
       * First do the up to 2 direct radial children.
       */
      if (it < t->base) {
	if (ir > 0) {
	  child_indices[n] = wavetreesphereface3d_sphereindices_to_index(t, it, 2*ir, id + 1);
	  n ++;
	}
	
	child_indices[n] = wavetreesphereface3d_sphereindices_to_index(t, it, 2*ir + 1, id + 1);
	n ++;
      }
    
      /*
       * Lastly do the 6 lateral children
       */
      for (i = 0; i < 3; i ++) {
	child_indices[n] = wavetreesphereface3d_sphereindices_to_index(t,
								       t->toffset[itd + 1] +
								       4*(it - t->toffset[itd]) + i,
								       2*ir, id + 1);
	n ++;
      }
      
      for (i = 0; i < 3; i ++) {
	child_indices[n] = wavetreesphereface3d_sphereindices_to_index(t,
								       t->toffset[itd + 1] +
								       4*(it - t->toffset[itd]) + i,
								       2*ir + 1, id + 1);
	n ++;
      }
    }
    
    /* printf("  Child Indices for %d (%d %d %d (%d)):\n", index, it, ir, id, t->rowstride); */
    /* for (i = 0; i < n; i ++) { */
    /*   printf("    %2d: %d\n", i, child_indices[i]); */
    /* } */

    return n;
  }

  ERROR("unreachable");
  return -1;
}


static int add_node(wavetreesphereface3d_t *t, int i, int d, double coeff)
{
  int j;
  int pc;
  int pcc;

  if (multiset_int_double_insert(t->S_v, i, d, coeff) < 0) {
    return -1;
  }
  
  j = wavetreesphereface3d_parent_index(t, i);

  pcc = wavetreesphereface3d_child_count(t, j, d - 1);
  if (pcc < 0) {
    return -1;
  }

  if (pcc == 1) {
    /* This was the parents first child so remove parent from S_d */
    multiset_int_remove(t->S_d, j, d - 1);
  }

  /* Add new node to S_d and remove from S_b */
  multiset_int_insert(t->S_d, i, d);
  multiset_int_remove(t->S_b, i, d);

  /* Add children of new node to S_b */
  if (d < t->degree) {
    pc = wavetreesphereface3d_get_child_indices(t, i, t->child_indices, MAX_CHILD_INDICES);
    if (pc < 0) {
      ERROR("failed to get child indices");
      return -1;
    }
    
    for (j = 0; j < pc; j ++) {
      if (t->child_indices[j] <= 0) {
	ERROR("got invalid child indices for node %d", i);
	for (j = 0; j < pc; j ++) {
	  ERROR("        : %d %d\n", j, t->child_indices[j]);
	}
	return -1;
      }
      multiset_int_insert(t->S_b, t->child_indices[j], d + 1);
    }
  }
  
  return 0;
}

static int remove_node(wavetreesphereface3d_t *t, int i, int d)
{
  int j;
  int pcc;
  int pc;
  
  if (multiset_int_double_remove(t->S_v, i, d) < 0) {
    ERROR("failed to remove index from S_v (index %d, depth %d)", i, d);
    return -1;
  }

  j = wavetreesphereface3d_parent_index(t, i);
  if (j < 0) {
    ERROR("failed to get parent index (index %d, depth %d)", i, d);
    return -1;
  }

  pcc = wavetreesphereface3d_child_count(t, j, d - 1);

  if (j > 0 && pcc == 0) {
    /* Add parent to S_d as it now has no children */
    multiset_int_insert(t->S_d, j, d - 1);
  }

  /* Remove node from S_d and add to S_b */
  multiset_int_remove(t->S_d, i, d);
  multiset_int_insert(t->S_b, i, d);

  /* Remove children from S_b */
  pc = wavetreesphereface3d_get_child_indices(t, i, t->child_indices, MAX_CHILD_INDICES);
  if (pc < 0) {
    ERROR("failed to get child indices");
    return -1;
  }
  
  for (j = 0; j < pc; j ++) {
    multiset_int_remove(t->S_b, t->child_indices[j], d + 1);
  }
  
  return 0;
}

int wavetreesphereface3d_child_count(wavetreesphereface3d_t *t, int index, int depth)
{
  int pc;
  int j;
  int cc;

  cc = 0;

  pc = wavetreesphereface3d_get_child_indices(t, index, t->child_indices, MAX_CHILD_INDICES);
  if (pc < 0) {
    return -1;
  }

  for (j = 0; j < pc; j ++) {
    if (multiset_int_double_is_element(t->S_v, t->child_indices[j], depth + 1)) {
      cc ++;
    }
  }
  
  return cc;
}


