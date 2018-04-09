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

#include "wavetreesphereface2d.h"

#include "multiset_int.h"
#include "multiset_int_double.h"

#include "slog.h"

typedef enum {
  UNDO_NONE = WT_PERTURB_NONE,
  UNDO_VALUE = WT_PERTURB_VALUE,
  UNDO_BIRTH = WT_PERTURB_BIRTH,
  UNDO_DEATH = WT_PERTURB_DEATH,
  UNDO_MOVE = WT_PERTURB_MOVE
} wavetreesphereface2d_undo_t;

struct _wavetreesphereface2d {

  manifold_t *manifold;
  double alpha;

  int degree_max;
  int base_triangles;
  
  int max_children;
  int *child_indices;

  multiset_int_double_t *S_v;
  multiset_int_t *S_b;
  multiset_int_t *S_d;

  wavetreesphereface2d_undo_t undo;
  int u_i;
  int u_j;
  int u_d;
  double u_v;

  chain_history_change_t last_step;
};

static int add_node(wavetreesphereface2d_t *t, int i, int d, double coeff);
static int remove_node(wavetreesphereface2d_t *t, int i, int d);

wavetreesphereface2d_t *
wavetreesphereface2d_create(manifold_t *manifold, double alpha)
{
  wavetreesphereface2d_t *r;

  r = malloc(sizeof(wavetreesphereface2d_t));
  if (r == NULL) {
    return NULL;
  }

  r->manifold = manifold;
  r->alpha = alpha;

  r->degree_max = manifold->degree + 1;
  r->base_triangles = manifold->ntrianglesatdepth(0);

  r->max_children = r->base_triangles;
  r->child_indices = malloc(sizeof(int) * r->max_children);
  if (r->child_indices == NULL) {
    ERROR("failed to allocate child indices");
    return NULL;
  }

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
  r->u_j = 0;
  r->u_d = 0;
  r->u_v = 0.0;

  memset(&(r->last_step), 0, sizeof(chain_history_change_t));
  r->last_step.header.type = CH_INITIALISE;
  
  return r;
}

void wavetreesphereface2d_destroy(wavetreesphereface2d_t *t)
{
  if (t != NULL) {
    multiset_int_double_destroy(t->S_v);
    multiset_int_destroy(t->S_d);
    multiset_int_destroy(t->S_b);

    free(t->child_indices);
      
    free(t);
  }
}

int
wavetreesphereface2d_save(const wavetreesphereface2d_t *t,
			  const char *filename)
{
  FILE *fp;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    ERROR("failed to save wavetree");
    return -1;
  }

  fprintf(fp, "%d %d\n", t->manifold->degree, t->base_triangles);
  fprintf(fp, "%.10g\n", t->alpha);

  if (multiset_int_double_write(t->S_v, fp) < 0) {
    ERROR("failed to write S_v");
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
wavetreesphereface2d_load(wavetreesphereface2d_t *t,
		const char *filename)
{
  FILE *fp;

  int degree;
  int base_triangles;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    ERROR("failed to open file");
    return -1;
  }

  if (fscanf(fp, "%d %d\n", &degree, &base_triangles) != 2) {
    ERROR("failed to read degree");
    return -1;
  }

  if (degree != t->manifold->degree) {
    ERROR("x-degree mismatch");
    return -1;
  }

  if (base_triangles != t->base_triangles) {
    ERROR("y-degree mismatch");
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

  fclose(fp);
    
  return 0;
}

static int encode_int(int v,
		      char *buffer,
		      int *offset,
		      int maxsize)
{
  int *p;
  
  if ((*offset) + sizeof(int) > maxsize) {
    return -1;
  }

  p = (int*)(buffer + (*offset));
  *p = v;
  (*offset) += sizeof(int);
  return 0;
}

static int encode_double(double v,
			 char *buffer,
			 int *offset,
			 int maxsize)
{
  double *p;
  
  if ((*offset) + sizeof(double) > maxsize) {
    return -1;
  }

  p = (double*)(buffer + (*offset));
  *p = v;
  (*offset) += sizeof(double);
  return 0;
}

#define ENCODEINT(i)							\
  if (encode_int(i, buffer, &offset, maxsize) < 0) {			\
    ERROR("buffer overflow %d %d (int)", offset, maxsize); \
    return -1;								\
  }

#define ENCODEDOUBLE(i)							\
  if (encode_double(i, buffer, &offset, maxsize) < 0) {			\
    ERROR("buffer overflow %d %d (double)", offset, maxsize); \
    return -1;								\
  }

int
wavetreesphereface2d_encode(wavetreesphereface2d_t *t,
			    char *buffer,
			    int maxsize)
{
  int offset;
  int d;
  int n;
  int i;

  int index;
  double value;
  
  offset = 0;

  /*
   * Header
   */
  
  ENCODEINT(t->manifold->degree);
  ENCODEINT(t->base_triangles);
  ENCODEDOUBLE(t->alpha);
  
  /*
   * No. elements
   */
  ENCODEINT(multiset_int_double_total_count(t->S_v));
  
  /*
   * Values
   */
  for (d = 0; d <= t->degree_max; d ++) {

    n = multiset_int_double_depth_count(t->S_v, d);
    for (i = 0; i < n; i ++) {

      if (multiset_int_double_nth_element(t->S_v, d, i, &index, &value) < 0) {
	ERROR("failed to get nth element");
	return -1;
      }

      /* printf("    encoding %d %d %f\n", d, index, value); */
      ENCODEINT(d);
      ENCODEINT(index);
      ENCODEDOUBLE(value);

    }
  }

  return offset;
}

static int decode_int(char *buffer, int *offset, int len, int *v)
{
  int *p;
  
  if ((*offset) + sizeof(int) > len) {
    return -1;
  }

  p = (int*)(buffer + (*offset));
  *v = *p;
  (*offset) += sizeof(int);
  return 0;
}

static int decode_double(char *buffer, int *offset, int len, double *v)
{
  double *p;
  
  if ((*offset) + sizeof(double) > len) {
    return -1;
  }

  p = (double*)(buffer + (*offset));
  *v = *p;
  (*offset) += sizeof(double);
  return 0;
}

#define DECODEINT(i)							\
  if (decode_int(buffer, &offset, len, i) < 0) {			\
    ERROR("buffer overflow %d %d (int)", offset, len); \
    return -1;								\
  }

#define DECODEDOUBLE(i)							\
  if (decode_double(buffer, &offset, len, i) < 0) {			\
    ERROR("buffer overflow %d %d (double)", offset, len); \
    return -1;								\
  }

int
wavetreesphereface2d_decode(wavetreesphereface2d_t *t,
			    char *buffer,
			    int len)
{
  int offset;
  int i;
  int n;

  int d;
  int index;
  double v;
  int iv;
  
  offset = 0;

  /*
   * Header
   */
  DECODEINT(&iv);
  if (iv != t->manifold->degree) {
    ERROR("degree  mismatch %d != %d", t->manifold->degree, iv);
    return -1;
  }
  DECODEINT(&iv);
  if (iv != t->base_triangles) {
    ERROR("base triangles mismatch %d != %d", t->base_triangles, iv);
    return -1;
  }
  DECODEDOUBLE(&v);
  if (v != t->alpha) {
    ERROR("alpha mismatch %f != %f", t->alpha, v);
    return -1;
  }

  /*
   * No. elements
   */
  DECODEINT(&n);
  
  /* Clear everything first before adding nodes from binary string */
  multiset_int_double_clear(t->S_v);
  multiset_int_clear(t->S_b);
  multiset_int_clear(t->S_d);

  /*
   * First Value we use for initialization
   */
  DECODEINT(&d);
  DECODEINT(&index);
  DECODEDOUBLE(&v);

  if (d != 0 && index != 0) {
    ERROR("unexpected initial node %d %d", d, index);
    return -1;
  }
  
  if (wavetreesphereface2d_initialize(t, v) < 0) {
    ERROR("failed to initialse root node");
    return -1;
  }
  
  /*
   * Remaining Values
   */
  for (i = 1; i < n; i ++) {

    DECODEINT(&d);
    DECODEINT(&index);
    DECODEDOUBLE(&v);

    if (add_node(t, index, d, v) < 0) {
      ERROR("failed to add node (%d %d %f)", index, d, v);
      return -1;
    }
  }

  return 0;
}

int 
wavetreesphereface2d_load_promote(wavetreesphereface2d_t *t,
				  const char *filename)
{
  FILE *fp;
  int i;
  int j;
  int coeffs;
  int idx;

  double coeff;

  int degree;
  int base_triangles;

  int d;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    ERROR("failed to open file");
    return -1;
  }

  if (fscanf(fp, "%d %d\n", &degree, &base_triangles) != 2) {
    ERROR("failed to read degree");
    return -1;
  }

  /* Promotion works for smaller or equal degree */
  if (degree > t->manifold->degree) {
    ERROR("degree mismatch %d > %d", 
	  degree, t->manifold->degree);
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

  if (fscanf(fp, "%d\n", &d) != 1) {
    ERROR("failed to read depth");
    return -1;
  }
  
  for (i = 0; i < d; i ++) {
    if (fscanf(fp, "%d %d\n", &j, &coeffs) != 2) {
      ERROR("failed to read depth no. coefficients");
      return -1;
    }

    if (j != i) {
      ERROR("depth index mismatch (%d %d)",
	    j, i);
      return -1;
    }

    for (j = 0; j < coeffs; j ++) {
      if (fscanf(fp, "%d %lf\n", &idx, &coeff) != 2) {
	ERROR("failed to read coefficient");
	return -1;
      }

      if (idx == 0) {
	if (i != 0) {
	  ERROR("index 0 not at depth 0");
	  return -1;
	}

	if (wavetreesphereface2d_initialize(t, coeff) < 0) {
	  ERROR("failed to initialize from 0 node");
	  return -1;
	}
      } else {
	if (add_node(t, idx, i, coeff) < 0) {
	  ERROR("failed to add node");
	  return -1;
	}
      }
    }
  }

  /*
   * Note that we don't read in S_b, S_d as we create these manually ourselves 
   * using the add_node function.
   */
  fclose(fp);
  return 0;
}

int wavetreesphereface2d_initialize(wavetreesphereface2d_t *t,
				    double dc)
{
  int d;
  int i;

  d = wavetreesphereface2d_depthofindex(t, 0);
  multiset_int_double_insert(t->S_v,
			     0,
			     d,
			     dc);

  for (i = 0; i < t->base_triangles; i ++) {
    multiset_int_insert(t->S_b,
			i + 1,
			d + 1);
  }
  
  if (multiset_int_total_count(t->S_b) != t->base_triangles) {
    ERROR("failed to initialize S_b");
    return -1;
  }

  if (multiset_int_double_total_count(t->S_v) != 1) {
    ERROR("failed to initialize S_v");
    return -1;
  }


  return 0;
}

double wavetreesphereface2d_dc(wavetreesphereface2d_t *t)
{
  double dc;

  if (multiset_int_double_get(t->S_v, 0, 0, &dc) < 0) {
    return -1.0;
  }

  return dc;
}

int
wavetreesphereface2d_get_size(wavetreesphereface2d_t *t)
{
  return t->manifold->ntotaltriangles + 1;
}

int
wavetreesphereface2d_get_ncoeff(wavetreesphereface2d_t *t)
{
  return t->manifold->ntotaltriangles + 1;
}

int wavetreesphereface2d_prunable_leaves(const wavetreesphereface2d_t *t)
{
  return multiset_int_total_count(t->S_d);
}

int wavetreesphereface2d_attachable_branches(const wavetreesphereface2d_t *t)
{
  return multiset_int_total_count(t->S_b);
}

int wavetreesphereface2d_coeff_count(const wavetreesphereface2d_t *t)
{
  return multiset_int_double_total_count(t->S_v);
}

int wavetreesphereface2d_map_to_array(const wavetreesphereface2d_t *t, double *a, int n)
{
  int d;
  int c;

  int i;
  int index;
  double value;
  double mean;

  /*
   * Root is mean of first level
   */
  if (multiset_int_double_nth_element(t->S_v, 0, 0, &index, &mean) < 0) {
    ERROR("failed to get root element");
    return -1;
  }

  if (index != 0) {
    ERROR("invalid index at root %d", index);
    return -1;
  }

  for (i = 0; i < t->base_triangles; i ++) {
    a[i] = mean;
  }

  /*
   * Depth 1 is mean + value
   */
  d = 1;
  c = multiset_int_double_depth_count(t->S_v, d);
  if (c > 0) {
    for (i = 0; i < c; i ++) {
      if (multiset_int_double_nth_element(t->S_v, d, i, &index, &value) < 0) {
	ERROR("failed to get nth element");
	return -1;
      }
      
      if (index < 0 || index >= n) {
	ERROR("index out of range %d (%d)", index, n);
	return -1;
      }
      
      a[index - 1] += value;
    }
  }

  /*
   * Rest follow
   */
  for (d = 2; d <= t->degree_max; d ++) {
    
    c = multiset_int_double_depth_count(t->S_v, d);
    if (c > 0) {
      
      for (i = 0; i < c; i ++) {
	
	if (multiset_int_double_nth_element(t->S_v, d, i, &index, &value) < 0) {
	  ERROR("failed to get nth element");
	  return -1;
	}
	
	if (index < 0 || index >= n) {
	  ERROR("index out of range %d (%d)", index, n);
	  return -1;
	}
	
	a[index - 1] = value;
      }
    }
  }

  return 0;
}

int
wavetreesphereface2d_propose_value(wavetreesphereface2d_t *t,
				   int i,
				   int d,
				   double value)
{
  if (t == NULL) {
    ERROR("null tree");
    return -1;
  }

  if (i < 0) {
    ERROR("invalid parameters");
    return -1;
  }

  if (multiset_int_double_is_element(t->S_v, i, d) == 0) {
    ERROR("node not active %d %d", 
	  i, 
	  multiset_int_double_total_count(t->S_v));
    return -1;
  }

  /*
   * Store the old value in undo information
   */
  t->undo = UNDO_VALUE;
  t->u_i = i;
  t->u_d = d;
  if (wavetreesphereface2d_get_coeff(t, i, &(t->u_v)) < 0) {
    ERROR("failed to get coefficient");
    return -1;
  }

  /*
   * Write the new value
   */
  if (multiset_int_double_set(t->S_v, i, d, value) < 0) {
    ERROR("failed to set new value");
    return -1;
  }

  /*
   * Record the proposal
   */
  t->last_step.header.type = CH_VALUE;
  t->last_step.header.accepted = 0;
  t->last_step.header.likelihood = 0.0;
  t->last_step.header.temperature = 0.0;
  t->last_step.header.hierarchical = 0.0;
  
  t->last_step.perturbation.value.node_depth = t->u_d;
  t->last_step.perturbation.value.node_id = t->u_i;
  t->last_step.perturbation.value.new_value = value;
  t->last_step.perturbation.value.old_value = t->u_v;

  return 0;
}

int
wavetreesphereface2d_propose_birth(wavetreesphereface2d_t *t,
				   int i,
				   int d,
				   double value)
{
  if (t == NULL) {
    ERROR("null tree");
    return -1;
  }

  if (i < 0) {
    ERROR("invalid parameters");
    return -1;
  }
  
  if (multiset_int_double_is_element(t->S_v, i, d)) {
    ERROR("non empty node %d %d", i, d);
    return -1;
  }

  if (wavetreesphereface2d_depthofindex(t, i) != d) {
    ERROR("depth mismatch, got %d, expected %d for index %d",
	  d, wavetreesphereface2d_depthofindex(t, i), i);
    return -1;
  }

  /*
   * Store the undo information
   */
  t->undo = UNDO_BIRTH;
  t->u_i = i;
  t->u_d = d;

  /*
   * Record the proposal
   */
  t->last_step.header.type = CH_BIRTH;
  t->last_step.header.accepted = 0;
  t->last_step.header.likelihood = 0.0;
  t->last_step.header.temperature = 0.0;
  t->last_step.header.hierarchical = 0.0;
  
  t->last_step.perturbation.birth.node_depth = t->u_d;
  t->last_step.perturbation.birth.node_id = t->u_i;
  t->last_step.perturbation.birth.new_value = value;

  return add_node(t, i, d, value);
}

int 
wavetreesphereface2d_propose_death(wavetreesphereface2d_t *t,
				   int i,
				   int d, 
				   double *old_value)
{
  if (t == NULL) {
    ERROR("null tree");
    return -1;
  }

  if (i < 0) {
    ERROR("invalid parameters");
    return -1;
  }

  if (multiset_int_double_is_element(t->S_v, i, d) == 0) {
    ERROR("empty node");
    return -1;
  }

  if (wavetreesphereface2d_get_coeff(t, i, &(t->u_v)) < 0) {
    ERROR("failed to get coefficient");
    return -1;
  }
  *old_value = t->u_v;
  t->undo = UNDO_DEATH;
  t->u_i = i;
  t->u_d = d;

  if (remove_node(t, i, d) < 0) {
    ERROR("failed to remove node");
    return -1;
  }

  /*
   * Record the proposal.
   */
  t->last_step.header.type = CH_DEATH;
  t->last_step.header.accepted = 0;
  t->last_step.header.likelihood = 0.0;
  t->last_step.header.temperature = 0.0;
  t->last_step.header.hierarchical = 0.0;
  
  t->last_step.perturbation.death.node_depth = t->u_d;
  t->last_step.perturbation.death.node_id = t->u_i;
  t->last_step.perturbation.death.old_value = t->u_v;

  return 0;
}

int 
wavetreesphereface2d_undo(wavetreesphereface2d_t *t)
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

  case UNDO_MOVE:
    if (remove_node(t, t->u_j, t->u_d) < 0) {
      ERROR("failed to remove node for undo move");
      return -1;
    }
    if (add_node(t, t->u_i, t->u_d, t->u_v) < 0) {
      ERROR("failed to add node for undo move");
      return -1;
    }
    break;

  default:
    ERROR("invalid undo type");
    return -1;
  }

  t->undo = UNDO_NONE;
  t->u_i = 0;
  t->u_j = 0;
  t->u_d = 0;
  t->u_v = 0.0;
    
  return 0;
}

int
wavetreesphereface2d_commit(wavetreesphereface2d_t *t)
{
  if (t->undo == UNDO_NONE) {
    ERROR("nothing to commit");
    return -1;
  }

  switch (t->undo) {
  case UNDO_VALUE:
    t->last_step.header.accepted = 1;
    break;

  case UNDO_BIRTH:
    t->last_step.header.accepted = 1;
    break;

  case UNDO_DEATH:
    t->last_step.header.accepted = 1;
    break;

  case UNDO_MOVE:
    t->last_step.header.accepted = 1;
    break;

  default:
    ERROR("invalid undo type for recording last move");
    return -1;
  }

  t->undo = UNDO_NONE;
  t->u_i = 0;
  t->u_j = 0;
  t->u_d = 0;
  t->u_v = 0.0;
  
  return 0;
}

int wavetreesphereface2d_valid(wavetreesphereface2d_t *t)
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

  for (d = 1; d <= t->degree_max; d ++) {
    n = multiset_int_double_depth_count(t->S_v, d);
    for (i = 0; i < n; i ++) {

      if (multiset_int_double_nth_element(t->S_v, d, i, &index, &value) < 0) {
	ERROR("failed to get nth element of S_v");
	return 0;
      }

      pi = wavetreesphereface2d_parent_index(t, index);
      if (!multiset_int_double_is_element(t->S_v, pi, d - 1)) {
	error_count ++;
	ERROR("parent %d not set for index %d", pi, index);
      }
    }
  }

  /*
   * Check S_b, elements in S_b must not be in S_v or S_d, parent of S_b must be in S_v
   */
  for (d = 1; d <= t->degree_max; d ++) {
    n = multiset_int_depth_count(t->S_b, d);
    for (i = 0; i < n; i ++) {

      if (multiset_int_nth_element(t->S_b, d, i, &index) < 0) {
	ERROR("failed to get nth element of S_b");
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

      pi = wavetreesphereface2d_parent_index(t, index);
      if (!multiset_int_double_is_element(t->S_v, pi, d - 1)) {
	error_count ++;
	ERROR("index %d parent %d (%d) doesn't exist in S_v", index, pi, d - 1);
      }
    }
  }

  /*
   * Check S_d, elements in S_d must be in S_v and not in S_b
   */
  for (d = 1; d <= t->degree_max; d ++) {
    n = multiset_int_depth_count(t->S_d, d);
    for (i = 0; i < n; i ++) {

      if (multiset_int_nth_element(t->S_d, d, i, &index) < 0) {
	ERROR("failed to get nth element of S_b");
	return 0;
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

/* Parent index of a given node */
int wavetreesphereface2d_parent_index(const wavetreesphereface2d_t *t, int c)
{
  int d;

  if (c == 0) {
    /* 0th has no parent */
    return -1;
  }

  d = wavetreesphereface2d_depthofindex(t, c);
  if (d == 1) {
    return 0;
  }

  return wavetreesphereface2d_depth_base(t, d - 1) + (c - wavetreesphereface2d_depth_base(t, d))/4;
}

int wavetreesphereface2d_child_count(const wavetreesphereface2d_t *t, int index, int depth)
{
  int j;
  int nchildren;
  int cc;

  cc = 0;

  if (wavetreesphereface2d_child_indices(t, index, depth, t->child_indices, &nchildren, t->max_children) < 0) {
    return -1;
  }

  for (j = 0; j < nchildren; j ++) {
    if (multiset_int_double_is_element(t->S_v, t->child_indices[j], depth + 1)) {
      cc ++;
    }
  }

  return cc;
}

int wavetreesphereface2d_depthofindex(const wavetreesphereface2d_t *t,
				      int i)
{
  int d;
  
  if (i == 0) {
    return 0;
  }

  d = 1;
  while (wavetreesphereface2d_depth_base(t, d) <= i) {
    d ++;
  }
  
  return d - 1;
}

int wavetreesphereface2d_max_child_count(wavetreesphereface2d_t *t)
{
  return t->max_children;
}

int
wavetreesphereface2d_child_indices(const wavetreesphereface2d_t *t,
				   int index,
				   int depth,
				   int *indices,
				   int *n,
				   int nmax)
{
  int i;
  
  int base;
  int child_base;
  
  if (depth == 0) {
    if (index != 0) {
      ERROR("error depth 0 and index is not 0");
      return -1;
    }

    if (nmax < t->base_triangles) {
      ERROR("error nmax too small %d < %d", nmax, t->base_triangles);
      return -1;
    }

    for (i = 0; i < t->base_triangles; i ++) {
      indices[i] = i + 1;
    }
    *n = t->base_triangles;

    return 0;
  } else {

    if (nmax < 3) {
      ERROR("error nmax too small %d < 3", nmax);
      return -1;
    }

    base = index - wavetreesphereface2d_depth_base(t, depth);

    child_base = wavetreesphereface2d_depth_base(t, depth + 1) + 4 * base;

    for (i = 0; i < 3; i ++) {
      indices[i] = child_base + i;
    }
    *n = 3;

    return 0;
  }

  ERROR("unreachable");
  return -1;
}

static int add_node(wavetreesphereface2d_t *t, int i, int d, double coeff)
{
  int j;
  int pcc;
  int nchildren;

  if (multiset_int_double_insert(t->S_v, i, d, coeff) < 0) {
    ERROR("failed to insert into S_v");
    return -1;
  }

  /*
   * Update parent CC
   */
  j = wavetreesphereface2d_parent_index(t, i);

  pcc = wavetreesphereface2d_child_count(t, j, d - 1);
  if (pcc < 0) {
    ERROR("failed to get pcc");
    return -1;
  }

  if (pcc == 1) {
    /* Remove parent from Sd, only needs to be done if this is its first child */
    multiset_int_remove(t->S_d, j, d - 1);
  }

  if (d == 0) {
    ERROR("0 depth");
    return -1;
  }
  multiset_int_insert(t->S_d, i, d);
  multiset_int_remove(t->S_b, i, d);

  if (wavetreesphereface2d_child_indices(t, i, d, t->child_indices, &nchildren, t->max_children) < 0) {
    return -1;
  }

  for (j = 0; j < nchildren; j ++) {
    multiset_int_insert(t->S_b, t->child_indices[j], d + 1);
  }
  
  /*printf("Add: %d %d\n", oset_int_count(t->S_b), oset_int_count(t->S_d)); */

  return 0;
}

static int remove_node(wavetreesphereface2d_t *t, int i, int d)
{
  int j;
  int pcc;
  int nchildren;
  
  if (d == 0) {
    ERROR("zero depth");
    return -1;
  }

  //printf("removed node: %d\n", i);

  if (multiset_int_double_remove(t->S_v, i, d) < 0) {
    ERROR("failed to remove from S_v");
    return -1;
  }

  /*
   * Update parent CC
   */
  j = wavetreesphereface2d_parent_index(t, i);

  pcc = wavetreesphereface2d_child_count(t, j, d - 1);

  /*
   * Update sets Sb and Sd
   */
  if (j > 0 && pcc == 0) {
    /* Add parent to S_d as it now has no children */
    multiset_int_insert(t->S_d, j, d - 1);
  }
  

  multiset_int_remove(t->S_d, i, d);
  multiset_int_insert(t->S_b, i, d);

  /* Remove children from S_b */
  if (wavetreesphereface2d_child_indices(t, i, d, t->child_indices, &nchildren, t->max_children) < 0) {
    return -1;
  }

  for (j = 0; j < nchildren; j ++) {
    multiset_int_remove(t->S_b, t->child_indices[j], d + 1);
  }
  
  /*  printf("Remove: %d %d\n", t->N_b, t->N_d); */

  return 0;
}

void wavetreesphereface2d_print_setinfo(wavetreesphereface2d_t *t)
{
  printf("Nb %d Nd %d CC %d\n", 
	 multiset_int_total_count(t->S_b), 
	 multiset_int_total_count(t->S_d), 
	 wavetreesphereface2d_coeff_count(t));
}

void wavetreesphereface2d_dump_sets(wavetreesphereface2d_t *t)
{
  printf("Sb:\n  ");
  multiset_int_dump(t->S_b);
  printf("Sd:\n  ");
  multiset_int_dump(t->S_d);
  printf("Sv:\n  ");
  multiset_int_double_dump(t->S_v);
}

void wavetreesphereface2d_dump_coeffs(wavetreesphereface2d_t *t)
{
  multiset_int_double_dump(t->S_v);
}

int wavetreesphereface2d_get_indices(wavetreesphereface2d_t *t, int *set, int *n)
{
  int d;
  int dn;
  int i;
  int j;
  double v;

  j = 0;

  for (d = 0; d <= t->degree_max; d ++) {

    dn = multiset_int_double_depth_count(t->S_v, d);
    if (dn == 0) {
      break;
    }

    for (i = 0; i < dn; i ++) {
      if (multiset_int_double_nth_element(t->S_v, d, i, &(set[j]), &v) < 0) {
	return -1;
      }
      j ++;
    }
  }

  *n = j;
  return 0;
}

int wavetreesphereface2d_maxdepth(const wavetreesphereface2d_t *t)
{
  return t->degree_max;
}

int wavetreesphereface2d_depth(wavetreesphereface2d_t *t)
{
  int d;
  int n;
  for (d = 0; d < t->degree_max; d ++) {
    n = multiset_int_double_depth_count(t->S_v, d);
    if (n == 0) {
      return d;
    }
  }

  return t->degree_max;
}

int wavetreesphereface2d_update_histogram(const wavetreesphereface2d_t *t, coefficient_histogram_t *hist)
{
  int d;
  int n;
  int i;

  int index;
  double value;
  
  for (d = 0; d < t->degree_max; d ++) {
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

int wavetreesphereface2d_depth_filter(void *wt, int subset, int index)
{
  wavetreesphereface2d_t *t = (wavetreesphereface2d_t*)wt;

  return wavetreesphereface2d_depthofindex(t, index) == subset;
}

const multiset_int_double_t *
wavetreesphereface2d_get_S_v(const wavetreesphereface2d_t *t)
{
  return t->S_v;
}

int
wavetreesphereface2d_set_from_S_v(wavetreesphereface2d_t *t,
				  const multiset_int_double_t *S_vp)
{
  int d;
  int n;
  int i;
  int index;
  double value;
  
  multiset_int_clear(t->S_b);
  multiset_int_clear(t->S_d);
  multiset_int_double_clear(t->S_v);

  if (multiset_int_double_get(S_vp, 0, 0, &value) < 0) {
    ERROR("input set doesn't have 0,0 element");
    return -1;
  }

  if (wavetreesphereface2d_initialize(t, value) < 0) {
    ERROR("failed to re-initialise tree");
    return -1;
  }
  
  for (d = 1; d <= t->degree_max; d ++) {
    n = multiset_int_double_depth_count(S_vp, d);
    if (n > 0) {

      for (i = 0; i < n; i ++) {

	if (multiset_int_double_nth_element(S_vp, d, i, &index, &value) < 0) {
	  ERROR("failed to get nth element");
	  return -1;
	}

	/*
	 * We have to use add_node to keep our S_d/S_b sets consistent
	 */
	if (add_node(t, index, d, value) < 0) {
	  ERROR("failed to add node");
	  return -1;
	}
      }
    }
  }

  return 0;
}

int
wavetreesphereface2d_get_last_perturbation(wavetreesphereface2d_t *t,
					   chain_history_change_t *step)
{
  memcpy(step, &t->last_step, sizeof(chain_history_change_t));

  return 0;
}

int
wavetreesphereface2d_set_invalid_perturbation(wavetreesphereface2d_t *t, wavetree_perturb_t p)
{
  memset(&(t->last_step), 0, sizeof(chain_history_change_t));
  t->last_step.header.type = p;

  return 0;
}

int
wavetreesphereface2d_perturb(wavetreesphereface2d_t *t,
			     wavetreesphereface2d_perturb_func_t f,
			     void *user,
			     double *prior_ratio)
{
  int d;
  int n;
  int i;
  int index;
  double value;
  double old_value;

  int pi;
  double pvalue;

  for (d = 0; d <= t->degree_max; d ++) {
    n = multiset_int_double_depth_count(t->S_v, d);
    for (i = 0; i < n; i ++) {

      if (multiset_int_double_nth_element(t->S_v, d, i, &index, &value) < 0) {
        ERROR("failed to get nth element");
        return -1;
      }

      /* Get parent coeff */
      if (index == 0) {
        pvalue = 0.0;
      } else {
        pi = wavetreesphereface2d_parent_index(t, index);
        if (wavetreesphereface2d_get_coeff(t, pi, &pvalue) < 0) {
          ERROR("failed to get parent coeff");
          return -1;
        }
      }
       
      old_value = value;
      if (f(user, index, 0, d, t->degree_max, pvalue, &value, prior_ratio) < 0) {
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

int
wavetreesphereface2d_get_coeff(const wavetreesphereface2d_t *t,
			       int i,
			       double *coeff)
{
  int d;

  d = wavetreesphereface2d_depthofindex(t, i);
  
  if (multiset_int_double_get(t->S_v, i, d, coeff) < 0) {
    return -1;
  }

  return 0;
}

/*
 * Functions for setting up a birth proposal
 */
int
wavetreesphereface2d_choose_birth_depth(const wavetreesphereface2d_t *t,
					double u,
					int maxdepth,
					int *depth,
					double *prob)
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

int
wavetreesphereface2d_reverse_birth_depth(const wavetreesphereface2d_t *t,
					 int depth,
					 int maxdepth,
					 double *prob)
{
  if (t == NULL) {
    return -1;
  }

  *prob = 1.0/((double)multiset_int_nonempty_count(t->S_d, maxdepth));
  return 0;
}

int
wavetreesphereface2d_choose_birth(const wavetreesphereface2d_t *t,
				  int depth,
				  double u,
				  int *coeff,
				  double *prob)
{
  int nindices;

  if (t == NULL) {
    return -1;
  }

  if (multiset_int_choose_index(t->S_b, depth, u, coeff, &nindices) < 0) {
    return -1;
  }

  *prob = 1.0/(double)(nindices);
  return 0;
}

int
wavetreesphereface2d_choose_birth_global(const wavetreesphereface2d_t *t,
					 double u,
					 int maxdepth,
					 int *depth,
					 int *coeff,
					 double *prob)
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

  return 0;
}



int
wavetreesphereface2d_reverse_birth(const wavetreesphereface2d_t *t,
				   int depth,
				   int coeff,
				   double *prob)
{
  if (t == NULL) {
    return -1;
  }
  
  *prob = 1.0/(double)(multiset_int_depth_count(t->S_d, depth));
  return 0;
}

int
wavetreesphereface2d_reverse_birth_global(const wavetreesphereface2d_t *t,
					  int maxdepth,
					  int depth,
					  int coeff,
					  double *prob)
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
int
wavetreesphereface2d_choose_death_depth(const wavetreesphereface2d_t *t,
					double u,
					int maxdepth,
					int *depth,
					double *prob)
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

int
wavetreesphereface2d_reverse_death_depth(const wavetreesphereface2d_t *t,
					 int depth,
					 int maxdepth,
					 double *prob)
{
  if (t == NULL) {
    return -1;
  }

  *prob = 1.0/((double)multiset_int_nonempty_count(t->S_b, maxdepth));
  return 0;
}

int
wavetreesphereface2d_choose_death(const wavetreesphereface2d_t *t,
				  int depth,
				  double u,
				  int *coeff,
				  double *prob)
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

int
wavetreesphereface2d_choose_death_global(const wavetreesphereface2d_t *t,
					 double u,
					 int maxdepth,
					 int *depth,
					 int *coeff,
					 double *prob)
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

int
wavetreesphereface2d_reverse_death(const wavetreesphereface2d_t *t,
				   int depth,
				   int coeff,
				   double *prob)
{
  if (t == NULL) {
    return -1;
  }
  
  *prob = 1.0/(double)(multiset_int_depth_count(t->S_b, depth));
  return 0;
}

int
wavetreesphereface2d_reverse_death_global(const wavetreesphereface2d_t *t,
					  int maxdepth,
					  int depth,
					  int coeff,
					  double *prob)
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
int
wavetreesphereface2d_choose_value_depth(const wavetreesphereface2d_t *t,
					double u,
					int maxdepth,
					int *depth,
					double *prob)
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

int
wavetreesphereface2d_choose_value(const wavetreesphereface2d_t *t,
				  int depth,
				  double u,
				  int *coeff,
				  double *prob)
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

int
wavetreesphereface2d_choose_value_global(const wavetreesphereface2d_t *t,
					 double u,
					 int maxdepth,
					 int *depth,
					 int *coeff,
					 double *prob)
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

int
wavetreesphereface2d_ncoeff_at_depth(int depth)
{
  if (depth < 0) {
    return -1;
  }
  
  if (depth == 0) {
    return 1;
  }

  if (depth == 1) {
    return 5;
  }

  if (depth == 2) {
    return 25;
  }

  return
    4 * (wavetreesphereface2d_ncoeff_at_depth(depth - 1) - wavetreesphereface2d_ncoeff_at_depth(1)) +
    wavetreesphereface2d_ncoeff_at_depth(depth - 1);
}
   

int
wavetreesphereface2d_ntriangle_at_depth(int depth)
{
  if (depth < 0) {
    return -1;
  }

  if (depth <= 2) {
    return 20;
  }

  return 4 * wavetreesphereface2d_ntriangle_at_depth(depth - 1);
}

int
wavetreesphereface2d_depth_base(const wavetreesphereface2d_t *t, int depth)
{
  if (depth < 0) {
    return -1;
  }
  
  if (depth == 0) {
    return 0;
  }

  if (depth == 1) {
    return 1;
  }

  return wavetreesphereface2d_depth_base(t, depth - 1) + t->manifold->ntrianglesatdepth(depth - 2);
}
