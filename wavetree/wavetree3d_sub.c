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
#include <math.h>

#include "wavetree3d_sub.h"

#include "multiset_int.h"
#include "multiset_int_double.h"

#include "slog.h"

typedef enum {
  UNDO_NONE = WT_PERTURB_NONE,
  UNDO_VALUE = WT_PERTURB_VALUE,
  UNDO_BIRTH = WT_PERTURB_BIRTH,
  UNDO_DEATH = WT_PERTURB_DEATH,
  UNDO_MOVE = WT_PERTURB_MOVE
} wavetree3d_sub_undo_t;

struct _wavetree3d_sub {


  int base_width;
  int base_height;
  int base_depth;
  int base_size;
  int *base_indices;

  int max_children;
  int *child_indices;
  
  int degree_max;
  int degree_min;

  int degree_width;
  int degree_height;
  int degree_depth;
  
  int width;
  int height;
  int depth;
  int size;

  double alpha;

  multiset_int_double_t *S_v;
  multiset_int_t *S_b;
  multiset_int_t *S_d;

  wavetree3d_sub_undo_t undo;
  int u_i;
  int u_d;
  double u_v;

  chain_history_change_t last_step;
};

static int add_node(wavetree3d_sub_t *t, int i, int d, double coeff);
static int remove_node(wavetree3d_sub_t *t, int i, int d);

wavetree3d_sub_t *wavetree3d_sub_create(int degree_width,
					int degree_height,
					int degree_depth,
					double alpha)
{
  wavetree3d_sub_t *r;
  int i;
  int j;
  int k;
  int l;

  r = malloc(sizeof(wavetree3d_sub_t));
  if (r == NULL) {
    return NULL;
  }

  r->degree_min = degree_width;
  if (degree_height < r->degree_min) {
    r->degree_min = degree_height;
  }
  if (degree_depth < r->degree_min) {
    r->degree_min = degree_depth;
  }

  r->degree_width = degree_width;
  r->degree_height = degree_height;
  r->degree_depth = degree_depth;

  r->width = 1 << degree_width;
  r->height = 1 << degree_height;
  r->depth = 1 << degree_depth;

  r->base_width = 1 << (degree_width - r->degree_min);
  r->base_height = 1 << (degree_height - r->degree_min);
  r->base_depth = 1 << (degree_depth - r->degree_min);
  r->base_size = r->base_width * r->base_height * r->base_depth;
  if (r->base_size <= 0) {
    ERROR("failed to compute base size");
    return NULL;
  }
  if (r->base_size == 1) {
    r->base_indices = NULL;
  } else {
    r->base_indices = malloc(sizeof(int) * r->base_size);
    if (r->base_indices == NULL) {
      ERROR("failed to allocate base indices");
      return NULL;
    }

    l = 0;
    for (k = 0; k < r->base_depth; k ++) {
      for (j = 0; j < r->base_height; j ++) {
	for (i = 0; i < r->base_width; i ++, l ++) {

	  r->base_indices[l] = i + j*r->width + k*r->width*r->height + 1;

	}
      }
    }
  }

  r->size = r->width * r->height * r->depth;
  
  if (r->base_size == 1) {

    r->degree_max = r->degree_width;
    r->max_children = 8;

  } else {
    
    r->degree_max = r->degree_min + 1;
    
    r->max_children = 8;
    if (r->base_size > r->max_children) {
      r->max_children = r->base_size;
    }
    
  }

  r->child_indices = malloc(sizeof(int) * r->max_children);
  if (r->child_indices == NULL) {
    ERROR("failed to allocate child indices");
    return NULL;
  }

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

  memset(&(r->last_step), 0, sizeof(chain_history_change_t));
  r->last_step.header.type = CH_INITIALISE;

  return r;
}

void wavetree3d_sub_destroy(wavetree3d_sub_t *t)
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
wavetree3d_sub_save(const wavetree3d_sub_t *t,
		const char *filename)
{
  FILE *fp;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    ERROR("failed to save wavetree");
    return -1;
  }

  fprintf(fp, "%d %d %d\n", t->degree_width, t->degree_height, t->degree_depth);
  fprintf(fp, "%d %d %d %d\n", t->width, t->height, t->depth, t->size);
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
wavetree3d_sub_load(wavetree3d_sub_t *t,
		    const char *filename)
{
  FILE *fp;

  int degree_width;
  int degree_height;
  int degree_depth;
  int width;
  int height;
  int depth;
  int size;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    ERROR("failed to open file");
    return -1;
  }

  if (fscanf(fp, "%d %d %d\n", &degree_width, &degree_height, &degree_depth) != 3) {
    ERROR("failed to read degree");
    return -1;
  }

  if (degree_width != t->degree_width ||
      degree_height != t->degree_height ||
      degree_depth != t->degree_depth) {
    ERROR("degree mismatch (%d %d %d) != (%d %d %d)",
	    degree_width, degree_height, degree_depth,
	    t->degree_width, t->degree_height, t->degree_depth);
    return -1;
  }

  if (fscanf(fp, "%d %d %d %d\n", &width, &height, &depth, &size) != 4) {
    ERROR("failed to read header");
    return -1;
  }

  if (width != t->width ||
      height != t->height ||
      depth != t->depth ||
      size != t->size) {
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
wavetree3d_sub_encode(wavetree3d_sub_t *t,
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
  /* printf("    encoding header: %d %d %d %f %d\n", */
  /* 	 t->degree_width, */
  /* 	 t->degree_height, */
  /* 	 t->degree_depth, */
  /* 	 t->alpha, */
  /* 	 multiset_int_double_total_count(t->S_v)); */
  
  ENCODEINT(t->degree_width);
  ENCODEINT(t->degree_height);
  ENCODEINT(t->degree_depth);
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
wavetree3d_sub_decode(wavetree3d_sub_t *t,
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
  if (iv != t->degree_width) {
    ERROR("degree width mismatch %d != %d", t->degree_width, iv);
    return -1;
  }
  DECODEINT(&iv);
  if (iv != t->degree_height) {
    ERROR("degree height mismatch %d != %d", t->degree_height, iv);
    return -1;
  }
  DECODEINT(&iv);
  if (iv != t->degree_depth) {
    ERROR("degree depth mismatch %d != %d", t->degree_depth, iv);
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
   * Values
   */  
  for (i = 0; i < n; i ++) {

    DECODEINT(&d);
    DECODEINT(&index);
    DECODEDOUBLE(&v);

    if (add_node(t, index, d, v) < 0) {
      ERROR("wavetree3d_sub_decode: failed to add node (%d %d %f)", index, d, v);
      return -1;
    }
  }

  return 0;
}

int
wavetree3d_sub_get_width(wavetree3d_sub_t *t)
{
  if (t != NULL) {
    return t->width;
  } 
  
  return -1;
}

int
wavetree3d_sub_get_height(wavetree3d_sub_t *t)
{
  if (t != NULL) {
    return t->height;
  } 
  
  return -1;
}

int
wavetree3d_sub_get_depth(wavetree3d_sub_t *t)
{
  if (t != NULL) {
    return t->depth;
  } 
  
  return -1;
}

int
wavetree3d_sub_get_size(wavetree3d_sub_t *t)
{
  if (t != NULL) {
    return t->size;
  }

  return -1;
}

int
wavetree3d_sub_get_ncoeff(wavetree3d_sub_t *t)
{
  if (t != NULL) {
    if (t->base_size > 1) {
      return t->size + 1;
    } else {
      return t->size;
    }
  }
	
  return -1;
}


int wavetree3d_sub_initialize(wavetree3d_sub_t *t,
			      double dc)
{
  int d;
  int i;

  /*
   * Clear all sets
   */
  multiset_int_double_clear(t->S_v);
  multiset_int_clear(t->S_b);
  multiset_int_clear(t->S_d);

  /*
   * Insert the dc value
   */
  d = wavetree3d_sub_depthofindex(t, 0);
  multiset_int_double_insert(t->S_v,
			     0,
			     d,
			     dc);

  /*
   * Add the potential birthing nodes
   */
  if (t->base_size == 1) {
    multiset_int_insert(t->S_b,
			wavetree3d_sub_UTR(t, 0),
			d + 1);
    multiset_int_insert(t->S_b,
			wavetree3d_sub_UBL(t, 0),
			d + 1);
    multiset_int_insert(t->S_b,
			wavetree3d_sub_UBR(t, 0),
			d + 1);
    multiset_int_insert(t->S_b,
			wavetree3d_sub_LTL(t, 0),
			d + 1);
    multiset_int_insert(t->S_b,
			wavetree3d_sub_LTR(t, 0),
			d + 1);
    multiset_int_insert(t->S_b,
			wavetree3d_sub_LBL(t, 0),
			d + 1);
    multiset_int_insert(t->S_b,
			wavetree3d_sub_LBR(t, 0),
			d + 1);
    
    if (multiset_int_double_total_count(t->S_v) != 1) {
      ERROR("failed to initialize S_v");
      return -1;
    }
    
    if (multiset_int_total_count(t->S_b) != 7) {
      ERROR("failed to initialize S_b");
      return -1;
    }
  } else {
    for (i = 0; i < t->base_size; i ++) {
      multiset_int_insert(t->S_b,
			  t->base_indices[i],
			  d + 1);
    }
    if (multiset_int_double_total_count(t->S_v) != 1) {
      ERROR("failed to initialize S_v");
      return -1;
    }
    
    if (multiset_int_total_count(t->S_b) != t->base_size) {
      ERROR("failed to initialize S_b");
      return -1;
    }
  }    

  return 0;
}

double wavetree3d_sub_dc(wavetree3d_sub_t *t)
{
  double dc;

  if (multiset_int_double_get(t->S_v, 0, 0, &dc) < 0) {
    return -1.0;
  }

  return dc;
}

int wavetree3d_sub_prunable_leaves(const wavetree3d_sub_t *t)
{
  return multiset_int_total_count(t->S_d);
}


int wavetree3d_sub_attachable_branches(const wavetree3d_sub_t *t)
{
  return multiset_int_total_count(t->S_b);
}

int wavetree3d_sub_coeff_count(const wavetree3d_sub_t *t)
{
  return multiset_int_double_total_count(t->S_v);
}

int wavetree3d_sub_map_to_array(wavetree3d_sub_t *t, 
				double *a, 
				int n)
{
  int d;
  int c;

  int i;
  int index;
  double value;
  double mean;

  if (t->base_size == 1) {

    for (d = 0; d <= t->degree_max; d ++) {
      
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
	  
	  a[index] = value;
	}
      }
    }
  } else {

    if (multiset_int_double_nth_element(t->S_v, 0, 0, &index, &mean) < 0) {
      ERROR("failed to get dc");
      return -1;
    }

    if (wavetree3d_sub_child_indices(t, 0, 0, t->child_indices, &c, t->max_children) < 0) {
      ERROR("failed to get 1st level children");
      return -1;
    }

    for (i = 0; i < c; i ++) {
      a[t->child_indices[i] - 1] = mean;
    }

    c = multiset_int_double_depth_count(t->S_v, 1);
    for (i = 0; i < c; i ++) {

      if (multiset_int_double_nth_element(t->S_v, 1, i, &index, &value) < 0) {
	ERROR("failed to get nth element");
	return -1;
      }

      index --;
	  
      if (index < 0 || index >= n) {
	ERROR("index out of range %d (%d)", index, n);
	return -1;
      }
	  
      a[index] += value;
    }

    for (d = 2; d <= t->degree_max; d ++) {
      
      c = multiset_int_double_depth_count(t->S_v, d);
      if (c > 0) {
	
	for (i = 0; i < c; i ++) {
	  
	  if (multiset_int_double_nth_element(t->S_v, d, i, &index, &value) < 0) {
	    ERROR("failed to get nth element");
	    return -1;
	  }

	  index --;
	  
	  if (index < 0 || index >= n) {
	    ERROR("index out of range %d (%d)", index, n);
	    return -1;
	  }
	  
	  a[index] = value;
	}
      }
    }
  }
  return 0;
}

int
wavetree3d_sub_propose_value(wavetree3d_sub_t *t,
			     int i,
			     int d,
			     double value)
{
  double old_value;

  if (t == NULL) {
    ERROR("null tree");
    return -1;
  }

  if (t->base_size == 1) {
    if (i < 0 || i >= t->size) {
      ERROR("invalid parameters %d %d %d", i, d, t->size);
      return -1;
    }
  } else {
    if (i < 0 || i > t->size) {
      ERROR("invalid parameters %d %d %d", i, d, t->size);
      return -1;
    }
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
wavetree3d_sub_propose_birth(wavetree3d_sub_t *t,
			     int i,
			     int d,
			     double value)
{
  if (t == NULL) {
    ERROR("null tree");
    return -1;
  }

  if (t->base_size == 1) {
    if (i < 0 || i >= t->size) {
      ERROR("invalid parameters %d %d %d", i, d, t->size);
      return -1;
    }
  } else {
    if (i < 0 || i > t->size) {
      ERROR("invalid parameters %d %d %d", i, d, t->size);
      return -1;
    }
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
wavetree3d_sub_propose_death(wavetree3d_sub_t *t,
			     int i,
			     int d,
			     double *old_value)
{
  if (t == NULL) {
    ERROR("null tree");
    return -1;
  }

  if (t->base_size == 1) { 
    if (i < 0 ||
	i >= t->size) {
      ERROR("invalid parameters");
      return -1;
    }
  } else {
    if (i < 0 ||
	i > t->size) {
      ERROR("invalid parameters");
      return -1;
    }
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
wavetree3d_sub_undo(wavetree3d_sub_t *t)
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
wavetree3d_sub_commit(wavetree3d_sub_t *t)
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

  default:
    ERROR("invalid undo type for recording last move");
    return -1;
  }
  
  t->undo = UNDO_NONE;
  t->u_i = 0;
  t->u_d = 0;
  t->u_v = 0.0;
  
  return 0;
}

int wavetree3d_sub_valid(wavetree3d_sub_t *t)
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

      pi = wavetree3d_sub_parent_index(t, index);
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

      pi = wavetree3d_sub_parent_index(t, index);
      if (!multiset_int_double_is_element(t->S_v, pi, d - 1)) {
	error_count ++;
	ERROR("index %d parent %d doesn't exist in S_v", index, pi);
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

void wavetree3d_sub_print_setinfo(wavetree3d_sub_t *t)
{
  printf("Nb %d Nd %d Nv %d\n", 
         multiset_int_total_count(t->S_b), 
         multiset_int_total_count(t->S_d), 
         multiset_int_double_total_count(t->S_v));
}

void wavetree3d_sub_dump_sets(wavetree3d_sub_t *t)
{
  printf("Sb:\n  ");
  multiset_int_dump(t->S_b);
  printf("Sd:\n  ");
  multiset_int_dump(t->S_d);
  printf("Sv:\n  ");
  multiset_int_double_dump(t->S_v);
}

void wavetree3d_sub_dump_coeffs(wavetree3d_sub_t *t)
{
  multiset_int_double_dump(t->S_v);
}

int wavetree3d_sub_depth(wavetree3d_sub_t *t) 
{
  int d;
  
  for (d = 1; d <= t->degree_max; d ++) {
    if (multiset_int_double_depth_count(t->S_v, d) == 0) {
      return d - 1;
    }
  }

  return -1;
}

int wavetree3d_sub_maxdepth(wavetree3d_sub_t *t)
{
  return t->degree_max;
}

int wavetree3d_sub_update_histogram(const wavetree3d_sub_t *t, coefficient_histogram_t *hist)
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
	ERROR("wavetree3d_sub_update_histogram: failed to get nth element");
	return -1;
      }

      coefficient_histogram_sample(hist, index, value);
    }
  }

  return 0;
}

int wavetree3d_sub_depth_filter(void *wt, int subset, int index)
{
  wavetree3d_sub_t *t = (wavetree3d_sub_t *)wt;

  return wavetree3d_sub_depthofindex(t, index) == subset;
}

/*
 * Functions for setting up a birth proposal
 */
int 
wavetree3d_sub_choose_birth_depth(const wavetree3d_sub_t *t, double u, int maxdepth, int *depth, double *prob)
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

int wavetree3d_sub_reverse_birth_depth(const wavetree3d_sub_t *t, int depth, int maxdepth, double *prob)
{
  if (t == NULL) {
    return -1;
  }

  *prob = 1.0/((double)multiset_int_nonempty_count(t->S_d, maxdepth));
  return 0;
}

int wavetree3d_sub_choose_birth(const wavetree3d_sub_t *t, int depth, double u, int *coeff, double *prob)
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

int wavetree3d_sub_choose_birth_global(const wavetree3d_sub_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob)
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
    if (multiset_int_choose_index_weighted(t->S_b, u, maxdepth, t->alpha, coeff, depth, prob) < 0) {
      return -1;
    }
  }

  return 0;
}

int wavetree3d_sub_reverse_birth(const wavetree3d_sub_t *t, int depth, int coeff, double *prob)
{
  if (t == NULL) {
    return -1;
  }
  
  *prob = 1.0/(double)(multiset_int_depth_count(t->S_d, depth));
  return 0;
}

int wavetree3d_sub_reverse_birth_global(const wavetree3d_sub_t *t, int maxdepth, int depth, int coeff, double *prob)
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
int wavetree3d_sub_choose_death_depth(const wavetree3d_sub_t *t, double u, int maxdepth, int *depth, double *prob)
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

int wavetree3d_sub_reverse_death_depth(const wavetree3d_sub_t *t, int depth, int maxdepth, double *prob)
{
  if (t == NULL) {
    return -1;
  }

  *prob = 1.0/((double)multiset_int_nonempty_count(t->S_b, maxdepth));
  return 0;
}

int wavetree3d_sub_choose_death(const wavetree3d_sub_t *t, int depth, double u, int *coeff, double *prob)
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


int wavetree3d_sub_choose_death_global(const wavetree3d_sub_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob)
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
    if (multiset_int_choose_index_weighted(t->S_d, u, maxdepth, t->alpha, coeff, depth, prob) < 0) {
      return -1;
    }
  }

  return 0;
}

int wavetree3d_sub_reverse_death(const wavetree3d_sub_t *t, int depth, int coeff, double *prob)
{
  if (t == NULL) {
    return -1;
  }
  
  *prob = 1.0/(double)(multiset_int_depth_count(t->S_b, depth));
  return 0;
}

int wavetree3d_sub_reverse_death_global(const wavetree3d_sub_t *t, int maxdepth, int depth, int coeff, double *prob)
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
int wavetree3d_sub_choose_value_depth(const wavetree3d_sub_t *t, double u, int maxdepth, int *depth, double *prob)
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

int wavetree3d_sub_choose_value(const wavetree3d_sub_t *t, int depth, double u, int *coeff, double *prob)
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


int wavetree3d_sub_choose_value_global(const wavetree3d_sub_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob)
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

const multiset_int_double_t *
wavetree3d_sub_get_S_v(const wavetree3d_sub_t *t)
{
  return t->S_v;
}

int
wavetree3d_sub_set_from_S_v(wavetree3d_sub_t *t,
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

  if (wavetree3d_sub_initialize(t, value) < 0) {
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
wavetree3d_sub_set_from_S_v_filtered(wavetree3d_sub_t *t,
				     const multiset_int_double_t *S_vp,
				     int maxdepth)
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

  if (wavetree3d_sub_initialize(t, value) < 0) {
    ERROR("failed to re-initialise tree");
    return -1;
  }
  
  for (d = 1; d <= maxdepth; d ++) {
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
wavetree3d_sub_get_last_perturbation(wavetree3d_sub_t *t,
				     chain_history_change_t *step)
{
  memcpy(step, &t->last_step, sizeof(chain_history_change_t));

  return 0;
}

int
wavetree3d_sub_set_invalid_perturbation(wavetree3d_sub_t *t, wavetree_perturb_t p)
{
  memset(&(t->last_step), 0, sizeof(chain_history_change_t));
  t->last_step.header.type = p;

  return 0;
}


int wavetree3d_sub_perturb(wavetree3d_sub_t *t,
			   wavetree3d_sub_perturb_func_t f,
			   void *user,
			   double *prior_ratio)
{
  int d;
  int n;
  int i;
  int index;
  double value;
  double old_value;

  int ii;
  int ij;
  int ik;

  int pi;
  double pvalue;

  for (d = 0; d < t->degree_max; d ++) {
    n = multiset_int_double_depth_count(t->S_v, d);
    for (i = 0; i < n; i ++) {

      if (multiset_int_double_nth_element(t->S_v, d, i, &index, &value) < 0) {
	ERROR("failed to get nth element");
	return -1;
      }

      if (wavetree3d_sub_3dindices(t, i, &ii, &ij, &ik) < 0) {
	ERROR("failed to get 3d indices");
	return -1;
      }

      /* Get parent coeff */
      if (index == 0) {
	pvalue = 0.0;
      } else {
	pi = wavetree3d_sub_parent_index(t, index);
	if (wavetree3d_sub_get_coeff(t, pi, d - 1, &pvalue) < 0) {
	  ERROR("failed to get parent coeff");
	  return -1;
	}
      }
       
      old_value = value;
      if (f(user, ii, ij, ik, d, t->degree_max, pvalue, &value, prior_ratio) < 0) {
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


int wavetree3d_sub_parent_index(const wavetree3d_sub_t *t, int c)
{
  int i;
  int j;
  int k;

  if (c == 0) {
    /* 0th has no parent */
    return -1;
  }

  if (t->base_size > 1) {

    if (wavetree3d_sub_3dindices(t, c, &i, &j, &k) < 0) {
      ERROR("failed to get 3d indices");
      return -1;
    }

    if (i < t->base_width &&
	j < t->base_height &&
	k < t->base_depth) {
      return 0;
    }

    if (i < (2*t->base_width) &&
	j < (2*t->base_height) &&
	k < (2*t->base_depth)) {

      i = i % t->base_width;
      j = j % t->base_height;
      k = k % t->base_depth;

      return wavetree3d_sub_from_3dindices(t, i, j, k);
    }

    /* divide by 2 method for levels 2+*/
    j /= 2;
    i /= 2;
    k /= 2;
    
    return wavetree3d_sub_from_3dindices(t, i, j, k);

  } else {

    if (wavetree3d_sub_3dindices(t, c, &i, &j, &k) < 0) {
      ERROR("failed to get 3d indices");
      return -1;
    }
    
    
    j /= 2;
    i /= 2;
    k /= 2;
    
    return wavetree3d_sub_from_3dindices(t, i, j, k);
  }
}

int wavetree3d_sub_3dindices(const wavetree3d_sub_t *t,
			     int i,
			     int *ii,
			     int *ij,
			     int *ik)
{
  int j;

  if (t->base_size == 1) {
    if (i < 0 ||
	i >= t->size) {
      return -1;
    }
  } else {
    if (i == 0) {
      *ii = -1;
      *ij = -1;
      *ik = -1;
      return 0;
    }
    
    if (i > t->size) {
      return -1;
    }

    i --;
  }
    
  j = i % (t->width * t->height);
  
  *ii = j % t->width;
  *ij = (j - (*ii)) / t->width;
  *ik = (i - j)/(t->width * t->height);

  return 0;
}

int wavetree3d_sub_from_3dindices(const wavetree3d_sub_t *t,
				  int ii,
				  int ij,
				  int ik)
{
  if (t == NULL) {
    return -1;
  }

  if (t->base_size == 1) {

    if (ii < 0 || ii >= t->width ||
	ij < 0 || ij >= t->height ||
	ik < 0 || ik >= t->depth) {
      ERROR("indices out of range");
      return -1;
    }
    
    return (ik * t->height + ij)*t->width + ii;

  } else {
    if (ii == (-1) &&
	ij == (-1) &&
	ik == (-1)) {
      return 0;
    }

    if (ii >= t->width ||
	ij >= t->height ||
	ik >= t->depth) {
      ERROR("indices out of range");
      return -1;
    }
    
    return (ik * t->height + ij)*t->width + ii + 1;
  }
}

int wavetree3d_sub_depthofindex(const wavetree3d_sub_t *t,
				int i)
{
  if (i == 0) {
    return 0;
  }

  return 1 + wavetree3d_sub_depthofindex(t, wavetree3d_sub_parent_index(t, i));
}

int wavetree3d_sub_get_coeff(const wavetree3d_sub_t *t,
			     int i,
			     int d,
			     double *coeff)
{
  if (multiset_int_double_get(t->S_v, i, d, coeff) < 0) {
    ERROR("failed to get coeff %d (depth %d)", i, d);
    return -1;
  }

  return 0;
}

static int r_generate_dyck_word(wavetree3d_sub_t *t, int index, char *buffer, int *len, int maxlength)
{
  int j;

  if (multiset_int_double_is_element(t->S_v, index, wavetree3d_sub_depthofindex(t, index)) == 0) {

    if (*len == maxlength) {
      return -1;
    }
    buffer[*len] = '.';
    (*len) ++;

    return 0;
  }

  if (*len == maxlength) {
    return -1;
  }
  buffer[*len] = '(';
  (*len) ++;

  
  j = wavetree3d_sub_UTL(t, index);
  if (j > 0) {
    if (r_generate_dyck_word(t, j, buffer, len, maxlength) < 0) {
      return -1;
    }
  }

  j = wavetree3d_sub_UTR(t, index);
  if (j > 0) {
    if (r_generate_dyck_word(t, j, buffer, len, maxlength) < 0) {
      return -1;
    }
  }

  j = wavetree3d_sub_UBL(t, index);
  if (j > 0) {
    if (r_generate_dyck_word(t, j, buffer, len, maxlength) < 0) {
      return -1;
    }
  }

  j = wavetree3d_sub_UBR(t, index);
  if (j > 0) {
    if (r_generate_dyck_word(t, j, buffer, len, maxlength) < 0) {
      return -1;
    }
  }

  j = wavetree3d_sub_LTL(t, index);
  if (j > 0) {
    if (r_generate_dyck_word(t, j, buffer, len, maxlength) < 0) {
      return -1;
    }
  }

  j = wavetree3d_sub_LTR(t, index);
  if (j > 0) {
    if (r_generate_dyck_word(t, j, buffer, len, maxlength) < 0) {
      return -1;
    }
  }

  j = wavetree3d_sub_LBL(t, index);
  if (j > 0) {
    if (r_generate_dyck_word(t, j, buffer, len, maxlength) < 0) {
      return -1;
    }
  }

  j = wavetree3d_sub_LBR(t, index);
  if (j > 0) {
    if (r_generate_dyck_word(t, j, buffer, len, maxlength) < 0) {
      return -1;
    }
  }

  if (*len == maxlength) {
    return -1;
  }
  buffer[*len] = ')';
  (*len) ++;

  return 0;
}


int wavetree3d_sub_generate_dyck_word(wavetree3d_sub_t *t, char *buffer, int maxlength)
{
  int len;

  len = 0;
  if (r_generate_dyck_word(t, 0, buffer, &len, maxlength) < 0) {
    return -1;
  }

  if (len < maxlength) {
    buffer[len] = '\0';
    return 0;
  }

  return -1;
}

static int db_open(uint64_t *binary, int *len)
{
  (*len) ++;
  return 0;
}

static int db_close(uint64_t *binary, int *len)
{
  (*binary) |= ((uint64_t)1) << (*len);
  (*len) ++;
  return 0;
}

static int db_leaf(uint64_t *binary, int *len)
{
  db_open(binary, len);
  db_close(binary, len);
  return 0;
}

static int r_generate_dyck_binary(wavetree3d_sub_t *t, int index, uint64_t *binary, int *len)
{
  int j;

  if (multiset_int_double_is_element(t->S_v, index, wavetree3d_sub_depthofindex(t, index)) == 0) {
    return db_leaf(binary, len);
  }

  db_open(binary, len);
  
  j = wavetree3d_sub_UTL(t, index);
  if (j > 0) {
    if (r_generate_dyck_binary(t, j, binary, len) < 0) {
      return -1;
    }
  } else {
    if (index != 0) {
      db_leaf(binary, len);
    }
  }

  j = wavetree3d_sub_UTR(t, index);
  if (j > 0) {
    if (r_generate_dyck_binary(t, j, binary, len) < 0) {
      return -1;
    }
  } else {
    db_leaf(binary, len);
  }

  j = wavetree3d_sub_UBL(t, index);
  if (j > 0) {
    if (r_generate_dyck_binary(t, j, binary, len) < 0) {
      return -1;
    }
  } else {
    db_leaf(binary, len);
  }


  j = wavetree3d_sub_UBR(t, index);
  if (j > 0) {
    if (r_generate_dyck_binary(t, j, binary, len) < 0) {
      return -1;
    }
  } else {
    db_leaf(binary, len);
  }
  
  j = wavetree3d_sub_LTL(t, index);
  if (j > 0) {
    if (r_generate_dyck_binary(t, j, binary, len) < 0) {
      return -1;
    }
  } else {
    if (index != 0) {
      db_leaf(binary, len);
    }
  }

  j = wavetree3d_sub_LTR(t, index);
  if (j > 0) {
    if (r_generate_dyck_binary(t, j, binary, len) < 0) {
      return -1;
    }
  } else {
    db_leaf(binary, len);
  }

  j = wavetree3d_sub_LBL(t, index);
  if (j > 0) {
    if (r_generate_dyck_binary(t, j, binary, len) < 0) {
      return -1;
    }
  } else {
    db_leaf(binary, len);
  }


  j = wavetree3d_sub_LBR(t, index);
  if (j > 0) {
    if (r_generate_dyck_binary(t, j, binary, len) < 0) {
      return -1;
    }
  } else {
    db_leaf(binary, len);
  }

  db_close(binary, len);

  return 0;
}

int wavetree3d_sub_generate_dyck_binary(wavetree3d_sub_t *t, uint64_t *binary)
{
  int len;

  *binary = 0;

  len = 0;
  if (r_generate_dyck_binary(t, 0, binary, &len) < 0) {
    return -1;
  }

  return 0;
}

static int add_node(wavetree3d_sub_t *t, int i, int d, double coeff)
{
  int j;
  int pcc;
  int nchildren;

  if (multiset_int_double_insert(t->S_v, i, d, coeff) < 0) {
    return -1;
  }
  
  j = wavetree3d_sub_parent_index(t, i);

  pcc = wavetree3d_sub_child_count(t, j, d - 1);
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
  if (wavetree3d_sub_child_indices(t, i, d, t->child_indices, &nchildren, t->max_children) < 0) {
    return -1;
  }

  for (j = 0; j < nchildren; j ++) {
    multiset_int_insert(t->S_b, t->child_indices[j], d + 1);
  }
  
  return 0;
}

static int remove_node(wavetree3d_sub_t *t, int i, int d)
{
  int j;
  int pcc;
  int nchildren;
  
  if (multiset_int_double_remove(t->S_v, i, d) < 0) {
    ERROR("failed to remove index from S_v (index %d, depth %d)", i, d);
    return -1;
  }

  j = wavetree3d_sub_parent_index(t, i);
  if (j < 0) {
    ERROR("failed to get parent index (index %d, depth %d)", i, d);
    return -1;
  }

  pcc = wavetree3d_sub_child_count(t, j, d - 1);

  if (j > 0 && pcc == 0) {
    /* Add parent to S_d as it now has no children */
    multiset_int_insert(t->S_d, j, d - 1);
  }

  /* Remove node from S_d and add to S_b */
  multiset_int_remove(t->S_d, i, d);
  multiset_int_insert(t->S_b, i, d);

  /* Remove children from S_b */
  if (wavetree3d_sub_child_indices(t, i, d, t->child_indices, &nchildren, t->max_children) < 0) {
    return -1;
  }

  for (j = 0; j < nchildren; j ++) {
    multiset_int_remove(t->S_b, t->child_indices[j], d + 1);
  }
  
  return 0;
}

int wavetree3d_sub_child_count(wavetree3d_sub_t *t, int index, int depth)
{
  int j;
  int nchildren;
  int cc;

  cc = 0;

  if (wavetree3d_sub_child_indices(t, index, depth, t->child_indices, &nchildren, t->max_children) < 0) {
    return -1;
  }

  for (j = 0; j < nchildren; j ++) {
    if (multiset_int_double_is_element(t->S_v, t->child_indices[j], depth + 1)) {
      cc ++;
    }
  }

  return cc;
}

int wavetree3d_sub_max_child_count(wavetree3d_sub_t *t)
{
  return t->max_children;
}

int wavetree3d_sub_child_indices(wavetree3d_sub_t *t, int index, int depth, int *indices, int *n, int nmax)
{
  int i;

  int ii;
  int jj;
  int kk;
  
  if (t->base_size == 1) {
    if (nmax < 8) {
      ERROR("nmax insufficient (8, %d)", nmax);
      return -1;
    }
    i = 0;
    indices[i] = wavetree3d_sub_UTL(t, index);
    if (indices[i] > 0) {
      i ++;
    }
    indices[i] = wavetree3d_sub_UTR(t, index);
    if (indices[i] > 0) {
      i ++;
    }
    indices[i] = wavetree3d_sub_UBL(t, index);
    if (indices[i] > 0) {
      i ++;
    }
    indices[i] = wavetree3d_sub_UBR(t, index);
    if (indices[i] > 0) {
      i ++;
    }
    indices[i] = wavetree3d_sub_LTL(t, index);
    if (indices[i] > 0) {
      i ++;
    }
    indices[i] = wavetree3d_sub_LTR(t, index);
    if (indices[i] > 0) {
      i ++;
    }
    indices[i] = wavetree3d_sub_LBL(t, index);
    if (indices[i] > 0) {
      i ++;
    }
    indices[i] = wavetree3d_sub_LBR(t, index);
    if (indices[i] > 0) {
      i ++;
    }

    *n = i;
    return 0;
    
  } else {

    if (depth == 0) {
      if (index != 0) {
	ERROR("invalid index for depth 0: %d", index);
	return -1;
      }

      if (nmax < t->base_size) {
	ERROR("array too small: %d < %d", nmax, t->base_size);
	return -1;
      }

      for (i = 0; i < t->base_size; i ++) {
	indices[i] = t->base_indices[i];
      }
      *n = t->base_size;

      return 0;
    }

    if (depth == 1) {

      if (wavetree3d_sub_3dindices(t,
				   index,
				   &ii, &jj, &kk) < 0) {
	ERROR("failed to get 3d indices");
	return -1;
      }

      i = 0;
      indices[i] = wavetree3d_sub_from_3dindices(t,
						 ii + t->base_width,
						 jj,
						 kk);
      if (indices[i] > 0) {
	i ++;
      }
      
      indices[i] = wavetree3d_sub_from_3dindices(t,
						 ii,
						 jj + t->base_height,
						 kk);
      if (indices[i] > 0) {
	i ++;
      }
      
      indices[i] = wavetree3d_sub_from_3dindices(t,
						 ii + t->base_width,
						 jj + t->base_height,
						 kk);
      if (indices[i] > 0) {
	i ++;
      }

      indices[i] = wavetree3d_sub_from_3dindices(t,
						 ii,
						 jj,
						 kk + t->base_depth);
      if (indices[i] > 0) {
	i ++;
      }
      indices[i] = wavetree3d_sub_from_3dindices(t,
						 ii + t->base_width,
						 jj,
						 kk + t->base_depth);
      if (indices[i] > 0) {
	i ++;
      }
      
      indices[i] = wavetree3d_sub_from_3dindices(t,
						 ii,
						 jj + t->base_height,
						 kk + t->base_depth);
      if (indices[i] > 0) {
	i ++;
      }
      
      indices[i] = wavetree3d_sub_from_3dindices(t,
						 ii + t->base_width,
						 jj + t->base_height,
						 kk + t->base_depth);
      if (indices[i] > 0) {
	i ++;
      }

      *n = i;
      return 0;

    } else {
      
      if (nmax < 8) {
	ERROR("nmax insufficient (8, %d)", nmax);
	return -1;
      }
      
      i = 0;
      indices[i] = wavetree3d_sub_UTL(t, index);
      if (indices[i] > 0) {
	i ++;
      }
      indices[i] = wavetree3d_sub_UTR(t, index);
      if (indices[i] > 0) {
	i ++;
      }
      indices[i] = wavetree3d_sub_UBL(t, index);
      if (indices[i] > 0) {
	i ++;
      }
      indices[i] = wavetree3d_sub_UBR(t, index);
      if (indices[i] > 0) {
	i ++;
      }
      indices[i] = wavetree3d_sub_LTL(t, index);
      if (indices[i] > 0) {
	i ++;
      }
      indices[i] = wavetree3d_sub_LTR(t, index);
      if (indices[i] > 0) {
	i ++;
      }
      indices[i] = wavetree3d_sub_LBL(t, index);
      if (indices[i] > 0) {
	i ++;
      }
      indices[i] = wavetree3d_sub_LBR(t, index);
      if (indices[i] > 0) {
	i ++;
      }

      *n = i;
      return 0;
    }

  }
    
  ERROR("unreachable");
  return -1;
}

#define WAVETREE3D_SUB_CHILD(t, i, di, dj, dk) \
  int ii, ij, ik; \
  \
  if (wavetree3d_sub_3dindices(t, i, &ii, &ij, &ik) < 0) { \
    return -1; \
  } \
  \
  ii = 2*ii + di;\
  ij = 2*ij + dj;\
  ik = 2*ik + dk;\
  \
  if (ii >= t->width || ij >= t->height || ik >= t->depth) {\
    return -1;\
  }\
  \
  return wavetree3d_sub_from_3dindices(t, ii, ij, ik)

int wavetree3d_sub_UTL(const wavetree3d_sub_t *t, int i)
{
  if (i == 0) {
    return -1;
  } else {
    WAVETREE3D_SUB_CHILD(t, i, 0, 0, 0);
  }
}

int wavetree3d_sub_UTR(const wavetree3d_sub_t *t, int i)
{
  WAVETREE3D_SUB_CHILD(t, i, 1, 0, 0);
}

int wavetree3d_sub_UBL(const wavetree3d_sub_t *t, int i)
{
  WAVETREE3D_SUB_CHILD(t, i, 0, 1, 0);
}

int wavetree3d_sub_UBR(const wavetree3d_sub_t *t, int i)
{
  WAVETREE3D_SUB_CHILD(t, i, 1, 1, 0);
}

int wavetree3d_sub_LTL(const wavetree3d_sub_t *t, int i)
{
  WAVETREE3D_SUB_CHILD(t, i, 0, 0, 1);
}

int wavetree3d_sub_LTR(const wavetree3d_sub_t *t, int i)
{
  WAVETREE3D_SUB_CHILD(t, i, 1, 0, 1);
}

int wavetree3d_sub_LBL(const wavetree3d_sub_t *t, int i)
{
  WAVETREE3D_SUB_CHILD(t, i, 0, 1, 1);
}

int wavetree3d_sub_LBR(const wavetree3d_sub_t *t, int i)
{
  WAVETREE3D_SUB_CHILD(t, i, 1, 1, 1);
}

double wavetree3d_sub_logpriorprobability(const wavetree3d_sub_t *t,
					  wavetree_pp_t *pp)
{
  double logprior = 0.0;
  int ii;
  int jj;
  int kk;
  int d;
  int c;
  int i;
  int index;
  double value;
  
  for (d = 0; d <= t->degree_max; d ++) {
      
    c = multiset_int_double_depth_count(t->S_v, d);
    if (c > 0) {
      
      for (i = 0; i < c; i ++) {
	
	if (multiset_int_double_nth_element(t->S_v, d, i, &index, &value) < 0) {
	  ERROR("failed to get nth element");
	  return -1;
	}

	if (wavetree3d_sub_3dindices(t, index, &ii, &jj, &kk) < 0) {
	  ERROR("failed to get 2d indices");
	  return -1;
	}
	
	logprior += log(wavetree_pp_prior_probability3d(pp,
							ii,
							jj,
							kk,
							d,
							t->degree_max,
							0.0,
							value));
      }
    }
  }

  return logprior;
}
