
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "wavetree2d_sub.h"

#include "multiset_int.h"
#include "multiset_int_double.h"

#include "slog.h"

typedef enum {
  UNDO_NONE = WT_PERTURB_NONE,
  UNDO_VALUE = WT_PERTURB_VALUE,
  UNDO_BIRTH = WT_PERTURB_BIRTH,
  UNDO_DEATH = WT_PERTURB_DEATH,
  UNDO_MOVE = WT_PERTURB_MOVE
} wavetree2d_sub_undo_t;

struct _wavetree2d_sub {

  int base_width;
  int base_height;
  int base_size;
  int *base_indices;

  int max_children;
  int *child_indices;
  
  int degree_max;
  int degree_min;
  
  int degree_width;
  int degree_height;

  int width;
  int height;
  int size;

  double alpha;

  multiset_int_double_t *S_v;
  multiset_int_t *S_b;
  multiset_int_t *S_d;

  wavetree2d_sub_undo_t undo;
  int u_i;
  int u_j;
  int u_d;
  double u_v;

  chain_history_change_t last_step;
};

static int add_node(wavetree2d_sub_t *t, int i, int d, double coeff);
static int remove_node(wavetree2d_sub_t *t, int i, int d);

wavetree2d_sub_t *wavetree2d_sub_create(int degree_width, int degree_height, double alpha)
{
  wavetree2d_sub_t *r;
  int i;
  int j;
  int l;
  
  r = malloc(sizeof(wavetree2d_sub_t));
  if (r == NULL) {
    return NULL;
  }

  r->degree_min = degree_width;
  r->degree_max = degree_width;
  if (degree_height < r->degree_min) {
    r->degree_min = degree_height;
  }
  if (degree_height > r->degree_max) {
    r->degree_max = degree_height;
  }

  r->degree_width = degree_width;
  r->degree_height = degree_height;
  r->width = 1 << degree_width;
  r->height = 1 << degree_height;
  r->size = r->width * r->height;

  r->base_width = 1 << (degree_width - r->degree_min);
  r->base_height = 1 << (degree_height - r->degree_min);
  r->base_size = r->base_width * r->base_height;
  if (r->base_size < 1) {
    ERROR("failed to compute base size");
    return NULL;
  }

  if (r->base_size == 1) {
    r->base_indices = malloc(sizeof(int) * 1);
    if (r->base_indices == NULL) {
      ERROR("failed to allocate base indices");
      return NULL;
    }
    r->base_indices[0] = 0;
    
  } else {
    r->base_indices = malloc(sizeof(int) * r->base_size);
    if (r->base_indices == NULL) {
      ERROR("failed to allocate base indices");
      return NULL;
    }

    l = 0;
    for (j = 0; j < r->base_height; j ++) {
      for (i = 0; i < r->base_width; i ++, l ++) {
	r->base_indices[l] = j*r->width + i + 1;
      }
    }
  }

  r->alpha = alpha;

  if (r->base_size == 1) {

    r->degree_max = r->degree_width;
    r->max_children = 4;

  } else {

    r->degree_max = r->degree_min + 1;

    r->max_children = 4;
    if (r->base_size > r->max_children) {
      r->max_children = r->base_size;
    }
    
  }

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

void wavetree2d_sub_destroy(wavetree2d_sub_t *t)
{
  if (t != NULL) {
    multiset_int_double_destroy(t->S_v);
    multiset_int_destroy(t->S_d);
    multiset_int_destroy(t->S_b);

    free(t->child_indices);
    free(t->base_indices);
      
    free(t);
  }
}

int
wavetree2d_sub_save(const wavetree2d_sub_t *t,
		const char *filename)
{
  FILE *fp;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    ERROR("failed to save wavetree");
    return -1;
  }

  fprintf(fp, "%d %d\n", t->degree_width, t->degree_height);
  fprintf(fp, "%d %d %d\n", t->width, t->height, t->size);
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
wavetree2d_sub_load(wavetree2d_sub_t *t,
		    const char *filename)
{
  FILE *fp;

  int degree_width;
  int degree_height;
  int width;
  int height;
  int size;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    ERROR("failed to open file");
    return -1;
  }

  if (fscanf(fp, "%d %d\n", &degree_width, &degree_height) != 2) {
    ERROR("failed to read degree");
    return -1;
  }

  if (degree_width != t->degree_width) {
    ERROR("x-degree mismatch");
    return -1;
  }

  if (degree_height != t->degree_height) {
    ERROR("y-degree mismatch");
    return -1;
  }

  if (fscanf(fp, "%d %d %d\n", &width, &height, &size) != 3) {
    ERROR("failed to read header");
    return -1;
  }

  if (width != t->width ||
      height != t->height ||
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
wavetree2d_sub_encode(wavetree2d_sub_t *t,
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
  
  ENCODEINT(t->degree_width);
  ENCODEINT(t->degree_height);
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
wavetree2d_sub_decode(wavetree2d_sub_t *t,
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

  /* First element should be root of tree */
  DECODEINT(&d);
  DECODEINT(&index);
  DECODEDOUBLE(&v);
  if (d != 0 || index != 0) {
    ERROR("first element is not root (%d, %d, %f)", index, d, v);
    return -1;
  }

  if (wavetree2d_sub_initialize(t, v) < 0) {
    ERROR("failed to initialise");
    return -1;
  }
  
  
  /*
   * Values
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
wavetree2d_sub_load_promote(wavetree2d_sub_t *t,
			    const char *filename)
{
  FILE *fp;
  int i;
  int j;
  int coeffs;
  int idx;

  double coeff;

  int degree_width;
  int degree_height;
  int width;
  int height;
  int size;

  int ii;
  int ij;
  int pidx;

  int d;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    ERROR("failed to open file");
    return -1;
  }

  if (fscanf(fp, "%d %d\n", &degree_width, &degree_height) != 2) {
    ERROR("failed to read degree");
    return -1;
  }

  /* Promotion works for smaller or equal degree */
  if (degree_width > t->degree_width) {
    ERROR("degree mismatch %d > %d", 
	  degree_width, t->degree_width);
    return -1;
  }

  if (fscanf(fp, "%d %d %d\n", &width, &height, &size) != 3) {
    ERROR("failed to read header");
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

	if (wavetree2d_sub_initialize(t, coeff) < 0) {
	  ERROR("failed to initialize from 0 node");
	  return -1;
	}
      } else {
	/*
	 * We need to convert the indexing between degrees. We do this by
	 * getting the 2D indices using the input file degree, and then
	 * transforming back using the current degree
	 */
	
	ii = idx % width;
	ij = (idx - ii) / width;

	pidx = ij * t->width + ii;

	if (add_node(t, pidx, i, coeff) < 0) {
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

int wavetree2d_sub_initialize(wavetree2d_sub_t *t,
			      double dc)
{
  int d;
  int i;

  multiset_int_double_clear(t->S_v);
  multiset_int_clear(t->S_b);
  multiset_int_clear(t->S_d);

  d = wavetree2d_sub_depthofindex(t, 0);
  multiset_int_double_insert(t->S_v,
			     0,
			     d,
			     dc);
  
  if (t->base_size == 1) {
    multiset_int_insert(t->S_b, 
			1, 
			d + 1);
    multiset_int_insert(t->S_b, 
			t->width, 
			d + 1);
    multiset_int_insert(t->S_b, 
			t->width + 1, 
			d + 1);
    
    if (multiset_int_total_count(t->S_b) != 3) {
      ERROR("failed to initialize S_b");
      return -1;
    }
  } else {
    for (i = 0; i < t->base_size; i ++) {
      multiset_int_insert(t->S_b,
			  t->base_indices[i],
			  d + 1);
    }
    if (multiset_int_total_count(t->S_b) != t->base_size) {
      ERROR("failed to initialize S_b");
      return -1;
    }
  }    

  if (multiset_int_double_total_count(t->S_v) != 1) {
    ERROR("failed to initialize S_v");
    return -1;
  }


  return 0;
}

double wavetree2d_sub_dc(wavetree2d_sub_t *t)
{
  double dc;

  if (multiset_int_double_get(t->S_v, 0, 0, &dc) < 0) {
    return -1.0;
  }

  return dc;
}

int
wavetree2d_sub_get_width(wavetree2d_sub_t *t)
{
  return t->width;
}

int
wavetree2d_sub_get_height(wavetree2d_sub_t *t)
{
  return t->height;
}

int
wavetree2d_sub_get_size(wavetree2d_sub_t *t)
{
  return t->size;
}

int
wavetree2d_sub_base_size(wavetree2d_sub_t *t)
{
  return t->base_size;
}

const int *
wavetree2d_sub_base_indices(wavetree2d_sub_t *t)
{
  return t->base_indices;
}

int
wavetree2d_sub_get_ncoeff(wavetree2d_sub_t *t)
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

int wavetree2d_sub_prunable_leaves(const wavetree2d_sub_t *t)
{
  return multiset_int_total_count(t->S_d);
}

int wavetree2d_sub_attachable_branches(const wavetree2d_sub_t *t)
{
  return multiset_int_total_count(t->S_b);
}

int wavetree2d_sub_coeff_count(const wavetree2d_sub_t *t)
{
  return multiset_int_double_total_count(t->S_v);
}

int wavetree2d_sub_map_to_array(const wavetree2d_sub_t *t, double *a, int n)
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

    if (wavetree2d_sub_child_indices(t, 0, 0, t->child_indices, &c, t->max_children) < 0) {
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

int wavetree2d_sub_map_impulse_to_array(const wavetree2d_sub_t *t,
					int coeff_index,
					double *a,
					int n)
{
  int c;
  int i;
  int ii;
  int ij;

  if (t->base_size == 1) {

    a[coeff_index] = 1.0;

  } else {

    if (coeff_index == 0) {

      if (wavetree2d_sub_child_indices(t, 0, 0, t->child_indices, &c, t->max_children) < 0) {
	ERROR("failed to get 1st level children");
	return -1;
      }

      for (i = 0; i < c; i ++) {
	a[t->child_indices[i] - 1] = 1.0;
      }

    } else {

      if (wavetree2d_sub_2dindices(t,
				   coeff_index,
				   &ii,
				   &ij) < 0) {
	ERROR("Failed to get 2d indices\n");
	return -1;
      }

      a[ij * t->width + ii] = 1.0;

    }
  }

  return 0;
}

int wavetree2d_sub_map_from_array(wavetree2d_sub_t *t, const double *a, int n)
{
  return wavetree2d_sub_create_from_array_with_threshold(t,
							 a,
							 n,
							 0.0);
}

int r_wavetree2d_sub_create_from_array_children(wavetree2d_sub_t *t,
						const double *a,
						int n,
						int index,
						int depth)
{
  int i;
  int offset;
  int child_indices[4];
  int c;

  offset = 0;
  if (t->base_size > 1) {
    offset = -1;
  }

  if (wavetree2d_sub_child_indices(t, index, depth, child_indices, &c, 4) < 0) {
    return -1;
  }

  for (i = 0; i < c; i ++) {

    if (add_node(t, child_indices[i], depth + 1, a[child_indices[i] + offset]) < 0) {
      return -1;
    }

    if (r_wavetree2d_sub_create_from_array_children(t, a, n, child_indices[i], depth + 1) < 0) {
      return -1;
    }
  }

  return 0;
}

int wavetree2d_sub_create_from_array_with_threshold(wavetree2d_sub_t *t,
						    const double *a,
						    int n,
						    double threshold)
{
  int i;
  int c;
  double mean;
  int *child_indices;
  int d;

  int index;
  double value;

  if (t->base_size == 1) {

    if (wavetree2d_sub_initialize(t, a[0]) < 0) {
      ERROR("Failed to reset wavetree");
      return -1;
    }

    if (r_wavetree2d_sub_create_from_array_children(t, a, n, 0, 0) < 0) {
      ERROR("failed to set children");
      return -1;
    }

  } else {

    child_indices = malloc(sizeof(int) * t->max_children);
    if (child_indices == NULL) {
      ERROR("Failed to allocate temporary array\n");
      return -1;
    }
    
    /*
     * Compute mean of first generation children and set root
     */
    if (wavetree2d_sub_child_indices(t, 0, 0, child_indices, &c, t->max_children) < 0) {
      ERROR("failed to get 1st level children");
      return -1;
    }

    mean = 0.0;
    for (i = 0; i < c; i ++) {
      mean = mean + a[child_indices[i] - 1];
    }
    mean = mean/(double)c;

    if (wavetree2d_sub_initialize(t, mean) < 0) {
      ERROR("Failed to reset wavetree");
      return -1;
    }

    /*
     * Set values for first generation children
     */
    for (i = 0; i < c; i ++) {

      if (add_node(t, child_indices[i], 1, a[child_indices[i] - 1] - mean) < 0) {
	ERROR("Failed to set %d 1st gen child: %d", i, child_indices[i]);
	return -1;
      }
      
    }

    /*
     * Do remaining children
     */
    for (i = 0; i < c; i ++) {
      if (r_wavetree2d_sub_create_from_array_children(t, a, n, child_indices[i], 1) < 0) {
	return -1;
      }
    }

    free(child_indices);

  }

  /*
   * Now prune low magnitude coefficients with no major child coefficients
   */
  if (threshold > 0.0) {
    for (d = t->degree_max; d > 1; d --) {

      c = multiset_int_double_depth_count(t->S_v, d);
      if (c > 0) {
	
	for (i = 0; i < c; i ++) {
	  
	  if (multiset_int_double_nth_element(t->S_v, d, i, &index, &value) < 0) {
	    ERROR("failed to get nth element");
	    return -1;
	  }

	  if (fabs(value) < threshold &&
	      wavetree2d_sub_child_count(t, index, d) == 0) {
	    if (remove_node(t, index, d) < 0) {
	      return -1;
	    }
	    i --;
	    c --;
	  }
	}
      }
    }
  }
  
  return 0;
}


int
wavetree2d_sub_propose_value(wavetree2d_sub_t *t,
int i,
int d,
double value)
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
  if (wavetree2d_sub_get_coeff(t, i, &(t->u_v)) < 0) {
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
wavetree2d_sub_propose_birth(wavetree2d_sub_t *t,
			     int i,
			     int d,
			     double value)
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
  
  if (multiset_int_double_is_element(t->S_v, i, d)) {
    ERROR("non empty node %d %d", i, d);
    return -1;
  }

  if (wavetree2d_sub_depthofindex(t, i) != d) {
    ERROR("depth mismatch, got %d, expected %d for index %d",
	  d, wavetree2d_sub_depthofindex(t, i), i);
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
wavetree2d_sub_propose_death(wavetree2d_sub_t *t,
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

  if (multiset_int_double_is_element(t->S_v, i, d) == 0) {
    ERROR("empty node");
    return -1;
  }

  if (wavetree2d_sub_get_coeff(t, i, &(t->u_v)) < 0) {
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
wavetree2d_sub_propose_move(wavetree2d_sub_t *t,
			    int i,
			    int d,
			    int new_i,
			    double new_value)
{
  if (t == NULL || 
      i < 0 ||
      i >= t->size) {
    ERROR("invalid parameters");
    return -1;
  }

  if (multiset_int_double_is_element(t->S_v, i, d) == 0) {
    ERROR("empty node");
    return -1;
  }
  
  if (multiset_int_double_is_element(t->S_v, new_i, d)) {
    ERROR("non-empty destination");
    return -1;
  }

  t->undo = UNDO_MOVE;
  t->u_i = i;
  t->u_j = new_i;
  t->u_d = d;
  if (wavetree2d_sub_get_coeff(t, i, &(t->u_v)) < 0) {
    ERROR("failed to get coefficient");
    return -1;
  }

  if (remove_node(t, i, d) < 0) {
    ERROR("failed to remove old node");
    return -1;
  }

  
  if (add_node(t, new_i, d, new_value) < 0) {
    ERROR("failed to add new node");
    return -1;
  }

  /*
   * Record the proposal
   */
  t->last_step.header.type = CH_MOVE;
  t->last_step.header.accepted = 0;
  t->last_step.header.likelihood = 0.0;
  t->last_step.header.temperature = 0.0;
  t->last_step.header.hierarchical = 0.0;
  
  t->last_step.perturbation.move.node_depth = t->u_d;
  t->last_step.perturbation.move.node_id = t->u_i;
  t->last_step.perturbation.move.new_node_id = t->u_j;
  t->last_step.perturbation.move.new_value = new_value;
  t->last_step.perturbation.move.old_value = t->u_v;

  return 0;
}

int 
wavetree2d_sub_undo(wavetree2d_sub_t *t)
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
wavetree2d_sub_commit(wavetree2d_sub_t *t)
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

int wavetree2d_sub_valid(wavetree2d_sub_t *t)
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

      pi = wavetree2d_sub_parent_index(t, index);
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

      pi = wavetree2d_sub_parent_index(t, index);
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
int wavetree2d_sub_parent_index(const wavetree2d_sub_t *t, int c)
{
  int i;
  int j;

  if (c == 0) {
    /* 0th has no parent */
    return -1;
  }

  if (t->base_size > 1) {

    if (wavetree2d_sub_2dindices(t, c, &i, &j) < 0) {
      ERROR("failed to get 2d indices");
      return -1;
    }

    /* Subtile indices have root as parent */
    if (i < t->base_width &&
	j < t->base_height) {
      /* Parent is 0 */
      return 0;
    }

    /* Direct descendents of subtile */
    if (i < (2*t->base_width) &&
	j < (2*t->base_height)) {
      i = i % t->base_width;
      j = j % t->base_height;

      return wavetree2d_sub_from_2dindices(t, i, j);
    }

    /* Rest follow the /2 pattern */
    i /= 2;
    j /= 2;

    return wavetree2d_sub_from_2dindices(t, i, j);

  } else {

    if (wavetree2d_sub_2dindices(t, c, &i, &j) < 0) {
      ERROR("failed to get 2d indices");
      return -1;
    }

    i /= 2;
    j /= 2;

    return wavetree2d_sub_from_2dindices(t, i, j);
  }
}

int wavetree2d_sub_child_count(const wavetree2d_sub_t *t, int index, int depth)
{
  int j;
  int nchildren;
  int cc;

  cc = 0;

  if (wavetree2d_sub_child_indices(t, index, depth, t->child_indices, &nchildren, t->max_children) < 0) {
    return -1;
  }

  for (j = 0; j < nchildren; j ++) {
    if (multiset_int_double_is_element(t->S_v, t->child_indices[j], depth + 1)) {
      cc ++;
    }
  }

  return cc;
}


int wavetree2d_sub_2dindices(const wavetree2d_sub_t *t,
			     int i,
			     int *ii,
			     int *ij)
{
  if (t->base_size == 1) {
    if (i >= t->size ||
	i < 0) {
      return -1;
    }
  } else {
    if (i == 0) {
      *ii = -1;
      *ij = -1;
      return 0;
    }

    if (i > t->size) {
      return -1;
    }

    i --;
  }
  
  *ii = i % t->width;
  *ij = (i - (*ii)) / t->width;
  
  return 0;
}

int wavetree2d_sub_from_2dindices(const wavetree2d_sub_t *t,
				  int ii,
				  int ij)
{
  if (t == NULL) {
    return -1;
  }

  if (t->base_size == 1) {
    if (ii < 0 ||
	ij < 0 ||
	ii >= t->width ||
	ij >= t->height) {
      return -1;
    }
    

    return ij*t->width + ii;

  } else {

    if (ii == (-1) &&
	ij == (-1)) {
      return 0;
    }

    if (ii < 0 ||
	ij < 0 ||
	ii >= t->width ||
	ij >= t->height) {
      return -1;
    }

    return ij*t->width + ii + 1;
  }    
}

int wavetree2d_sub_depthofindex(const wavetree2d_sub_t *t,
				int i)
{
  if (i == 0) {
    return 0;
  }

  return 1 + wavetree2d_sub_depthofindex(t, wavetree2d_sub_parent_index(t, i));
}

int wavetree2d_sub_max_child_count(wavetree2d_sub_t *t)
{
  return t->max_children;
}

int wavetree2d_sub_child_indices(const wavetree2d_sub_t *t, int index, int depth, int *indices, int *n, int nmax)
{
  int i;

  int ii;
  int jj;
  
  if (t->base_size == 1) {
    if (nmax < 4) {
      ERROR("nmax insufficient (4, %d)", nmax);
      return -1;
    }
    i = 0;
    indices[i] = wavetree2d_sub_TL(t, index);
    if (indices[i] > 0) {
      i ++;
    }
    indices[i] = wavetree2d_sub_TR(t, index);
    if (indices[i] > 0) {
      i ++;
    }
    indices[i] = wavetree2d_sub_BL(t, index);
    if (indices[i] > 0) {
      i ++;
    }
    indices[i] = wavetree2d_sub_BR(t, index);
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

      if (wavetree2d_sub_2dindices(t,
				   index,
				   &ii, &jj) < 0) {
	ERROR("failed to get 2d indices");
	return -1;
      }

      i = 0;
      indices[i] = wavetree2d_sub_from_2dindices(t,
						 ii + t->base_width,
						 jj);
      if (indices[i] > 0) {
	i ++;
      }
      
      indices[i] = wavetree2d_sub_from_2dindices(t,
						 ii,
						 jj + t->base_height);
      if (indices[i] > 0) {
	i ++;
      }
      
      indices[i] = wavetree2d_sub_from_2dindices(t,
						 ii + t->base_width,
						 jj + t->base_height);
      if (indices[i] > 0) {
	i ++;
      }

      *n = i;
      return 0;

    } else {
      
      if (nmax < 4) {
	ERROR("nmax insufficient (4, %d)", nmax);
	return -1;
      }
      
      i = 0;
      indices[i] = wavetree2d_sub_TL(t, index);
      if (indices[i] > 0) {
	i ++;
      }
      indices[i] = wavetree2d_sub_TR(t, index);
      if (indices[i] > 0) {
	i ++;
      }
      indices[i] = wavetree2d_sub_BL(t, index);
      if (indices[i] > 0) {
	i ++;
      }
      indices[i] = wavetree2d_sub_BR(t, index);
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

int wavetree2d_sub_TL(const wavetree2d_sub_t *t, int i)
{
  int ii;
  int ij;

  if (i == 0) {
    return -1;
  }

  if (wavetree2d_sub_2dindices(t, i, &ii, &ij) < 0) {
    return -1;
  }

  ii = 2*ii;
  ij = 2*ij;

  if (ii >= t->width || ij >= t->height) {
    return -1;
  }

  return wavetree2d_sub_from_2dindices(t, ii, ij);
}

int wavetree2d_sub_TR(const wavetree2d_sub_t *t, int i)
{
  int ii;
  int ij;

  if (wavetree2d_sub_2dindices(t, i, &ii, &ij) < 0) {
    return -1;
  }

  ii = 2*ii + 1;
  ij = 2*ij;

  if (ii >= t->width || ij >= t->height) {
    return -1;
  }

  return wavetree2d_sub_from_2dindices(t, ii, ij);
}

int wavetree2d_sub_BL(const wavetree2d_sub_t *t, int i)
{
  int ii;
  int ij;

  if (wavetree2d_sub_2dindices(t, i, &ii, &ij) < 0) {
    return -1;
  }

  ii = 2*ii;
  ij = 2*ij + 1;

  if (ii >= t->width || ij >= t->height) {
    return -1;
  }

  return wavetree2d_sub_from_2dindices(t, ii, ij);
}

int wavetree2d_sub_BR(const wavetree2d_sub_t *t, int i)
{
  int ii;
  int ij;

  if (wavetree2d_sub_2dindices(t, i, &ii, &ij) < 0) {
    return -1;
  }

  ii = 2*ii + 1;
  ij = 2*ij + 1;

  if (ii >= t->width || ij >= t->height) {
    return -1;
  }

  return wavetree2d_sub_from_2dindices(t, ii, ij);
}

static int add_node(wavetree2d_sub_t *t, int i, int d, double coeff)
{
  int j;
  int pcc;
  int nchildren;

  if (multiset_int_double_insert(t->S_v, i, d, coeff) < 0) {
    ERROR("error: failed to insert into S_v");
    return -1;
  }

  /*
   * Update parent CC
   */
  j = wavetree2d_sub_parent_index(t, i);

  pcc = wavetree2d_sub_child_count(t, j, d - 1);
  if (pcc < 0) {
    ERROR("failed to get pcc (%d, %d: %d)", j, d - 1, i);
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

  if (wavetree2d_sub_child_indices(t, i, d, t->child_indices, &nchildren, t->max_children) < 0) {
    return -1;
  }

  for (j = 0; j < nchildren; j ++) {
    multiset_int_insert(t->S_b, t->child_indices[j], d + 1);
  }
  
  /*printf("Add: %d %d\n", oset_int_count(t->S_b), oset_int_count(t->S_d)); */

  return 0;
}

static int remove_node(wavetree2d_sub_t *t, int i, int d)
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
  j = wavetree2d_sub_parent_index(t, i);

  pcc = wavetree2d_sub_child_count(t, j, d - 1);

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
  if (wavetree2d_sub_child_indices(t, i, d, t->child_indices, &nchildren, t->max_children) < 0) {
    return -1;
  }

  for (j = 0; j < nchildren; j ++) {
    multiset_int_remove(t->S_b, t->child_indices[j], d + 1);
  }
  
  /*  printf("Remove: %d %d\n", t->N_b, t->N_d); */

  return 0;
}

void wavetree2d_sub_print_setinfo(wavetree2d_sub_t *t)
{
  printf("Nb %d Nd %d CC %d\n", 
	 multiset_int_total_count(t->S_b), 
	 multiset_int_total_count(t->S_d), 
	 wavetree2d_sub_coeff_count(t));
}

void wavetree2d_sub_dump_sets(wavetree2d_sub_t *t)
{
  printf("Sb:\n  ");
  multiset_int_dump(t->S_b);
  printf("Sd:\n  ");
  multiset_int_dump(t->S_d);
  printf("Sv:\n  ");
  multiset_int_double_dump(t->S_v);
}

void wavetree2d_sub_dump_coeffs(wavetree2d_sub_t *t)
{
  multiset_int_double_dump(t->S_v);
}

static int r_generate_dyck_word(wavetree2d_sub_t *t, int index, char *buffer, int *len, int maxlength)
{
  int j;

  if (multiset_int_double_is_element(t->S_v,
				     index,
				     wavetree2d_sub_depthofindex(t, index)) == 0) {

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

  
  j = wavetree2d_sub_TL(t, index);
  if (j > 0) {
    if (r_generate_dyck_word(t, j, buffer, len, maxlength) < 0) {
      return -1;
    }
  }

  j = wavetree2d_sub_TR(t, index);
  if (j > 0) {
    if (r_generate_dyck_word(t, j, buffer, len, maxlength) < 0) {
      return -1;
    }
  }

  j = wavetree2d_sub_BL(t, index);
  if (j > 0) {
    if (r_generate_dyck_word(t, j, buffer, len, maxlength) < 0) {
      return -1;
    }
  }

  j = wavetree2d_sub_BR(t, index);
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

int wavetree2d_sub_generate_dyck_word(wavetree2d_sub_t *t, char *buffer, int maxlength)
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

static int r_generate_dyck_binary(wavetree2d_sub_t *t, int index, uint64_t *binary, int *len)
{
  int j;

  if (multiset_int_double_is_element(t->S_v,
				     index,
				     wavetree2d_sub_depthofindex(t, index)) == 0) {
    return db_leaf(binary, len);
  }

  db_open(binary, len);
  
  j = wavetree2d_sub_TL(t, index);
  if (j > 0) {
    if (r_generate_dyck_binary(t, j, binary, len) < 0) {
      return -1;
    }
  } else {
    if (index != 0) {
      db_leaf(binary, len);
    }
  }

  j = wavetree2d_sub_TR(t, index);
  if (j > 0) {
    if (r_generate_dyck_binary(t, j, binary, len) < 0) {
      return -1;
    }
  } else {
    db_leaf(binary, len);
  }

  j = wavetree2d_sub_BL(t, index);
  if (j > 0) {
    if (r_generate_dyck_binary(t, j, binary, len) < 0) {
      return -1;
    }
  } else {
    db_leaf(binary, len);
  }


  j = wavetree2d_sub_BR(t, index);
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

int wavetree2d_sub_generate_dyck_binary(wavetree2d_sub_t *t, uint64_t *binary)
{
  int len;

  *binary = 0;

  len = 0;
  if (r_generate_dyck_binary(t, 0, binary, &len) < 0) {
    return -1;
  }

  return 0;
}

int wavetree2d_sub_get_indices(wavetree2d_sub_t *t, int *set, int *n)
{
  int d;
  int dn;
  int i;
  int j;
  int madegree_width;
  double v;

  j = 0;

  madegree_width = t->degree_width;
  if (t->degree_height > madegree_width) {
    madegree_width = t->degree_height;
  }
  
  for (d = 0; d <= madegree_width; d ++) {

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

int wavetree2d_sub_maxdepth(const wavetree2d_sub_t *t)
{
  return t->degree_max;
}

int wavetree2d_sub_depth(wavetree2d_sub_t *t)
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

int wavetree2d_sub_update_histogram(const wavetree2d_sub_t *t, coefficient_histogram_t *hist)
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

int wavetree2d_sub_depth_filter(void *wt, int subset, int index)
{
  wavetree2d_sub_t *t = (wavetree2d_sub_t*)wt;

  return wavetree2d_sub_depthofindex(t, index) == subset;
}

const multiset_int_double_t *
wavetree2d_sub_get_S_v(const wavetree2d_sub_t *t)
{
  return t->S_v;
}

int
wavetree2d_sub_set_from_S_v(wavetree2d_sub_t *t,
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

  if (wavetree2d_sub_initialize(t, value) < 0) {
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
wavetree2d_sub_get_model(const wavetree2d_sub_t *t,
			 int nmax,
			 int *indices,
			 double *values,
			 int *_n)
{
  int d;
  int n;
  int i;
  int index;
  double value;
  
  int k;

  k = 0;
  
  for (d = 0; d <= t->degree_max; d ++) {

    n = multiset_int_double_depth_count(t->S_v, d);
    for (i = 0; i < n; i ++) {

      if (multiset_int_double_nth_element(t->S_v, d, i, &index, &value) < 0) {
	ERROR("failed to get nth element");
	return -1;
      }

      if (k >= nmax) {
	ERROR("too many coefficients.");
	return -1;
      }
      
      indices[k] = index;
      values[k] = value;

      k ++;
      
    }
  }

  *_n = k;
  return 0;
}

int
wavetree2d_sub_set_model(wavetree2d_sub_t *t,
			 int *indices,
			 double *values,
			 int n)
{
  int i;
  int d;

  if (n <= 0) {
    ERROR("invalid number of coefficients (%d)", n);
    return -1;
  }
  
  multiset_int_double_clear(t->S_v);
  multiset_int_clear(t->S_b);
  multiset_int_clear(t->S_d);

  /*
   * First entry should be root element
   */
  if (indices[0] != 0) {
    ERROR("first element should be root (%d)", indices[0]);
    return -1;
  }

  if (wavetree2d_sub_initialize(t, values[0]) < 0) {
    ERROR("failed to initialize");
    return -1;
  }
  
  /*
   * Values
   */  
  for (i = 1; i < n; i ++) {

    d = wavetree2d_sub_depthofindex(t, indices[i]);
    if (d < 0) {
      ERROR("failed to get depth of index %d", indices[i]);
      return -1;
    }

    if (add_node(t, indices[i], d, values[i]) < 0) {
      ERROR("failed to add node (%d %d %f)", indices[i], d, values[i]);
      return -1;
    }
  }

  return 0;
  
}



int
wavetree2d_sub_get_last_perturbation(wavetree2d_sub_t *t,
				 chain_history_change_t *step)
{
  memcpy(step, &t->last_step, sizeof(chain_history_change_t));

  return 0;
}

int
wavetree2d_sub_set_invalid_perturbation(wavetree2d_sub_t *t, wavetree_perturb_t p)
{
  memset(&(t->last_step), 0, sizeof(chain_history_change_t));
  t->last_step.header.type = p;

  return 0;
}

int wavetree2d_sub_perturb(wavetree2d_sub_t *t,
			   wavetree2d_sub_perturb_func_t f,
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

  int pi;
  double pvalue;

  for (d = 0; d < t->degree_max; d ++) {
    n = multiset_int_double_depth_count(t->S_v, d);
    for (i = 0; i < n; i ++) {

      if (multiset_int_double_nth_element(t->S_v, d, i, &index, &value) < 0) {
        ERROR("failed to get nth element");
        return -1;
      }

      if (wavetree2d_sub_2dindices(t, i, &ii, &ij) < 0) {
        ERROR("failed to get 2d indices");
        return -1;
      }

      /* Get parent coeff */
      if (index == 0) {
        pvalue = 0.0;
      } else {
        pi = wavetree2d_sub_parent_index(t, index);
        if (wavetree2d_sub_get_coeff(t, pi, &pvalue) < 0) {
          ERROR("failed to get parent coeff");
          return -1;
        }
      }
       
      old_value = value;
      if (f(user, ii, ij, d, t->degree_max, pvalue, &value, prior_ratio) < 0) {
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

int wavetree2d_sub_get_coeff(const wavetree2d_sub_t *t,
			     int i,
			     double *coeff)
{
  int d;

  d = wavetree2d_sub_depthofindex(t, i);
  
  if (multiset_int_double_get(t->S_v, i, d, coeff) < 0) {
    return -1;
  }

  return 0;
}

/*
 * Functions for setting up a birth proposal
 */

int wavetree2d_sub_choose_birth_depth(const wavetree2d_sub_t *t, double u, int maxdepth, int *depth, double *prob)
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

int wavetree2d_sub_reverse_birth_depth(const wavetree2d_sub_t *t, int depth, int maxdepth, double *prob)
{
  if (t == NULL) {
    return -1;
  }

  *prob = 1.0/((double)multiset_int_nonempty_count(t->S_d, maxdepth));
  return 0;
}

int wavetree2d_sub_choose_birth(const wavetree2d_sub_t *t, int depth, double u, int *coeff, double *prob)
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

int wavetree2d_sub_choose_birth_global(const wavetree2d_sub_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob)
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



int wavetree2d_sub_reverse_birth(const wavetree2d_sub_t *t, int depth, int coeff, double *prob)
{
  if (t == NULL) {
    return -1;
  }
  
  *prob = 1.0/(double)(multiset_int_depth_count(t->S_d, depth));
  return 0;
}

int wavetree2d_sub_reverse_birth_global(const wavetree2d_sub_t *t, int maxdepth, int depth, int coeff, double *prob)
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
int wavetree2d_sub_choose_death_depth(const wavetree2d_sub_t *t, double u, int maxdepth, int *depth, double *prob)
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

int wavetree2d_sub_reverse_death_depth(const wavetree2d_sub_t *t, int depth, int maxdepth, double *prob)
{
  if (t == NULL) {
    return -1;
  }

  *prob = 1.0/((double)multiset_int_nonempty_count(t->S_b, maxdepth));
  return 0;
}

int wavetree2d_sub_choose_death(const wavetree2d_sub_t *t, int depth, double u, int *coeff, double *prob)
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

int wavetree2d_sub_choose_death_global(const wavetree2d_sub_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob)
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

int wavetree2d_sub_reverse_death(const wavetree2d_sub_t *t, int depth, int coeff, double *prob)
{
  if (t == NULL) {
    return -1;
  }
  
  *prob = 1.0/(double)(multiset_int_depth_count(t->S_b, depth));
  return 0;
}

int wavetree2d_sub_reverse_death_global(const wavetree2d_sub_t *t, int maxdepth, int depth, int coeff, double *prob)
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
int wavetree2d_sub_choose_value_depth(const wavetree2d_sub_t *t, double u, int maxdepth, int *depth, double *prob)
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

int wavetree2d_sub_choose_value(const wavetree2d_sub_t *t, int depth, double u, int *coeff, double *prob)
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

int wavetree2d_sub_choose_value_global(const wavetree2d_sub_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob)
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

/*
 * Move proposal code
 */
int wavetree2d_sub_choose_move_depth(const wavetree2d_sub_t *t, double u, int maxdepth, int *depth, double *prob)
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

int wavetree2d_sub_choose_move(const wavetree2d_sub_t *t, int depth, double u, int *coeff, double *prob)
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

int wavetree2d_sub_choose_move_global(const wavetree2d_sub_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob)
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

static int check_move_sibling(const wavetree2d_sub_t *t, int depth, int coeff, int di, int dj, int *sibling)
{
  int scoeff;
  int pi;

  int ii;
  int jj;

  if (wavetree2d_sub_2dindices(t, coeff, &ii, &jj) < 0) {
    ERROR("check_move_sibling: failed to get 2d indices");
    return 0;
  }

  ii += di;
  jj += dj;

  if (ii < 0 || ii >= t->width ||
      jj < 0 || jj >= t->width) {
    return 0;
  }

  scoeff = wavetree2d_sub_from_2dindices(t, ii, jj);
  if (scoeff < 0) {
    ERROR("check_move_sibling: failed to get sibling index");
    return 0;
  }

  /*
   * Must be same depth
   */
  if (wavetree2d_sub_depthofindex(t, scoeff) != depth) {
    return 0;
  }

  /*
   * Must not already be set
   */
  if (multiset_int_double_is_element(t->S_v, scoeff, depth)) {
    return 0;
  }

  /* 
   * Parent index must be currently set.
   */
  pi = wavetree2d_sub_parent_index(t, scoeff);
  if (multiset_int_double_is_element(t->S_v, pi, depth - 1) == 0) {
    return 0;
  }

  /*
   * scoeff is a valid sibling for a move proposal. 
   */
  *sibling = scoeff;
  return -1;
}

int
wavetree2d_sub_move_available_siblings(const wavetree2d_sub_t *t, int depth, int coeff, int *siblings, int *nsibling)
{
  int c;

  c = 0;

  if (check_move_sibling(t, depth, coeff, -1, -1, &siblings[c])) {
    c ++;
  }
  if (check_move_sibling(t, depth, coeff, 0, -1, &siblings[c])) {
    c ++;
  }
  if (check_move_sibling(t, depth, coeff, 1, -1, &siblings[c])) {
    c ++;
  }

  if (check_move_sibling(t, depth, coeff, -1, 0, &siblings[c])) {
    c ++;
  }
  if (check_move_sibling(t, depth, coeff, 1, 0, &siblings[c])) {
    c ++;
  }

  if (check_move_sibling(t, depth, coeff, -1, 1, &siblings[c])) {
    c ++;
  }
  if (check_move_sibling(t, depth, coeff, 0, 1, &siblings[c])) {
    c ++;
  }
  if (check_move_sibling(t, depth, coeff, 1, 1, &siblings[c])) {
    c ++;
  }

  *nsibling = c;

  return 0;
}

int wavetree2d_sub_choose_move_sibling(const wavetree2d_sub_t *t, double u, int depth, int coeff, int *sibling, double *prob)
{
  int valid_siblings[8];
  int nvalid;
  int si;

  if (wavetree2d_sub_move_available_siblings(t, depth, coeff, valid_siblings, &nvalid) < 0) {
    return -1;
  }

  if (nvalid == 0) {
    return -1;
  }

  si = (int)(u * (double)nvalid);
  *sibling = valid_siblings[si];
  *prob = 1.0/((double)nvalid);

  return 0;
}

int wavetree2d_sub_reverse_choose_move_sibling(const wavetree2d_sub_t *t, int depth, int coeff, int sibling, double *prob)
{
  int valid_siblings[8];
  int nvalid;
  int si;
  int i;

  if (wavetree2d_sub_move_available_siblings(t, depth, sibling, valid_siblings, &nvalid) < 0) {
    return -1;
  }

  if (nvalid == 0) {
    /*
     * This means that only the original coeff is the possible choice 
     */
    *prob = 1.0;
    return 0;
  }

  if (nvalid == 8) {
    /*
     * This shouldn't happen.
     */
    ERROR("too many valid siblings");
    return -1;
  }

  /*
   * Sanity check to make sure we're not double counting.
   */
  si = -1;
  for (i = 0; i < nvalid; i ++) {
    if (valid_siblings[i] == coeff) {
      si = i;
      break;
    }
  }

  if (si >= 0) {
    ERROR("found coeff in list of valid siblings");

    ERROR("  sibling: %d", sibling);
    ERROR("    coeff: %d", coeff);
    for (i = 0; i < nvalid; i ++) {
      ERROR("  %d %d", i, valid_siblings[i]);
    }
    return -1;
  }

  *prob = 1.0/((double)(nvalid + 1));

  return 0;
}


double wavetree2d_sub_logpriorprobability(const wavetree2d_sub_t *t,
					  wavetree_pp_t *pp)
{
  double logprior = 0.0;
  int ii;
  int jj;
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

	if (wavetree2d_sub_2dindices(t, index, &ii, &jj) < 0) {
	  ERROR("failed to get 2d indices");
	  return -1;
	}
	
	logprior += log(wavetree_pp_prior_probability2d(pp,
							ii,
							jj,
							d,
							t->degree_max,
							0.0,
							value));
      }
    }
  }

  return logprior;
}
