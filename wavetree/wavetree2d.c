
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "wavetree2d.h"

#include "multiset_int.h"

#include "slog.h"

typedef enum {
  UNDO_NONE = WT_PERTURB_NONE,
  UNDO_VALUE = WT_PERTURB_VALUE,
  UNDO_BIRTH = WT_PERTURB_BIRTH,
  UNDO_DEATH = WT_PERTURB_DEATH,
  UNDO_MOVE = WT_PERTURB_MOVE
} wavetree2d_undo_t;

struct _wavetree2d {

  int degree_max;
  int degree_width;
  int degree_height;

  int width;
  int height;
  int size;

  double alpha;

  multiset_int_double_t *S_v;
  multiset_int_t *S_b;
  multiset_int_t *S_d;

  wavetree2d_undo_t undo;
  int u_i;
  int u_j;
  int u_d;
  double u_v;

  chain_history_change_t last_step;
};

static int depthOfIndex(const wavetree2d_t *t,
			int i);

static int add_node(wavetree2d_t *t, int i, int d, double coeff);
static int remove_node(wavetree2d_t *t, int i, int d);

wavetree2d_t *wavetree2d_create(int degree_width, int degree_height, double alpha)
{
  wavetree2d_t *r;

  r = malloc(sizeof(wavetree2d_t));
  if (r == NULL) {
    return NULL;
  }

  r->degree_width = degree_width;
  r->degree_height = degree_height;
  r->width = 1 << degree_width;
  r->height = 1 << degree_height;
  r->size = r->width * r->height;

  r->degree_max = r->degree_width;
  if (r->degree_height > r->degree_max) {
    r->degree_max = r->degree_height;
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
  r->u_j = 0;
  r->u_d = 0;
  r->u_v = 0.0;

  memset(&(r->last_step), 0, sizeof(chain_history_change_t));
  r->last_step.header.type = CH_INITIALISE;
  
  return r;
}

void wavetree2d_destroy(wavetree2d_t *t)
{
  if (t != NULL) {
    multiset_int_double_destroy(t->S_v);
    multiset_int_destroy(t->S_d);
    multiset_int_destroy(t->S_b);
    free(t);
  }
}

int
wavetree2d_save(const wavetree2d_t *t,
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
wavetree2d_load(wavetree2d_t *t,
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
wavetree2d_encode(wavetree2d_t *t,
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
wavetree2d_decode(wavetree2d_t *t,
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

  /*
   * Values
   */  
  for (i = 0; i < n; i ++) {

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
wavetree2d_load_promote(wavetree2d_t *t,
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

	if (wavetree2d_initialize(t, coeff) < 0) {
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

int wavetree2d_initialize(wavetree2d_t *t,
			  double dc)
{
  multiset_int_double_insert(t->S_v,
			     0,
			     depthOfIndex(t, 0),
			     dc);
  
  multiset_int_insert(t->S_b, 
		      1, 
		      depthOfIndex(t, 1));
  multiset_int_insert(t->S_b, 
		      t->width, 
		      depthOfIndex(t, t->width));
  multiset_int_insert(t->S_b, 
		      t->width + 1, 
		      depthOfIndex(t, t->width + 1));

  if (multiset_int_double_total_count(t->S_v) != 1) {
    ERROR("failed to initialize S_v");
    return -1;
  }

  if (multiset_int_total_count(t->S_b) != 3) {
    ERROR("failed to initialize S_b");
    return -1;
  }

  return 0;
}

double wavetree2d_dc(wavetree2d_t *t)
{
  double dc;

  if (multiset_int_double_get(t->S_v, 0, 0, &dc) < 0) {
    return -1.0;
  }

  return dc;
}

int
wavetree2d_get_width(wavetree2d_t *t)
{
  return t->width;
}

int
wavetree2d_get_height(wavetree2d_t *t)
{
  return t->height;
}

int
wavetree2d_get_size(wavetree2d_t *t)
{
  return t->size;
}

int wavetree2d_prunable_leaves(const wavetree2d_t *t)
{
  return multiset_int_total_count(t->S_d);
}

int wavetree2d_attachable_branches(const wavetree2d_t *t)
{
  return multiset_int_total_count(t->S_b);
}

int wavetree2d_coeff_count(const wavetree2d_t *t)
{
  return multiset_int_double_total_count(t->S_v);
}

int wavetree2d_map_to_array(const wavetree2d_t *t, double *a, int n)
{
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

        if (index < 0 || index >= n) {
          ERROR("index out of range %d (%d)", index, n);
          return -1;
        }

        a[index] = value;
      }
    }
  }

  return 0;
}

int wavetree2d_map_from_array(wavetree2d_t *t, const double *a, int n)
{
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

        if (index < 0 || index >= n) {
          ERROR("index out of range %d (%d)", index, n);
          return -1;
        }

        if (multiset_int_double_set(t->S_v, index, d, a[index]) < 0) {
	  ERROR("failed to set value");
	  return -1;
	}
      }
    }
  }

  return 0;
}

int
wavetree2d_propose_value(wavetree2d_t *t,
			 int i,
			 int d,
			 double value)
{
  if (t == NULL || 
      i < 0 ||
      i >= t->size) {
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
  if (wavetree2d_get_coeff(t, i, &(t->u_v)) < 0) {
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
wavetree2d_propose_birth(wavetree2d_t *t,
			 int i,
			 int d,
			 double value)
{
  if (t == NULL || 
      i < 0 ||
      i >= t->size) {
    ERROR("invalid parameters");
    return -1;
  }

  if (multiset_int_double_is_element(t->S_v, i, d)) {
    ERROR("non empty node %d %d", i, d);
    return -1;
  }

  if (depthOfIndex(t, i) != d) {
    ERROR("depth mismatch, got %d, expected %d for index %d",
	  d, depthOfIndex(t, i), i);
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
wavetree2d_propose_death(wavetree2d_t *t,
			 int i,
			 int d, 
			 double *old_value)
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

  if (wavetree2d_get_coeff(t, i, &(t->u_v)) < 0) {
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
wavetree2d_propose_move(wavetree2d_t *t,
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
  if (wavetree2d_get_coeff(t, i, &(t->u_v)) < 0) {
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
wavetree2d_undo(wavetree2d_t *t)
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
wavetree2d_commit(wavetree2d_t *t)
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

int wavetree2d_valid(wavetree2d_t *t)
{
  return -1;
}

/* Parent index of a given node */
int wavetree2d_parent_index(const wavetree2d_t *t, int c)
{
  int i;
  int j;

  if (c == 0) {
    /* 0th has no parent */
    return -1;
  }

  i = c%t->width;
  j = (c - i)/t->width;

  if (j % 2 == 1) {
    j --;
  }
  
  if (i % 2 == 1) {
    i --;
  }

  j /= 2;
  i /= 2;

  return j*t->width + i;
}

int wavetree2d_child_count(const wavetree2d_t *t, int index, int depth)
{
  int j;
  int cc;

  cc = 0;

  j = wavetree2d_TL(t, index);
  if (j >= 0) {
    if (multiset_int_double_is_element(t->S_v, j, depth + 1)) {
      cc ++;
    }
  }

  j = wavetree2d_TR(t, index);
  if (j >= 0) {
    if (multiset_int_double_is_element(t->S_v, j, depth + 1)) {
      cc ++;
    }
  }

  j = wavetree2d_BL(t, index);
  if (j >= 0) {
    if (multiset_int_double_is_element(t->S_v, j, depth + 1)) {
      cc ++;
    }
  }

  j = wavetree2d_BR(t, index);
  if (j >= 0) {
    if (multiset_int_double_is_element(t->S_v, j, depth + 1)) {
      cc ++;
    }
  }

  return cc;
}


int wavetree2d_2dindices(const wavetree2d_t *t,
			 int i,
			 int *ii,
			 int *ij)
{
  if (i >= t->size ||
      i < 0) {
    return -1;
  }

  *ii = i % t->width;
  *ij = (i - (*ii)) / t->width;

  return 0;
}

int wavetree2d_from_2dindices(const wavetree2d_t *t,
			     int ii,
			     int ij)
{
  if (t == NULL ||
      ii < 0 ||
      ij < 0 ||
      ii >= t->width ||
      ij >= t->height) {
    return -1;
  }

  return ij*t->width + ii;
}

int wavetree2d_depthofindex(const wavetree2d_t *t,
			    int i)
{
  return depthOfIndex(t, i);
}

static int depthOfIndex(const wavetree2d_t *t,
			int i)
{
  int ii;
  int jj;
  int d;

  if (wavetree2d_2dindices(t, i, &ii, &jj) < 0) {
    return -1;
  }

  if (ii < jj) {
    ii = jj;
  }

  d = 0;
  while (ii > 0) {
    d ++;
    ii >>= 1;
  }

  return d;
}

int wavetree2d_TL(const wavetree2d_t *t, int i)
{
  int ii;
  int ij;
  int j;

  if (i == 0) {
    return -1;
  }

  if (wavetree2d_2dindices(t, i, &ii, &ij) < 0) {
    return -1;
  }

  ii = 2*ii;
  ij = 2*ij;

  if (ii >= t->width || ij >= t->height) {
    return -1;
  }

  j = wavetree2d_from_2dindices(t, ii, ij);
  if (j <= i) {
    ERROR("error in %d out %d", i, j);
    return -1;
  }
  return ij*t->width + ii;
}

int wavetree2d_TR(const wavetree2d_t *t, int i)
{
  int ii;
  int ij;
  int j;

  if (wavetree2d_2dindices(t, i, &ii, &ij) < 0) {
    return -1;
  }

  ii = 2*ii + 1;
  ij = 2*ij;

  if (ii >= t->width || ij >= t->height) {
    return -1;
  }

  j = wavetree2d_from_2dindices(t, ii, ij);
  if (j <= i) {
    ERROR("error in %d out %d", i, j);
    return -1;
  }
  return j;
}

int wavetree2d_BL(const wavetree2d_t *t, int i)
{
  int ii;
  int ij;
  int j;

  if (wavetree2d_2dindices(t, i, &ii, &ij) < 0) {
    return -1;
  }

  ii = 2*ii;
  ij = 2*ij + 1;

  if (ii >= t->width || ij >= t->height) {
    return -1;
  }

  j = wavetree2d_from_2dindices(t, ii, ij);
  if (j <= i) {
    ERROR("error in %d out %d", i, j);
    return -1;
  }
  return j;
}

int wavetree2d_BR(const wavetree2d_t *t, int i)
{
  int ii;
  int ij;
  int j;

  if (wavetree2d_2dindices(t, i, &ii, &ij) < 0) {
    return -1;
  }

  ii = 2*ii + 1;
  ij = 2*ij + 1;

  if (ii >= t->width || ij >= t->height) {
    return -1;
  }

  j = wavetree2d_from_2dindices(t, ii, ij);
  if (j <= i) {
    ERROR("error in %d out %d", i, j);
    return -1;
  }
  return j;
}

static int add_node(wavetree2d_t *t, int i, int d, double coeff)
{
  int j;
  int pcc;

  if (multiset_int_double_insert(t->S_v, i, d, coeff) < 0) {
    ERROR("failed to insert into S_v");
    return -1;
  }

  /*
   * Update parent CC
   */
  j = wavetree2d_parent_index(t, i);

  pcc = wavetree2d_child_count(t, j, d - 1);
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
  
  j = wavetree2d_TL(t,i);
  if (j > 0) {
    multiset_int_insert(t->S_b, j, d + 1);
  }

  j = wavetree2d_TR(t,i);
  if (j > 0) {
    multiset_int_insert(t->S_b, j, d + 1);
  }

  j = wavetree2d_BL(t,i);
  if (j > 0) {
    multiset_int_insert(t->S_b, j, d + 1);
  }

  j = wavetree2d_BR(t,i);
  if (j > 0) {
    multiset_int_insert(t->S_b, j, d + 1);
  }

  /*printf("Add: %d %d\n", oset_int_count(t->S_b), oset_int_count(t->S_d)); */

  return 0;
}

static int remove_node(wavetree2d_t *t, int i, int d)
{
  int j;
  int pcc;
  
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
  j = wavetree2d_parent_index(t, i);

  pcc = wavetree2d_child_count(t, j, d - 1);

  /*
   * Update sets Sb and Sd
   */
  if (j > 0 && pcc == 0) {
    /* Add parent to S_d as it now has no children */
    multiset_int_insert(t->S_d, j, d - 1);
  }
  

  multiset_int_remove(t->S_d, i, d);
  multiset_int_insert(t->S_b, i, d);

  j = wavetree2d_TL(t,i);
  if (j > 0) {
    multiset_int_remove(t->S_b, j, d + 1);
  }

  j = wavetree2d_TR(t,i);
  if (j > 0) {
    multiset_int_remove(t->S_b, j, d + 1);
  }

  j = wavetree2d_BL(t,i);
  if (j > 0) {
    multiset_int_remove(t->S_b, j, d + 1);
  }

  j = wavetree2d_BR(t,i);
  if (j > 0) {
    multiset_int_remove(t->S_b, j, d + 1);
  }
  
  /*  printf("Remove: %d %d\n", t->N_b, t->N_d); */

  return 0;
}

void wavetree2d_print_setinfo(wavetree2d_t *t)
{
  printf("Nb %d Nd %d CC %d\n", 
	 multiset_int_total_count(t->S_b), 
	 multiset_int_total_count(t->S_d), 
	 wavetree2d_coeff_count(t));
}

void wavetree2d_dump_sets(wavetree2d_t *t)
{
  printf("Sb:\n  ");
  multiset_int_dump(t->S_b);
  printf("Sd:\n  ");
  multiset_int_dump(t->S_d);
  printf("Sv:\n  ");
  multiset_int_double_dump(t->S_v);
}

void wavetree2d_dump_coeffs(wavetree2d_t *t)
{
  multiset_int_double_dump(t->S_v);
}

static int r_generate_dyck_word(wavetree2d_t *t, int index, char *buffer, int *len, int maxlength)
{
  int j;

  if (multiset_int_double_is_element(t->S_v,
				     index,
				     wavetree2d_depthofindex(t, index)) == 0) {

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

  
  j = wavetree2d_TL(t, index);
  if (j > 0) {
    if (r_generate_dyck_word(t, j, buffer, len, maxlength) < 0) {
      return -1;
    }
  }

  j = wavetree2d_TR(t, index);
  if (j > 0) {
    if (r_generate_dyck_word(t, j, buffer, len, maxlength) < 0) {
      return -1;
    }
  }

  j = wavetree2d_BL(t, index);
  if (j > 0) {
    if (r_generate_dyck_word(t, j, buffer, len, maxlength) < 0) {
      return -1;
    }
  }

  j = wavetree2d_BR(t, index);
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

int wavetree2d_generate_dyck_word(wavetree2d_t *t, char *buffer, int maxlength)
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

static int r_generate_dyck_binary(wavetree2d_t *t, int index, uint64_t *binary, int *len)
{
  int j;

  if (multiset_int_double_is_element(t->S_v,
				     index,
				     wavetree2d_depthofindex(t, index)) == 0) {
    return db_leaf(binary, len);
  }

  db_open(binary, len);
  
  j = wavetree2d_TL(t, index);
  if (j > 0) {
    if (r_generate_dyck_binary(t, j, binary, len) < 0) {
      return -1;
    }
  } else {
    if (index != 0) {
      db_leaf(binary, len);
    }
  }

  j = wavetree2d_TR(t, index);
  if (j > 0) {
    if (r_generate_dyck_binary(t, j, binary, len) < 0) {
      return -1;
    }
  } else {
    db_leaf(binary, len);
  }

  j = wavetree2d_BL(t, index);
  if (j > 0) {
    if (r_generate_dyck_binary(t, j, binary, len) < 0) {
      return -1;
    }
  } else {
    db_leaf(binary, len);
  }


  j = wavetree2d_BR(t, index);
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

int wavetree2d_generate_dyck_binary(wavetree2d_t *t, uint64_t *binary)
{
  int len;

  *binary = 0;

  len = 0;
  if (r_generate_dyck_binary(t, 0, binary, &len) < 0) {
    return -1;
  }

  return 0;
}

int wavetree2d_get_indices(wavetree2d_t *t, int *set, int *n)
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

int wavetree2d_maxdepth(const wavetree2d_t *t)
{
  return t->degree_max;
}

int wavetree2d_depth(wavetree2d_t *t)
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

int wavetree2d_update_histogram(const wavetree2d_t *t, coefficient_histogram_t *hist)
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

int wavetree2d_depth_filter(void *wt, int subset, int index)
{
  wavetree2d_t *t = (wavetree2d_t*)wt;

  return wavetree2d_depthofindex(t, index) == subset;
}

const multiset_int_double_t *
wavetree2d_get_S_v(const wavetree2d_t *t)
{
  return t->S_v;
}

int
wavetree2d_set_from_S_v(wavetree2d_t *t,
			multiset_int_double_t *S_vp)
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

  if (wavetree2d_initialize(t, value) < 0) {
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
wavetree2d_get_last_perturbation(wavetree2d_t *t,
				 chain_history_change_t *step)
{
  memcpy(step, &t->last_step, sizeof(chain_history_change_t));

  return 0;
}

int
wavetree2d_set_invalid_perturbation(wavetree2d_t *t, wavetree_perturb_t p)
{
  memset(&(t->last_step), 0, sizeof(chain_history_change_t));
  t->last_step.header.type = p;

  return 0;
}

int wavetree2d_perturb(wavetree2d_t *t,
		       wavetree2d_perturb_func_t f,
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

      if (wavetree2d_2dindices(t, i, &ii, &ij) < 0) {
        ERROR("failed to get 3d indices");
        return -1;
      }

      /* Get parent coeff */
      if (index == 0) {
        pvalue = 0.0;
      } else {
        pi = wavetree2d_parent_index(t, index);
        if (wavetree2d_get_coeff(t, pi, &pvalue) < 0) {
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

int wavetree2d_get_coeff(const wavetree2d_t *t,
			 int i,
			 double *coeff)
{
  int d;

  d = wavetree2d_depthofindex(t, i);
  
  if (multiset_int_double_get(t->S_v, i, d, coeff) < 0) {
    return -1;
  }

  return 0;
}

/*
 * Functions for setting up a birth proposal
 */

int wavetree2d_choose_birth_depth(const wavetree2d_t *t, double u, int maxdepth, int *depth, double *prob)
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

int wavetree2d_reverse_birth_depth(const wavetree2d_t *t, int depth, int maxdepth, double *prob)
{
  if (t == NULL) {
    return -1;
  }

  *prob = 1.0/((double)multiset_int_nonempty_count(t->S_d, maxdepth));
  return 0;
}

int wavetree2d_choose_birth(const wavetree2d_t *t, int depth, double u, int *coeff, double *prob)
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

int wavetree2d_choose_birth_global(const wavetree2d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob)
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



int wavetree2d_reverse_birth(const wavetree2d_t *t, int depth, int coeff, double *prob)
{
  if (t == NULL) {
    return -1;
  }
  
  *prob = 1.0/(double)(multiset_int_depth_count(t->S_d, depth));
  return 0;
}

int wavetree2d_reverse_birth_global(const wavetree2d_t *t, int maxdepth, int depth, int coeff, double *prob)
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
int wavetree2d_choose_death_depth(const wavetree2d_t *t, double u, int maxdepth, int *depth, double *prob)
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

int wavetree2d_reverse_death_depth(const wavetree2d_t *t, int depth, int maxdepth, double *prob)
{
  if (t == NULL) {
    return -1;
  }

  *prob = 1.0/((double)multiset_int_nonempty_count(t->S_b, maxdepth));
  return 0;
}

int wavetree2d_choose_death(const wavetree2d_t *t, int depth, double u, int *coeff, double *prob)
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

int wavetree2d_choose_death_global(const wavetree2d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob)
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

int wavetree2d_reverse_death(const wavetree2d_t *t, int depth, int coeff, double *prob)
{
  if (t == NULL) {
    return -1;
  }
  
  *prob = 1.0/(double)(multiset_int_depth_count(t->S_b, depth));
  return 0;
}

int wavetree2d_reverse_death_global(const wavetree2d_t *t, int maxdepth, int depth, int coeff, double *prob)
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
int wavetree2d_choose_value_depth(const wavetree2d_t *t, double u, int maxdepth, int *depth, double *prob)
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

int wavetree2d_choose_value(const wavetree2d_t *t, int depth, double u, int *coeff, double *prob)
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

int wavetree2d_choose_value_global(const wavetree2d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob)
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
int wavetree2d_choose_move_depth(const wavetree2d_t *t, double u, int maxdepth, int *depth, double *prob)
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

int wavetree2d_choose_move(const wavetree2d_t *t, int depth, double u, int *coeff, double *prob)
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

int wavetree2d_choose_move_global(const wavetree2d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob)
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

static int check_move_sibling(const wavetree2d_t *t, int depth, int coeff, int di, int dj, int *sibling)
{
  int scoeff;
  int pi;

  int ii;
  int jj;

  if (wavetree2d_2dindices(t, coeff, &ii, &jj) < 0) {
    ERROR("check_move_sibling: failed to get 2d indices");
    return 0;
  }

  ii += di;
  jj += dj;

  if (ii < 0 || ii >= t->width ||
      jj < 0 || jj >= t->width) {
    return 0;
  }

  scoeff = wavetree2d_from_2dindices(t, ii, jj);
  if (scoeff < 0) {
    ERROR("check_move_sibling: failed to get sibling index");
    return 0;
  }

  /*
   * Must be same depth
   */
  if (wavetree2d_depthofindex(t, scoeff) != depth) {
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
  pi = wavetree2d_parent_index(t, scoeff);
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
wavetree2d_move_available_siblings(const wavetree2d_t *t, int depth, int coeff, int *siblings, int *nsibling)
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

int wavetree2d_choose_move_sibling(const wavetree2d_t *t, double u, int depth, int coeff, int *sibling, double *prob)
{
  int valid_siblings[8];
  int nvalid;
  int si;

  if (wavetree2d_move_available_siblings(t, depth, coeff, valid_siblings, &nvalid) < 0) {
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

int wavetree2d_reverse_choose_move_sibling(const wavetree2d_t *t, int depth, int coeff, int sibling, double *prob)
{
  int valid_siblings[8];
  int nvalid;
  int si;
  int i;

  if (wavetree2d_move_available_siblings(t, depth, sibling, valid_siblings, &nvalid) < 0) {
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


