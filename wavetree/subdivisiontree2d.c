
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_sf_trig.h>

#include "subdivisiontree2d.h"

#include "multiset_int.h"

#include "slog.h"

typedef enum {
  UNDO_NONE = 0,
  UNDO_VALUE,
  UNDO_BIRTH,
  UNDO_DEATH
} subdivisiontree2d_undo_t;

typedef double (*basis_function_t)(double x, double y, double c, double overlap);

struct _subdivisiontree2d {

  int degree;
  int width;
  int size;
  int totalcoefficients;

  double alpha;
  double overlap;

  int *mask;
  int *cc;
  double *coeff;

  multiset_int_t *S_v;
  multiset_int_t *S_b;
  multiset_int_t *S_d;

  subdivisiontree2d_undo_t undo;
  int u_d;
  int u_i;
  double u_v;

  basis_function_t basis;
};

static int TL(const subdivisiontree2d_t *t, int i);
static int TR(const subdivisiontree2d_t *t, int i);
static int BL(const subdivisiontree2d_t *t, int i);
static int BR(const subdivisiontree2d_t *t, int i);

static int add_node(subdivisiontree2d_t *t, int i, int d, double coeff);
static int remove_node(subdivisiontree2d_t *t, int i, int d);

subdivisiontree2d_t *
subdivisiontree2d_create(int max_depth,
			 double alpha,
			 double overlap,
			 subdivision_basis_t basis)
{
  subdivisiontree2d_t *r;

  if (max_depth < 1) {
    ERROR("invalid depth %d", max_depth);
    return NULL;
  }

  r = malloc(sizeof(subdivisiontree2d_t));
  if (r == NULL) {
    return NULL;
  }

  r->degree = max_depth;
  r->width = 1 << max_depth;
  r->size = r->width * r->width;
  r->totalcoefficients = subdivisiontree2d_total_coefficients(max_depth);

  r->alpha = alpha;
  r->overlap = overlap;

  r->mask = (int*)malloc(sizeof(int) * r->totalcoefficients);
  if (r->mask == NULL) {
    ERROR("failed to allocate mask");
    return NULL;
  }

  r->cc = (int*)malloc(sizeof(int) * r->totalcoefficients);
  if (r->cc == NULL) {
    ERROR("failed to allocate cc");
    return NULL;
  }

  r->coeff = (double*)malloc(sizeof(double) * r->totalcoefficients);
  if (r->coeff == NULL) {
    ERROR("failed to allocate coeffs");
    return NULL;
  }
  
  r->S_v = multiset_int_create();
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

  /* Setting to zero initializes most things correctly. */
  memset(r->mask, 0, sizeof(int) * r->totalcoefficients);
  memset(r->cc, 0, sizeof(int) * r->totalcoefficients);
  memset(r->coeff, 0, sizeof(double) * r->totalcoefficients);

  r->undo = UNDO_NONE;
  r->u_d = -1;
  r->u_i = -1;
  r->u_v = 0.0;

  switch (basis) {
  case SUBDIVISION_BASIS_CONSTANT:
    r->basis = basis_constant;
    break;
    
  case SUBDIVISION_BASIS_PYRAMID:
    r->basis = basis_pyramid;
    break;
    
  case SUBDIVISION_BASIS_LANCZOS:
    r->basis = basis_lanczos;
    break;

  default:
    ERROR("invalid basis");
    return NULL;
  }

  return r;
}

void
subdivisiontree2d_destroy(subdivisiontree2d_t *t)
{
  if (t != NULL) {
    multiset_int_destroy(t->S_v);
    multiset_int_destroy(t->S_d);
    multiset_int_destroy(t->S_b);
    free(t->coeff);
    free(t->cc);
    free(t->mask);
    free(t);
  }
}

int subdivisiontree2d_initialize(subdivisiontree2d_t *t,
				 double dc)
{
  t->mask[0] = 1;
  t->coeff[0] = dc;

  multiset_int_insert(t->S_v, 0, 0);

  multiset_int_insert(t->S_b, TL(t, 0), 1);
  multiset_int_insert(t->S_b, TR(t, 0), 1);
  multiset_int_insert(t->S_b, BL(t, 0), 1);
  multiset_int_insert(t->S_b, BR(t, 0), 1);

  if (multiset_int_total_count(t->S_v) != 1) {
    ERROR("failed to initialize S_v");
    return -1;
  }

  if (multiset_int_total_count(t->S_b) != 4) {
    ERROR("failed to initialize S_b");
    return -1;
  }

  return 0;
}

int
subdivisiontree2d_save(const subdivisiontree2d_t *t,
		       const char *filename)
{
  FILE *fp;
  int i;
  int j;
  int coeffs;
  int idx;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    ERROR("failed to save subdivisiontree");
    return -1;
  }

  fprintf(fp, "%d\n", t->degree);
  fprintf(fp, "%d %d %d\n", t->width, t->size, t->totalcoefficients);
  fprintf(fp, "%.10g\n", t->alpha);

  for (i = 0; i <= t->degree; i ++) {
    coeffs = multiset_int_depth_count(t->S_v, i);
    fprintf(fp, "%d %d\n", i, coeffs);

    for (j = 0; j < coeffs; j ++) {
      
      if (multiset_int_nth_element(t->S_v, i, j, &idx) < 0) {
	ERROR("failed to get nth element");
	return -1;
      }
      
      if (t->mask[idx] != 1) {
	ERROR("element in S_v not in mask");
	return -1;
      }

      fprintf(fp, "%d %d %.10f\n", idx, (int)t->cc[idx], t->coeff[idx]);
    }
  }

  if (multiset_int_write(t->S_b, fp) < 0) {
    ERROR("failed to write set b");
    return -1;
  }

  if (multiset_int_write(t->S_d, fp) < 0) {
    ERROR("failed to write set d");
    return -1;
  }

  fclose(fp);

  return 0;
}

int 
subdivisiontree2d_load(subdivisiontree2d_t *t,
		       const char *filename)
{
  FILE *fp;
  int i;
  int j;
  int coeffs;
  int idx;

  int cc;
  double coeff;

  int degree;
  int width;
  int size;
  int total;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    ERROR("failed to open file");
    return -1;
  }

  if (fscanf(fp, "%d\n", &degree) != 1) {
    ERROR("failed to read degree");
    return -1;
  }

  if (fscanf(fp, "%d %d %d\n", &width, &size, &total) != 3) {
    ERROR("failed to read header");
    return -1;
  }

  if (degree != t->degree ||
      width != t->width ||
      size != t->size ||
      total != t->totalcoefficients) {
    ERROR("size mismatch (%d, %d), (%d %d), (%d %d), (%d %d)",
	  degree, t->degree,
	  width, t->width,
	  size, t->size,
	  total, t->totalcoefficients);
    return -1;
  }

  if (fscanf(fp, "%lf\n", &(t->alpha)) != 1) {
    ERROR("failed to read alpha");
    return -1;
  }

  /* Clear everything first */
  memset(t->mask, 0, sizeof(int) * t->totalcoefficients);
  memset(t->cc, 0, sizeof(int) * t->totalcoefficients);
  memset(t->coeff, 0, sizeof(double) * t->totalcoefficients);
  multiset_int_clear(t->S_v);
  multiset_int_clear(t->S_b);
  multiset_int_clear(t->S_d);

  for (i = 0; i <= t->degree; i ++) {
    if (fscanf(fp, "%d %d\n", &j, &coeffs) != 2) {
      ERROR("failed to read depth no. coefficients");
      return -1;
    }

    if (j != i) {
      ERROR("depth index mismatch");
      return -1;
    }

    for (j = 0; j < coeffs; j ++) {
      if (fscanf(fp, "%d %d %lf\n", &idx, &cc, &coeff) != 3) {
	ERROR("failed to read coefficient");
	return -1;
      }

      t->mask[idx] = 1;
      t->cc[idx] = cc;
      t->coeff[idx] = coeff;

      if (multiset_int_insert(t->S_v, i, idx) < 0) {
	ERROR("failed to insert coefficient");
	return -1;
      }
    }
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


int 
subdivisiontree2d_load_promote(subdivisiontree2d_t *t,
			       const char *filename)
{
  FILE *fp;
  int i;
  int j;
  int coeffs;
  int idx;

  int cc;
  double coeff;

  int degree;
  int width;
  int size;
  int total;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    ERROR("failed to open file");
    return -1;
  }

  if (fscanf(fp, "%d\n", &degree) != 1) {
    ERROR("failed to read degree");
    return -1;
  }

  /* Promotion works for smaller or equal degree */
  if (degree > t->degree) {
    ERROR("degree mismatch %d > %d", 
	  degree, t->degree);
    return -1;
  }

  if (fscanf(fp, "%d %d %d\n", &width, &size, &total) != 3) {
    ERROR("failed to read header");
    return -1;
  }

  if (fscanf(fp, "%lf\n", &(t->alpha)) != 1) {
    ERROR("failed to read alpha");
    return -1;
  }

  /* Clear everything first */
  memset(t->mask, 0, sizeof(char) * t->size);
  memset(t->cc, 0, sizeof(short) * t->size);
  memset(t->coeff, 0, sizeof(double) * t->size);

  multiset_int_clear(t->S_v);
  multiset_int_clear(t->S_b);
  multiset_int_clear(t->S_d);

  for (i = 0; i <= degree; i ++) {
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
      if (fscanf(fp, "%d %d %lf\n", &idx, &cc, &coeff) != 3) {
	ERROR("failed to read coefficient");
	return -1;
      }

      if (idx == 0) {
	if (i != 0) {
	  ERROR("index 0 not at depth 0");
	  return -1;
	}

	if (subdivisiontree2d_initialize(t, coeff) < 0) {
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



static int r_coeff_count(const subdivisiontree2d_t *t, int i)
{
  if (i < 0 || i >= t->totalcoefficients) {
    return 0;
  }

  if (t->mask[i] == 0) {
    return 0;
  }

  return 1 + 
    r_coeff_count(t, TL(t, i)) + 
    r_coeff_count(t, TR(t, i)) + 
    r_coeff_count(t, BL(t, i)) +
    r_coeff_count(t, BR(t, i));
}

int subdivisiontree2d_coeff_count(const subdivisiontree2d_t *t)
{
  return t->mask[0] + 
    r_coeff_count(t, TL(t, 0)) +
    r_coeff_count(t, TR(t, 0)) + 
    r_coeff_count(t, BL(t, 0)) + 
    r_coeff_count(t, BR(t, 0));
}

int subdivisiontree2d_coeff_count_scan(const subdivisiontree2d_t *t)
{
  int i;
  int cc;

  cc = 0;
  for (i = 0; i < t->totalcoefficients; i ++) {
    cc += t->mask[i];
  }

  return cc;
}

static int r_map_to_array(const subdivisiontree2d_t *t,
			  int i,
			  int xi, 
			  int xj,
			  int yi, 
			  int yj,
			  double *a,
			  int n)
{
  int x;
  int y;
  int cx;
  int cy;

  double dx;
  double x0;
  double y0;

  double u;
  double v;

  int dij;
  int dt;

  int oxi, oxj;
  int oyi, oyj;

  if (i < 0 || t->mask[i] == 0) {
    return 0;
  }

  if (i == 0) {
    /*
     * Set
     */
    dx = 2.0/((double)(xj - xi + 1));
    x0 = -1.0 + dx/2.0;
    y0 = x0;

    for (y = yi, v = y0; y <= yj; y ++, v += dx) {
      for (x = xi, u = x0; x <= xj; x ++, u += dx) {
	//a[y*t->width + x] = t->coeff[i];
	a[y*t->width + x] = t->basis(u, v, t->coeff[i], t->overlap);
      }
    }

  } else {
    /*
     * Add
     */
    if (t->coeff[i] != 0.0) {
      dx = 2.0/((double)(xj - xi + 1));
      x0 = -1.0 + dx/2.0;
      y0 = x0;

      if (t->overlap > 0.0) {
	
	if (t->overlap <= 0.5) {
	  dij = (xj - xi + 1)/2;
	} else {
	  dij = (xj - xi + 1);
	}
	
	/*
	 * Adjust xi
	 */
	dt = xi - dij;
	if (dt < 0) {
	  dt = 0;
	}
	x0 -= dx*((double)(xi - dt));
	oxi = dt;
	
	/* 
	 * Adjust xj 
	 */
	dt = xj + dij;
	if (dt >= t->width) {
	  dt = t->width - 1;
	}
	oxj = dt;
	
	/*
	 * Adjust yi
	 */
	dt = yi - dij;
	if (dt < 0) {
	  dt = 0;
	}
	y0 -= dx*((double)(yi - dt));
	oyi = dt;
	
	/* 
	 * Adjust yj 
	 */
	dt = yj + dij;
	if (dt >= t->width) {
	  dt = t->width - 1;
	}
	oyj = dt;
      } else {
	oxi = xi;
	oxj = xj;
	oyi = yi;
	oyj = yj;
      }

      for (y = oyi, v = y0; y <= oyj; y ++, v += dx) {
	for (x = oxi, u = x0; x <= oxj; x ++, u += dx) {
	  a[y*t->width + x] += t->basis(u, v, t->coeff[i], t->overlap);
	}
      }
    }
  }

  cx = (xi + xj)/2;
  cy = (yi + yj)/2;

  if (r_map_to_array(t, TL(t, i), xi, cx, yi, cy, a, n) < 0) {
    return -1;
  }
  if (r_map_to_array(t, TR(t, i), cx + 1, xj, yi, cy, a, n) < 0) {
    return -1;
  }
  if (r_map_to_array(t, BL(t, i), xi, cx, cy + 1, yj, a, n) < 0) {
    return -1;
  }
  if (r_map_to_array(t, BR(t, i), cx + 1, xj, cy + 1, yj, a, n) < 0) {
    return -1;
  }

  return 0;
}

int subdivisiontree2d_map_to_array(const subdivisiontree2d_t *t, 
				   double *a, 
				   int n)
{
  if (n != t->size) {
    ERROR("subdivisiontree2d_map_to_array: size mismatch");
    return -1;
  }
  
  return r_map_to_array(t, 0, 0, t->width - 1, 0, t->width -1, a, n);
}

static int paint_coefficient(const subdivisiontree2d_t *t,
			     double *a,
			     int n,
			     int index,
			     int depth,
			     double value)
{
  int ii, ij, rs;
  int pw, ph;
  int xi, xj;
  int yi, yj;

  double dx;
  double x0, y0;

  int dij;
  int dt;

  double u, v;
  int x, y;

  if (subdivisiontree2d_2dindices(t, index, &ii, &ij, &rs) < 0) {
    return -1;
  }

  pw = t->width/rs;
  ph = t->width/rs;
  xi = pw * ii;
  xj = xi + pw - 1;
  yi = ph * ij;
  yj = yi + ph - 1;

  dx = 2.0/((double)(pw));

  x0 = -1.0 + dx/2.0;
  y0 = x0;
      
  if (t->overlap > 0.0) {
    
    if (t->overlap <= 0.5) {
      dij = (xj - xi + 1)/2;
    } else {
      dij = (xj - xi + 1);
    }
    
    /*
     * Adjust xi
     */
    dt = xi - dij;
    if (dt < 0) {
      dt = 0;
    }
    x0 -= dx*((double)(xi - dt));
    xi = dt;
    
    /* 
     * Adjust xj 
     */
    dt = xj + dij;
    if (dt >= t->width) {
      dt = t->width - 1;
    }
    xj = dt;
    
    /*
     * Adjust yi
     */
    dt = yi - dij;
    if (dt < 0) {
      dt = 0;
    }
    y0 -= dx*((double)(yi - dt));
    yi = dt;
    
    /* 
     * Adjust yj 
     */
    dt = yj + dij;
    if (dt >= t->width) {
      dt = t->width - 1;
    }
    yj = dt;
  }

  
  for (y = yi, v = y0; y <= yj; y ++, v += dx) {
    for (x = xi, u = x0; x <= xj; x ++, u += dx) {
      a[y*t->width + x] += t->basis(u, v, value, t->overlap);
    }
  }

  return 0;
}

int subdivisiontree2d_propose_to_array(const subdivisiontree2d_t *t,
				       double *a,
				       int n)
{
  switch (t->undo) {
  case UNDO_BIRTH:
    /* New index in undo information, new value in actual coefficient */
    if (paint_coefficient(t, a, n, t->u_i, t->u_d, t->coeff[t->u_i]) < 0) {
      ERROR("failed to paint birth value");
      return -1;
    }
    break;

  case UNDO_DEATH:
    /* Removed index in undo information, previous value in undo information */
    if (paint_coefficient(t, a, n, t->u_i, t->u_d, -t->u_v) < 0) {
      ERROR("failed to paint death value");
      return -1;
    }
    break;

  case UNDO_VALUE:
    /* Index in undo information, old value in undo information, new value in actual coefficient */
    if (paint_coefficient(t, a, n, t->u_i, t->u_d, t->coeff[t->u_i] - t->u_v) < 0) {
      ERROR("failed to paint death value");
      return -1;
    }
    break;

  default:
    ERROR("invalid undo state");
    return -1;
  }
  return 0;
}

int subdivisiontree2d_revert_to_array(const subdivisiontree2d_t *t,
				      double *a,
				      int n)
{
  switch (t->undo) {
  case UNDO_BIRTH:
    /* New index in undo information, new value in actual coefficient */
    if (paint_coefficient(t, a, n, t->u_i, t->u_d, -t->coeff[t->u_i]) < 0) {
      ERROR("failed to paint birth value");
      return -1;
    }
    break;

  case UNDO_DEATH:
    /* Removed index in undo information, previous value in undo information */
    if (paint_coefficient(t, a, n, t->u_i, t->u_d, t->u_v) < 0) {
      ERROR("failed to paint death value");
      return -1;
    }
    break;

  case UNDO_VALUE:
    /* Index in undo information, old value in undo information, new value in actual coefficient */
    if (paint_coefficient(t, a, n, t->u_i, t->u_d, t->u_v - t->coeff[t->u_i]) < 0) {
      ERROR("failed to paint death value");
      return -1;
    }
    break;

  default:
    ERROR("invalid undo state");
    return -1;
  }
  return 0;
}

int
subdivisiontree2d_propose_value(subdivisiontree2d_t *t,
				int i,
				int d,
				double value)
{
  if (t == NULL || 
      i < 0 ||
      i >= t->totalcoefficients) {
    ERROR("invalid parameters");
    return -1;
  }

  if (t->mask[i] != 1) {
    ERROR("node not active");
    return -1;
  }

  /*
   * Store the old value in undo information
   */
  t->undo = UNDO_VALUE;
  t->u_d = d;
  t->u_i = i;
  t->u_v = t->coeff[i];

  /*
   * Write the new value
   */
  t->coeff[i] = value;

  return 0;
}

int
subdivisiontree2d_propose_birth(subdivisiontree2d_t *t,
				int i,
				int d,
				double value)
{
  if (t == NULL || 
      i < 0 ||
      i >= t->totalcoefficients) {
    ERROR("invalid parameters %p %i %d %f",
	  t, i, d, value);
    return -1;
  }

  if (d <= 0 || d > t->degree) {
    ERROR("invalid depth");
    return -1;
  }

  if (t->mask[i] != 0) {
    ERROR("non empty node %d %d", i, d);
    return -1;
  }

  /*
   * Store the undo information
   */
  t->undo = UNDO_BIRTH;
  t->u_d = d;
  t->u_i = i;

  return add_node(t, i, d, value);
}

int 
subdivisiontree2d_propose_death(subdivisiontree2d_t *t,
				int i,
				int d,
				double *old_value)
{
  if (t == NULL || 
      i < 0 ||
      i >= t->totalcoefficients) {
    ERROR("invalid parameters");
    return -1;
  }

  if (t->mask[i] != 1) {
    ERROR("empty node");
    return -1;
  }

  *old_value = t->coeff[i];
  t->undo = UNDO_DEATH;
  t->u_d = d;
  t->u_i = i;
  t->u_v = t->coeff[i];

  return remove_node(t, i, d);
}



int
subdivisiontree2d_commit(subdivisiontree2d_t *t)
{
  if (t->undo == UNDO_NONE) {
    ERROR("nothing to commit");
    return -1;
  }

  t->undo = UNDO_NONE;
  t->u_i = -1;
  t->u_d = -1;
  t->u_v = 0.0;
  
  return 0;
}

int 
subdivisiontree2d_undo(subdivisiontree2d_t *t)
{
  switch (t->undo) {
  case UNDO_NONE:
    ERROR("nothing to undo");
    return -1;

  case UNDO_VALUE:
    t->coeff[t->u_i] = t->u_v;
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
  t->u_i = -1;
  t->u_d = -1;
  t->u_v = 0.0;
  
  return 0;
}

int 
subdivisiontree2d_get_width(subdivisiontree2d_t *t)
{
  if (t != NULL) {
    return t->width;
  }
  
  return -1;
}

int
subdivisiontree2d_get_size(subdivisiontree2d_t *t)
{
  if (t != NULL) {
    return t->size;
  }

  return -1;
}

int 
subdivisiontree2d_update_histogram(const subdivisiontree2d_t *t, coefficient_histogram_t *hist)
{
  int i;

  for (i = 0; i < t->totalcoefficients; i ++) {
    if (t->mask[i]) {
      if (coefficient_histogram_sample(hist, i, t->coeff[i]) < 0) {
	return -1;
      }
    }
  }

  return 0;
}

int 
subdivisiontree2d_depth_filter(void *wt, int subset, int index)
{
  subdivisiontree2d_t *t = (subdivisiontree2d_t*)wt;

  return subdivisiontree2d_depthofindex(t, index) == subset;
}

/*
 * Functions for setting up a birth proposal
 */
int 
subdivisiontree2d_choose_birth_depth(const subdivisiontree2d_t *t, double u, int maxdepth, int *depth, double *prob)
{
  int ndepths;

  if (t == NULL) {
    return -1;
  }

  if (multiset_int_choose_depth(t->S_b, u, maxdepth, depth, &ndepths) < 0) {
    return -1;
  }

  *prob = 1.0/((double)ndepths);
  return 0;}

int
subdivisiontree2d_reverse_birth_depth(const subdivisiontree2d_t *t, int depth, int maxdepth, double *prob)
{
  if (t == NULL) {
    return -1;
  }

  *prob = 1.0/((double)multiset_int_nonempty_count(t->S_d, maxdepth));
  return 0;}

int 
subdivisiontree2d_choose_birth(const subdivisiontree2d_t *t, int depth, double u, int *coeff, double *prob)
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
subdivisiontree2d_choose_birth_global(const subdivisiontree2d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob)
{
  int nindices;

  if (t == NULL) {
    return -1;
  }

  if (multiset_int_choose_index_globally(t->S_b, u, maxdepth, coeff, depth, &nindices) < 0) {
    return -1;
  }

  *prob = 1.0/(double)(nindices);
  return 0;
}

int 
subdivisiontree2d_reverse_birth(const subdivisiontree2d_t *t, int depth, int i, double *prob)
{
  if (t == NULL) {
    return -1;
  }
  
  *prob = 1.0/(double)(multiset_int_depth_count(t->S_d, depth));
  return 0;
}

int 
subdivisiontree2d_reverse_birth_global(const subdivisiontree2d_t *t, int maxdepth, int depth, int coeff, double *prob)
{
  if (t == NULL) {
    return -1;
  }

  *prob = 1.0/(double)(multiset_int_restricted_total_count(t->S_d, maxdepth));
  return 0;
}

/*
 * Functions for setting up a death proposal
 */
int subdivisiontree2d_choose_death_depth(const subdivisiontree2d_t *t, double u, int maxdepth, int *depth, double *prob)
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

int subdivisiontree2d_reverse_death_depth(const subdivisiontree2d_t *t, int depth, int maxdepth, double *prob)
{
  if (t == NULL) {
    return -1;
  }

  *prob = 1.0/((double)multiset_int_nonempty_count(t->S_b, maxdepth));
  return 0;
}

int 
subdivisiontree2d_choose_death(const subdivisiontree2d_t *t, int depth, double u, int *coeff, double *prob)
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
subdivisiontree2d_choose_death_global(const subdivisiontree2d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob)
{
  int nindices;

  if (t == NULL) {
    return -1;
  }

  if (multiset_int_choose_index_globally(t->S_d, u, maxdepth, coeff, depth, &nindices) < 0) {
    return -1;
  }

  *prob = 1.0/(double)(nindices);
  return 0;
}

int 
subdivisiontree2d_reverse_death(const subdivisiontree2d_t *t, int depth, int i, double *prob)
{
  if (t == NULL) {
    return -1;
  }
  
  *prob = 1.0/(double)(multiset_int_depth_count(t->S_b, depth));
  return 0;
}

int 
subdivisiontree2d_reverse_death_global(const subdivisiontree2d_t *t, int maxdepth, int depth, int coeff, double *prob)
{
  if (t == NULL) {
    return -1;
  }

  *prob = 1.0/(double)(multiset_int_restricted_total_count(t->S_b, maxdepth));
  return 0;
}

/*
 * Functions for setting up a value proposal
 */
int subdivisiontree2d_choose_value_depth(const subdivisiontree2d_t *t, double u, int maxdepth, int *depth, double *prob)
{
  int ndepths;

  if (t == NULL) {
    return -1;
  }

  if (multiset_int_choose_depth(t->S_v, u, maxdepth, depth, &ndepths) < 0) {
    return -1;
  }

  *prob = 1.0/((double)ndepths);
  return 0;
}

int 
subdivisiontree2d_choose_value(const subdivisiontree2d_t *t, int depth, double u, int *coeff, double *prob)
{
  int nindices;

  if (t == NULL) {
    return -1;
  }

  if (u < 0.0 || u > 1.0) {
    ERROR("subdivisiontree2d_choose_value: invalid random %f", u);
    return -1;
  }

  if (multiset_int_choose_index(t->S_v, depth, u, coeff, &nindices) < 0) {
    return -1;
  }

  if (nindices == 0) {
    ERROR("subdivisiontree2d_choose_value: no coefficients at depth %d", depth);
    return -1;
  }

  *prob = 1.0/(double)(nindices);
  return 0;
}

int
subdivisiontree2d_choose_value_global(const subdivisiontree2d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob)
{
  int nindices;

  if (t == NULL) {
    return -1;
  }

  if (multiset_int_choose_index_globally(t->S_v, u, maxdepth, coeff, depth, &nindices) < 0) {
    return -1;
  }

  *prob = 1.0/(double)(nindices);
  return 0;
}

static int r_subdivisiontree2d_perturb(const subdivisiontree2d_t *t, 
				       int depth,
				       double pcoeff,
				       int i, 
				       subdivisiontree2d_perturb_func_t f, 
				       void *user,
				       double *prior_ratio)
{
  if (i < 0) {
    return 0;
  }

  if (t->mask[i] == 0) {
    return 0;
  }

  if (f(user, i, 0, depth, t->degree, pcoeff, &(t->coeff[i]), prior_ratio) < 0) {
    return -1;
  }

  if (r_subdivisiontree2d_perturb(t, depth + 1, t->coeff[i], TL(t, i), f, user, prior_ratio) < 0) {
    return -1;
  }
                                                       
  if (r_subdivisiontree2d_perturb(t, depth + 1, t->coeff[i], TR(t, i), f, user, prior_ratio) < 0) {
    return -1;
  }

  if (r_subdivisiontree2d_perturb(t, depth + 1, t->coeff[i], BL(t, i), f, user, prior_ratio) < 0) {
    return -1;
  }
  
  if (r_subdivisiontree2d_perturb(t, depth + 1, t->coeff[i], BR(t, i), f, user, prior_ratio) < 0) {
    return -1;
  }

  return 0;
}

int subdivisiontree2d_perturb(subdivisiontree2d_t *t,
			      subdivisiontree2d_perturb_func_t f,
			      void *user, 
			      double *prior_ratio)
{
  return r_subdivisiontree2d_perturb(t, 0, 0.0, 0, f, user, prior_ratio);
}



int subdivisiontree2d_parent_index(subdivisiontree2d_t *t, int c)
{
  int ii;
  int ij;
  int rs;

  if (c == 0) {
    return -1;
  }

  if (subdivisiontree2d_2dindices(t, c, &ii, &ij, &rs) < 0) {
    return -1;
  }

  return subdivisiontree2d_from_2dindices(t, ii/2, ij/2, rs/2);
}

double 
subdivisiontree2d_dc(subdivisiontree2d_t *t)
{
  if (t != NULL) {
    return t->coeff[0];
  }

  return 0.0;
}

int subdivisiontree2d_get_coeff(const subdivisiontree2d_t *t,
				int i,
				double *coeff)
{
  if (t != NULL && i < t->totalcoefficients) {
    *coeff = t->coeff[i];
    return 0;
  }

  return -1;
}

int subdivisiontree2d_valid(subdivisiontree2d_t *t)
{
  int error_count;
  int cc, cc_scan, cc_set;
  int local_cc;
  int i;
  int j;
  int d;

  error_count = 0;

  /*
   * Coefficient count test 
   */
  cc = subdivisiontree2d_coeff_count(t);
  cc_scan = subdivisiontree2d_coeff_count_scan(t);
  cc_set = multiset_int_total_count(t->S_v);
  if (cc != cc_scan ||
      cc != cc_set) {
    ERROR("coeff counts don't match %d != %d != %d",
	  cc, cc_scan, cc_set);
    error_count ++;
  }

  /*
   * Test S_b, S_d and cc are correct
   */
  for (i = 0; i < t->totalcoefficients; i ++) {


    if (t->mask[i]) {

      d = subdivisiontree2d_depthofindex(t, i);

      if (!multiset_int_is_element(t->S_v, i, d)) {
	ERROR(
		"subdivisiontree2d_valid: %d mask set but not in S_v idx %d dep %d",
		i,
		i,
		d);
	error_count ++;
      }

      local_cc = 0;

      j = TL(t, i);
      if (j > 0) { 
	if (t->mask[j]) {
	  local_cc ++;
	} else {
	  if (!multiset_int_is_element(t->S_b, j, d + 1)) {
	    ERROR(
		    "subdivisiontree2d_valid: %d (TL child of %d) not in S_b", j, i);
	    error_count ++;
	  }
	}
      }

      j = TR(t, i);
      if (j > 0) { 
	if (t->mask[j]) {
	  local_cc ++;
	} else {
	  if (!multiset_int_is_element(t->S_b, j, d + 1)) {
	    ERROR(
		    "subdivisiontree2d_valid: %d (TR child of %d) not in S_b", j, i);
	    error_count ++;
	  }
	}
      }
      
      j = BL(t, i);
      if (j > 0) { 
	if (t->mask[j]) {
	  local_cc ++;
	} else {
	  if (!multiset_int_is_element(t->S_b, j, d + 1)) {
	    ERROR(
		    "subdivisiontree2d_valid: %d (BL child of %d) not in S_b", j, i);
	    error_count ++;
	  }
	}
      }

      j = BR(t, i);
      if (j > 0) { 
	if (t->mask[j]) {
	  local_cc ++;
	} else {
	  if (!multiset_int_is_element(t->S_b, j, d + 1)) {
	    ERROR("%d (BR child of %d) not in S_b", j, i);
	    error_count ++;
	  }
	}
      }

      if (local_cc != t->cc[i]) {
	ERROR("cc incorrect (%d != %d)", local_cc, t->cc[i]);
	error_count ++;
      }

      if (local_cc == 0 && i > 0) {
	if (!multiset_int_is_element(t->S_d, i, d)) {
	  ERROR("childless node not in S_d (%d)", 
		i);
	  error_count ++;
	}
      }
    }
  }

  if (error_count > 0) {
    subdivisiontree2d_print(t);
    return 0;
  }

  return -1;
}

/*
 * Exposed internal functions
 */

int
subdivisiontree2d_total_coefficients(int max_depth)
{
  int l;

  if (max_depth < 0) {
    return -1;
  } else if (max_depth == 0) {
    return 1;
  } else {
    l = subdivisiontree2d_total_coefficients(max_depth - 1);
    if (l < 0) {
      return -1;
    }

    return l + (1 << (2 * max_depth));
  }
}

int
subdivisiontree2d_depthofindex(const subdivisiontree2d_t *t,
			       int index)
{
  if (index == 0) {
    return 0;
  }

  return 1 + subdivisiontree2d_depthofindex(t, (index-1)/4);
}

int 
subdivisiontree2d_2dindices(const subdivisiontree2d_t *t,
			    int index,
			    int *ii,
			    int *ij,
			    int *rs)
{
  int d;
  int s;
  int r;
  int tc;
  int ti;

  s = 0;
  r = 1;
  for (d = 0; d <= t->degree; d ++) {
    tc = subdivisiontree2d_total_coefficients(d);

    if (index < tc) {
      ti = index - s;

      *ii = ti % r;
      *ij = ti / r;
      *rs = r;

      return 0;
    }

    s = tc;
    r *= 2;
  }

  return -1;
}

int
subdivisiontree2d_rs_from_depth(int depth)
{
  return 1 << (depth - 1);
}

int
subdivisiontree2d_from_2dindices(const subdivisiontree2d_t *t,
				 int ii,
				 int ij,
				 int rs)
{
  int o;
  int s;
  int i;

  s = rs/2;
  o = 0;
  while (s > 0) {
    o += s*s;
    s /= 2;
  }

  i = o + ij * rs + ii;
  if (i >= t->totalcoefficients) {
    return -1;
  }
  return i;
}

double basis_constant(double x, double y, double c, double overlap)
{
  return c;
}

double basis_pyramid(double x, double y, double c, double overlap)
{
  double rx;

  rx = 1.0 + 2.0*overlap;

  x = fabs(x);
  if (x > rx) {
    return 0.0;
  }

  y = fabs(y);
  if (y > rx) {
    return 0.0;
  }

  if (x > y) {
    return c*(1.0 - x/rx);
  } else {
    return c*(1.0 - y/rx);
  }
  //return c*(1.0 - x/rx)*(1.0 - y/rx);
}

double basis_lanczos(double x, double y, double c, double overlap)
{
  double a;

  a = overlap*2.0 + 1.0;

  return
    c * 
    gsl_sf_sinc(x)*gsl_sf_sinc(x/a) * 
    gsl_sf_sinc(y)*gsl_sf_sinc(y/a);
}

int
subdivisiontree2d_TL(const subdivisiontree2d_t *t,
		     int i)
{
  return TL(t, i);
}

int
subdivisiontree2d_TR(const subdivisiontree2d_t *t,
		     int i)
{
  return TR(t, i);
}

int
subdivisiontree2d_BL(const subdivisiontree2d_t *t,
		     int i)
{
  return BL(t, i);
}

int
subdivisiontree2d_BR(const subdivisiontree2d_t *t,
		     int i)
{
  return BR(t, i);
}

int
subdivisiontree2d_print(const subdivisiontree2d_t *t)
{
  int i;

  for (i = 0; i < t->totalcoefficients; i ++) {
    if (t->mask[i]) {
      printf("%d: %f %d\n", i, t->coeff[i], t->cc[i]);
    }
  }
  
  printf("S_b: ");
  multiset_int_dump(t->S_b);

  printf("S_d: ");
  multiset_int_dump(t->S_d);

  return 0;
}

/*
 * Internal functions
 */
static int TL(const subdivisiontree2d_t *t, int i)
{
  int ii;
  int ij;
  int rs;

  if (subdivisiontree2d_2dindices(t,
				  i,
				  &ii,
				  &ij,
				  &rs) < 0) {
    return -1;
  }

  return subdivisiontree2d_from_2dindices(t,
					  2*ii,
					  2*ij,
					  2*rs);
}

static int TR(const subdivisiontree2d_t *t, int i)
{
  int ii;
  int ij;
  int rs;

  if (subdivisiontree2d_2dindices(t,
				  i,
				  &ii,
				  &ij,
				  &rs) < 0) {
    return -1;
  }

  return subdivisiontree2d_from_2dindices(t,
					  2*ii + 1,
					  2*ij,
					  2*rs);
}

static int BL(const subdivisiontree2d_t *t, int i)
{
  int ii;
  int ij;
  int rs;

  if (subdivisiontree2d_2dindices(t,
				  i,
				  &ii,
				  &ij,
				  &rs) < 0) {
    return -1;
  }

  return subdivisiontree2d_from_2dindices(t,
					  2*ii,
					  2*ij + 1,
					  2*rs);
}

static int BR(const subdivisiontree2d_t *t, int i)
{
  int ii;
  int ij;
  int rs;


  if (subdivisiontree2d_2dindices(t,
				  i,
				  &ii,
				  &ij,
				  &rs) < 0) {
    return -1;
  }

  return subdivisiontree2d_from_2dindices(t,
					  2*ii + 1,
					  2*ij + 1,
					  2*rs);
}



static int add_node(subdivisiontree2d_t *t, int i, int d, double coeff)
{
  int j;

  if (i < 0 || i >= t->totalcoefficients) {
    ERROR("index out of range %d (%d)",
	  i, 
	  t->totalcoefficients);
    return -1;
  }

  if (t->mask[i] != 0) {
    ERROR("node not empty");
    return -1;
  }

  t->mask[i] = 1;
  t->coeff[i] = coeff;
  if (multiset_int_insert(t->S_v, i, d) != 1) {
    ERROR("failed to add to S_v");
  }

  /*
   * Update parent CC
   */
  j = subdivisiontree2d_parent_index(t, i);
  /* printf("add_node: %d (%d)\n", i, j); */

  if (j >= 0) {
    t->cc[j] ++;
    /* printf("Incrementing parent %d from %d\n", j, i); */
  } else {
    printf("add_node: got parent %d for node %d\n", j, i);
  }

  /*
   * Update sets Sb and Sd
   */
  if (t->cc[j] == 1) {
    /* Remove parent from Sd, only needs to be done if this is its first child */
    multiset_int_remove(t->S_d, j, d - 1);
  }

  multiset_int_insert(t->S_d, i, d);
  multiset_int_remove(t->S_b, i, d);
  
  j = TL(t,i);
  if (j > 0) {
    multiset_int_insert(t->S_b, j, d + 1);
  }

  j = TR(t,i);
  if (j > 0) {
    multiset_int_insert(t->S_b, j, d + 1);
  }

  j = BL(t,i);
  if (j > 0) {
    multiset_int_insert(t->S_b, j, d + 1);
  }

  j = BR(t,i);
  if (j > 0) {
    multiset_int_insert(t->S_b, j, d + 1);
  }

  /*printf("Add: %d %d\n", multiset_int_count(t->S_b), multiset_int_count(t->S_d)); */

  return 0;
}


static int remove_node(subdivisiontree2d_t *t, int i, int d)
{
  int j;

  if (t->mask[i] != 1) {
    ERROR("node not full");
    return -1;
  }

  //printf("removed node: %d\n", i);

  t->mask[i] = 0;
  t->coeff[i] = 0.0;
  multiset_int_remove(t->S_v, i, d);

  /*
   * Update parent CC
   */
  j = subdivisiontree2d_parent_index(t, i);
  if (j >= 0) {
    t->cc[j] --;
  }
  
  /*
   * Update sets Sb and Sd
   */
  if (t->cc[j] == 0 && j > 0) {
    /* Add parent to Sd, now has no children. */
    multiset_int_insert(t->S_d, j, d - 1);
  }

  multiset_int_remove(t->S_d, i, d);
  multiset_int_insert(t->S_b, i, d);

  j = TL(t,i);
  if (j > 0) {
    multiset_int_remove(t->S_b, j, d + 1);
  }

  j = TR(t,i);
  if (j > 0) {
    multiset_int_remove(t->S_b, j, d + 1);
  }

  j = BL(t,i);
  if (j > 0) {
    multiset_int_remove(t->S_b, j, d + 1);
  }

  j = BR(t,i);
  if (j > 0) {
    multiset_int_remove(t->S_b, j, d + 1);
  }
  
  /*  printf("Remove: %d %d\n", t->N_b, t->N_d); */

  return 0;
}

