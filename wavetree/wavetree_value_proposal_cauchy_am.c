
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "wavetree_value_proposal.h"

#include "slog.h"


/*
 * Depth dependent Cauchy proposal
 */

struct cauchy_am {
  coefficient_histogram_t *histogram;
  
  gsl_rng *rng;

  int errork;

  double *std;    /* [histogram->ncoeff] */

  double epsilon; /* Min std deviation */
  double A;       /* Max std deviation */
  double tau;     /* Target alpha */
};

static int cauchy_am_init(void *user);

static int cauchy_am_perturb(wavetree_prior_t *prior,
			     void *user,
			     int i,
			     int j,
			     int k,
			     int level,
			     int maxlevel,
			     double parent_coeff,
			     double temperature,
			     double *coeff,
			     double *prior_ratio);

static int cauchy_am_errors(void *user,
			    int *prior_errors);

static int cauchy_am_destroy(void *user);

static int cauchy_am_save(void *user, const char *filename);

static int cauchy_am_restore(void *user, const char *filename);

wavetree_value_t *
wavetree_value_create_cauchy_am(double std0,
				double epsilon,
				double A,
				double tau,
				unsigned long int seed,
				coefficient_histogram_t *histogram)
{
  wavetree_value_t *w;
  struct cauchy_am *g;
  int i;

  if (histogram == NULL) {
    ERROR("NULL histogram");
    return NULL;
  }

  w = malloc(sizeof(wavetree_value_t));
  if (w == NULL) {
    ERROR("failed to allocate structure");
    return NULL;
  }

  w->init = cauchy_am_init;
  w->perturb = cauchy_am_perturb;
  w->errors = cauchy_am_errors;
  w->destroy = cauchy_am_destroy;
  w->save = cauchy_am_save;
  w->restore = cauchy_am_restore;

  g = malloc(sizeof(struct cauchy_am));
  if (g == NULL) {
    ERROR("failed to allocate cauchy structure");
    return NULL;
  }

  g->rng = gsl_rng_alloc(gsl_rng_taus);
  if (g->rng == NULL) {
    ERROR("failed to create rng");
    return NULL;
  }

  gsl_rng_set(g->rng, seed);
  
  g->errork = 0;

  g->std = malloc(sizeof(double) * histogram->ncoeff);
  if (g->std == NULL) {
    ERROR("failed to allocate std array");
    return NULL;
  }

  for (i = 0; i < histogram->ncoeff; i ++) {
    g->std[i] = std0;
  }

  g->histogram = histogram;

  
  g->epsilon = epsilon;
  g->A = A;
  g->tau = tau;

  w->user = g;

  return w;
}

static int cauchy_am_init(void *user)
{
  struct cauchy_am *g = (struct cauchy_am *)user;

  g->errork = 0;

  return 0;
}

static int cauchy_am_perturb(wavetree_prior_t *prior,
			     void *user,
			     int i,
			     int j,
			     int k,
			     int level,
			     int maxlevel,
			     double parent_coeff,
			     double temperature,
			     double *coeff,
			     double *prior_ratio)
{
  struct cauchy_am *g = (struct cauchy_am *)user;
  double new_coeff;
  double p;
  double pp;

  double alpha;
  int index;

  index = coefficient_histogram_coord_to_index(g->histogram, i, j, k, level);
  if (index < 0) {
    ERROR("failed to convert coords to index");
    return -1;
  }
  
  alpha = g->histogram->valpha_mean[index];

  /*
   * alpha is in log space
   */
  if (alpha != 0.0) {
    /*
     * Update the std dev based on the previous alpha and parameters
     */
    /* if (index == 0) { */
    /*   printf("%10.6f (%10.6f) %10.6f -> %10.6f (%f %f)\n", exp(alpha), alpha, g->tau, g->std[index], g->epsilon, g->A); */
    /* } */
    
    if (exp(alpha) < g->tau) {
      g->std[index] *= (1.0 - 0.05);
    } else {
      g->std[index] *= (1.0 + 0.05);
    }
    /* g->std[index] = g->std[index] + 0.01*(exp(alpha) - g->tau); */
    
    if (g->std[index] < g->epsilon) {
      g->std[index] = g->epsilon;
    }
    
    if (g->std[index] > g->A) {
      g->std[index] = g->A;
    }
  }

  /* if (index == 0) { */
  /*   printf("am: %f %f (%f)\n", g->std[index], alpha, g->tau); */
  /* } */

  new_coeff = (*coeff) + gsl_ran_cauchy(g->rng, g->std[index]);

  if (prior->valid(prior->user,
		   i,
		   j,
		   k, 
		   level,
		   maxlevel,
		   parent_coeff,
		   new_coeff)) {
  } else {
    g->errork ++;
  }

  p = prior->prob(prior->user,
		  i, j, k, level, maxlevel, parent_coeff, *coeff);
  
  pp = prior->prob(prior->user,
		   i, j, k, level, maxlevel, parent_coeff, new_coeff);

  *prior_ratio *= pp/p;

  *coeff = new_coeff;

  return 0;
}

static int cauchy_am_errors(void *user,
			    int *prior_errors)
{
  struct cauchy_am *g = (struct cauchy_am *)user;

  *prior_errors = g->errork;
  
  return 0;
}

static int cauchy_am_destroy(void *user)
{
  struct cauchy_am *g = (struct cauchy_am *)user;
  
  gsl_rng_free(g->rng);

  free(g->std);
  free(g);

  return 0;
}

static int cauchy_am_save(void *user, const char *filename)
{
  struct cauchy_am *g = (struct cauchy_am *)user;
  FILE *fp;
  int n;
  
  fp = fopen(filename, "w");
  if (fp == NULL) {
    ERROR("failed to create file");
    return -1;
  }

  if (fwrite(&g->histogram->ncoeff, sizeof(int), 1, fp) != 1) {
    ERROR("failed to write ncoeff");
    return -1;
  }

  n = fwrite(g->std, sizeof(double), g->histogram->ncoeff, fp);
  if (n != g->histogram->ncoeff) {
    ERROR("failed to write std %d %d", n, g->histogram->ncoeff);
    return -1;
  }

  fclose(fp);

  return 0;
}

static int cauchy_am_restore(void *user, const char *filename)
{
  struct cauchy_am *g = (struct cauchy_am *)user;
  FILE *fp;
  int ncoeff;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    ERROR("failed to open file");
    return -1;
  }

  if (fread(&ncoeff, sizeof(int), 1, fp) != 1) {
    ERROR("failed to read ncoeff");
    return -1;
  }

  if (ncoeff != g->histogram->ncoeff) {
    ERROR("ncoeff mismatch %d != %d", ncoeff, g->histogram->ncoeff);
    return -1;
  }
  
  if (fread(g->std, sizeof(double), ncoeff, fp) != ncoeff) {
    ERROR("failed to read std");
    return -1;
  }

  fclose(fp);

  return 0;
}
  

