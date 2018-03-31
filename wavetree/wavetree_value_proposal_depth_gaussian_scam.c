
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "wavetree_value_proposal.h"

#include "slog.h"

/*
 * Depth dependent Gaussian proposal
 */

struct depth_gaussian_scam {
  coefficient_histogram_t *histogram;
  
  gsl_rng *rng;

  int errork;

  int ndepths;
  double *std;

  double epsilon;
  double s;
  int threshold;
};

static int depth_gaussian_scam_init(void *user);

static int depth_gaussian_scam_perturb(wavetree_prior_t *prior,
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

static int depth_gaussian_scam_errors(void *user,
				  int *prior_errors);

static int depth_gaussian_scam_destroy(void *user);

wavetree_value_t *
wavetree_value_create_depth_gaussian_scam(int ndepths,
					  double *std,
					  double epsilon,
					  double s,
					  int threshold,
					  unsigned long int seed,
					  coefficient_histogram_t *histogram)
{
  wavetree_value_t *w;
  struct depth_gaussian_scam *g;
  int i;

  if (ndepths <= 0) {
    ERROR("invalid parameters");
    return NULL;
  }

  if (histogram == NULL) {
    ERROR("NULL histogram");
    return NULL;
  }

  w = malloc(sizeof(wavetree_value_t));
  if (w == NULL) {
    ERROR("failed to allocate structure\n");
    return NULL;
  }

  w->init = depth_gaussian_scam_init;
  w->perturb = depth_gaussian_scam_perturb;
  w->errors = depth_gaussian_scam_errors;
  w->destroy = depth_gaussian_scam_destroy;
  w->save = NULL;
  w->restore = NULL;

  g = malloc(sizeof(struct depth_gaussian_scam));
  if (g == NULL) {
    ERROR("failed to allocate gaussian structure");
    return NULL;
  }

  g->rng = gsl_rng_alloc(gsl_rng_taus);
  if (g->rng == NULL) {
    ERROR("failed to create rng");
    return NULL;
  }

  gsl_rng_set(g->rng, seed);
  
  g->errork = 0;

  g->std = malloc(sizeof(double) * ndepths);
  if (g->std == NULL) {
    ERROR("failed to allocate restore array");
    return NULL;
  }

  g->ndepths = ndepths;
  for (i = 0; i < ndepths; i ++) {
    g->std[i] = std[i];
  }

  g->histogram = histogram;

  g->epsilon = epsilon;
  g->s = s;
  g->threshold = threshold;

  w->user = g;

  return w;
}

static int depth_gaussian_scam_init(void *user)
{
  struct depth_gaussian_scam *g = (struct depth_gaussian_scam *)user;

  g->errork = 0;

  return 0;
}

static int depth_gaussian_scam_perturb(wavetree_prior_t *prior,
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
  struct depth_gaussian_scam *g = (struct depth_gaussian_scam *)user;
  int d;
  double new_coeff;
  double p;
  double pp;

  double mean;
  double std;
  int n;
  int index;

  d = level;
  if (d >= g->ndepths) {
    d = g->ndepths - 1;
  }

  index = coefficient_histogram_coord_to_index(g->histogram, i, j, k, level);
  if (index < 0) {
    ERROR("failed to convert coords to index");
    return -1;
  }

  n = coefficient_histogram_get_coefficient_mean_std(g->histogram,
						     index,
						     &mean,
						     &std);
  if (n < g->threshold) {
    /* Use default initial gaussian width during "burnin period" */
    new_coeff = (*coeff) + gsl_ran_gaussian_ziggurat(g->rng, g->std[d]);
  } else {
    /* Use recorded coefficient variance for gaussian proposal width */
    new_coeff = (*coeff) + gsl_ran_gaussian_ziggurat(g->rng, sqrt(g->s * (std*std + g->epsilon)));
  }

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

static int depth_gaussian_scam_errors(void *user,
				 int *prior_errors)
{
  struct depth_gaussian_scam *g = (struct depth_gaussian_scam *)user;

  *prior_errors = g->errork;
  
  return 0;
}

static int depth_gaussian_scam_destroy(void *user)
{
  struct depth_gaussian_scam *g = (struct depth_gaussian_scam *)user;
  
  gsl_rng_free(g->rng);

  free(g->std);
  free(g);

  return 0;
}