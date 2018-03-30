

#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "wavetree_value_proposal.h"

#include "slog.h"

/*
 * Depth dependent Cauchy proposal
 */

struct depth_cauchy {
  gsl_rng *rng;

  int errork;

  int ndepths;
  double *std;
};

static int depth_cauchy_init(void *user);

static int depth_cauchy_perturb(wavetree_prior_t *prior,
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

static int depth_cauchy_errors(void *user,
			       int *prior_errors);

static int depth_cauchy_destroy(void *user);

wavetree_value_t *
wavetree_value_create_depth_cauchy(int ndepths,
				   double *std,
				   unsigned long int seed)
{
  wavetree_value_t *w;
  struct depth_cauchy *g;
  int i;

  if (ndepths <= 0) {
    ERROR("invalid parameters");
    return NULL;
  }

  w = malloc(sizeof(wavetree_value_t));
  if (w == NULL) {
    ERROR("failed to allocate structure");
    return NULL;
  }

  w->init = depth_cauchy_init;
  w->perturb = depth_cauchy_perturb;
  w->errors = depth_cauchy_errors;
  w->destroy = depth_cauchy_destroy;
  w->save = NULL;
  w->restore = NULL;

  g = malloc(sizeof(struct depth_cauchy));
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

  g->std = malloc(sizeof(double) * ndepths);
  if (g->std == NULL) {
    ERROR("failed to allocate restore array");
    return NULL;
  }

  g->ndepths = ndepths;
  for (i = 0; i < ndepths; i ++) {
    g->std[i] = std[i];
  }

  w->user = g;

  return w;
}

static int depth_cauchy_init(void *user)
{
  struct depth_cauchy *g = (struct depth_cauchy *)user;

  g->errork = 0;

  return 0;
}

static int depth_cauchy_perturb(wavetree_prior_t *prior,
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
  struct depth_cauchy *g = (struct depth_cauchy *)user;
  int d;
  double new_coeff;
  double p;
  double pp;

  d = level;
  if (d >= g->ndepths) {
    d = g->ndepths - 1;
  }

  new_coeff = (*coeff) + gsl_ran_cauchy(g->rng, g->std[d]);

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

static int depth_cauchy_errors(void *user,
			       int *prior_errors)
{
  struct depth_cauchy *g = (struct depth_cauchy *)user;

  *prior_errors = g->errork;
  
  return 0;
}

static int depth_cauchy_destroy(void *user)
{
  struct depth_cauchy *g = (struct depth_cauchy *)user;
  
  gsl_rng_free(g->rng);

  free(g->std);
  free(g);

  return 0;
}
