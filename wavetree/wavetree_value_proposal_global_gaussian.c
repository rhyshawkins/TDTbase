#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "wavetree_value_proposal.h"

#include "slog.h"

struct global_gaussian {
  gsl_rng *rng;

  int errork;

  double std;
};

static int global_gaussian_init(void *user);

static int global_gaussian_perturb(wavetree_prior_t *prior,
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

static int global_gaussian_errors(void *user,
				  int *prior_errors);

static int global_gaussian_destroy(void *user);
	       
wavetree_value_t *
wavetree_value_create_global_gaussian(double std,
				      unsigned long int seed)
{
  wavetree_value_t *w;
  struct global_gaussian *g;

  w = malloc(sizeof(wavetree_value_t));
  if (w == NULL) {
    ERROR("failed to allocate structure");
    return NULL;
  }

  w->init = global_gaussian_init;
  w->perturb = global_gaussian_perturb;
  w->errors = global_gaussian_errors;
  w->destroy = global_gaussian_destroy;
  w->save = NULL;
  w->restore = NULL;

  g = malloc(sizeof(struct global_gaussian));
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

  g->std = std;

  w->user = g;

  return w;
}

static int global_gaussian_init(void *user)
{
  struct global_gaussian *g = (struct global_gaussian *)user;

  g->errork = 0;

  return 0;
}

static int global_gaussian_perturb(wavetree_prior_t *prior,
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
  struct global_gaussian *g = (struct global_gaussian *)user;

  double new_coeff = (*coeff) + sqrt(temperature) * gsl_ran_gaussian_ziggurat(g->rng, g->std);

  double p;
  double pp;

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

static int global_gaussian_errors(void *user,
				  int *prior_errors)
{
  struct global_gaussian *g = (struct global_gaussian *)user;

  *prior_errors = g->errork;
  
  return 0;
}

static int global_gaussian_destroy(void *user)
{
  struct global_gaussian *g = (struct global_gaussian *)user;

  gsl_rng_free(g->rng);
  free(g);
  
  return 0;
}

