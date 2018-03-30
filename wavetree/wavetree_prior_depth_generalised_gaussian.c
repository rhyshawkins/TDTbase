
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "wavetree_prior.h"

/*
 * Coefficients generalised gaussian at each depth
 */

struct depth_generalised_gaussian {
  gsl_rng *rng;

  int ndepths;
  double *va;
  double beta;
};

static int 
range_depth_generalised_gaussian(void *user,
				 int i,
				 int j,
				 int k,
				 int level,
				 int maxlevel,
				 double parent_coeff,
				 double *vmin,
				 double *vmax)
{
  struct depth_generalised_gaussian *s; 
  int d;

  s = (struct depth_generalised_gaussian *)user;

  d = level;
  if (d >= s->ndepths) {
    d = s->ndepths - 1;
  }

  /*
   * Just use 3 standard deviations for width
   */
  *vmax = s->va[d] * 3.0;
  *vmin = -(*vmax);

  return 0;
}

static int 
sample_depth_generalised_gaussian(void *user,
				  int i,
				  int j,
				  int k,
				  int level,
				  int maxlevel,
				  double parent_coeff,
				  double *coeff)
{
  struct depth_generalised_gaussian *s; 
  int d;

  s = (struct depth_generalised_gaussian *)user;

  d = level;
  if (d >= s->ndepths) {
    d = s->ndepths - 1;
  }
  *coeff = gsl_ran_exppow(s->rng, s->va[d], s->beta);
  
  return 0;
}

static double
prob_depth_generalised_gaussian(void *user,
		   int i,
		   int j,
		   int k,
		   int level,
		   int maxlevel,
		   double parent_coeff,
		   double coeff)
{
  struct depth_generalised_gaussian *s; 
  int d;
  double p;
  s = (struct depth_generalised_gaussian *)user;

  d = level;
  if (d >= s->ndepths) {
    d = s->ndepths - 1;
  }

  p = gsl_ran_exppow_pdf(coeff, s->va[d], s->beta);
  /* printf("  p = %g (%d %f %f %f)\n", p, d, coeff, s->va[d], s->beta); */
  return p;
}

static int 
valid_depth_generalised_gaussian(void *user,
		    int i,
		    int j,
		    int k,
		    int level,
		    int maxlevel,
		    double parent_coeff,
		    double coeff)
{
  /* All values are valid */
  return -1;
}

static int
setscale_depth_generalised_gaussian(void *user,
				    double newscale,
				    double *oldscale)
{
  return -1;
}

static int 
destroy_depth_generalised_gaussian(void *user)
{
  struct depth_generalised_gaussian *s;

  s = (struct depth_generalised_gaussian *)user;

  gsl_rng_free(s->rng);
  free(s->va);

  free(s);

  return 0;
}

wavetree_prior_t*
wavetree_prior_create_depth_generalised_gaussian(int ndepths,
						 double *va,
						 double beta,
						 unsigned long int seed)
{
  wavetree_prior_t *w;
  struct depth_generalised_gaussian *s;
  int i;

  w = malloc(sizeof(wavetree_prior_t));
  if (w == NULL) {
    return NULL;
  }

  s = malloc(sizeof(struct depth_generalised_gaussian));
  if (s == NULL) {
    return NULL;
  }

  s->rng = gsl_rng_alloc(gsl_rng_taus);
  if (s->rng == NULL) {
    return NULL;
  }

  gsl_rng_set(s->rng, seed);

  s->ndepths = ndepths;

  s->va = malloc(sizeof(double) * ndepths);
  if (s->va == NULL) {
    return NULL;
  }

  for (i = 0; i < ndepths; i ++) {
    s->va[i] = va[i];
  }
  s->beta = beta;

  w->user = s;
  w->range = range_depth_generalised_gaussian;
  w->sample = sample_depth_generalised_gaussian;
  w->prob = prob_depth_generalised_gaussian;
  w->valid = valid_depth_generalised_gaussian;
  w->setscale = setscale_depth_generalised_gaussian;
  w->destroy = destroy_depth_generalised_gaussian;

  return w;
}


