
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "wavetree_prior.h"

struct globally_laplacian {
  gsl_rng *rng;

  double b;
};

static int
range_globally_laplacian(void *user,
			 int i,
			 int j,
			 int k,
			 int level,
			 int maxlevel,
			 double parent_coeff,
			 double *vmin,
			 double *vmax);

static int 
sample_globally_laplacian(void *user,
			  int i,
			  int j,
			  int k,
			  int level,
			  int maxlevel,
			  double parent_coeff,
			  double *coeff);

static double
prob_globally_laplacian(void *user,
			int i,
			int j,
			int k,
			int level,
			int maxlevel,
			double parent_coeff,
			double coeff);

static int 
valid_globally_laplacian(void *user,
			 int i,
			 int j,
			 int k,
			 int level,
			 int maxlevel,
			 double parent_coeff,
			 double coeff);

static int
setscale_globally_laplacian(void *user,
			    double scale,
			    double *oldscale);

static int
destroy_globally_laplacian(void *user);

wavetree_prior_t *
wavetree_prior_create_globally_laplacian(double b,
					 unsigned long int seed)
{
  wavetree_prior_t *w;
  struct globally_laplacian *s;

  w = malloc(sizeof(wavetree_prior_t));
  if (w == NULL) {
    return NULL;
  }

  s = malloc(sizeof(struct globally_laplacian));
  if (s == NULL) {
    return NULL;
  }

  s->rng = gsl_rng_alloc(gsl_rng_taus);
  if (s->rng == NULL) {
    return NULL;
  }

  gsl_rng_set(s->rng, seed);

  s->b = b;

  w->user = s;
  w->range = range_globally_laplacian;
  w->sample = sample_globally_laplacian;
  w->prob = prob_globally_laplacian;
  w->valid = valid_globally_laplacian;
  w->setscale = setscale_globally_laplacian;
  w->destroy = destroy_globally_laplacian;

  return w;
}

void
wavetree_prior_free_globally_laplacian(wavetree_prior_t *w)
{
  struct globally_laplacian *s;

  if (w != NULL) {
    s = (struct globally_laplacian *)w->user;
    gsl_rng_free(s->rng);
    free(s);
    free(w);
  }
}


static int
range_globally_laplacian(void *user,
			 int i,
			 int j,
			 int k,
			 int level,
			 int maxlevel,
			 double parent_coeff,
			 double *vmin,
			 double *vmax)
{
  struct globally_laplacian *s; 

  s = (struct globally_laplacian *)user;

  /*
   * Just use 3 standard deviations for width
   */
  *vmax = s->b * 3.0;
  *vmin = -(*vmax);

  return 0;
}

static int 
sample_globally_laplacian(void *user,
			  int i,
			  int j,
			  int k,
			  int level,
			  int maxlevel,
			  double parent_coeff,
			  double *coeff)
{
  struct globally_laplacian *s; 

  s = (struct globally_laplacian *)user;

  *coeff = gsl_ran_laplace(s->rng, s->b);
  
  return 0;
}

static double
prob_globally_laplacian(void *user,
			int i,
			int j,
			int k,
			int level,
			int maxlevel,
			double parent_coeff,
			double coeff)
{
  struct globally_laplacian *s; 

  s = (struct globally_laplacian *)user;

  return gsl_ran_laplace_pdf(coeff, s->b);
}

static int 
valid_globally_laplacian(void *user,
			 int i,
			 int j,
			 int k,
			 int level,
			 int maxlevel,
			 double parent_coeff,
			 double coeff)
{
  return -1;
}

static int
setscale_globally_laplacian(void *user,
			    double scale,
			    double *oldscale)
{
  struct globally_laplacian *s; 

  s = (struct globally_laplacian *)user;

  if (oldscale != NULL) {
    *oldscale = s->b;
  }

  if (scale > 0.0) {
    s->b = scale;

    return 0;
  }

  return -1;
}

static int
destroy_globally_laplacian(void *user)
{
  struct globally_laplacian *s; 

  s = (struct globally_laplacian *)user;

  gsl_rng_free(s->rng);

  free(s);

  return 0;
}
