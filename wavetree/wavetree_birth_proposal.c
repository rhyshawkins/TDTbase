
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "wavetree_birth_proposal.h"

#include "slog.h"

/*
 * Local user arguments
 */
void
wavetree_birth_destroy(wavetree_bd_t *p)
{
  p->destroy(p->user);
  free(p);
}

/*
 * Local proposal functions
 */
static int birth_from_prior(void *user,
			    int i,
			    int j,
			    int k,
			    int level,
			    int maxlevel,
			    double parent_coeff,
			    double *coeff,
			    double *prob);

static int death_from_prior(void *user,
			    int i,
			    int j,
			    int k,
			    int level,
			    int maxlevel,
			    double parent_coeff,
			    double coeff,
			    double *prob);

static int destroy_from_prior(void *user);

/*
 * Local creation functions
 */

wavetree_bd_t *
wavetree_birth_create_birth_from_prior(wavetree_prior_t *p)
{
  wavetree_bd_t *w;

  w = malloc(sizeof(wavetree_bd_t));
  if (w == NULL) {
    ERROR("failed to allocate structure");
    return NULL;
  }
  
  w->user = p;

  w->birth = birth_from_prior;
  w->death = death_from_prior;
  w->destroy = destroy_from_prior;

  return w;
}

static int 
birth_from_prior(void *user,
		 int i,
		 int j,
		 int k,
		 int level,
		 int maxlevel,
		 double parent_coeff,
		 double *coeff,
		 double *prob)
{
  wavetree_prior_t *wp;

  wp = (wavetree_prior_t*)user;
  if (wp->sample(wp->user, i, j, k, level, maxlevel, parent_coeff, coeff) < 0) {
    return -1;
  }

  *prob = wp->prob(wp->user, i, j, k, level, maxlevel, parent_coeff, *coeff);

  return 0;
}

static int 
death_from_prior(void *user,
		 int i,
		 int j,
		 int k,
		 int level,
		 int maxlevel,
		 double parent_coeff,
		 double coeff,
		 double *prob)
{
  wavetree_prior_t *wp;

  wp = (wavetree_prior_t*)user;

  *prob = wp->prob(wp->user, i, j, k, level, maxlevel, parent_coeff, coeff);
  return 0;
}

static int
destroy_from_prior(void *user)
{
  return 0;
}
		   

/*
 * Gaussian proposal
 */

struct bd_gaussian {
  gsl_rng *rng;

  double sigma;
};

static int birth_gaussian(void *user,
			  int i,
			  int j,
			  int k,
			  int level,
			  int maxlevel,
			  double parent_coeff,
			  double *coeff,
			  double *prob);

static int death_gaussian(void *user,
			  int i,
			  int j,
			  int k,
			  int level,
			  int maxlevel,
			  double parent_coeff,
			  double coeff,
			  double *prob);

static int destroy_gaussian(void *user);

wavetree_bd_t *
wavetree_birth_create_gaussian(wavetree_prior_t *p,
			       double std,
			       unsigned long int seed)
{
  wavetree_bd_t *w;
  struct bd_gaussian *s;

  w = malloc(sizeof(wavetree_bd_t));
  if (w == NULL) {
    ERROR("failed to allocate structure");
    return NULL;
  }

  s = malloc(sizeof(struct bd_gaussian));
  if (s == NULL) {
    ERROR("failed to allocate user structure");
    return NULL;
  }

  s->rng = gsl_rng_alloc(gsl_rng_taus);
  if (s->rng == NULL) {
    ERROR("failed to create rng");
    return NULL;
  }

  gsl_rng_set(s->rng, seed);
  
  s->sigma = std;
  
  w->user = s;
  w->birth = birth_gaussian;
  w->death = death_gaussian;
  w->destroy = destroy_gaussian;

  return w;
}

void
wavetree_birth_free_gaussian(wavetree_bd_t *w)
{
  struct bd_gaussian *s;
  if (w != NULL) {

    s = (struct bd_gaussian *)w->user;
    gsl_rng_free(s->rng);
    free(s);

    free(w);
  }
}

static int 
birth_gaussian(void *user,
	       int i,
	       int j,
	       int k,
	       int level,
	       int maxlevel,
	       double parent_coeff,
	       double *coeff,
	       double *prob)
{
  struct bd_gaussian *s;
  
  s = (struct bd_gaussian *)user;

  *coeff = gsl_ran_gaussian_ziggurat(s->rng, s->sigma);
  *prob = gsl_ran_gaussian_pdf(*coeff, s->sigma);

  return 0;
}

static int 
death_gaussian(void *user,
	       int i,
	       int j,
	       int k,
	       int level,
	       int maxlevel,
	       double parent_coeff,
	       double coeff,
	       double *prob)
{
  struct bd_gaussian *s;

  s = (struct bd_gaussian *)user;

  *prob = gsl_ran_gaussian_pdf(coeff, s->sigma);
  return 0;
}

static int destroy_gaussian(void *user)
{
  struct bd_gaussian *s;

  s = (struct bd_gaussian *)user;
  
  gsl_rng_free(s->rng);
  free(s);

  return 0;
}

/*
 * Gaussian based on depth proposal
 */

struct bd_depth_gaussian {
  gsl_rng *rng;

  int ndepths;
  double *sigma;
};

static int birth_depth_gaussian(void *user,
				int i,
				int j,
				int k,
				int level,
				int maxlevel,
				double parent_coeff,
				double *coeff,
				double *prob);

static int death_depth_gaussian(void *user,
				int i,
				int j,
				int k,
				int level,
				int maxlevel,
				double parent_coeff,
				double coeff,
				double *prob);

static int destroy_depth_gaussian(void *user);

wavetree_bd_t *
wavetree_birth_create_depth_gaussian(wavetree_prior_t *p,
				     int ndepths,
				     double *std,
				     unsigned long int seed)
{
  wavetree_bd_t *w;
  struct bd_depth_gaussian *s;
  int i;

  w = malloc(sizeof(wavetree_bd_t));
  if (w == NULL) {
    ERROR("failed to allocate structure");
    return NULL;
  }

  s = malloc(sizeof(struct bd_depth_gaussian));
  if (s == NULL) {
    ERROR("failed to allocate user structure");
    return NULL;
  }

  s->rng = gsl_rng_alloc(gsl_rng_taus);
  if (s->rng == NULL) {
    ERROR("failed to create rng");
    return NULL;
  }

  gsl_rng_set(s->rng, seed);
  
  s->sigma = malloc(sizeof(double) * ndepths);
  if (s->sigma == NULL) {
    ERROR("failed to create array");
    return NULL;
  }
  
  for (i = 0; i < ndepths; i ++) {
    s->sigma[i] = std[i];
  }
  s->ndepths = ndepths;
  
  w->user = s;
  w->birth = birth_depth_gaussian;
  w->death = death_depth_gaussian;
  w->destroy = destroy_depth_gaussian;

  return w;
}

static int 
birth_depth_gaussian(void *user,
	       int i,
	       int j,
	       int k,
	       int level,
	       int maxlevel,
	       double parent_coeff,
	       double *coeff,
	       double *prob)
{
  struct bd_depth_gaussian *s;
  int d;

  s = (struct bd_depth_gaussian *)user;

  d = level;
  if (d >= s->ndepths) {
    d = s->ndepths - 1;
  }

  *coeff = gsl_ran_gaussian_ziggurat(s->rng, s->sigma[d]);
  *prob = gsl_ran_gaussian_pdf(*coeff, s->sigma[d]);

  return 0;
}

static int 
death_depth_gaussian(void *user,
	       int i,
	       int j,
	       int k,
	       int level,
	       int maxlevel,
	       double parent_coeff,
	       double coeff,
	       double *prob)
{
  struct bd_depth_gaussian *s;
  int d;

  s = (struct bd_depth_gaussian *)user;

  d = level;
  if (d >= s->ndepths) {
    d = s->ndepths - 1;
  }

  *prob = gsl_ran_gaussian_pdf(coeff, s->sigma[d]);
  return 0;
}

static int destroy_depth_gaussian(void *user)
{
  struct bd_depth_gaussian *s;

  s = (struct bd_depth_gaussian *)user;
  
  gsl_rng_free(s->rng);
  free(s->sigma);
  free(s);
  
  return 0;
}


