//
//    Wavetree Library : A library for performed trans-dimensional tree inversion,
//    See
//
//      R Hawkins and M Sambridge, "Geophysical imaging using trans-dimensional trees",
//      Geophysical Journal International, 2015, 203:2, 972 - 1000,
//      https://doi.org/10.1093/gji/ggv326
//    
//    Copyright (C) 2014 - 2018 Rhys Hawkins
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "wavetree_prior.h"

struct globally_uniform {
  gsl_rng *rng;

  double vmin;
  double vmax;
};

static int
range_globally_uniform(void *user,
		       int i,
		       int j,
		       int k,
		       int level,
		       int maxlevel,
		       double parent_coeff,
		       double *vmin,
		       double *vmax);

static int 
sample_globally_uniform(void *user,
			int i,
			int j,
			int k,
			int level,
			int maxlevel,
			double parent_coeff,
			double *coeff);

static double
prob_globally_uniform(void *user,
		      int i,
		      int j,
		      int k,
		      int level,
		      int maxlevel,
		      double parent_coeff,
		      double coeff);

static int 
valid_globally_uniform(void *user,
		       int i,
		       int j,
		       int k,
		       int level,
		       int maxlevel,
		       double parent_coeff,
		       double coeff);

static int
setscale_globally_uniform(void *user,
			  double scale,
			  double *oldcale);

static int
destroy_globally_uniform(void *user);

wavetree_prior_t *
wavetree_prior_create_globally_uniform(double vmin,
				       double vmax,
				       unsigned long int seed)
{
  wavetree_prior_t *w;
  struct globally_uniform *s;

  w = malloc(sizeof(wavetree_prior_t));
  if (w == NULL) {
    return NULL;
  }

  s = malloc(sizeof(struct globally_uniform));
  if (s == NULL) {
    return NULL;
  }

  s->rng = gsl_rng_alloc(gsl_rng_taus);
  if (s->rng == NULL) {
    return NULL;
  }

  gsl_rng_set(s->rng, seed);

  s->vmin = vmin;
  s->vmax = vmax;

  w->user = s;
  w->range = range_globally_uniform;
  w->sample = sample_globally_uniform;
  w->prob = prob_globally_uniform;
  w->valid = valid_globally_uniform;
  w->setscale = setscale_globally_uniform;
  w->destroy = destroy_globally_uniform;

  return w;
}

void
wavetree_prior_free_globally_uniform(wavetree_prior_t *w)
{
  struct globally_uniform *s;

  if (w != NULL) {
    s = (struct globally_uniform *)w->user;
    gsl_rng_free(s->rng);
    free(s);
    free(w);
  }
}


static int
range_globally_uniform(void *user,
		       int i,
		       int j,
		       int k,
		       int level,
		       int maxlevel,
		       double parent_coeff,
		       double *vmin,
		       double *vmax)
{
  struct globally_uniform *s; 

  s = (struct globally_uniform *)user;

  *vmin = s->vmin;
  *vmax = s->vmax;

  return 0;
}

static int 
sample_globally_uniform(void *user,
			int i,
			int j,
			int k,
			int level,
			int maxlevel,
			double parent_coeff,
			double *coeff)
{
  struct globally_uniform *s; 

  s = (struct globally_uniform *)user;

  *coeff = s->vmin + (s->vmax - s->vmin)*gsl_rng_uniform(s->rng);
  
  return 0;
}

static double
prob_globally_uniform(void *user,
		      int i,
		      int j,
		      int k,
		      int level,
		      int maxlevel,
		      double parent_coeff,
		      double coeff)
{
  struct globally_uniform *s; 

  s = (struct globally_uniform *)user;

  if (coeff >= s->vmin && coeff <= s->vmax) {
    return 1.0/(s->vmax - s->vmin);
  }

  return 0.0;
}

static int 
valid_globally_uniform(void *user,
		       int i,
		       int j,
		       int k,
		       int level,
		       int maxlevel,
		       double parent_coeff,
		       double coeff)
{
  struct globally_uniform *s; 

  s = (struct globally_uniform *)user;

  if (coeff >= s->vmin && coeff <= s->vmax) {
    return -1;
  }
  
  return 0;
}

static int
setscale_globally_uniform(void *user,
			  double scale,
			  double *oldcale)
{
  return -1;
}

static int
destroy_globally_uniform(void *user)
{
  struct globally_uniform *s; 

  s = (struct globally_uniform *)user;

  gsl_rng_free(s->rng);

  free(s);

  return 0;
}
