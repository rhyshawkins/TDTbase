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

/*
 * Depth uniform
 */

struct depth_uniform {
  gsl_rng *rng;

  int ndepths;
  double *vmin;
  double *vmax;
};

static int
range_depth_uniform(void *user,
		    int i,
		    int j,
		    int k,
		    int level,
		    int maxlevel,
		    double parent_coeff,
		    double *vmin,
		    double *vmax);
		    
static int 
sample_depth_uniform(void *user,
		     int i,
		     int j,
		     int k,
		     int level,
		     int maxlevel,
		     double parent_coeff,
		     double *coeff);

static double
prob_depth_uniform(void *user,
		   int i,
		   int j,
		   int k,
		   int level,
		   int maxlevel,
		   double parent_coeff,
		   double coeff);

static int 
valid_depth_uniform(void *user,
		    int i,
		    int j,
		    int k,
		    int level,
		    int maxlevel,
		    double parent_coeff,
		    double coeff);

static int
setscale_depth_uniform(void *user,
		       double newscale,
		       double *oldscale);


static int
destroy_depth_uniform(void *user);


wavetree_prior_t*
wavetree_prior_create_depth_uniform(int ndepths,
				    double *vmin,
				    double *vmax,
				    unsigned long int seed)
{
  wavetree_prior_t *w;
  struct depth_uniform *s;
  int i;
  
  w = malloc(sizeof(wavetree_prior_t));
  if (w == NULL) {
    return NULL;
  }

  s = malloc(sizeof(struct depth_uniform));
  if (s == NULL) {
    return NULL;
  }

  s->rng = gsl_rng_alloc(gsl_rng_taus);
  if (s->rng == NULL) {
    return NULL;
  }

  gsl_rng_set(s->rng, seed);

  s->ndepths = ndepths;

  s->vmin = malloc(sizeof(double) * ndepths);
  if (s->vmin == NULL) {
    return NULL;
  }

  s->vmax = malloc(sizeof(double) * ndepths);
  if (s->vmax == NULL) {
    return NULL;
  }

  for (i = 0; i < ndepths; i ++) {
    s->vmin[i] = vmin[i];
    s->vmax[i] = vmax[i];
  }
  
  w->user = s;
  w->range = range_depth_uniform;
  w->sample = sample_depth_uniform;
  w->prob = prob_depth_uniform;
  w->valid = valid_depth_uniform;
  w->setscale = setscale_depth_uniform;
  w->destroy = destroy_depth_uniform;

  return w;
}

static int
range_depth_uniform(void *user,
		    int i,
		    int j,
		    int k,
		    int level,
		    int maxlevel,
		    double parent_coeff,
		    double *vmin,
		    double *vmax)
{
  struct depth_uniform *s; 
  int d;

  s = (struct depth_uniform *)user;

  d = level;
  if (d >= s->ndepths) {
    d = s->ndepths - 1;
  }
  
  *vmin = s->vmin[d];
  *vmax = s->vmax[d];

  return 0;
}

static int 
sample_depth_uniform(void *user,
		     int i,
		     int j,
		     int k,
		     int level,
		     int maxlevel,
		     double parent_coeff,
		     double *coeff)
{
  struct depth_uniform *s; 
  int d;

  s = (struct depth_uniform *)user;

  d = level;
  if (d >= s->ndepths) {
    d = s->ndepths - 1;
  }

  *coeff = s->vmin[d] + (s->vmax[d] - s->vmin[d])*gsl_rng_uniform(s->rng);
  
  return 0;
}

static double
prob_depth_uniform(void *user,
		   int i,
		   int j,
		   int k,
		   int level,
		   int maxlevel,
		   double parent_coeff,
		   double coeff)
{
  struct depth_uniform *s; 
  int d;

  s = (struct depth_uniform *)user;

  d = level;
  if (d >= s->ndepths) {
    d = s->ndepths - 1;
  }

  if (coeff >= s->vmin[d] && coeff <= s->vmax[d]) {
    return 1.0/(s->vmax[d] - s->vmin[d]);
  }
  
  return 0.0;
}


static int 
valid_depth_uniform(void *user,
		    int i,
		    int j,
		    int k,
		    int level,
		    int maxlevel,
		    double parent_coeff,
		    double coeff)
{
  struct depth_uniform *s; 
  int d;

  s = (struct depth_uniform *)user;

  d = level;
  if (d >= s->ndepths) {
    d = s->ndepths - 1;
  }

  if (coeff >= s->vmin[d] && coeff <= s->vmax[d]) {
    return -1;
  } else {
    return 0;
  }
}

static int
setscale_depth_uniform(void *user,
		       double newscale,
		       double *oldscale)
{
  return -1;
}

static int
destroy_depth_uniform(void *user)
{
  struct depth_uniform *s;

  s = (struct depth_uniform *)user;
  gsl_rng_free(s->rng);
  free(s->vmin);
  free(s->vmax);
  
  free(s);

  return 0;
}

