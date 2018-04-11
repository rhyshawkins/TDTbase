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
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "wavetree_value_proposal.h"

/*
 * Depth dependent Gaussian proposal
 */

struct depth_gaussian_am {
  coefficient_histogram_t *histogram;
  
  gsl_rng *rng;

  int errork;

  int ndepths;
  double *std;

  double epsilon;
  double s;
  int threshold;
};

static int depth_gaussian_am_init(void *user);

static int depth_gaussian_am_perturb(wavetree_prior_t *prior,
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

static int depth_gaussian_am_errors(void *user,
				  int *prior_errors);

static int depth_gaussian_am_destroy(void *user);

wavetree_value_t *
wavetree_value_create_depth_gaussian_am(int ndepths,
					double *std,
					double epsilon,
					double s,
					int threshold,
					unsigned long int seed,
					coefficient_histogram_t *histogram)
{
  wavetree_value_t *w;
  struct depth_gaussian_am *g;
  int i;

  if (ndepths <= 0) {
    fprintf(stderr, "wavetree_value_create_depth_gaussian_am: invalid parameters\n");
    return NULL;
  }

  if (histogram == NULL) {
    fprintf(stderr, "wavetree_value_create_depth_gaussian_am: NULL histogram\n");
    return NULL;
  }

  w = malloc(sizeof(wavetree_value_t));
  if (w == NULL) {
    fprintf(stderr, 
	    "wavetree_value_create_depth_gaussian_am: "
	    "failed to allocate structure\n");
    return NULL;
  }

  w->init = depth_gaussian_am_init;
  w->perturb = depth_gaussian_am_perturb;
  w->errors = depth_gaussian_am_errors;
  w->destroy = depth_gaussian_am_destroy;

  g = malloc(sizeof(struct depth_gaussian_am));
  if (g == NULL) {
    fprintf(stderr, 
	    "wavetree_value_create_depth_gaussian_am: "
	    "failed to allocate gaussian structure\n");
    return NULL;
  }

  g->rng = gsl_rng_alloc(gsl_rng_taus);
  if (g->rng == NULL) {
    fprintf(stderr, 
	    "wavetree_value_create_depth_gaussian_am: "
	    "failed to create rng\n");
    return NULL;
  }

  gsl_rng_set(g->rng, seed);
  
  g->errork = 0;

  g->std = malloc(sizeof(double) * ndepths);
  if (g->std == NULL) {
    fprintf(stderr, 
	    "wavetree_value_greate_depth_gaussian_am: "
	    "failed to allocate restore array\n");
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

static int depth_gaussian_am_init(void *user)
{
  struct depth_gaussian_am *g = (struct depth_gaussian_am *)user;

  g->errork = 0;

  return 0;
}

static int depth_gaussian_am_perturb(wavetree_prior_t *prior,
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
  struct depth_gaussian_am *g = (struct depth_gaussian_am *)user;
  int d;
  double new_coeff;
  double p;
  double pp;

  double mean;
  double std;
  int n;

  d = level;
  if (d >= g->ndepths) {
    d = g->ndepths - 1;
  }

  n = coefficient_histogram_get_coefficient_mean_std(g->histogram,
						     i,
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

static int depth_gaussian_am_errors(void *user,
				 int *prior_errors)
{
  struct depth_gaussian_am *g = (struct depth_gaussian_am *)user;

  *prior_errors = g->errork;
  
  return 0;
}

static int depth_gaussian_am_destroy(void *user)
{
  struct depth_gaussian_am *g = (struct depth_gaussian_am *)user;
  
  gsl_rng_free(g->rng);

  free(g->std);
  free(g);

  return 0;
}
