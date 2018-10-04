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

#include "wavetreepp.h"

#include "wavetree_value_proposal.h"
#include "wavetree_birth_proposal.h"

#include "slog.h"

wavetree_pp_t *
wavetree_pp_load(const char *filename, unsigned long int seed, coefficient_histogram_t *histogram)
{
  FILE *fp;

  wavetree_prior_t *prior;
  wavetree_bd_t *bd;
  wavetree_value_t *value;

  wavetree_pp_t *pp;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    return NULL;
  }

  prior = wavetree_pp_load_prior(fp, seed);
  if (prior == NULL) {
    return NULL;
  }

  bd = wavetree_pp_load_bd(fp, seed, prior);
  if (bd == NULL) {
    return NULL;
  }

  value = wavetree_pp_load_value(fp, seed, histogram);
  if (value == NULL) {
    return NULL;
  }

  pp = malloc(sizeof(wavetree_pp_t));
  if (pp == NULL) {
    return NULL;
  }

  pp->prior = prior;
  pp->bd = bd;
  pp->value = value;

  fclose(fp);

  return pp;
}

wavetree_prior_t *
wavetree_pp_load_prior(FILE *fp, unsigned long int seed)
{
  char buffer[256];
  double vmin[16];
  double vmax[16];
  int ndepths;
  double beta;
  int i;

  wavetree_prior_t *r;

  if (fp == NULL) {
    return NULL;
  }

  if (fscanf(fp, "%s\n", buffer) != 1) {
    ERROR("failed to read prior name");
    return NULL;
  }

  if (strcmp(buffer, "uniform") == 0) {

    if (fscanf(fp, "%lf %lf\n", &(vmin[0]), &(vmax[0])) != 2) {
      ERROR("failed to read uniform prior bounds");
      return NULL;
    }

    r = wavetree_prior_create_globally_uniform(vmin[0], vmax[0], seed);
    if (r == NULL) {
      ERROR("failed to create Globally Uniform prior");
      return NULL;
    }

  } else if (strcmp(buffer, "laplace") == 0) {

    if (fscanf(fp, "%lf\n", &beta) != 1) {
      ERROR("faield to read laplace prior width");
      return NULL;
    }

    r = wavetree_prior_create_globally_laplacian(beta, seed);
    if (r == NULL) {
      ERROR("failed to create Globally Laplacian prior");
      return NULL;
    }
    
  } else if (strcmp(buffer, "depthuniform") == 0) {

    if (fscanf(fp, "%d\n", &ndepths) != 1) {
      ERROR("failed to read n depths for Depth Uniform");
      return NULL;
    }

    if (ndepths < 0 || ndepths > 16) {
      ERROR("invalid ndepths %d", ndepths);
      return NULL;
    }

    for (i = 0; i < ndepths; i ++) {
      if (fscanf(fp, "%lf %lf\n", &(vmin[i]), &(vmax[i])) != 2) {
	ERROR("failed to read std dev %d", i);
	return NULL;
      }
    }	
    
    r = wavetree_prior_create_depth_uniform(ndepths,
					    vmin,
					    vmax,
					    seed);
    if (r == NULL) {
      ERROR("failed to create Depth Uniform prior");
      return NULL;
    }

  /* } else if (strcmp(buffer, "depthlaplacian") == 0) { */

  /*   if (fscanf(fp, "%d\n", &ndepths) != 1) { */
  /*     ERROR("failed to read n depths for Depth Laplacian"); */
  /*     return NULL; */
  /*   } */

  /*   if (ndepths < 0 || ndepths > 16) { */
  /*     ERROR("invalid ndepths %d", ndepths); */
  /*     return NULL; */
  /*   } */

  /*   for (i = 0; i < ndepths; i ++) { */
  /*     if (fscanf(fp, "%lf\n", &(vmin[i])) != 1) { */
  /* 	ERROR("failed to read std dev %d", i); */
  /* 	return NULL; */
  /*     } */
  /*   } */

  /*   r = wavetree_prior_create_depth_laplacian(ndepths, vmin, seed); */

  /*   if (r == NULL) { */
  /*     ERROR("failed to create Depth Uniform prior"); */
  /*     return NULL; */
  /*   } */
    
  } else if (strcmp(buffer, "depthgeneralisedgaussian") == 0) {

    if (fscanf(fp, "%lf\n", &beta) != 1) {
      ERROR("failed to read prior beta");
      return NULL;
    }

    if (fscanf(fp, "%d\n", &ndepths) != 1) {
      ERROR("failed to read n depths for Generalise Gaussian");
      return NULL;
    }

    if (ndepths < 0 || ndepths > 16) {
      ERROR("invalid ndepths %d", ndepths);
      return NULL;
    }

    for (i = 0; i < ndepths; i ++) {
      if (fscanf(fp, "%lf\n", &(vmin[i])) != 1) {
	ERROR("failed to read std dev %d", i);
	return NULL;
      }
    }	

    r = wavetree_prior_create_depth_generalised_gaussian(ndepths, vmin, beta, seed);
    if (r == NULL) {
      ERROR("failed to create generalised gaussian");
      return NULL;
    }
  } else {
    ERROR("invalid bd proprosal name: %s", buffer);
    return NULL;
  }
 

  return r;
}

wavetree_bd_t *
wavetree_pp_load_bd(FILE *fp, unsigned long int seed, wavetree_prior_t *prior)
{
  char buffer[256];
  double std[16];
  int ndepths;
  int i;

  wavetree_bd_t *r;

  if (fp == NULL) {
    return NULL;
  }

  if (fscanf(fp, "%s\n", buffer) != 1) {
    ERROR("failed to read bd proposal name");
    return NULL;
  }

  if (strcmp(buffer, "priorbirth") == 0) {

    r = wavetree_birth_create_birth_from_prior(prior);
    if (r == NULL) {
      ERROR("failed to create birth from prior");
      return NULL;
    }

  } else if (strcmp(buffer, "gaussianbirth") == 0) {

    if (fscanf(fp, "%lf\n", &std[0]) != 1) {
      ERROR("failed to read std deviation");
      return NULL;
    }

    r = wavetree_birth_create_gaussian(prior,
				       std[0],
				       seed);
    if (r == NULL) {
      ERROR("failed to create birth from gaussian");
      return NULL;
    }

  } else if (strcmp(buffer, "depthgaussianbirth") == 0) {

    if (fscanf(fp, "%d\n", &ndepths) != 1) {
      ERROR("failed to read n depths for depth gaussian");
      return NULL;
    }

    if (ndepths < 0 || ndepths > 16) {
      ERROR("invalid ndepths %d", ndepths);
      return NULL;
    }

    for (i = 0; i < ndepths; i ++) {
      if (fscanf(fp, "%lf\n", &(std[i])) != 1) {
	ERROR("failed to read std dev %d", i);
	return NULL;
      }
    }	

    r = wavetree_birth_create_depth_gaussian(prior,
					     ndepths,
					     std,
					     seed);
    if (r == NULL) {
      ERROR("failed to create depth gaussian");
      return NULL;
    }

  } else {
    ERROR("invalid bd proprosal name: %s", buffer);
    return NULL;
  }

  return r;
  
}

wavetree_value_t *
wavetree_pp_load_value(FILE *fp, unsigned long int seed, coefficient_histogram_t *histogram)
{
  char buffer[256];
  double std[16];
  int ndepths;
  int i;
  double epsilon;
  double tau;
  double A;
  double s;
  int threshold;

  wavetree_value_t *r;

  if (fp == NULL) {
    return NULL;
  }

  if (fscanf(fp, "%s\n", buffer) != 1) {
    ERROR("failed to read value proposal name");
    return NULL;
  }

  if (strcmp(buffer, "gaussianperturb") == 0) {

    if (fscanf(fp, "%lf\n", &std[0]) != 1) {
      ERROR("failed to read std deviation");
      return NULL;
    }

    r = wavetree_value_create_global_gaussian(std[0],
					      seed);
    if (r == NULL) {
      ERROR("failed to create birth from gaussian");
      return NULL;
    }

  } else if (strcmp(buffer, "depthgaussianperturb") == 0) {

    if (fscanf(fp, "%d\n", &ndepths) != 1) {
      ERROR("failed to read n depths for depth gaussian");
      return NULL;
    }

    if (ndepths < 0 || ndepths > 16) {
      ERROR("invalid ndepths %d", ndepths);
      return NULL;
    }

    for (i = 0; i < ndepths; i ++) {
      if (fscanf(fp, "%lf\n", &(std[i])) != 1) {
	ERROR("failed to read std dev %d", i);
	return NULL;
      }
    }	

    r = wavetree_value_create_depth_gaussian(ndepths,
					     std,
					     seed);
    if (r == NULL) {
      ERROR("failed to create depth gaussian");
      return NULL;
    }

  } else if (strcmp(buffer, "depthgaussianscamperturb") == 0) {

    if (fscanf(fp, "%d\n", &ndepths) != 1) {
      ERROR("failed to read n depths for depth gaussian scam");
      return NULL;
    }

    if (ndepths < 0 || ndepths > 16) {
      ERROR("invalid ndepths %d", ndepths);
      return NULL;
    }

    for (i = 0; i < ndepths; i ++) {
      if (fscanf(fp, "%lf\n", &(std[i])) != 1) {
	ERROR("failed to read std dev %d", i);
	return NULL;
      }
    }	

    if (fscanf(fp, "%lf\n", &epsilon) != 1) {
      ERROR("failed to read epsilon");
      return NULL;
    }

    if (fscanf(fp, "%lf\n", &s) != 1) {
      ERROR("failed to read s (scam)");
      return NULL;
    }

    if (fscanf(fp, "%d\n", &threshold) != 1) {
      ERROR("failed to read threshold (scam)");
      return NULL;
    }
     
    r = wavetree_value_create_depth_gaussian_scam(ndepths,
						  std,
						  epsilon,
						  s,
						  threshold,
						  seed,
						  histogram);

    if (r == NULL) {
      ERROR("failed to create depth gaussian scam");
      return NULL;
    }

  } else if (strcmp(buffer, "gaussianamperturb") == 0) {

    if (fscanf(fp, "%lf\n", &(std[0])) != 1) {
      ERROR("failed to read std0 for gaussian am");
      return NULL;
    }
    if (fscanf(fp, "%lf\n", &epsilon) != 1) {
      ERROR("failed to read epsilon for gaussian am");
      return NULL;
    }
    if (fscanf(fp, "%lf\n", &A) != 1) {
      ERROR("failed to read A for gaussian am");
      return NULL;
    }
    if (fscanf(fp, "%lf\n", &tau) != 1) {
      ERROR("failed to read tau for gaussian am");
      return NULL;
    }

    r = wavetree_value_create_gaussian_am(std[0],
					  epsilon,
					  A,
					  tau,
					  seed,
					  histogram);
    if (r == NULL) {
      ERROR("failed to create gaussian am");
      return NULL;
    }

  } else if (strcmp(buffer, "depthcauchyperturb") == 0) {

    if (fscanf(fp, "%d\n", &ndepths) != 1) {
      ERROR("failed to read n depths for depth cauchy");
      return NULL;
    }

    if (ndepths < 0 || ndepths > 16) {
      ERROR("invalid ndepths %d", ndepths);
      return NULL;
    }

    for (i = 0; i < ndepths; i ++) {
      if (fscanf(fp, "%lf\n", &(std[i])) != 1) {
	ERROR("failed to read std dev %d", i);
	return NULL;
      }
    }	

    r = wavetree_value_create_depth_cauchy(ndepths,
					     std,
					     seed);
    if (r == NULL) {
      ERROR("failed to create depth cauchy");
      return NULL;
    }

  } else if (strcmp(buffer, "cauchyamperturb") == 0) {

    if (fscanf(fp, "%lf\n", &(std[0])) != 1) {
      ERROR("failed to read std0 for cauchy am");
      return NULL;
    }
    
    if (fscanf(fp, "%lf\n", &epsilon) != 1) {
      ERROR("failed to read epsilon for cauchy am");
      return NULL;
    }
    if (fscanf(fp, "%lf\n", &A) != 1) {
      ERROR("failed to read A for cauchy am");
      return NULL;
    }
    if (fscanf(fp, "%lf\n", &tau) != 1) {
      ERROR("failed to read tau for cauchy am");
      return NULL;
    }

    r = wavetree_value_create_cauchy_am(std[0],
					epsilon,
					A,
					tau,
					seed,
					histogram);
    if (r == NULL) {
      ERROR("failed to create cauchy am");
      return NULL;
    }

  } else {
    ERROR("invalid value proprosal name: %s", buffer);
    return NULL;
  }

  return r;
  
}


/*
 * Helper functions
 */
int
wavetree_pp_value_init(wavetree_pp_t *w)
{
  if (w->value->init(w->value->user) < 0) {
    return -1;
  }

  return 0;
}

int
wavetree_pp_propose_value2d(wavetree_pp_t *w, 
			    int i, int j, 
			    int level, int maxlevel, double parent_coeff,
			    double temperature,
			    double *coeff,
			    double *prior_ratio)
{
  if (w->value->perturb(w->prior, 
			w->value->user,
			i,
			j,
			0,
			level,
			maxlevel,
			parent_coeff,
			temperature,
			coeff,
			prior_ratio) < 0) {
    return -1;
  }
  
  return 0;
}

int 
wavetree_pp_value_error_count(wavetree_pp_t *w)
{
  int prior_errors;
  if (w->value->errors(w->value->user, &prior_errors) < 0) {
    return -1;
  }

  return prior_errors;
}

double
wavetree_pp_prior_probability2d(wavetree_pp_t *w,
				int i,
				int j,
				int depth,
				int maxdepth,
				double parent_coeff,
				double coeff)
{
  return w->prior->prob(w->prior->user,
			i,
			j,
			0,
			depth,
			maxdepth,
			parent_coeff,
			coeff);
}

int
wavetree_pp_prior_range2d(wavetree_pp_t *w,
			  int i,
			  int j,
			  int depth,
			  int maxdepth,
			  double parent_coeff,
			  double *vmin,
			  double *vmax)
{
  return w->prior->range(w->prior->user,
			 i,
			 j,
			 0,
			 depth,
			 maxdepth,
			 parent_coeff,
			 vmin,
			 vmax);
}

int
wavetree_pp_birth2d(wavetree_pp_t *w,
		    int i,
		    int j,
		    int depth,
		    int maxdepth,
		    double parent_coeff,
		    double *coeff,
		    double *prob,
		    int *valid)
{
  if (w->bd->birth(w->bd->user,
		   i,
		   j,
		   0,
		   depth,
		   maxdepth,
		   parent_coeff,
		   coeff,
		   prob) < 0) {
    ERROR("failed to do birth proposal");
    return -1;
  }

  *valid = w->prior->valid(w->prior->user,
			   i,
			   j,
			   0,
			   depth,
			   maxdepth,
			   parent_coeff,
			   *coeff);

  return 0;
}

int
wavetree_pp_death2d(wavetree_pp_t *w,
		    int i,
		    int j,
		    int depth,
		    int maxdepth,
		    double parent_coeff,
		    double coeff,
		    double *prob)
{
  if (w->bd->death(w->bd->user,
		   i, 
		   j,
		   0,
		   depth,
		   maxdepth,
		   parent_coeff,
		   coeff,
		   prob) < 0) {
    ERROR("failed to do death proposal");
    return -1;
  }

  return 0;
}

int
wavetree_pp_propose_value3d(wavetree_pp_t *w, 
			    int i, int j, int k,
			    int level, int maxlevel, double parent_coeff,
			    double temperature,
			    double *coeff,
			    double *prior_ratio)
{
  if (w->value->perturb(w->prior, 
			w->value->user,
			i,
			j,
			k,
			level,
			maxlevel,
			parent_coeff,
			temperature,
			coeff,
			prior_ratio) < 0) {
    return -1;
  }
  
  return 0;
}

double
wavetree_pp_prior_probability3d(wavetree_pp_t *w,
				int i,
				int j,
				int k,
				int depth,
				int maxdepth,
				double parent_coeff,
				double coeff)
{
  return w->prior->prob(w->prior->user,
			i,
			j,
			k,
			depth,
			maxdepth,
			parent_coeff,
			coeff);
}

int
wavetree_pp_prior_range3d(wavetree_pp_t *w,
			  int i,
			  int j,
			  int k,
			  int depth,
			  int maxdepth,
			  double parent_coeff,
			  double *vmin,
			  double *vmax)
{
  return w->prior->range(w->prior->user,
			 i,
			 j,
			 k,
			 depth,
			 maxdepth,
			 parent_coeff,
			 vmin,
			 vmax);
}


int
wavetree_pp_birth3d(wavetree_pp_t *w,
		    int i,
		    int j,
		    int k,
		    int depth,
		    int maxdepth,
		    double parent_coeff,
		    double *coeff,
		    double *prob,
		    int *valid)
{
  if (w->bd->birth(w->bd->user,
		   i,
		   j,
		   k,
		   depth,
		   maxdepth,
		   parent_coeff,
		   coeff,
		   prob) < 0) {
    ERROR("failed to do birth proposal");
    return -1;
  }

  *valid = w->prior->valid(w->prior->user,
			   i,
			   j,
			   k,
			   depth,
			   maxdepth,
			   parent_coeff,
			   *coeff);

  return 0;
}


int
wavetree_pp_death3d(wavetree_pp_t *w,
		    int i,
		    int j,
		    int k,
		    int depth,
		    int maxdepth,
		    double parent_coeff,
		    double coeff,
		    double *prob)
{
  if (w->bd->death(w->bd->user,
		   i, 
		   j,
		   k,
		   depth,
		   maxdepth,
		   parent_coeff,
		   coeff,
		   prob) < 0) {
    ERROR("failed to do death proposal");
    return -1;
  }

  return 0;
}

int
wavetree_pp_prior_sample3d(wavetree_pp_t *w,
			   int i,
			   int j,
			   int k,
			   int depth,
			   int maxdepth,
			   double parent_coeff,
			   double *coeff)
{
  if (w->prior->sample(w->prior->user,
		       i, j, k,
		       depth,
		       maxdepth,
		       parent_coeff,
		       coeff) < 0) {
    ERROR("failed to sample prior");
    return -1;
  }

  return 0;
}

int
wavetree_pp_setscale(wavetree_pp_t *w,
		     double newscale,
		     double *oldscale)
{
  return w->prior->setscale(w->prior->user,
			    newscale,
			    oldscale);
}

/*
 * General cleanup function
 */
void
wavetree_pp_destroy(wavetree_pp_t *p)
{
  if (p != NULL) {

    p->prior->destroy(p->prior->user);
    free(p->prior);

    p->bd->destroy(p->bd->user);
    free(p->bd);

    p->value->destroy(p->value->user);
    free(p->value);

    free(p);
  }
}

/*
 * Creation functions
 */
wavetree_pp_t *
wavetree_pp_create_globally_uniform(double vmin,
				    double vmax,
				    double vstd,
				    int kmax,
				    unsigned long int seed)
{
  wavetree_pp_t *w;

  w = malloc(sizeof(wavetree_pp_t));
  if (w == NULL) {
    return NULL;
  }

  w->prior = wavetree_prior_create_globally_uniform(vmin, vmax, seed);
  if (w->prior == NULL) {
    ERROR("failed to create prior");
    return NULL;
  }

  w->bd = wavetree_birth_create_birth_from_prior(w->prior);
  if (w->bd == NULL) {
    ERROR("failed to create birth proposal");
    return NULL;
  }

  w->value = wavetree_value_create_global_gaussian(vstd, seed);
  if (w->value == NULL) {
    ERROR("failed to create value proposal");
    return NULL;
  }

  return w;
}

/*
 * Globally uniform prior on coefficients, gaussian birth proposal, gaussian move proposal.
 */
wavetree_pp_t *
wavetree_pp_create_globally_uniform_with_gaussian_proposal(double vmin,
							   double vmax,
							   double vstd,
							   double bd_std,
							   double move_std,
							   int kmax,
							   unsigned long int seed)
{
  wavetree_pp_t *w;

  w = malloc(sizeof(wavetree_pp_t));
  if (w == NULL) {
    return NULL;
  }

  w->prior = wavetree_prior_create_globally_uniform(vmin, vmax, seed);
  if (w->prior == NULL) {
    ERROR("failed to create prior");
    return NULL;
  }

  w->bd = wavetree_birth_create_gaussian(w->prior, bd_std, seed);
  if (w->bd == NULL) {
    ERROR("failed to create birth proposal");
    NULL;
  }

  w->value = wavetree_value_create_global_gaussian(vstd, seed);
  if (w->value == NULL) {
    ERROR("failed to create value proposal");
    return NULL;
  }

  return w;
}

wavetree_pp_t *
wavetree_pp_create_depth_dependent_uniform(double *vmin,
					   double *vmax,
					   double *vstd,
					   int ndepths,
					   int kmax,
					   double move_std,
					   unsigned long int seed)
{
  wavetree_pp_t *w;

  w = malloc(sizeof(wavetree_pp_t));
  if (w == NULL) {
    return NULL;
  }

  w->prior = wavetree_prior_create_depth_uniform(ndepths, vmin, vmax, seed);
  if (w->prior == NULL) {
    ERROR("failed to create prior");
    return NULL;
  }

  w->bd = wavetree_birth_create_birth_from_prior(w->prior);
  if (w->bd == NULL) {
    ERROR("failed to create birth proposal");
    return NULL;
  }

  w->value = wavetree_value_create_depth_gaussian(ndepths, vstd, seed);
  if (w->value == NULL) {
    ERROR("failed to create value proposal");
    return NULL;
  }

  return w;
}

wavetree_pp_t *
wavetree_pp_create_depth_dependent_uniform_with_gaussian_proposal(double *vmin,
								  double *vmax,
								  double *vstd,
								  double *bd_std,
								  int ndepths,
								  int kmax,
								  double move_std,
								  unsigned long int seed)
{
  wavetree_pp_t *w;

  w = malloc(sizeof(wavetree_pp_t));
  if (w == NULL) {
    return NULL;
  }

  w->prior = wavetree_prior_create_depth_uniform(ndepths, vmin, vmax, seed);
  if (w->prior == NULL) {
    ERROR("failed to create prior");
    return NULL;
  }

  w->bd = wavetree_birth_create_depth_gaussian(w->prior, ndepths, bd_std, seed);
  if (w->bd == NULL) {
    ERROR("failed to create birth proposal");
    NULL;
  }

  w->value = wavetree_value_create_depth_gaussian(ndepths, vstd, seed);
  if (w->value == NULL) {
    ERROR("failed to create value proposal");
    return NULL;
  }

  return w;
}

wavetree_pp_t *
wavetree_pp_create_depth_dependent_generalised_gaussian(double *va,
							double beta,
							double *vstd,
							int ndepths,
							int kmax,
							double move_std,
							unsigned long int seed)
{
  wavetree_pp_t *w;

  w = malloc(sizeof(wavetree_pp_t));
  if (w == NULL) {
    return NULL;
  }

  w->prior = wavetree_prior_create_depth_generalised_gaussian(ndepths, va, beta, seed);
  if (w->prior == NULL) {
    ERROR("failed to create prior");
    return NULL;
  }

  w->bd = wavetree_birth_create_birth_from_prior(w->prior);
  if (w->bd == NULL) {
    ERROR("failed to create birth proposal");
    return NULL;
  }

  w->value = wavetree_value_create_depth_gaussian(ndepths, vstd, seed);
  if (w->value == NULL) {
    ERROR("failed to create value proposal");
    return NULL;
  }

  return w;
}

/*
 * Depth dependent generalised gaussian, gaussian birth
 */
wavetree_pp_t *
wavetree_pp_create_depth_dependent_generalised_gaussian_with_gaussian_proposal(double *va,
									       double beta,
									       double *vstd,
									       double *bd_std,
									       int ndepths,
									       int kmax,
									       double move_std,
									       unsigned long int seed)
{
  wavetree_pp_t *w;

  w = malloc(sizeof(wavetree_pp_t));
  if (w == NULL) {
    return NULL;
  }

  w->prior = wavetree_prior_create_depth_generalised_gaussian(ndepths, va, beta, seed);
  if (w->prior == NULL) {
    ERROR("failed to create prior");
    return NULL;
  }

  w->bd = wavetree_birth_create_depth_gaussian(w->prior, ndepths, bd_std, seed);
  if (w->bd == NULL) {
    ERROR("failed to create birth proposal");
    NULL;
  }

  w->value = wavetree_value_create_depth_gaussian(ndepths, vstd, seed);
  if (w->value == NULL) {
    ERROR("failed to create value proposal");
    return NULL;
  }

  return w;
}
