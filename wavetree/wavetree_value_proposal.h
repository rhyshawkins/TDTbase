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

#ifndef wavetree_value_proposal_h
#define wavetree_value_proposal_h

#include "wavetree_prior.h"
#include "coefficient_histogram.h"

typedef int (*wavetree_pp_value_initialize_t)(void *user);

typedef int (*wavetree_pp_value_proposal_t)(wavetree_prior_t *prior,
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

typedef int (*wavetree_pp_value_errors_t)(void *user,
					  int *prior_errors);
					    
typedef int (*wavetree_pp_value_restore_t)(void *user,
					   double *coeff);

typedef int (*wavetree_pp_save_t)(void *user,
				  const char *filename);

typedef int (*wavetree_pp_restore_t)(void *user,
				     const char *filename);

struct _wavetree_pp_value {
  void *user;
  
  wavetree_pp_value_initialize_t init;  /* Resets error count */
  wavetree_pp_value_proposal_t perturb; /* Modify coefficients */
  wavetree_pp_value_errors_t errors;    /* Returns count of errors */
  wavetree_pp_destroy_t destroy;

  wavetree_pp_save_t save;              /* Save state to file */
  wavetree_pp_restore_t restore;        /* Restore state from file */
};

typedef struct _wavetree_pp_value wavetree_value_t;

/*
 * General destroy routine.
 */
void
wavetree_value_destroy(wavetree_value_t *p);

/*
 * Global gaussian proposal
 */
wavetree_value_t *
wavetree_value_create_global_gaussian(double std,
				      unsigned long int seed);

/*
 * Depth dependent gaussian proposal
 */
wavetree_value_t *
wavetree_value_create_depth_gaussian(int ndepths,
				     double *std,
				     unsigned long int seed);

wavetree_value_t *
wavetree_value_create_depth_cauchy(int ndepths,
				   double *std,
				   unsigned long int seed);

/*
 * SCAM adaptive proposal (Haario 2005) with initial depth dependent gaussian proposal.
 */
wavetree_value_t *
wavetree_value_create_depth_gaussian_scam(int ndepths,
					  double *std,
					  double epsilon,
					  double s,
					  int threshold,
					  unsigned long int seed,
					  coefficient_histogram_t *histogram);

wavetree_value_t *
wavetree_value_create_gaussian_am(double std0,
				  double epsilon,
				  double A,
				  double tau,
				  unsigned long int seed,
				  coefficient_histogram_t *histogram);

wavetree_value_t *
wavetree_value_create_cauchy_am(double std0,
				double epsilon,
				double A,
				double tau,
				unsigned long int seed,
				coefficient_histogram_t *histogram);


#endif /* wavetree_value_proposal_h */
