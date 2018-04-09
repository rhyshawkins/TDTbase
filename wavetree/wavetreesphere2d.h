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


#ifndef wavetreesphere2d_h
#define wavetreesphere2d_h

#include <stdint.h>

#include "coefficient_histogram.h"

typedef struct _wavetreesphere2d wavetreesphere2d_t;

wavetreesphere2d_t *
wavetreesphere2d_create(int max_depth);

void 
wavetreesphere2d_destroy(wavetreesphere2d_t *t);

int
wavetreesphere2d_save(const wavetreesphere2d_t *t,
		      const char *filename);

int 
wavetreesphere2d_load(wavetreesphere2d_t *t,
		      const char *filename);

int 
wavetreesphere2d_load_promote(wavetreesphere2d_t *t,
			      const char *filename);

int
wavetreesphere2d_get_width(wavetreesphere2d_t *t);

int
wavetreesphere2d_get_size(wavetreesphere2d_t *t);

int 
wavetreesphere2d_initialize(wavetreesphere2d_t *t,
			    double dc);

double 
wavetreesphere2d_dc(wavetreesphere2d_t *t);

int 
wavetreesphere2d_coeff_count(const wavetreesphere2d_t *t);

int 
wavetreesphere2d_map_to_array(const wavetreesphere2d_t *t, 
			      double *a, 
			      int n);

int
wavetreesphere2d_propose_value(wavetreesphere2d_t *t,
			       int i,
			       int d,
			       double value);

int
wavetreesphere2d_propose_birth(wavetreesphere2d_t *t,
			       int i,
			       int d,
			       double value);

int 
wavetreesphere2d_propose_death(wavetreesphere2d_t *t,
			       int i,
			       int d,
			       double *old_value);

int 
wavetreesphere2d_undo(wavetreesphere2d_t *t);

int
wavetreesphere2d_commit(wavetreesphere2d_t *t);

int 
wavetreesphere2d_valid(wavetreesphere2d_t *t);

void 
wavetreesphere2d_dump_sets(wavetreesphere2d_t *t);

void 
wavetreesphere2d_dump_coeffs(wavetreesphere2d_t *t);

int 
wavetreesphere2d_get_indices(wavetreesphere2d_t *t, int *set, int *n);

int 
wavetreesphere2d_depth(wavetreesphere2d_t *t);

int 
wavetreesphere2d_maxdepth(wavetreesphere2d_t *t);

int 
wavetreesphere2d_update_histogram(const wavetreesphere2d_t *t, coefficient_histogram_t *hist);

int 
wavetreesphere2d_depth_filter(void *wt, int subset, int index);

/*
 * Functions for setting up a birth proposal
 */
int 
wavetreesphere2d_choose_birth_depth(const wavetreesphere2d_t *t, 
				    double u, 
				    int maxdepth, 
				    int *depth, 
				    double *prob);

int 
wavetreesphere2d_reverse_birth_depth(const wavetreesphere2d_t *t, 
				     int depth, 
				     int maxdepth, 
				     double *prob);

int 
wavetreesphere2d_choose_birth(const wavetreesphere2d_t *t, 
			      int depth, 
			      double u, 
			      int *coeff, 
			      double *prob);

int 
wavetreesphere2d_choose_birth_global(const wavetreesphere2d_t *t, 
				     double u, 
				     int maxdepth, 
				     int *depth, 
				     int *coeff, 
				     double *prob);

int 
wavetreesphere2d_reverse_birth(const wavetreesphere2d_t *t, 
			       int depth, 
			       int coeff, 
			       double *prob);

int 
wavetreesphere2d_reverse_birth_global(const wavetreesphere2d_t *t, 
				      int maxdepth, 
				      int depth, 
				      int coeff, 
				      double *prob);

/*
 * Functions for setting up a death proposal
 */
int 
wavetreesphere2d_choose_death_depth(const wavetreesphere2d_t *t, 
				    double u, 
				    int maxdepth, 
				    int *depth, 
				    double *prob);

int 
wavetreesphere2d_reverse_death_depth(const wavetreesphere2d_t *t, 
				     int depth, 
				     int maxdepth, 
				     double *prob);

int 
wavetreesphere2d_choose_death(const wavetreesphere2d_t *t, 
			      int depth, 
			      double u, 
			      int *coeff, 
			      double *prob);

int 
wavetreesphere2d_choose_death_global(const wavetreesphere2d_t *t, 
				     double u, 
				     int maxdepth, 
				     int *depth, 
				     int *coeff, 
				     double *prob);

int 
wavetreesphere2d_reverse_death(const wavetreesphere2d_t *t, 
			       int depth, 
			       int coeff, 
			       double *prob);

int 
wavetreesphere2d_reverse_death_global(const wavetreesphere2d_t *t, 
				      int maxdepth, 
				      int depth, 
				      int coeff, 
				      double *prob);

/*
 * Functions for setting up a value proposal
 */
int wavetreesphere2d_choose_value_depth(const wavetreesphere2d_t *t, 
					double u, 
					int maxdepth, 
					int *depth, 
					double *prob);

int wavetreesphere2d_choose_value(const wavetreesphere2d_t *t, 
				  int depth, 
				  double u, 
				  int *coeff, 
				  double *prob);

int wavetreesphere2d_choose_value_global(const wavetreesphere2d_t *t, 
					 double u, 
					 int maxdepth, 
					 int *depth, 
					 int *coeff, 
					 double *prob);

typedef int (*wavetreesphere2d_perturb_func_t)(void *user, 
					       int i, 
					       int j, 
					       int level,
					       int maxlevel,
					       double parent_coeff,
					       double *coeff,
					       double *prior_ratio);

int wavetreesphere2d_perturb(wavetreesphere2d_t *t,
			     wavetreesphere2d_perturb_func_t f,
			     void *user,
			     double *prior_ratio);

int wavetreesphere2d_parent_index(wavetreesphere2d_t *t, int c);

int wavetreesphere2d_2dindices(const wavetreesphere2d_t *t,
			       int i,
			       int *ii,
			       int *ij);

int wavetreesphere2d_from_2dindices(const wavetreesphere2d_t *t,
				    int ii,
				    int ij);

int wavetreesphere2d_depthofindex(const wavetreesphere2d_t *t,
				  int i);

int wavetreesphere2d_get_coeff(const wavetreesphere2d_t *t,
			       int i,
			       double *coeff);

/*
 * Internal functions exposed for testing
 */

int wavetreesphere2d_TL(const wavetreesphere2d_t *t, int i);
int wavetreesphere2d_TR(const wavetreesphere2d_t *t, int i);
int wavetreesphere2d_BL(const wavetreesphere2d_t *t, int i);
int wavetreesphere2d_BR(const wavetreesphere2d_t *t, int i);


/*
 * Old functions
 */


#endif /* wavetreesphere2d_h */
			  
