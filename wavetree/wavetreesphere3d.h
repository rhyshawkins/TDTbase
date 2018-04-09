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

#ifndef wavetreesphere3d_h
#define wavetreesphere3d_h

#include <stdint.h>

#include "manifold.h"
#include "coefficient_histogram.h"

typedef struct _wavetreesphere3d wavetreesphere3d_t;

wavetreesphere3d_t *wavetreesphere3d_create(manifold_t *m,
					    int max_depth, 
					    double alpha);

void wavetreesphere3d_destroy(wavetreesphere3d_t *t);

int
wavetreesphere3d_save(const wavetreesphere3d_t *t,
		      const char *filename);

int 
wavetreesphere3d_load(wavetreesphere3d_t *t,
		      const char *filename);

int 
wavetreesphere3d_load_promote(wavetreesphere3d_t *t,
			      const char *filename);

int
wavetreesphere3d_get_width(wavetreesphere3d_t *t);

int
wavetreesphere3d_get_height(wavetreesphere3d_t *t);

int
wavetreesphere3d_get_depth(wavetreesphere3d_t *t);

int
wavetreesphere3d_get_image_size(wavetreesphere3d_t *t);

int
wavetreesphere3d_get_coeff_size(wavetreesphere3d_t *t);

int wavetreesphere3d_initialize(wavetreesphere3d_t *t,
				double dc);

double 
wavetreesphere3d_dc(wavetreesphere3d_t *t);

int 
wavetreesphere3d_prunable_leaves(const wavetreesphere3d_t *t);

int 
wavetreesphere3d_attachable_branches(const wavetreesphere3d_t *t);

int 
wavetreesphere3d_coeff_count(const wavetreesphere3d_t *t);

int 
wavetreesphere3d_map_to_array(wavetreesphere3d_t *t, 
			      double *a,
			      int n);

int
wavetreesphere3d_propose_value(wavetreesphere3d_t *t,
			       int i,
			       int d,
			       double value);

int
wavetreesphere3d_propose_birth(wavetreesphere3d_t *t,
			       int i,
			       int d,
			       double value);

int 
wavetreesphere3d_propose_death(wavetreesphere3d_t *t,
			       int i,
			       int d,
			       double *old_value);

int 
wavetreesphere3d_undo(wavetreesphere3d_t *t);

int
wavetreesphere3d_commit(wavetreesphere3d_t *t);

int 
wavetreesphere3d_valid(wavetreesphere3d_t *t);

void 
wavetreesphere3d_print_setinfo(wavetreesphere3d_t *t);

void 
wavetreesphere3d_dump_sets(wavetreesphere3d_t *t);

void 
wavetreesphere3d_dump_coeffs(wavetreesphere3d_t *t);

int 
wavetreesphere3d_depth(wavetreesphere3d_t *t);

int 
wavetreesphere3d_maxdepth(wavetreesphere3d_t *t);

int 
wavetreesphere3d_update_histogram(const wavetreesphere3d_t *t, coefficient_histogram_t *hist);

int 
wavetreesphere3d_depth_filter(void *wt, int subset, int index);

/*
 * Functions for setting up a birth proposal
 */
int 
wavetreesphere3d_choose_birth_depth(const wavetreesphere3d_t *t, double u, int maxdepth, int *depth, double *prob);

int 
wavetreesphere3d_reverse_birth_depth(const wavetreesphere3d_t *t, int depth, int maxdepth, double *prob);

int 
wavetreesphere3d_choose_birth(const wavetreesphere3d_t *t, int depth, double u, int *coeff, double *prob);

int 
wavetreesphere3d_choose_birth_global(const wavetreesphere3d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

int 
wavetreesphere3d_reverse_birth(const wavetreesphere3d_t *t, int depth, int coeff, double *prob);

int 
wavetreesphere3d_reverse_birth_global(const wavetreesphere3d_t *t, int maxdepth, int depth, int coeff, double *prob);

/*
 * Functions for setting up a death proposal
 */
int
wavetreesphere3d_choose_death_depth(const wavetreesphere3d_t *t, double u, int maxdepth, int *depth, double *prob);

int 
wavetreesphere3d_reverse_death_depth(const wavetreesphere3d_t *t, int depth, int maxdepth, double *prob);

int 
wavetreesphere3d_choose_death(const wavetreesphere3d_t *t, int depth, double u, int *coeff, double *prob);

int 
wavetreesphere3d_choose_death_global(const wavetreesphere3d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

int 
wavetreesphere3d_reverse_death(const wavetreesphere3d_t *t, int depth, int coeff, double *prob);

int 
wavetreesphere3d_reverse_death_global(const wavetreesphere3d_t *t, int maxdepth, int depth, int coeff, double *prob);

/*
 * Functions for setting up a value proposal
 */
int 
wavetreesphere3d_choose_value_depth(const wavetreesphere3d_t *t, double u, int maxdepth, int *depth, double *prob);

int 
wavetreesphere3d_choose_value(const wavetreesphere3d_t *t, int depth, double u, int *coeff, double *prob);

int 
wavetreesphere3d_choose_value_global(const wavetreesphere3d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

typedef int (*wavetreesphere3d_perturb_func_t)(void *user, 
					       int i, 
					       int j, 
					       int k,
					       int level,
					       int maxlevel,
					       double parent_coeff,
					       double *coeff,
					       double *prior_ratio);
int 
wavetreesphere3d_perturb(wavetreesphere3d_t *t,
			 wavetreesphere3d_perturb_func_t f,
			 void *user,
			 double *prior_ratio);

int 
wavetreesphere3d_parent_index(wavetreesphere3d_t *t, int c);

int 
wavetreesphere3d_sphere2dindices(wavetreesphere3d_t *t,
				 int i,
				 int *ii,
				 int *ik);

int 
wavetreesphere3d_from_sphere2dindices(wavetreesphere3d_t *t,
				      int ii,
				      int ik);

int 
wavetreesphere3d_depthofindex(wavetreesphere3d_t *t,
			      int i);

int 
wavetreesphere3d_get_coeff(wavetreesphere3d_t *t,
			   int i,
			   int d,
			   double *coeff);

int
wavetreesphere3d_get_child_indices(wavetreesphere3d_t *t,
				   int index,
				   int *indices,
				   int maxindices);

int
wavetreesphere3d_child_count(wavetreesphere3d_t *t, int index, int depth);

int
wavetreesphere3d_index_to_offset(wavetreesphere3d_t *t,
				 int depth,
				 int index);

int
wavetreesphere3d_offset_to_index(wavetreesphere3d_t *t,
				 int offset,
				 int *depth,
				 int *index);

#endif /* wavetreesphere3d_h */
			  
