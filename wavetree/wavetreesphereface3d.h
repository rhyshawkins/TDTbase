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

#ifndef wavetreesphereface3d_h
#define wavetreesphereface3d_h

#include <stdint.h>

#include "manifold.h"
#include "coefficient_histogram.h"

typedef struct _wavetreesphereface3d wavetreesphereface3d_t;

wavetreesphereface3d_t *wavetreesphereface3d_create(manifold_t *m,
						    int max_depth, 
						    double alpha);

void wavetreesphereface3d_destroy(wavetreesphereface3d_t *t);

int
wavetreesphereface3d_save(const wavetreesphereface3d_t *t,
			  const char *filename);

int 
wavetreesphereface3d_load(wavetreesphereface3d_t *t,
			  const char *filename);

int 
wavetreesphereface3d_load_promote(wavetreesphereface3d_t *t,
				  const char *filename);


int
wavetreesphereface3d_get_rowstride(wavetreesphereface3d_t *t);

int
wavetreesphereface3d_get_width(wavetreesphereface3d_t *t);

int
wavetreesphereface3d_get_height(wavetreesphereface3d_t *t);

int
wavetreesphereface3d_get_depth(wavetreesphereface3d_t *t);

int
wavetreesphereface3d_get_image_size(wavetreesphereface3d_t *t);

int
wavetreesphereface3d_get_coeff_size(wavetreesphereface3d_t *t);

int
wavetreesphereface3d_initialize(wavetreesphereface3d_t *t,
				double dc);

double 
wavetreesphereface3d_dc(wavetreesphereface3d_t *t);

int 
wavetreesphereface3d_prunable_leaves(const wavetreesphereface3d_t *t);

int 
wavetreesphereface3d_attachable_branches(const wavetreesphereface3d_t *t);

int 
wavetreesphereface3d_coeff_count(const wavetreesphereface3d_t *t);

int 
wavetreesphereface3d_map_to_array(wavetreesphereface3d_t *t, 
				  double *a,
				  int n);

int
wavetreesphereface3d_propose_value(wavetreesphereface3d_t *t,
				   int i,
				   int d,
				   double value);

int
wavetreesphereface3d_propose_birth(wavetreesphereface3d_t *t,
				   int i,
				   int d,
				   double value);

int 
wavetreesphereface3d_propose_death(wavetreesphereface3d_t *t,
				   int i,
				   int d,
				   double *old_value);

int 
wavetreesphereface3d_undo(wavetreesphereface3d_t *t);

int
wavetreesphereface3d_commit(wavetreesphereface3d_t *t);

int 
wavetreesphereface3d_valid(wavetreesphereface3d_t *t);

void 
wavetreesphereface3d_print_setinfo(wavetreesphereface3d_t *t);

void 
wavetreesphereface3d_dump_sets(wavetreesphereface3d_t *t);

void 
wavetreesphereface3d_dump_coeffs(wavetreesphereface3d_t *t);

int 
wavetreesphereface3d_depth(wavetreesphereface3d_t *t);

int 
wavetreesphereface3d_maxdepth(wavetreesphereface3d_t *t);

int 
wavetreesphereface3d_update_histogram(const wavetreesphereface3d_t *t, coefficient_histogram_t *hist);

int 
wavetreesphereface3d_depth_filter(void *wt, int subset, int index);

/*
 * Functions for setting up a birth proposal
 */
int 
wavetreesphereface3d_choose_birth_depth(const wavetreesphereface3d_t *t, double u, int maxdepth, int *depth, double *prob);

int 
wavetreesphereface3d_reverse_birth_depth(const wavetreesphereface3d_t *t, int depth, int maxdepth, double *prob);

int 
wavetreesphereface3d_choose_birth(const wavetreesphereface3d_t *t, int depth, double u, int *coeff, double *prob);

int 
wavetreesphereface3d_choose_birth_global(const wavetreesphereface3d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

int 
wavetreesphereface3d_reverse_birth(const wavetreesphereface3d_t *t, int depth, int coeff, double *prob);

int 
wavetreesphereface3d_reverse_birth_global(const wavetreesphereface3d_t *t, int maxdepth, int depth, int coeff, double *prob);

/*
 * Functions for setting up a death proposal
 */
int
wavetreesphereface3d_choose_death_depth(const wavetreesphereface3d_t *t, double u, int maxdepth, int *depth, double *prob);

int 
wavetreesphereface3d_reverse_death_depth(const wavetreesphereface3d_t *t, int depth, int maxdepth, double *prob);

int 
wavetreesphereface3d_choose_death(const wavetreesphereface3d_t *t, int depth, double u, int *coeff, double *prob);

int 
wavetreesphereface3d_choose_death_global(const wavetreesphereface3d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

int 
wavetreesphereface3d_reverse_death(const wavetreesphereface3d_t *t, int depth, int coeff, double *prob);

int 
wavetreesphereface3d_reverse_death_global(const wavetreesphereface3d_t *t, int maxdepth, int depth, int coeff, double *prob);

/*
 * Functions for setting up a value proposal
 */
int 
wavetreesphereface3d_choose_value_depth(const wavetreesphereface3d_t *t, double u, int maxdepth, int *depth, double *prob);

int 
wavetreesphereface3d_choose_value(const wavetreesphereface3d_t *t, int depth, double u, int *coeff, double *prob);

int 
wavetreesphereface3d_choose_value_global(const wavetreesphereface3d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

typedef int (*wavetreesphereface3d_perturb_func_t)(void *user, 
						   int i, 
						   int j, 
						   int k,
						   int level,
						   int maxlevel,
						   double parent_coeff,
						   double *coeff,
						   double *prior_ratio);
int 
wavetreesphereface3d_perturb(wavetreesphereface3d_t *t,
			     wavetreesphereface3d_perturb_func_t f,
			     void *user,
			     double *prior_ratio);

int 
wavetreesphereface3d_parent_index(wavetreesphereface3d_t *t, int c);

int
wavetreesphereface3d_index_to_separatesphereindices(wavetreesphereface3d_t *wt,
						    int i,
						    int *t,
						    int *r,
						    int *td,
						    int *rd,
						    int *d);

int 
wavetreesphereface3d_index_to_sphereindices(wavetreesphereface3d_t *wt,
					    int i,
					    int *t,
					    int *r,
					    int *d);

int 
wavetreesphereface3d_sphereindices_to_index(wavetreesphereface3d_t *wt,
					    int t,
					    int r,
					    int d);

int 
wavetreesphereface3d_depthofindex(wavetreesphereface3d_t *wt,
				  int i);

int 
wavetreesphereface3d_get_coeff(wavetreesphereface3d_t *t,
			       int i,
			       int d,
			       double *coeff);

int
wavetreesphereface3d_get_child_indices(wavetreesphereface3d_t *t,
				       int index,
				       int *indices,
				       int maxindices);

int
wavetreesphereface3d_child_count(wavetreesphereface3d_t *t, int index, int depth);

int
wavetreesphereface3d_index_to_offset(wavetreesphereface3d_t *t,
				     int depth,
				     int index);

int
wavetreesphereface3d_offset_to_index(wavetreesphereface3d_t *t,
				     int offset,
				     int *depth,
				     int *index);

#endif /* wavetreesphereface3d_h */
			  
