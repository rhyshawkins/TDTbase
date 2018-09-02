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

#ifndef wavetree2d_sub_h
#define wavetree2d_sub_h

#include <stdint.h>

#include "coefficient_histogram.h"
#include "multiset_int_double.h"
#include "chain_history.h"
#include "wavetree.h"
#include "wavetreepp.h"

typedef struct _wavetree2d_sub wavetree2d_sub_t;

wavetree2d_sub_t *wavetree2d_sub_create(int xdegree,
				int ydegree,
				double alpha);

void wavetree2d_sub_destroy(wavetree2d_sub_t *t);

int
wavetree2d_sub_save(const wavetree2d_sub_t *t,
		    const char *filename);

int 
wavetree2d_sub_load(wavetree2d_sub_t *t,
		    const char *filename);

int 
wavetree2d_sub_load_promote(wavetree2d_sub_t *t,
			    const char *filename);

int
wavetree2d_sub_encode(wavetree2d_sub_t *t,
		      char *buffer,
		      int max_size);

int
wavetree2d_sub_decode(wavetree2d_sub_t *t,
		      char *buffer,
		      int buffer_len);


int
wavetree2d_sub_get_width(wavetree2d_sub_t *t);

int
wavetree2d_sub_get_height(wavetree2d_sub_t *t);

int
wavetree2d_sub_get_size(wavetree2d_sub_t *t);

int
wavetree2d_sub_base_size(wavetree2d_sub_t *t);

const int *
wavetree2d_sub_base_indices(wavetree2d_sub_t *t);

int
wavetree2d_sub_get_ncoeff(wavetree2d_sub_t *t);

int wavetree2d_sub_initialize(wavetree2d_sub_t *t,
			      double dc);

double wavetree2d_sub_dc(wavetree2d_sub_t *t);

int wavetree2d_sub_prunable_leaves(const wavetree2d_sub_t *t);

int wavetree2d_sub_attachable_branches(const wavetree2d_sub_t *t);

int wavetree2d_sub_coeff_count(const wavetree2d_sub_t *t);

int wavetree2d_sub_coeff_count_recurse(const wavetree2d_sub_t *t);

int wavetree2d_sub_coeff_count_scan(const wavetree2d_sub_t *t);

int wavetree2d_sub_coeff_count_set(const wavetree2d_sub_t *t);

int wavetree2d_sub_map_to_array(const wavetree2d_sub_t *t, 
				double *a, 
				int n);

int wavetree2d_sub_map_impulse_to_array(const wavetree2d_sub_t *t,
					int coeff_index,
					double *a,
					int n);

int wavetree2d_sub_map_from_array(wavetree2d_sub_t *t, 
				  const double *a, 
				  int n);

int wavetree2d_sub_create_from_array_with_threshold(wavetree2d_sub_t *t,
						    const double *a,
						    int n,
						    double threshold);

int
wavetree2d_sub_propose_value(wavetree2d_sub_t *t,
			 int i,
			 int d,
			 double value);

int
wavetree2d_sub_propose_birth(wavetree2d_sub_t *t,
			 int i,
			 int d,
			 double value);

int 
wavetree2d_sub_propose_death(wavetree2d_sub_t *t,
			 int i,
			 int d,
			 double *old_value);

int 
wavetree2d_sub_propose_move(wavetree2d_sub_t *t,
			int i,
			int d,
			int new_i,
			double new_value);

int 
wavetree2d_sub_undo(wavetree2d_sub_t *t);

int
wavetree2d_sub_commit(wavetree2d_sub_t *t);

typedef int (*wavetree2d_sub_birth_func_t)(void *user, 
				       int i,
				       int j,
				       int level,
				       int maxlevel,
				       double parent_coeff, 
				       double *coeff);
int wavetree2d_sub_valid(wavetree2d_sub_t *t);

void wavetree2d_sub_print_mask(wavetree2d_sub_t *t, int depth);

void wavetree2d_sub_print_setinfo(wavetree2d_sub_t *t);

void wavetree2d_sub_dump_sets(wavetree2d_sub_t *t);

void wavetree2d_sub_dump_coeffs(wavetree2d_sub_t *t);

int wavetree2d_sub_generate_dyck_word(wavetree2d_sub_t *t, char *buffer, int maxlength);

int wavetree2d_sub_generate_dyck_binary(wavetree2d_sub_t *t, uint64_t *binary);

int wavetree2d_sub_get_indices(wavetree2d_sub_t *t, int *set, int *n);

int wavetree2d_sub_depth(wavetree2d_sub_t *t);

int wavetree2d_sub_maxdepth(const wavetree2d_sub_t *t);

int wavetree2d_sub_update_histogram(const wavetree2d_sub_t *t, coefficient_histogram_t *hist);

int wavetree2d_sub_depth_filter(void *wt, int subset, int index);

const multiset_int_double_t *
wavetree2d_sub_get_S_v(const wavetree2d_sub_t *t);

int
wavetree2d_sub_set_from_S_v(wavetree2d_sub_t *t,
			    const multiset_int_double_t *S_vp);

int
wavetree2d_sub_get_model(const wavetree2d_sub_t *t,
			 int nmax,
			 int *indices,
			 double *values,
			 int *n);

int
wavetree2d_sub_set_model(wavetree2d_sub_t *t,
			 int *indices,
			 double *values,
			 int n);

int
wavetree2d_sub_get_last_perturbation(wavetree2d_sub_t *t,
				 chain_history_change_t *step);

int
wavetree2d_sub_set_invalid_perturbation(wavetree2d_sub_t *t, wavetree_perturb_t p);
			    

/*
 * Functions for setting up a birth proposal
 */
int wavetree2d_sub_choose_birth_depth(const wavetree2d_sub_t *t, double u, int maxdepth, int *depth, double *prob);

int wavetree2d_sub_reverse_birth_depth(const wavetree2d_sub_t *t, int depth, int maxdepth, double *prob);

int wavetree2d_sub_choose_birth(const wavetree2d_sub_t *t, int depth, double u, int *coeff, double *prob);

int wavetree2d_sub_choose_birth_global(const wavetree2d_sub_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

int wavetree2d_sub_reverse_birth(const wavetree2d_sub_t *t, int depth, int coeff, double *prob);

int wavetree2d_sub_reverse_birth_global(const wavetree2d_sub_t *t, int maxdepth, int depth, int coeff, double *prob);

/*
 * Functions for setting up a death proposal
 */
int wavetree2d_sub_choose_death_depth(const wavetree2d_sub_t *t, double u, int maxdepth, int *depth, double *prob);

int wavetree2d_sub_reverse_death_depth(const wavetree2d_sub_t *t, int depth, int maxdepth, double *prob);

int wavetree2d_sub_choose_death(const wavetree2d_sub_t *t, int depth, double u, int *coeff, double *prob);

int wavetree2d_sub_choose_death_global(const wavetree2d_sub_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

int wavetree2d_sub_reverse_death(const wavetree2d_sub_t *t, int depth, int coeff, double *prob);

int wavetree2d_sub_reverse_death_global(const wavetree2d_sub_t *t, int maxdepth, int depth, int coeff, double *prob);

/*
 * Functions for setting up a value proposal
 */
int wavetree2d_sub_choose_value_depth(const wavetree2d_sub_t *t, double u, int maxdepth, int *depth, double *prob);

int wavetree2d_sub_choose_value(const wavetree2d_sub_t *t, int depth, double u, int *coeff, double *prob);

int wavetree2d_sub_choose_value_global(const wavetree2d_sub_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

/*
 * Functions for performing a move proposal
 */
int wavetree2d_sub_choose_move_depth(const wavetree2d_sub_t *t, double u, int maxdepth, int *depth, double *prob);

int wavetree2d_sub_choose_move(const wavetree2d_sub_t *t, int depth, double u, int *coeff, double *prob);

int wavetree2d_sub_choose_move_global(const wavetree2d_sub_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

int wavetree2d_sub_move_available_siblings(const wavetree2d_sub_t *t, int depth, int coeff, int *siblings, int *nsibling);

int wavetree2d_sub_choose_move_sibling(const wavetree2d_sub_t *t, double u, int depth, int coeff, int *sibling, double *prob);

int wavetree2d_sub_reverse_choose_move_sibling(const wavetree2d_sub_t *t, int depth, int coeff, int sibling, double *prob);

/*
 *
 */

typedef int (*wavetree2d_sub_perturb_func_t)(void *user, 
					 int i, 
					 int j, 
					 int level,
					 int maxlevel,
					 double parent_coeff,
					 double *coeff,
					 double *prior_ratio);
int wavetree2d_sub_perturb(wavetree2d_sub_t *t,
			   wavetree2d_sub_perturb_func_t f,
			   void *user,
			   double *prior_ratio);

int wavetree2d_sub_parent_index(const wavetree2d_sub_t *t, int c);

int wavetree2d_sub_2dindices(const wavetree2d_sub_t *t,
			     int i,
			     int *ii,
			     int *ij);

int wavetree2d_sub_from_2dindices(const wavetree2d_sub_t *t,
				  int ii,
				  int ij);

int wavetree2d_sub_depthofindex(const wavetree2d_sub_t *t,
				int i);

int wavetree2d_sub_get_coeff(const wavetree2d_sub_t *t,
			     int i,
			     double *coeff);

double wavetree2d_sub_logpriorprobability(const wavetree2d_sub_t *t,
					  wavetree_pp_t *pp);
					  
double wavetree2d_sub_mean_abs_deviation(const wavetree2d_sub_t *t);

/*
 * Internal functions exposed for testing
 */

int wavetree2d_sub_max_child_count(wavetree2d_sub_t *t);

int wavetree2d_sub_child_count(const wavetree2d_sub_t *t, int index, int depth);

int wavetree2d_sub_child_indices(const wavetree2d_sub_t *t, int index, int depth, int *indices, int *n, int nmax);

int wavetree2d_sub_TL(const wavetree2d_sub_t *t, int i);
int wavetree2d_sub_TR(const wavetree2d_sub_t *t, int i);
int wavetree2d_sub_BL(const wavetree2d_sub_t *t, int i);
int wavetree2d_sub_BR(const wavetree2d_sub_t *t, int i);

#endif /* wavetree2d_sub_h */
			  
