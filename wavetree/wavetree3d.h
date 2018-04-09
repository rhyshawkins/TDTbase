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

#ifndef wavetree3d_h
#define wavetree3d_h

#include <stdint.h>

#include "coefficient_histogram.h"
#include "multiset_int_double.h"
#include "chain_history.h"
#include "wavetree.h"

typedef struct _wavetree3d wavetree3d_t;

wavetree3d_t *wavetree3d_create(int h_width,
				int h_height,
				int h_depth,
				double alpha);

void wavetree3d_destroy(wavetree3d_t *t);

int
wavetree3d_save(const wavetree3d_t *t,
		const char *filename);

int 
wavetree3d_load(wavetree3d_t *t,
		const char *filename);

int
wavetree3d_encode(wavetree3d_t *t,
		  char *buffer,
		  int maxsize);

int
wavetree3d_decode(wavetree3d_t *t,
		  char *buffer,
		  int len);

int 
wavetree3d_load_promote(wavetree3d_t *t,
			const char *filename);

int
wavetree3d_get_width(wavetree3d_t *t);

int
wavetree3d_get_height(wavetree3d_t *t);

int
wavetree3d_get_depth(wavetree3d_t *t);

int
wavetree3d_get_size(wavetree3d_t *t);

int wavetree3d_initialize(wavetree3d_t *t,
			  double dc);

double wavetree3d_dc(wavetree3d_t *t);

int wavetree3d_prunable_leaves(const wavetree3d_t *t);

int wavetree3d_attachable_branches(const wavetree3d_t *t);

int wavetree3d_coeff_count(const wavetree3d_t *t);

int wavetree3d_map_to_array(const wavetree3d_t *t, 
			    double *a, 
			    int n);

int
wavetree3d_propose_value(wavetree3d_t *t,
			 int i,
			 int d,
			 double value);

int
wavetree3d_propose_birth(wavetree3d_t *t,
			 int i,
			 int d,
			 double value);

int 
wavetree3d_propose_death(wavetree3d_t *t,
			 int i,
			 int d,
			 double *old_value);

int 
wavetree3d_undo(wavetree3d_t *t);

int
wavetree3d_commit(wavetree3d_t *t);

int wavetree3d_valid(wavetree3d_t *t);

void wavetree3d_print_setinfo(wavetree3d_t *t);

void wavetree3d_dump_sets(wavetree3d_t *t);

void wavetree3d_dump_coeffs(wavetree3d_t *t);

/* int wavetree3d_generate_dyck_word(wavetree3d_t *t, char *buffer, int maxlength); */

/* int wavetree3d_generate_dyck_binary(wavetree3d_t *t, uint64_t *binary); */

int wavetree3d_depth(wavetree3d_t *t);

int wavetree3d_maxdepth(wavetree3d_t *t);

int wavetree3d_update_histogram(const wavetree3d_t *t, coefficient_histogram_t *hist);

int wavetree3d_depth_filter(void *wt, int subset, int index);

const multiset_int_double_t *
wavetree3d_get_S_v(const wavetree3d_t *t);

int
wavetree3d_set_from_S_v(wavetree3d_t *t,
			const multiset_int_double_t *S_v);

int
wavetree3d_get_last_perturbation(wavetree3d_t *t,
				 chain_history_change_t *step);

int
wavetree3d_set_invalid_perturbation(wavetree3d_t *t, wavetree_perturb_t p);
			    
/*
 * Functions for setting up a birth proposal
 */
int wavetree3d_choose_birth_depth(const wavetree3d_t *t, double u, int maxdepth, int *depth, double *prob);

int wavetree3d_reverse_birth_depth(const wavetree3d_t *t, int depth, int maxdepth, double *prob);

int wavetree3d_choose_birth(const wavetree3d_t *t, int depth, double u, int *coeff, double *prob);

int wavetree3d_choose_birth_global(const wavetree3d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

int wavetree3d_reverse_birth(const wavetree3d_t *t, int depth, int coeff, double *prob);

int wavetree3d_reverse_birth_global(const wavetree3d_t *t, int maxdepth, int depth, int coeff, double *prob);

/*
 * Functions for setting up a death proposal
 */
int wavetree3d_choose_death_depth(const wavetree3d_t *t, double u, int maxdepth, int *depth, double *prob);

int wavetree3d_reverse_death_depth(const wavetree3d_t *t, int depth, int maxdepth, double *prob);

int wavetree3d_choose_death(const wavetree3d_t *t, int depth, double u, int *coeff, double *prob);

int wavetree3d_choose_death_global(const wavetree3d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

int wavetree3d_reverse_death(const wavetree3d_t *t, int depth, int coeff, double *prob);

int wavetree3d_reverse_death_global(const wavetree3d_t *t, int maxdepth, int depth, int coeff, double *prob);

/*
 * Functions for setting up a value proposal
 */
int wavetree3d_choose_value_depth(const wavetree3d_t *t, double u, int maxdepth, int *depth, double *prob);

int wavetree3d_choose_value(const wavetree3d_t *t, int depth, double u, int *coeff, double *prob);

int wavetree3d_choose_value_global(const wavetree3d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

typedef int (*wavetree3d_perturb_func_t)(void *user, 
					 int i, 
					 int j, 
					 int k,
					 int level,
					 int maxlevel,
					 double parent_coeff,
					 double *coeff,
					 double *prior_ratio);
int wavetree3d_perturb(wavetree3d_t *t,
		       wavetree3d_perturb_func_t f,
		       void *user,
		       double *prior_ratio);

int wavetree3d_parent_index(wavetree3d_t *t, int c);

int wavetree3d_3dindices(const wavetree3d_t *t,
			 int i,
			 int *ii,
			 int *ij,
			 int *ik);

int wavetree3d_from_3dindices(const wavetree3d_t *t,
			      int ii,
			      int ij,
			      int ik);

int wavetree3d_depthofindex(const wavetree3d_t *t,
			    int i);

int wavetree3d_get_coeff(const wavetree3d_t *t,
			 int i,
			 int d,
			 double *coeff);

int wavetree3d_generate_dyck_word(wavetree3d_t *t, char *buffer, int maxlength);

int wavetree3d_generate_dyck_binary(wavetree3d_t *t, uint64_t *binary);


/*
 * Internal functions exposed for testing
 */

int wavetree3d_child_count(wavetree3d_t *t, int index, int depth);

int wavetree3d_UTL(const wavetree3d_t *t, int i);
int wavetree3d_UTR(const wavetree3d_t *t, int i);
int wavetree3d_UBL(const wavetree3d_t *t, int i);
int wavetree3d_UBR(const wavetree3d_t *t, int i);
int wavetree3d_LTL(const wavetree3d_t *t, int i);
int wavetree3d_LTR(const wavetree3d_t *t, int i);
int wavetree3d_LBL(const wavetree3d_t *t, int i);
int wavetree3d_LBR(const wavetree3d_t *t, int i);


#endif /* wavetree3d_h */
			  
