#ifndef wavetree2d_h
#define wavetree2d_h

#include <stdint.h>

#include "coefficient_histogram.h"
#include "multiset_int_double.h"
#include "chain_history.h"
#include "wavetree.h"

typedef struct _wavetree2d wavetree2d_t;

wavetree2d_t *wavetree2d_create(int xdegree,
				int ydegree,
				double alpha);

void wavetree2d_destroy(wavetree2d_t *t);

int
wavetree2d_save(const wavetree2d_t *t,
		const char *filename);

int 
wavetree2d_load(wavetree2d_t *t,
		const char *filename);

int 
wavetree2d_load_promote(wavetree2d_t *t,
			const char *filename);

int
wavetree2d_get_width(wavetree2d_t *t);

int
wavetree2d_get_height(wavetree2d_t *t);

int
wavetree2d_get_size(wavetree2d_t *t);

int wavetree2d_initialize(wavetree2d_t *t,
			  double dc);

double wavetree2d_dc(wavetree2d_t *t);

int wavetree2d_prunable_leaves(const wavetree2d_t *t);

int wavetree2d_attachable_branches(const wavetree2d_t *t);

int wavetree2d_coeff_count(const wavetree2d_t *t);

int wavetree2d_coeff_count_recurse(const wavetree2d_t *t);

int wavetree2d_coeff_count_scan(const wavetree2d_t *t);

int wavetree2d_coeff_count_set(const wavetree2d_t *t);

int wavetree2d_map_to_array(const wavetree2d_t *t, 
			    double *a, 
			    int n);

int wavetree2d_map_from_array(wavetree2d_t *t, 
			      const double *a, 
			      int n);

int
wavetree2d_propose_value(wavetree2d_t *t,
			 int i,
			 int d,
			 double value);

int
wavetree2d_propose_birth(wavetree2d_t *t,
			 int i,
			 int d,
			 double value);

int 
wavetree2d_propose_death(wavetree2d_t *t,
			 int i,
			 int d,
			 double *old_value);

int 
wavetree2d_propose_move(wavetree2d_t *t,
			int i,
			int d,
			int new_i,
			double new_value);

int 
wavetree2d_undo(wavetree2d_t *t);

int
wavetree2d_commit(wavetree2d_t *t);

typedef int (*wavetree2d_birth_func_t)(void *user, 
				       int i,
				       int j,
				       int level,
				       int maxlevel,
				       double parent_coeff, 
				       double *coeff);
int wavetree2d_valid(wavetree2d_t *t);

void wavetree2d_print_mask(wavetree2d_t *t, int depth);

void wavetree2d_print_setinfo(wavetree2d_t *t);

void wavetree2d_dump_sets(wavetree2d_t *t);

void wavetree2d_dump_coeffs(wavetree2d_t *t);

int wavetree2d_generate_dyck_word(wavetree2d_t *t, char *buffer, int maxlength);

int wavetree2d_generate_dyck_binary(wavetree2d_t *t, uint64_t *binary);

int wavetree2d_get_indices(wavetree2d_t *t, int *set, int *n);

int wavetree2d_depth(wavetree2d_t *t);

int wavetree2d_maxdepth(const wavetree2d_t *t);

int wavetree2d_update_histogram(const wavetree2d_t *t, coefficient_histogram_t *hist);

int wavetree2d_depth_filter(void *wt, int subset, int index);

const multiset_int_double_t *
wavetree2d_get_S_v(const wavetree2d_t *t);

int
wavetree2d_set_from_S_v(wavetree2d_t *t,
			multiset_int_double_t *S_vp);

int
wavetree2d_get_last_perturbation(wavetree2d_t *t,
				 chain_history_change_t *step);

int
wavetree2d_set_invalid_perturbation(wavetree2d_t *t, wavetree_perturb_t p);
			    

/*
 * Functions for setting up a birth proposal
 */
int wavetree2d_choose_birth_depth(const wavetree2d_t *t, double u, int maxdepth, int *depth, double *prob);

int wavetree2d_reverse_birth_depth(const wavetree2d_t *t, int depth, int maxdepth, double *prob);

int wavetree2d_choose_birth(const wavetree2d_t *t, int depth, double u, int *coeff, double *prob);

int wavetree2d_choose_birth_global(const wavetree2d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

int wavetree2d_reverse_birth(const wavetree2d_t *t, int depth, int coeff, double *prob);

int wavetree2d_reverse_birth_global(const wavetree2d_t *t, int maxdepth, int depth, int coeff, double *prob);

/*
 * Functions for setting up a death proposal
 */
int wavetree2d_choose_death_depth(const wavetree2d_t *t, double u, int maxdepth, int *depth, double *prob);

int wavetree2d_reverse_death_depth(const wavetree2d_t *t, int depth, int maxdepth, double *prob);

int wavetree2d_choose_death(const wavetree2d_t *t, int depth, double u, int *coeff, double *prob);

int wavetree2d_choose_death_global(const wavetree2d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

int wavetree2d_reverse_death(const wavetree2d_t *t, int depth, int coeff, double *prob);

int wavetree2d_reverse_death_global(const wavetree2d_t *t, int maxdepth, int depth, int coeff, double *prob);

/*
 * Functions for setting up a value proposal
 */
int wavetree2d_choose_value_depth(const wavetree2d_t *t, double u, int maxdepth, int *depth, double *prob);

int wavetree2d_choose_value(const wavetree2d_t *t, int depth, double u, int *coeff, double *prob);

int wavetree2d_choose_value_global(const wavetree2d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

/*
 * Functions for performing a move proposal
 */
int wavetree2d_choose_move_depth(const wavetree2d_t *t, double u, int maxdepth, int *depth, double *prob);

int wavetree2d_choose_move(const wavetree2d_t *t, int depth, double u, int *coeff, double *prob);

int wavetree2d_choose_move_global(const wavetree2d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

int wavetree2d_move_available_siblings(const wavetree2d_t *t, int depth, int coeff, int *siblings, int *nsibling);

int wavetree2d_choose_move_sibling(const wavetree2d_t *t, double u, int depth, int coeff, int *sibling, double *prob);

int wavetree2d_reverse_choose_move_sibling(const wavetree2d_t *t, int depth, int coeff, int sibling, double *prob);

/*
 *
 */

typedef int (*wavetree2d_perturb_func_t)(void *user, 
					 int i, 
					 int j, 
					 int level,
					 int maxlevel,
					 double parent_coeff,
					 double *coeff,
					 double *prior_ratio);
int wavetree2d_perturb(wavetree2d_t *t,
		       wavetree2d_perturb_func_t f,
		       void *user,
		       double *prior_ratio);

int wavetree2d_parent_index(const wavetree2d_t *t, int c);

int wavetree2d_child_count(const wavetree2d_t *t, int index, int depth);

int wavetree2d_2dindices(const wavetree2d_t *t,
			 int i,
			 int *ii,
			 int *ij);

int wavetree2d_from_2dindices(const wavetree2d_t *t,
			     int ii,
			     int ij);

int wavetree2d_depthofindex(const wavetree2d_t *t,
			    int i);

int wavetree2d_get_coeff(const wavetree2d_t *t,
			 int i,
			 double *coeff);

/*
 * Internal functions exposed for testing
 */

int wavetree2d_TL(const wavetree2d_t *t, int i);
int wavetree2d_TR(const wavetree2d_t *t, int i);
int wavetree2d_BL(const wavetree2d_t *t, int i);
int wavetree2d_BR(const wavetree2d_t *t, int i);


/*
 * Old functions
 */

#if 0
int wavetree2d_birth(wavetree2d_t *t, 
		     int i, 
		     wavetree2d_birth_func_t f, 
		     void *user);

int wavetree2d_death(wavetree2d_t *t, 
		     int i, 
		     double *coeff_delta);

int wavetree2d_clone(wavetree2d_t *t, 
		     const wavetree2d_t *src);
#endif

#endif /* wavetree2d_h */
			  
