#ifndef wavetreesphereface2d_h
#define wavetreesphereface2d_h

#include <stdint.h>

#include "manifold.h"
#include "coefficient_histogram.h"
#include "multiset_int_double.h"
#include "chain_history.h"
#include "wavetree.h"

typedef struct _wavetreesphereface2d wavetreesphereface2d_t;

wavetreesphereface2d_t *
wavetreesphereface2d_create(manifold_t *manifold,
			    double alpha);

void
wavetreesphereface2d_destroy(wavetreesphereface2d_t *t);

int
wavetreesphereface2d_save(const wavetreesphereface2d_t *t,
			  const char *filename);

int 
wavetreesphereface2d_load(wavetreesphereface2d_t *t,
			  const char *filename);

int 
wavetreesphereface2d_load_promote(wavetreesphereface2d_t *t,
				  const char *filename);

int
wavetreesphereface2d_encode(wavetreesphereface2d_t *t,
			    char *buffer,
			    int maxsize);

int
wavetreesphereface2d_decode(wavetreesphereface2d_t *t,
			    char *buffer,
			    int len);

int
wavetreesphereface2d_get_size(wavetreesphereface2d_t *t);

int
wavetreesphereface2d_get_ncoeff(wavetreesphereface2d_t *t);

int
wavetreesphereface2d_initialize(wavetreesphereface2d_t *t,
				double dc);

double
wavetreesphereface2d_dc(wavetreesphereface2d_t *t);

int
wavetreesphereface2d_prunable_leaves(const wavetreesphereface2d_t *t);

int
wavetreesphereface2d_attachable_branches(const wavetreesphereface2d_t *t);

int
wavetreesphereface2d_coeff_count(const wavetreesphereface2d_t *t);

int
wavetreesphereface2d_coeff_count_recurse(const wavetreesphereface2d_t *t);

int
wavetreesphereface2d_coeff_count_scan(const wavetreesphereface2d_t *t);

int
wavetreesphereface2d_coeff_count_set(const wavetreesphereface2d_t *t);

int
wavetreesphereface2d_map_to_array(const wavetreesphereface2d_t *t, 
				  double *a, 
				  int n);

int
wavetreesphereface2d_map_from_array(wavetreesphereface2d_t *t, 
				    const double *a, 
				    int n);

int
wavetreesphereface2d_propose_value(wavetreesphereface2d_t *t,
				   int i,
				   int d,
				   double value);

int
wavetreesphereface2d_propose_birth(wavetreesphereface2d_t *t,
				   int i,
				   int d,
				   double value);

int 
wavetreesphereface2d_propose_death(wavetreesphereface2d_t *t,
				   int i,
				   int d,
				   double *old_value);

int 
wavetreesphereface2d_propose_move(wavetreesphereface2d_t *t,
				  int i,
				  int d,
				  int new_i,
				  double new_value);

int 
wavetreesphereface2d_undo(wavetreesphereface2d_t *t);

int
wavetreesphereface2d_commit(wavetreesphereface2d_t *t);

typedef int (*wavetreesphereface2d_birth_func_t)(void *user, 
						 int i,
						 int j,
						 int level,
						 int maxlevel,
						 double parent_coeff, 
						 double *coeff);
int
wavetreesphereface2d_valid(wavetreesphereface2d_t *t);

void
wavetreesphereface2d_print_mask(wavetreesphereface2d_t *t, int depth);

void
wavetreesphereface2d_print_setinfo(wavetreesphereface2d_t *t);

void
wavetreesphereface2d_dump_sets(wavetreesphereface2d_t *t);

void
wavetreesphereface2d_dump_coeffs(wavetreesphereface2d_t *t);

int
wavetreesphereface2d_get_indices(wavetreesphereface2d_t *t, int *set, int *n);

int
wavetreesphereface2d_depth(wavetreesphereface2d_t *t);

int
wavetreesphereface2d_maxdepth(const wavetreesphereface2d_t *t);

int
wavetreesphereface2d_update_histogram(const wavetreesphereface2d_t *t, coefficient_histogram_t *hist);

int
wavetreesphereface2d_depth_filter(void *wt, int subset, int index);

const multiset_int_double_t *
wavetreesphereface2d_get_S_v(const wavetreesphereface2d_t *t);

int
wavetreesphereface2d_set_from_S_v(wavetreesphereface2d_t *t,
				  const multiset_int_double_t *S_vp);

int
wavetreesphereface2d_get_last_perturbation(wavetreesphereface2d_t *t,
					   chain_history_change_t *step);

int
wavetreesphereface2d_set_invalid_perturbation(wavetreesphereface2d_t *t, wavetree_perturb_t p);


/*
 * Functions for setting up a birth proposal
 */
int
wavetreesphereface2d_choose_birth_depth(const wavetreesphereface2d_t *t,
					double u,
					int maxdepth,
					int *depth,
					double *prob);

int
wavetreesphereface2d_reverse_birth_depth(const wavetreesphereface2d_t *t,
					 int depth,
					 int maxdepth,
					 double *prob);

int
wavetreesphereface2d_choose_birth(const wavetreesphereface2d_t *t,
				  int depth,
				  double u,
				  int *coeff,
				  double *prob);

int
wavetreesphereface2d_choose_birth_global(const wavetreesphereface2d_t *t,
					 double u,
					 int maxdepth,
					 int *depth,
					 int *coeff,
					 double *prob);

int
wavetreesphereface2d_reverse_birth(const wavetreesphereface2d_t *t,
				   int depth,
				   int coeff,
				   double *prob);

int
wavetreesphereface2d_reverse_birth_global(const wavetreesphereface2d_t *t,
					  int maxdepth,
					  int depth,
					  int coeff,
					  double *prob);

/*
 * Functions for setting up a death proposal
 */
int
wavetreesphereface2d_choose_death_depth(const wavetreesphereface2d_t *t,
					double u,
					int maxdepth,
					int *depth,
					double *prob);

int
wavetreesphereface2d_reverse_death_depth(const wavetreesphereface2d_t *t,
					 int depth,
					 int maxdepth,
					 double *prob);

int
wavetreesphereface2d_choose_death(const wavetreesphereface2d_t *t,
				  int depth,
				  double u,
				  int *coeff,
				  double *prob);

int
wavetreesphereface2d_choose_death_global(const wavetreesphereface2d_t *t,
					 double u,
					 int maxdepth,
					 int *depth,
					 int *coeff,
					 double *prob);

int
wavetreesphereface2d_reverse_death(const wavetreesphereface2d_t *t,
				   int depth,
				   int coeff,
				   double *prob);

int
wavetreesphereface2d_reverse_death_global(const wavetreesphereface2d_t *t,
					  int maxdepth,
					  int depth,
					  int coeff,
					  double *prob);

/*
 * Functions for setting up a value proposal
 */
int
wavetreesphereface2d_choose_value_depth(const wavetreesphereface2d_t *t,
					double u,
					int maxdepth,
					int *depth,
					double *prob);

int
wavetreesphereface2d_choose_value(const wavetreesphereface2d_t *t,
				  int depth,
				  double u,
				  int *coeff,
				  double *prob);

int
wavetreesphereface2d_choose_value_global(const wavetreesphereface2d_t *t,
					 double u,
					 int maxdepth,
					 int *depth,
					 int *coeff,
					 double *prob);

/*
 * Functions for performing a move proposal
 */
int
wavetreesphereface2d_choose_move_depth(const wavetreesphereface2d_t *t,
				       double u,
				       int maxdepth,
				       int *depth,
				       double *prob);

int
wavetreesphereface2d_choose_move(const wavetreesphereface2d_t *t,
				 int depth,
				 double u,
				 int *coeff,
				 double *prob);

int
wavetreesphereface2d_choose_move_global(const wavetreesphereface2d_t *t,
					double u,
					int maxdepth,
					int *depth,
					int *coeff,
					double *prob);

int
wavetreesphereface2d_move_available_siblings(const wavetreesphereface2d_t *t,
					     int depth,
					     int coeff,
					     int *siblings,
					     int *nsibling);

int
wavetreesphereface2d_choose_move_sibling(const wavetreesphereface2d_t *t,
					 double u,
					 int depth,
					 int coeff,
					 int *sibling,
					 double *prob);

int
wavetreesphereface2d_reverse_choose_move_sibling(const wavetreesphereface2d_t *t,
						 int depth,
						 int coeff,
						 int sibling,
						 double *prob);

/*
 *
 */
typedef
int (*wavetreesphereface2d_perturb_func_t)(void *user, 
					   int i, 
					   int j, 
					   int level,
					   int maxlevel,
					   double parent_coeff,
					   double *coeff,
					   double *prior_ratio);

int
wavetreesphereface2d_perturb(wavetreesphereface2d_t *t,
			     wavetreesphereface2d_perturb_func_t f,
			     void *user,
			     double *prior_ratio);

int
wavetreesphereface2d_parent_index(const wavetreesphereface2d_t *t,
				  int c);

int
wavetreesphereface2d_2dindices(const wavetreesphereface2d_t *t,
			       int i,
			       int *ii,
			       int *ij);

int
wavetreesphereface2d_from_2dindices(const wavetreesphereface2d_t *t,
				    int ii,
				    int ij);

int
wavetreesphereface2d_depthofindex(const wavetreesphereface2d_t *t,
				  int i);

int
wavetreesphereface2d_get_coeff(const wavetreesphereface2d_t *t,
			       int i,
			       double *coeff);

/*
 * Internal functions exposed for testing
 */

int
wavetreesphereface2d_max_child_count(wavetreesphereface2d_t *t);

int
wavetreesphereface2d_child_count(const wavetreesphereface2d_t *t,
				 int index,
				 int depth);

int
wavetreesphereface2d_child_indices(const wavetreesphereface2d_t *t,
				   int index,
				   int depth,
				   int *indices,
				   int *n,
				   int nmax);

int
wavetreesphereface2d_ncoeff_at_depth(int depth);

int
wavetreesphereface2d_ntriangle_at_depth(int depth);

int
wavetreesphereface2d_depth_base(const wavetreesphereface2d_t *t, int depth);

#endif /* wavetreesphereface2d_h */
			  
