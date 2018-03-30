#ifndef wavetree3d_sub_h
#define wavetree3d_sub_h

#include <stdint.h>

#include "coefficient_histogram.h"
#include "multiset_int_double.h"
#include "chain_history.h"
#include "wavetree.h"
#include "wavetreepp.h"

typedef struct _wavetree3d_sub wavetree3d_sub_t;

wavetree3d_sub_t *wavetree3d_sub_create(int h_width,
					int h_height,
					int h_depth,
					double alpha);

void wavetree3d_sub_destroy(wavetree3d_sub_t *t);

int
wavetree3d_sub_save(const wavetree3d_sub_t *t,
		    const char *filename);

int 
wavetree3d_sub_load(wavetree3d_sub_t *t,
		    const char *filename);

int
wavetree3d_sub_encode(wavetree3d_sub_t *t,
		      char *buffer,
		      int maxsize);

int
wavetree3d_sub_decode(wavetree3d_sub_t *t,
		  char *buffer,
		  int len);

int 
wavetree3d_sub_load_promote(wavetree3d_sub_t *t,
			const char *filename);

int
wavetree3d_sub_get_width(wavetree3d_sub_t *t);

int
wavetree3d_sub_get_height(wavetree3d_sub_t *t);

int
wavetree3d_sub_get_depth(wavetree3d_sub_t *t);

int
wavetree3d_sub_get_size(wavetree3d_sub_t *t);

int
wavetree3d_sub_get_ncoeff(wavetree3d_sub_t *t);

int wavetree3d_sub_initialize(wavetree3d_sub_t *t,
			  double dc);

double wavetree3d_sub_dc(wavetree3d_sub_t *t);

int wavetree3d_sub_prunable_leaves(const wavetree3d_sub_t *t);

int wavetree3d_sub_attachable_branches(const wavetree3d_sub_t *t);

int wavetree3d_sub_coeff_count(const wavetree3d_sub_t *t);

int wavetree3d_sub_map_to_array(wavetree3d_sub_t *t, 
				double *a, 
				int n);

int
wavetree3d_sub_propose_value(wavetree3d_sub_t *t,
			     int i,
			     int d,
			     double value);

int
wavetree3d_sub_propose_birth(wavetree3d_sub_t *t,
			 int i,
			 int d,
			 double value);

int 
wavetree3d_sub_propose_death(wavetree3d_sub_t *t,
			 int i,
			 int d,
			 double *old_value);

int 
wavetree3d_sub_undo(wavetree3d_sub_t *t);

int
wavetree3d_sub_commit(wavetree3d_sub_t *t);

int wavetree3d_sub_valid(wavetree3d_sub_t *t);

void wavetree3d_sub_print_setinfo(wavetree3d_sub_t *t);

void wavetree3d_sub_dump_sets(wavetree3d_sub_t *t);

void wavetree3d_sub_dump_coeffs(wavetree3d_sub_t *t);

/* int wavetree3d_sub_generate_dyck_word(wavetree3d_sub_t *t, char *buffer, int maxlength); */

/* int wavetree3d_sub_generate_dyck_binary(wavetree3d_sub_t *t, uint64_t *binary); */

int wavetree3d_sub_depth(wavetree3d_sub_t *t);

int wavetree3d_sub_maxdepth(wavetree3d_sub_t *t);

int wavetree3d_sub_update_histogram(const wavetree3d_sub_t *t, coefficient_histogram_t *hist);

int wavetree3d_sub_depth_filter(void *wt, int subset, int index);

const multiset_int_double_t *
wavetree3d_sub_get_S_v(const wavetree3d_sub_t *t);

int
wavetree3d_sub_set_from_S_v(wavetree3d_sub_t *t,
			    const multiset_int_double_t *S_vp);

int
wavetree3d_sub_set_from_S_v_filtered(wavetree3d_sub_t *t,
				     const multiset_int_double_t *S_vp,
				     int maxdepth);

int
wavetree3d_sub_get_last_perturbation(wavetree3d_sub_t *t,
				 chain_history_change_t *step);

int
wavetree3d_sub_set_invalid_perturbation(wavetree3d_sub_t *t, wavetree_perturb_t p);
			    
/*
 * Functions for setting up a birth proposal
 */
int wavetree3d_sub_choose_birth_depth(const wavetree3d_sub_t *t, double u, int maxdepth, int *depth, double *prob);

int wavetree3d_sub_reverse_birth_depth(const wavetree3d_sub_t *t, int depth, int maxdepth, double *prob);

int wavetree3d_sub_choose_birth(const wavetree3d_sub_t *t, int depth, double u, int *coeff, double *prob);

int wavetree3d_sub_choose_birth_global(const wavetree3d_sub_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

int wavetree3d_sub_reverse_birth(const wavetree3d_sub_t *t, int depth, int coeff, double *prob);

int wavetree3d_sub_reverse_birth_global(const wavetree3d_sub_t *t, int maxdepth, int depth, int coeff, double *prob);

/*
 * Functions for setting up a death proposal
 */
int wavetree3d_sub_choose_death_depth(const wavetree3d_sub_t *t, double u, int maxdepth, int *depth, double *prob);

int wavetree3d_sub_reverse_death_depth(const wavetree3d_sub_t *t, int depth, int maxdepth, double *prob);

int wavetree3d_sub_choose_death(const wavetree3d_sub_t *t, int depth, double u, int *coeff, double *prob);

int wavetree3d_sub_choose_death_global(const wavetree3d_sub_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

int wavetree3d_sub_reverse_death(const wavetree3d_sub_t *t, int depth, int coeff, double *prob);

int wavetree3d_sub_reverse_death_global(const wavetree3d_sub_t *t, int maxdepth, int depth, int coeff, double *prob);

/*
 * Functions for setting up a value proposal
 */
int wavetree3d_sub_choose_value_depth(const wavetree3d_sub_t *t, double u, int maxdepth, int *depth, double *prob);

int wavetree3d_sub_choose_value(const wavetree3d_sub_t *t, int depth, double u, int *coeff, double *prob);

int wavetree3d_sub_choose_value_global(const wavetree3d_sub_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

typedef int (*wavetree3d_sub_perturb_func_t)(void *user, 
					 int i, 
					 int j, 
					 int k,
					 int level,
					 int maxlevel,
					 double parent_coeff,
					 double *coeff,
					 double *prior_ratio);
int wavetree3d_sub_perturb(wavetree3d_sub_t *t,
		       wavetree3d_sub_perturb_func_t f,
		       void *user,
		       double *prior_ratio);

int wavetree3d_sub_parent_index(const wavetree3d_sub_t *t, int c);

int wavetree3d_sub_3dindices(const wavetree3d_sub_t *t,
			     int i,
			     int *ii,
			     int *ij,
			     int *ik);

int wavetree3d_sub_from_3dindices(const wavetree3d_sub_t *t,
				  int ii,
				  int ij,
				  int ik);

int wavetree3d_sub_depthofindex(const wavetree3d_sub_t *t,
				int i);

int wavetree3d_sub_get_coeff(const wavetree3d_sub_t *t,
			     int i,
			     int d,
			     double *coeff);

int wavetree3d_sub_generate_dyck_word(wavetree3d_sub_t *t, char *buffer, int maxlength);

int wavetree3d_sub_generate_dyck_binary(wavetree3d_sub_t *t, uint64_t *binary);

double wavetree3d_sub_logpriorprobability(const wavetree3d_sub_t *t,
					  wavetree_pp_t *pp);

/*
 * Internal functions exposed for testing
 */

int wavetree3d_sub_max_child_count(wavetree3d_sub_t *t);

int wavetree3d_sub_child_count(wavetree3d_sub_t *t, int index, int depth);

int wavetree3d_sub_child_indices(wavetree3d_sub_t *t, int index, int depth, int *indices, int *n, int nmax);

int wavetree3d_sub_UTL(const wavetree3d_sub_t *t, int i);
int wavetree3d_sub_UTR(const wavetree3d_sub_t *t, int i);
int wavetree3d_sub_UBL(const wavetree3d_sub_t *t, int i);
int wavetree3d_sub_UBR(const wavetree3d_sub_t *t, int i);
int wavetree3d_sub_LTL(const wavetree3d_sub_t *t, int i);
int wavetree3d_sub_LTR(const wavetree3d_sub_t *t, int i);
int wavetree3d_sub_LBL(const wavetree3d_sub_t *t, int i);
int wavetree3d_sub_LBR(const wavetree3d_sub_t *t, int i);


#endif /* wavetree3d_sub_h */
			  
