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
