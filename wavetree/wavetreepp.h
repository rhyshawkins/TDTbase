#ifndef wavetree_pp_h
#define wavetree_pp_h

#include <stdio.h>

#include "wavetree_prior.h"
#include "wavetree_value_proposal.h"
#include "wavetree_birth_proposal.h"
#include "coefficient_histogram.h"

typedef int (*wavetree_pp_move_proposal_t)(void *user,
					   int i1,
					   int j1,
					   int k1,
					   int level1,
					   int i2,
					   int j2, 
					   int k2,
					   int level2,
					   int maxlevel,
					   double coeff,
					   double *new_coeff,
					   double *prob);

struct _wavetree_move {
  void *user;
  wavetree_pp_move_proposal_t move;
};
typedef struct _wavetree_move wavetree_move_t;

struct _wavetree_pp {
  
  wavetree_prior_t *prior;
  wavetree_bd_t *bd;
  wavetree_value_t *value;
  wavetree_move_t move;

};
typedef struct _wavetree_pp wavetree_pp_t;

/*
 * Load a proposal from a file
 */
wavetree_pp_t *
wavetree_pp_load(const char *filename, unsigned long int seed, coefficient_histogram_t *histogram);

wavetree_prior_t *
wavetree_pp_load_prior(FILE *fp, unsigned long int seed);

wavetree_bd_t *
wavetree_pp_load_bd(FILE *fp, unsigned long int seed, wavetree_prior_t *prior);

wavetree_value_t *
wavetree_pp_load_value(FILE *fp, unsigned long int seed, coefficient_histogram_t *histogram);

/*
 * Helper functions
 */

int
wavetree_pp_value_init(wavetree_pp_t *w);

int
wavetree_pp_propose_value2d(wavetree_pp_t *w, 
			    int i, int j,
			    int level, int maxlevel, double parent_coeff,
			    double temperature,
			    double *coeff,
			    double *prior_ratio);

int 
wavetree_pp_value_error_count(wavetree_pp_t *w);

double
wavetree_pp_prior_probability2d(wavetree_pp_t *w,
				int i,
				int j,
				int depth,
				int maxdepth,
				double parent_coeff,
				double coeff);

int
wavetree_pp_prior_range2d(wavetree_pp_t *w,
			  int i,
			  int j,
			  int depth,
			  int maxdepth,
			  double parent_coeff,
			  double *vmin,
			  double *vmax);

int
wavetree_pp_birth2d(wavetree_pp_t *w,
		    int i,
		    int j,
		    int depth,
		    int maxdepth,
		    double parent_coeff,
		    double *coeff,
		    double *prob,
		    int *valid);

int
wavetree_pp_death2d(wavetree_pp_t *w,
		    int i,
		    int j,
		    int depth,
		    int maxdepth,
		    double parent_coeff,
		    double coeff,
		    double *prob);

int
wavetree_pp_propose_value3d(wavetree_pp_t *w, 
			    int i, int j, int k,
			    int level, int maxlevel, double parent_coeff,
			    double temperature,
			    double *coeff,
			    double *prior_ratio);

double
wavetree_pp_prior_probability3d(wavetree_pp_t *w,
				int i,
				int j,
				int k,
				int depth,
				int maxdepth,
				double parent_coeff,
				double coeff);

int
wavetree_pp_prior_range3d(wavetree_pp_t *w,
			  int i,
			  int j,
			  int k,
			  int depth,
			  int maxdepth,
			  double parent_coeff,
			  double *vmin,
			  double *vmax);

int
wavetree_pp_birth3d(wavetree_pp_t *w,
		    int i,
		    int j,
		    int k,
		    int depth,
		    int maxdepth,
		    double parent_coeff,
		    double *coeff,
		    double *prob,
		    int *valid);

int
wavetree_pp_death3d(wavetree_pp_t *w,
		    int i,
		    int j,
		    int k,
		    int depth,
		    int maxdepth,
		    double parent_coeff,
		    double coeff,
		    double *prob);

int
wavetree_pp_prior_sample3d(wavetree_pp_t *w,
			   int i,
			   int j,
			   int k,
			   int depth,
			   int maxdepth,
			   double parent_coeff,
			   double *coeff);

int
wavetree_pp_setscale(wavetree_pp_t *w,
		     double newscale,
		     double *oldscale);

/*
 * General cleanup routine.
 */
void
wavetree_pp_destroy(wavetree_pp_t *p);

/*
 * Globally uniform prior on coefficients, birth from prior, move with prior.
 */
wavetree_pp_t *
wavetree_pp_create_globally_uniform(double vmin,
				    double vmax,
				    double vstd,
				    int kmax,
				    unsigned long int seed);

/*
 * Globally uniform prior on coefficients, gaussian birth proposal, gaussian move proposal.
 */
wavetree_pp_t *
wavetree_pp_create_globally_uniform_with_gaussian_proposal(double vmin,
							   double vmax,
							   double vstd,
							   double bd_std,
							   double move_std,
							   int kmax,
							   unsigned long int seed);

void
wavetree_pp_free_globally_uniform_with_gaussian_proposal(wavetree_pp_t *w);

/*
 * Depth dependant prior on coefficients, birth from prior
 */
wavetree_pp_t *
wavetree_pp_create_depth_dependent_uniform(double *vmin,
					   double *vmax,
					   double *vstd,
					   int ndepths,
					   int kmax,
					   double move_std,
					   unsigned long int seed);

/*
 * Depth dependant prior on coefficients, depth dependant gaussian proposal
 */
wavetree_pp_t *
wavetree_pp_create_depth_dependent_uniform_with_gaussian_proposal(double *vmin,
								  double *vmax,
								  double *vstd,
								  double *bd_std,
								  int ndepths,
								  int kmax,
								  double move_std,
								  unsigned long int seed);

/*
 * Depth dependent generalised gaussian, birth from prior
 */
wavetree_pp_t *
wavetree_pp_create_depth_dependent_generalised_gaussian(double *va,
							double beta,
							double *vstd,
							int ndepths,
							int kmax,
							double move_std,
							unsigned long int seed);

/*
 * Depth dependent generalised gaussian, gaussian birth
 */
wavetree_pp_t *
wavetree_pp_create_depth_dependent_generalised_gaussian_with_gaussian_proposal(double *va,
									       double beta,
									       double *vstd,
									       double *bd_std,
									       int ndepths,
									       int kmax,
									       double move_std,
									       unsigned long int seed);



#endif /* wavetree_pp_h */
