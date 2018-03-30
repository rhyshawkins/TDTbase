#ifndef wavetree_prior_h
#define wavetree_prior_h

typedef int (*wavetree_pp_priorrange_t)(void *user,
				       int i,
				       int j,
				       int k,
				       int level,
				       int maxlevel,
				       double parent_coeff,
				       double *vmin,
				       double *vmax);

typedef int (*wavetree_pp_priorsample_t)(void *user,
					 int i,
					 int j,
					 int k,
					 int level,
					 int maxlevel,
					 double parent_coeff,
					 double *coeff);

typedef double (*wavetree_pp_priorprobability_t)(void *user,
						 int i,
						 int j,
						 int k,
						 int level,
						 int maxlevel,
						 double parent_coeff,
						 double coeff);

typedef int (*wavetree_pp_valid_t)(void *user,
				   int i,
				   int j,
				   int k,
				   int level,
				   int maxlevel,
				   double parent_coeff,
				   double coeff);

typedef int (*wavetree_pp_setscale_t)(void *user,
				      double scale,
				      double *oldscale);

typedef int (*wavetree_pp_destroy_t)(void *user);

struct _wavetree_prior {
  void *user;
  wavetree_pp_priorrange_t range;
  wavetree_pp_priorsample_t sample;
  wavetree_pp_priorprobability_t prob;
  wavetree_pp_valid_t valid;
  wavetree_pp_setscale_t setscale;
  wavetree_pp_destroy_t destroy;
};
typedef struct _wavetree_prior wavetree_prior_t;


/*
 * General cleanup routine.
 */
void
wavetree_prior_destroy(wavetree_prior_t *p);


/*
 * Coefficients uniform across entire transform space
 */

wavetree_prior_t *
wavetree_prior_create_globally_uniform(double vmin,
				       double vmax,
				       unsigned long int seed);

/*
 * Coefficients Laplacian across entire transform space
 */

wavetree_prior_t *
wavetree_prior_create_globally_laplacian(double beta,
					 unsigned long int seed);

/*
 * Coefficients uniform at each resolution, for depths greater than max depth
 * we use the last depth prior.
 */
wavetree_prior_t*
wavetree_prior_create_depth_uniform(int ndepths,
				    double *vmin,
				    double *vmax,
				    unsigned long int seed);


/*
 * Coefficients generalised gaussian at each depth
 */
wavetree_prior_t*
wavetree_prior_create_depth_generalised_gaussian(int ndepths,
						 double *va,
						 double beta,
						 unsigned long int seed);

#endif /* wavetree_prior_h */
