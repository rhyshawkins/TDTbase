#ifndef wavetree_birth_proposal_h
#define wavetree_birth_proposal_h

#include "wavetree_prior.h"

typedef int (*wavetree_pp_birth_proposal_t)(void *user,
					    int i,
					    int j,
					    int k,
					    int level,
					    int maxlevel,
					    double parent_coeff,
					    double *coeff,
					    double *prob);

typedef int (*wavetree_pp_death_proposal_t)(void *user,
					    int i,
					    int j,
					    int k,
					    int level,
					    int maxlevel,
					    double parent_coeff,
					    double coeff,
					    double *prob);

struct _wavetree_bd {
  void *user;
  wavetree_pp_birth_proposal_t birth;
  wavetree_pp_death_proposal_t death;
  wavetree_pp_destroy_t destroy;
};
typedef struct _wavetree_bd wavetree_bd_t;

/*
 * General destroy function
 */
void
wavetree_birth_destroy(wavetree_bd_t *p);

wavetree_bd_t *
wavetree_birth_create_birth_from_prior(wavetree_prior_t *p);

/*
 * Gaussian birth proposal uniform for all coefficients
 */
wavetree_bd_t *
wavetree_birth_create_gaussian(wavetree_prior_t *p,
			       double std,
			       unsigned long int seed);

void
wavetree_birth_free_gaussian(wavetree_bd_t *w);

/*
 * Gaussian birth proposal with std deviation depending on depth
 */
wavetree_bd_t *
wavetree_birth_create_depth_gaussian(wavetree_prior_t *p,
				     int ndepths,
				     double *std,
				     unsigned long int seed);

void
wavetree_birth_free_depth_gaussian(wavetree_bd_t *w);


#endif /* wavetree_birth_proposal_h */
