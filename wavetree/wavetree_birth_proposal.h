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
