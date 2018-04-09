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

#ifndef subdivisiontree2d_h
#define subdivisiontree2d_h

#include "coefficient_histogram.h"

typedef struct _subdivisiontree2d subdivisiontree2d_t;

typedef enum {
  SUBDIVISION_BASIS_CONSTANT,
  SUBDIVISION_BASIS_PYRAMID,
  SUBDIVISION_BASIS_LANCZOS
} subdivision_basis_t;

subdivisiontree2d_t *
subdivisiontree2d_create(int max_depth,
			 double alpha,
			 double overlap,
			 subdivision_basis_t basis);

void
subdivisiontree2d_destroy(subdivisiontree2d_t *t);

int
subdivisiontree2d_save(const subdivisiontree2d_t *t,
		       const char *filename);

int 
subdivisiontree2d_load(subdivisiontree2d_t *t,
		       const char *filename);

int 
subdivisiontree2d_load_promote(subdivisiontree2d_t *t,
			       const char *filename);

int
subdivisiontree2d_get_width(subdivisiontree2d_t *t);

int
subdivisiontree2d_get_size(subdivisiontree2d_t *t);

int 
subdivisiontree2d_initialize(subdivisiontree2d_t *t,
			     double dc);

int subdivisiontree2d_coeff_count(const subdivisiontree2d_t *t);

int subdivisiontree2d_coeff_count_scan(const subdivisiontree2d_t *t);

int subdivisiontree2d_map_to_array(const subdivisiontree2d_t *t, 
				   double *a, 
				   int n);

int subdivisiontree2d_propose_to_array(const subdivisiontree2d_t *t,
				       double *a,
				       int n);

int subdivisiontree2d_revert_to_array(const subdivisiontree2d_t *t,
				      double *a,
				      int n);

int
subdivisiontree2d_propose_value(subdivisiontree2d_t *t,
				int i,
				int d,
				double value);

int
subdivisiontree2d_propose_birth(subdivisiontree2d_t *t,
				int i,
				int d,
				double value);

int 
subdivisiontree2d_propose_death(subdivisiontree2d_t *t,
				int i,
				int d,
				double *old_value);

int
subdivisiontree2d_commit(subdivisiontree2d_t *t);

int 
subdivisiontree2d_undo(subdivisiontree2d_t *t);

int
subdivisiontree2d_update_histogram(const subdivisiontree2d_t *t, coefficient_histogram_t *hist);

int 
subdivisiontree2d_depth_filter(void *wt, int subset, int index);

/*
 * Functions for setting up a birth proposal
 */
int 
subdivisiontree2d_choose_birth_depth(const subdivisiontree2d_t *t, double u, int maxdepth, int *depth, double *prob);

int
subdivisiontree2d_reverse_birth_depth(const subdivisiontree2d_t *t, int depth, int maxdepth, double *prob);

int 
subdivisiontree2d_choose_birth(const subdivisiontree2d_t *t, int depth, double u, int *coeff, double *prob);

int
subdivisiontree2d_choose_birth_global(const subdivisiontree2d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

int 
subdivisiontree2d_reverse_birth(const subdivisiontree2d_t *t, int depth, int coeff, double *prob);

int 
subdivisiontree2d_reverse_birth_global(const subdivisiontree2d_t *t, int maxdepth, int depth, int coeff, double *prob);

/*
 * Functions for setting up a death proposal
 */
int subdivisiontree2d_choose_death_depth(const subdivisiontree2d_t *t, double u, int maxdepth, int *depth, double *prob);

int subdivisiontree2d_reverse_death_depth(const subdivisiontree2d_t *t, int depth, int maxdepth, double *prob);

int 
subdivisiontree2d_choose_death(const subdivisiontree2d_t *t, int depth, double u, int *coeff, double *prob);

int
subdivisiontree2d_choose_death_global(const subdivisiontree2d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);

int 
subdivisiontree2d_reverse_death(const subdivisiontree2d_t *t, int depth, int coeff, double *prob);

int 
subdivisiontree2d_reverse_death_global(const subdivisiontree2d_t *t, int maxdepth, int depth, int coeff, double *prob);

/*
 * Functions for setting up a value proposal
 */

int 
subdivisiontree2d_choose_value_depth(const subdivisiontree2d_t *t, double u, int maxdepth, int *depth, double *prob);

int 
subdivisiontree2d_choose_value(const subdivisiontree2d_t *t, int depth, double u, int *idx, double *value);

int
subdivisiontree2d_choose_value_global(const subdivisiontree2d_t *t, double u, int maxdepth, int *depth, int *coeff, double *prob);


typedef int (*subdivisiontree2d_perturb_func_t)(void *user, 
						int i, 
						int j,
						int level,
						int maxlevel,
						double parent_coeff,
						double *coeff,
						double *prior_ratio);

int subdivisiontree2d_perturb(subdivisiontree2d_t *t,
			      subdivisiontree2d_perturb_func_t f,
			      void *user,
			      double *prior_ratio);

int 
subdivisiontree2d_parent_index(subdivisiontree2d_t *t, int c);

double 
subdivisiontree2d_dc(subdivisiontree2d_t *t);

int 
subdivisiontree2d_get_coeff(const subdivisiontree2d_t *t,
			    int i,
			    double *coeff);

int 
subdivisiontree2d_valid(subdivisiontree2d_t *t);

/*
 * Internal functions exposed for testing.
 */
int
subdivisiontree2d_total_coefficients(int max_depth);

int
subdivisiontree2d_depthofindex(const subdivisiontree2d_t *t,
			       int index);

int 
subdivisiontree2d_2dindices(const subdivisiontree2d_t *t,
			    int index,
			    int *ii,
			    int *ij,
			    int *rs);

int
subdivisiontree2d_rs_from_depth(int depth);

int
subdivisiontree2d_from_2dindices(const subdivisiontree2d_t *t,
				 int ii,
				 int ij,
				 int rs);

int
subdivisiontree2d_TL(const subdivisiontree2d_t *t,
		     int i);

int
subdivisiontree2d_TR(const subdivisiontree2d_t *t,
		     int i);

int
subdivisiontree2d_BL(const subdivisiontree2d_t *t,
		     int i);

int
subdivisiontree2d_BR(const subdivisiontree2d_t *t,
		     int i);

int
subdivisiontree2d_print(const subdivisiontree2d_t *t);

/*
 * Internal functions exposed for testing
 */
double basis_constant(double x, double y, double c, double overlap);
double basis_pyramid(double x, double y, double c, double overlap);
double basis_lanczos(double x, double y, double c, double overlap);

/*
 * Old functions
 */
#if 0
int
subdivisiontree2d_simplify_top_down(subdivisiontree2d_t *t);

int
subdivisiontree2d_simplify_bottom_up(subdivisiontree2d_t *t);

#endif

#endif /* subdivisiontree2d_h */
