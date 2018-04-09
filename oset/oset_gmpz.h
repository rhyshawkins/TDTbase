//
//    Ordered set library, used for maintaining arbitrary tree based models
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
#ifndef oset_gmpz_h
#define oset_gmpz_h

#include <gmp.h>

typedef struct _oset_gmpz oset_gmpz_t;

oset_gmpz_t *
oset_gmpz_create(void);

void
oset_gmpz_delete(oset_gmpz_t *s);

int
oset_gmpz_insert(oset_gmpz_t *s, mpz_t i, int depth);

int 
oset_gmpz_remove(oset_gmpz_t *s, mpz_t i);

int
oset_gmpz_count(const oset_gmpz_t *s);

int
oset_gmpz_nth_element(const oset_gmpz_t *s, int n, mpz_t e, int *depth);

int 
oset_gmpz_is_element(const oset_gmpz_t *s, mpz_t i);

int
oset_gmpz_clone(oset_gmpz_t *dst, const oset_gmpz_t *src);

int 
oset_gmpz_intersection(oset_gmpz_t *dst, const oset_gmpz_t *a, const oset_gmpz_t *b);

int
oset_gmpz_inorder(const oset_gmpz_t *s);

int 
oset_gmpz_weighted_choose(const oset_gmpz_t *s, 
			  double alpha, 
			  double u, 
			  double *sum_weights, 
			  double *weight);

double
oset_gmpz_weighted_sum(const oset_gmpz_t *s, 
		       double alpha);

double
oset_gmpz_weight(const oset_gmpz_t *s, 
		 double alpha, 
		 mpz_t i);

double
oset_gmpz_inverse_weighted_sum(const oset_gmpz_t *s, 
			       double alpha);

double
oset_gmpz_inverse_weight(const oset_gmpz_t *s, 
			 double alpha, 
			 mpz_t i);

int 
oset_gmpz_inverse_weighted_choose(const oset_gmpz_t *s, 
				  double alpha, 
				  double u, 
				  double *sum_weights, 
				  double *weight);

void
oset_gmpz_dump(const oset_gmpz_t *s);

#endif /* oset_gmpz_h */

