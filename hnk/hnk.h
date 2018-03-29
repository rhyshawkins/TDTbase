//
//    HNK Library : A library for computing combinations of arrangements of
//    general trees for the Trans-dimensional Tree algorithm. See
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
#ifndef hnk_h
#define hnk_h
				
#include <gmp.h>

typedef struct _hnk hnk_t;

typedef int (*hnk_compute_maxk_at_h_t)(hnk_t *t,
				       int h,
				       int nroot);


typedef int (*hnk_compute_hnk_t)(hnk_t *t,
				 int h,
				 int k,
				 mpz_t hnk);

hnk_t *
hnk_create(int maxh,
	   int maxk,
	   int maxsplit,
	   hnk_compute_maxk_at_h_t cmaxk,
	   hnk_compute_hnk_t chnk,
	   hnk_t *subtree);

hnk_t *
hnk_create_aggregate(int maxh,
		     int maxk,
		     int maxsplit,
		     int nsubtrees,
		     ...);

hnk_t *
hnk_create_aggregate_empty(int maxh,
			   int maxk,
			   int maxsplit,
			   int nsubtrees);

int
hnk_aggregate_set_subtree(hnk_t *t,
			  int i,
			  hnk_t *subtree);

int hnk_aggregate_maxk_at_h(hnk_t *t,
			    int h,
			    int nroot);

int hnk_aggregate_hnk(hnk_t *t,
		      int h,
		      int k,
		      mpz_t hnk);

int hnk_save(hnk_t *t,
	     const char *filename);

int hnk_restore(hnk_t *t,
		const char *filename);

void
hnk_destroy(hnk_t *t);

int
hnk_get_maxk_at_h(hnk_t *t, 
		  int h);

int
hnk_get_maxk_at_h_storage(hnk_t *t, 
			  int h);

int
hnk_get_hnk(hnk_t *t, 
	    int h, 
	    int k, 
	    mpz_t hnk);

int
hnk_is_hnk_memoized(hnk_t *t,
		    int h,
		    int k);

int
hnk_highest_memoized_k(hnk_t *t,
		       int h);

/*
 * Returns ratio of hnk(k+1)/hnk(k) for given k
 */
int
hnk_get_kplus1_ratio(hnk_t *t, 
		     int h, 
		     int k, 
		     double *ratio);

int
hnk_get_subtree_maxk_at_h(hnk_t *t, 
			  int h);

int
hnk_get_subtree_hnk(hnk_t *t, 
		    int h, 
		    int k, 
		    mpz_t hnk);

hnk_t *
hnk_get_subtree(hnk_t *t);

int
hnk_save_text(hnk_t *t, 
	      const char *filename);
	   
hnk_t *
hnk_load_text(const char *filename, 
	      hnk_compute_maxk_at_h_t cmaxk,
	      hnk_compute_hnk_t chnk,
	      hnk_t *subtree);
  

int
hnk_general_split(hnk_t *t,
		  hnk_t *subtree,
		  int h,
		  int k,
		  int nsplit,
		  mpz_t count);

int
hnk_aggregate_split(hnk_t *t,
		    hnk_t **subtree,
		    int nsubtree,
		    int index,
		    int h,
		    int k,
		    int nsplit,
		    mpz_t count);

int
hnk_naggregatesplits(int nsubtree);

#endif /* */
