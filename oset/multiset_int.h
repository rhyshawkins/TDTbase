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
#ifndef multiset_int_h
#define multiset_int_h

#include <stdio.h>

typedef struct _multiset_int multiset_int_t;

multiset_int_t *
multiset_int_create(void);

void
multiset_int_destroy(multiset_int_t *s);

int
multiset_int_clear(multiset_int_t *s);

int
multiset_int_insert(multiset_int_t *s, int index, int depth);

int
multiset_int_remove(multiset_int_t *s, int index, int depth);

int
multiset_int_total_count(const multiset_int_t *s);

int 
multiset_int_restricted_total_count(const multiset_int_t *s, int maxdepth);

int 
multiset_int_depth_count(const multiset_int_t *s, int depth);

int
multiset_int_nonempty_count(const multiset_int_t *s, int maxdepth);

int 
multiset_int_is_element(const multiset_int_t *s, int index, int depth);

int
multiset_int_nth_element(const multiset_int_t *s, int depth, int i, int *index);

int 
multiset_int_choose_depth(multiset_int_t *s,
			  double u,
			  int maxdepth,
			  int *depth,
			  int *ndepths);

int
multiset_int_choose_index(multiset_int_t *s,
			  int depth,
			  double u,
			  int *index,
			  int *nindices);

int
multiset_int_choose_index_globally(multiset_int_t *s,
				   double u,
				   int maxdepth,
				   int *index,
				   int *ndepths,
				   int *nindices);

int 
multiset_int_choose_index_weighted(multiset_int_t *s,
				   double u,
				   int maxdepth,
				   double depthweight,
				   int *index,
				   int *depth,
				   double *prob);

int
multiset_int_reverse_choose_index_weighted(multiset_int_t *s,
					   int maxdepth,
					   double depthweight,
					   int index,
					   int depth,
					   double *prob);

void
multiset_int_dump(const multiset_int_t *s);

int
multiset_int_write(const multiset_int_t *s, FILE *fp);

int 
multiset_int_read(multiset_int_t *s, FILE *fp);

#endif /* multiset_int_h */
