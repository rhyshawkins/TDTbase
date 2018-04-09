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
#ifndef multiset_int_double_h
#define multiset_int_double_h

#include <stdio.h>

typedef struct _multiset_int_double multiset_int_double_t;

multiset_int_double_t *
multiset_int_double_create(void);

void
multiset_int_double_destroy(multiset_int_double_t *s);

int
multiset_int_double_clear(multiset_int_double_t *s);

int
multiset_int_double_clone(multiset_int_double_t *dest,
			  const multiset_int_double_t *src);

int
multiset_int_double_insert(multiset_int_double_t *s, int index, int depth, double value);

int
multiset_int_double_get(const multiset_int_double_t *s, int index, int depth, double *value);

int
multiset_int_double_set(multiset_int_double_t *s, int index, int depth, double value);

int
multiset_int_double_remove(multiset_int_double_t *s, int index, int depth);

int
multiset_int_double_total_count(const multiset_int_double_t *s);

int 
multiset_int_double_restricted_total_count(const multiset_int_double_t *s, int maxdepth);

int 
multiset_int_double_depth_count(const multiset_int_double_t *s, int depth);

int
multiset_int_double_nonempty_count(const multiset_int_double_t *s, int maxdepth);

int 
multiset_int_double_is_element(const multiset_int_double_t *s, int index, int depth);

int
multiset_int_double_nth_element(const multiset_int_double_t *s, int depth, int i, int *index, double *value);

int 
multiset_int_double_choose_depth(const multiset_int_double_t *s,
				 double u,
				 int maxdepth,
				 int *depth,
				 int *ndepths);

int
multiset_int_double_choose_index(const multiset_int_double_t *s,
				 int depth,
				 double u,
				 int *index,
				 int *nindices);

int
multiset_int_double_choose_index_globally(const multiset_int_double_t *s,
					  double u,
					  int maxdepth,
					  int *index,
					  int *ndepths,
					  int *nindices);

int 
multiset_int_double_choose_index_weighted(const multiset_int_double_t *s,
					  double u,
					  int maxdepth,
					  double depthweight,
					  int *index,
					  int *depth,
					  double *prob);

int
multiset_int_double_reverse_choose_index_weighted(const multiset_int_double_t *s,
						  int maxdepth,
						  double depthweight,
						  int index,
						  int depth,
						  double *prob);

void
multiset_int_double_dump(const multiset_int_double_t *s);

int
multiset_int_double_write(const multiset_int_double_t *s, FILE *fp);

int 
multiset_int_double_read(multiset_int_double_t *s, FILE *fp);

typedef size_t (*multiset_int_double_write_t)(const void *p, size_t size, size_t nmemb, void *fp);
typedef size_t (*multiset_int_double_read_t)(void *p, size_t size, size_t nmemb, void *fp);

int
multiset_int_double_write_binary(multiset_int_double_t *s,
				 multiset_int_double_write_t write_function,
				 void *fp);

int
multiset_int_double_read_binary(multiset_int_double_t *s,
				multiset_int_double_read_t read_function,
				void *fp);


#endif /* multiset_int_double_h */
