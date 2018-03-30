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
multiset_int_total_count(multiset_int_t *s);

int 
multiset_int_restricted_total_count(multiset_int_t *s, int maxdepth);

int 
multiset_int_depth_count(multiset_int_t *s, int depth);

int
multiset_int_nonempty_count(multiset_int_t *s, int maxdepth);

int 
multiset_int_is_element(multiset_int_t *s, int index, int depth);

int
multiset_int_nth_element(multiset_int_t *s, int depth, int i, int *index);

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
