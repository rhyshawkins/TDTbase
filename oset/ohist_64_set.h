#ifndef ohist_64_h
#define ohist_64_h

#include <stdint.h>

typedef struct _ohist_64_set ohist_64_set_t;

ohist_64_set_t *
ohist_64_set_create(int maxk);

void
ohist_64_set_destroy(ohist_64_set_t *s);

int
oset_clear(ohist_64_set_t *s);

int
ohist_64_set_insert(ohist_64_set_t *s, 
		    uint64_t idx, 
		    int k, 
		    const int *set,
		    int incr,
		    int *pset,
		    int *matching);

int
ohist_64_set_nelements(const ohist_64_set_t *s, int k);

int 
ohist_64_set_nth_element(const ohist_64_set_t *s, 
			 int k,
			 int n,
			 uint64_t *idx,
			 int *count);

void
ohist_64_set_dump(const ohist_64_set_t *s);

#endif /* ohist_64_set_h */

