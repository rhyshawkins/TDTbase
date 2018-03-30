#ifndef ohist_64_h
#define ohist_64_h

#include <stdint.h>

typedef struct _ohist_64 ohist_64_t;

ohist_64_t *
ohist_64_create(void);

void
ohist_64_destroy(ohist_64_t *s);

int
ohist_64_clear(ohist_64_t *s);

int
ohist_64_insert(ohist_64_t *s, uint64_t idx, int k, int incr);

int
ohist_64_nelements(const ohist_64_t *s, int k);

int 
ohist_64_nth_element(const ohist_64_t *s, 
		     int k,
		     int n,
		     uint64_t *idx,
		     int *count);

void
ohist_64_dump(const ohist_64_t *s);

#endif /* ohist_64_h */

