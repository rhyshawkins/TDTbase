#ifndef ttree_int_h
#define ttree_int_h

typedef struct _ttree_int ttree_int_t;

ttree_int_t *
ttree_int_create(int maxk);

void
ttree_int_destroy(ttree_int_t *t);

int
ttree_int_insert(ttree_int_t *t, int k, const char *string, int incr);

int
ttree_int_get(ttree_int_t *t, int k, const char *string, int *count);

typedef int (*ttree_int_iterate_t)(void *user, const char *string, int count);
int
ttree_int_iterate(ttree_int_t *t, int k, ttree_int_iterate_t iterator, void *user);

#endif /* ttree_int_h */
