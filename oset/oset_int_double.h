#ifndef oset_int_double_h
#define oset_int_double_h

typedef struct _oset_int_double oset_int_double_t;

typedef enum {
  OSET_INT_DOUBLE_IGNORE = 0,
  OSET_INT_DOUBLE_OVERWRITE,
  OSET_INT_DOUBLE_SUM
} oset_int_double_multiple_insertion_action_t;

oset_int_double_t *
oset_int_double_create(void);

void
oset_int_double_destroy(oset_int_double_t *s);

int
oset_int_double_clear(oset_int_double_t *s);

int
oset_int_double_insert(oset_int_double_t *s,
		       int i,
		       double value,
		       oset_int_double_multiple_insertion_action_t insert_action);

int 
oset_int_double_remove(oset_int_double_t *s, int i);

int
oset_int_double_count(const oset_int_double_t *s);

int 
oset_int_double_nth_element(const oset_int_double_t *s, 
			    int n,
			    int *idx,
			    double *value);
int 
oset_int_is_element(const oset_int_double_t *s, int i);

void
oset_int_dump(const oset_int_double_t *s);

#endif /* oset_int_double_h */

