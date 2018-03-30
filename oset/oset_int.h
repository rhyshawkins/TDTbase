#ifndef oset_int_h
#define oset_int_h

typedef struct _oset_int oset_int_t;

oset_int_t *
oset_int_create(void);

void
oset_int_destroy(oset_int_t *s);

int
oset_int_clear(oset_int_t *s);

int
oset_int_insert(oset_int_t *s, int i, int depth);

int 
oset_int_remove(oset_int_t *s, int i);

int
oset_int_count(const oset_int_t *s);

int 
oset_int_nth_element(const oset_int_t *s, 
		     int n,
		     int *idx,
		     int *depth);

int 
oset_int_is_element(const oset_int_t *s, 
		    int i);

int
oset_int_clone(oset_int_t *dst, 
	       const oset_int_t *src);

int 
oset_int_intersection(oset_int_t *dst, 
		      const oset_int_t *a, 
		      const oset_int_t *b);

int
oset_int_inorder(const oset_int_t *s);

int 
oset_int_weighted_choose(const oset_int_t *s, 
			 double alpha,
			 double u, 
			 int maxdepth,
			 double *sum_weights,
			 double *weight);

double
oset_int_weighted_sum(const oset_int_t *s,
		      double alpha,
		      int maxdepth);

double
oset_int_weight(const oset_int_t *s, double alpha, int i);

double
oset_int_inverse_weighted_sum(const oset_int_t *s,
			      double alpha,
			      int maxdepth);

double
oset_int_inverse_weight(const oset_int_t *s, 
			double alpha,
			int i);

int 
oset_int_inverse_weighted_choose(const oset_int_t *s, 
				 double alpha,
				 double u,
				 int maxweight,
				 double *sum_weights, 
				 double *weight);

void
oset_int_dump(const oset_int_t *s);

#endif /* oset_int_h */

