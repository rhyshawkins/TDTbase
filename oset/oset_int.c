
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "oset_int.h"

#include "slog.h"

#define USE_BSEARCH

static const int INCREMENT = 1024;

struct _oset_int {
  int size;

  int n;
  int *s;
  int *depth;
};

static int oset_int_expand(oset_int_t *s);

static int oset_int_find_exact(const oset_int_t *s, 
			       int target,
			       int start,
			       int end);

static int oset_int_find_insertion_index(const oset_int_t *s,
					 int idx,
					 int start,
					 int end);

oset_int_t *
oset_int_create(void)
{
  oset_int_t *t;

  t = (oset_int_t *)malloc(sizeof(oset_int_t));
  if (t == NULL) {
    return NULL;
  }

  t->size = INCREMENT;
  t->n = 0;
  t->s = (int *)malloc(sizeof(int) * INCREMENT);
  if (t->s == NULL) {
    return NULL;
  }
  t->depth = (int *)malloc(sizeof(int) * INCREMENT);
  if (t->depth == NULL) {
    return NULL;
  }
  
  memset(t->s, 0, sizeof(int) * INCREMENT);
  memset(t->depth, 0, sizeof(int) * INCREMENT);

  return t;
}

void
oset_int_destroy(oset_int_t *s)
{
  if (s != NULL) {
    free(s->depth);
    free(s->s);
    free(s);
  }
}

int
oset_int_clear(oset_int_t *s)
{
  if (s == NULL) {
    return -1;
  }

  s->n = 0;
  memset(s->s, 0, sizeof(int) * s->size);
  memset(s->depth, 0, sizeof(int) * s->size);

  return 0;
}

int
oset_int_insert(oset_int_t *s, int i, int depth)
{
  int j;
  int ii;

  if (s->n == s->size) {
    if (oset_int_expand(s) < 0) {
      return -1;
    }
  }

#ifdef USE_BSEARCH

  ii = oset_int_find_insertion_index(s,
				     i,
				     0,
				     s->n - 1);

  if (ii < 0) {
    /* 
     * Element already present
     */
    return 0;
  }

  for (j = s->n; j > ii; j --) {
    s->s[j] = s->s[j - 1];
    s->depth[j] = s->depth[j - 1];
  }
  s->s[ii] = i;
  s->depth[ii] = depth;
  s->n ++;
  return 1;

#else

  if (s->n == 0 || i < s->s[0]) {
    /* Start add */
    for (j = s->n; j > 0; j --) {
      s->s[j] = s->s[j - 1];
      s->depth[j] = s->depth[j - 1];
    }
    s->s[0] = i;
    s->depth[0] = depth;
    s->n ++;

    return 1;
    
  } else if (i > s->s[s->n - 1]) {

    /* End add */
    s->s[s->n] = i;
    s->depth[s->n] = depth;
    s->n ++;

    return 1;

  } else {

    ii = -1;
    for (j = 0; j < s->n; j ++) {
      if (s->s[j] == i) {
	return 0;
      } else if (s->s[j] > i) {
	ii = j;
	break;
      }
    }

    if (ii >= 0) {
      for (j = s->n; j > ii; j --) {
	s->s[j] = s->s[j - 1];
	s->depth[j] = s->depth[j - 1];
      }
      s->s[ii] = i;
      s->depth[ii] = depth;
      s->n ++;

      return 1;
    }
  }

  ERROR("error inserting (%d %d)", i, s->n);
  oset_int_dump(s);
  return 0;

#endif
}

int 
oset_int_remove(oset_int_t *s, int i)
{
  int j;
  int di;

#ifdef USE_BSEARCH

  di = oset_int_find_exact(s, 
			   i,
			   0,
			   s->n - 1);

  if (di >= 0) { 
    for (j = di; j < (s->n - 1); j ++) {
      s->s[j] = s->s[j + 1];
      s->depth[j] = s->depth[j + 1];
    }
    s->n --;
    if (s->n < 0) {
      ERROR("error removing: %d, %d", i, s->n);
      return -1;
    }
    return 1;
  }

#else

  di = -1;
  for (j = 0; j < s->n; j ++) {
    if (s->s[j] == i) {
      di = j;
      break;
    } 

    if (s->s[j] > i) {
      return 0;
    }
  }
#endif

  return 0;
}

int
oset_int_count(const oset_int_t *s)
{
  return s->n;
}

int 
oset_int_nth_element(const oset_int_t *s, int n, int *idx, int *depth)
{
  if (n < 0 || n >= s->n) {
    return -1;
  }

  if (idx != NULL) {
    *idx = s->s[n];
  }
  if (depth != NULL) {
    *depth = s->depth[n];
  }

  return s->s[n];
}

int 
oset_int_is_element(const oset_int_t *s, int i)
{

#ifdef USE_BSEARCH

  return (oset_int_find_exact(s, 
			      i,
			      0, 
			      s->n - 1) >= 0);

#else

  for (j = 0; j < s->n; j ++) {
    if (s->s[j] == i) {
      return -1;
    }

    if (s->s[j] > i) {
      return 0;
    }
  }

  return 0;

#endif 
}

int
oset_int_clone(oset_int_t *dst, const oset_int_t *src)
{
  int i;

  if (src->n > dst->size) {
    /*
     * If necessary, increase the destination size to match the source. 
     */
    free(dst->s);
    dst->s = (int *)malloc(sizeof(int) * src->size);
    if (dst->s == NULL) {
      return -1;
    }
    dst->size = src->size;
  }

  dst->n = src->n;
  for (i = 0; i < src->n; i ++) {
    dst->s[i] = src->s[i];
    dst->depth[i] = src->depth[i];
  }

  return 0;
}

int 
oset_int_intersection(oset_int_t *dst, const oset_int_t *a, const oset_int_t *b)
{
  int ai;
  int bi;
  
  dst->n = 0;
  ai = 0;
  bi = 0;

  while (ai < a->n && bi < b->n) {

    if (a->s[ai] == b->s[bi]) {

      if (dst->n == dst->size) {
	if (oset_int_expand(dst) < 0) {
	  return -1;
	}
      }

      dst->s[dst->n] = a->s[ai];
      dst->depth[dst->n] = a->depth[ai];
      dst->n ++;

      ai ++;
      bi ++;
    } else if (a->s[ai] < b->s[bi]) {
      ai ++;
    } else {
      bi ++;
    }
  }

  return dst->n;
}

static int oset_int_expand(oset_int_t *s)
{
  int new_size;
  int *new_s;
  int *new_depth;
  int i;

  new_size = s->size + INCREMENT;
  new_s = (int *)malloc(sizeof(int) * new_size);
  if (new_s == NULL) {
    return -1;
  }

  new_depth = (int *)malloc(sizeof(int) * new_size);
  if (new_depth == NULL) {
    return -1;
  }

  for (i = 0; i < s->n; i ++) {
    new_s[i] = s->s[i];
    new_depth[i] = s->depth[i];
  }

  free(s->depth);
  free(s->s);
  s->s = new_s;
  s->depth = new_depth;
  s->size = new_size;

  return 0;
}

int
oset_int_inorder(const oset_int_t *s) 
{
  int i;

  for (i = 1; i < s->n; i ++) {
    if (s->s[i - 1] >= s->s[i]) {
      return 0;
    }
  }
  return -1;
}

int
oset_int_choose(const oset_int_t *s, double u, double *sum_weights, double *weight)
{
  *sum_weights = 1.0;
  return (int)(u * (double)s->n);
}

int 
oset_int_weighted_choose(const oset_int_t *s, 
			 double alpha,
			 double u, 
			 int maxdepth,
			 double *sum_weights, 
			 double *weight)
{
  int sum;
  int i;
  double w;
  double t;

  sum = oset_int_weighted_sum(s, alpha, maxdepth);

  t = 0.0;
  for (i = 0; i < s->n; i ++) {
    if (s->depth[i] <= maxdepth) {
      w = pow((double)s->depth[i], alpha);
      t += w/sum;
      if (t > u) {
	*weight = w;
	*sum_weights = (double)sum;
	return i;
      }
    }
  }

  return -1;
}

int 
oset_int_inverse_weighted_choose(const oset_int_t *s, 
				 double alpha, 
				 double u, 
				 int maxdepth,
				 double *sum_weights, 
				 double *weight)
{
  double sum;
  int i;
  double t;
  double w;

  sum = oset_int_inverse_weighted_sum(s, alpha, maxdepth);

  t = 0.0;
  for (i = 0; i < s->n; i ++) {
    if (s->depth[i] <= maxdepth) {
      w = pow((double)s->depth[i], -alpha);
      t += w/sum;
      if (t > u) {
	*weight = w;
	*sum_weights = (double)sum;

	/* printf(" %d %d: %d\n", (int)(u * (double)s->n), i, s->s[i]); */
	/* oset_int_dump(s); */
	return i;
      }
    }
  }

  return -1;
}

double
oset_int_weighted_sum(const oset_int_t *s, double alpha, int maxdepth)
{
  int i;
  double sum;

  sum = 0.0;
  for (i = 0; i < s->n; i ++) {
    if (s->depth[i] <= maxdepth) {
      sum += pow((double)s->depth[i], alpha);
    }
  }
  
  return sum;
}

double
oset_int_weight(const oset_int_t *s, double alpha, int i)
{
  int j;

  for (j = 0; j < s->n; j ++) {
    if (s->s[j] == i) {
      return pow((double)s->depth[j], alpha);
    }
  }

  return 0.0;
}

double
oset_int_inverse_weighted_sum(const oset_int_t *s, double alpha, int maxdepth)
{
  int i;
  double sum;

  sum = 0.0;
  for (i = 0; i < s->n; i ++) {
    if (s->depth[i] <= maxdepth) {
      sum += pow((double)s->depth[i], -alpha);
    }
  }
  
  return sum;
}

double
oset_int_inverse_weight(const oset_int_t *s, double alpha, int i)
{
  int j;

  for (j = 0; j < s->n; j ++) {
    if (s->s[j] == i) {
      return pow((double)s->depth[j], -alpha);
    }
  }

  return 0.0;
}

void
oset_int_dump(const oset_int_t *s)
{
  int i;

  for (i = 0; i < s->n; i ++) {
    printf("(%d,%d) ", s->s[i], (int)s->depth[i]);
  }
  printf("\n");
}

static int oset_int_find_exact(const oset_int_t *s, 
			       int target,
			       int start,
			       int end)
{
  /* Return index of element in s->s[] array that matches target */
  int c;

  if (start > end) {
    return -1;
  } else if (start == end) {
    
    if (s->s[start] == target) {
      return start;
    } else {
      return -1;
    }
  }

  c = (start + end)/2;

  if (s->s[c] < target) {
    return oset_int_find_exact(s, target, c + 1, end);
  } else if (s->s[c] > target) {
    return oset_int_find_exact(s, target, start, c - 1);
  } else {
    return c;
  }
}

static int oset_int_find_insertion_index(const oset_int_t *s,
					 int idx,
					 int start,
					 int end)
{
  int c;

  if (start > end) {
    return start;
  }

  if (start == end) {
    if (idx == s->s[start]) {
      return -1;
    } else if (idx < s->s[start]) {
      return start;
    } else {
      return start + 1;
    }
  }

  c = (start + end)/2;

  if (s->s[c] < idx) {
    return oset_int_find_insertion_index(s, idx, c + 1, end);
  } else if (s->s[c] > idx) {
    return oset_int_find_insertion_index(s, idx, start, c - 1);
  } else {
    return -1;
  }
}
