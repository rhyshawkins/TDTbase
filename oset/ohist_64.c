
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ohist_64.h"

#include "slog.h"

#define USE_BSEARCH

static const int K_INCREMENT = 32;
static const int INCREMENT = 1024;

struct _ohist_64 {
  int k_size;

  int *set_n;
  int *set_size;
  uint64_t **s;
  int **count;
};

static int ohist_64_expand(ohist_64_t *s, int k);

static int ohist_64_find_insertion_index(const ohist_64_t *s,
					 uint64_t idx,
					 int k,
					 int start,
					 int end);

ohist_64_t *
ohist_64_create(void)
{
  ohist_64_t *t;
  int i;

  t = (ohist_64_t *)malloc(sizeof(ohist_64_t));
  if (t == NULL) {
    return NULL;
  }

  t->k_size = K_INCREMENT;

  t->set_n = (int*)malloc(sizeof(int) * t->k_size);
  if (t->set_n == NULL) {
    return NULL;
  }
  memset(t->set_n, 0, sizeof(int) * t->k_size);

  t->set_size = (int*)malloc(sizeof(int) * t->k_size);
  if (t->set_size == NULL) {
    return NULL;
  }
  for (i = 0; i < t->k_size; i ++) {
    t->set_size[i] = INCREMENT;
  }

  t->s = (uint64_t **)malloc(sizeof(uint64_t*) * t->k_size);
  if (t->s == NULL) {
    return NULL;
  }
  t->count = (int **)malloc(sizeof(int*) * t->k_size);
  if (t->count == NULL) {
    return NULL;
  }

  for (i = 0; i < t->k_size; i ++) {
    t->s[i] = (uint64_t*)malloc(sizeof(uint64_t) * INCREMENT);
    if (t->s[i] == NULL) {
      return NULL;
    }
    memset(t->s[i], 0, sizeof(uint64_t) * INCREMENT);

    t->count[i] = (int*)malloc(sizeof(int) * INCREMENT);
    if (t->count[i] == NULL) {
      return NULL;
    }
    memset(t->count[i], 0, sizeof(int) * INCREMENT);
  }
      
  return t;
}

void
ohist_64_destroy(ohist_64_t *s)
{
  int i;

  if (s != NULL) {
    for (i = 0; i < s->k_size; i ++) {
      free(s->count[i]);
      free(s->s[i]);
    }
    free(s->count);
    free(s->s);
    free(s->set_n);
    free(s->set_size);

    free(s);
  }
}

int
ohist_64_clear(ohist_64_t *s)
{
  int i;

  if (s == NULL) {
    return -1;
  }

  for (i = 0; i < s->k_size; i ++) {
    s->set_n[i] = 0;
    memset(s->s[i], 0, sizeof(uint64_t) * s->set_size[i]);
    memset(s->count[i], 0, sizeof(int) * s->set_size[i]);
  }

  return 0;
}

int
ohist_64_insert(ohist_64_t *s, uint64_t idx, int k, int incr)
{
  int j;
  int ii;

  if (k < 0 || k >= s->k_size) {
    ERROR("k out of range %d", k);
    return -1;
  }

  if (s->set_n[k] == s->set_size[k]) {
    if (ohist_64_expand(s, k) < 0) {
      ERROR("failed to expand arrays");
      return -1;
    }
  }

#ifdef USE_BSEARCH

  ii = ohist_64_find_insertion_index(s,
				     idx,
				     k,
				     0,
				     s->set_n[k] - 1);

  if (ii < 0) {
    /*
     * Error
     */
    ERROR("ohist_64_insert: failed to find insertion point\n");
    return -1;
  }

  if (s->s[k][ii] == idx) {
    /* 
     * Element already present
     */
    s->count[k][ii] += incr;
    return 0;

  } else {
    for (j = s->set_n[k]; j > ii; j --) {
      s->s[k][j] = s->s[k][j - 1];
      s->count[k][j] = s->count[k][j - 1];
    }
    s->s[k][ii] = idx;
    s->count[k][ii] = incr;
    s->set_n[k] ++;
    return 1;
  }
#else

  if (s->n == 0 || i < s->s[0]) {
    /* Start add */
    for (j = s->n; j > 0; j --) {
      s->s[j] = s->s[j - 1];
      s->count[j] = s->count[j - 1];
    }
    s->s[0] = i;
    s->count[0] = count;
    s->n ++;

    return 1;
    
  } else if (i > s->s[s->n - 1]) {

    /* End add */
    s->s[s->n] = i;
    s->count[s->n] = count;
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
	s->count[j] = s->count[j - 1];
      }
      s->s[ii] = i;
      s->count[ii] = count;
      s->n ++;

      return 1;
    }
  }

  ERROR("error inserting (%d %d)", i, s->n);
  ohist_64_dump(s);
  return 0;

#endif
}


int
ohist_64_nelements(const ohist_64_t *s, int k)
{
  if (k < 0 || k >= s->k_size) {
    return -1;
  }

  return s->set_n[k];
}

int 
ohist_64_nth_element(const ohist_64_t *s, int k, int n, uint64_t *idx, int *count)
{
  if (k < 0 || k >= s->k_size) {
    return -1;
  }

  if (n < 0 || n >= s->set_n[k]) {
    return -1;
  }

  if (idx != NULL) {
    *idx = s->s[k][n];
  }
  if (count != NULL) {
    *count = s->count[k][n];
  }

  return 0;
}


static int ohist_64_expand(ohist_64_t *s, int k)
{
  int new_size;
  uint64_t *new_s;
  int *new_count;
  int i;

  if (k < 0 || k >= s->k_size) {
    return -1;
  }

  new_size = s->set_size[k] + INCREMENT;
  new_s = (uint64_t *)malloc(sizeof(uint64_t) * new_size);
  if (new_s == NULL) {
    return -1;
  }

  new_count = (int *)malloc(sizeof(int) * new_size);
  if (new_count == NULL) {
    return -1;
  }

  for (i = 0; i < s->set_n[k]; i ++) {
    new_s[i] = s->s[k][i];
    new_count[i] = s->count[k][i];
  }

  free(s->count[k]);
  free(s->s[k]);
  s->s[k] = new_s;
  s->count[k] = new_count;
  s->set_size[k] = new_size;

  return 0;
}


void
ohist_64_dump(const ohist_64_t *s)
{
  int i;
  int k;

  for (k = 0; k < s->k_size; k ++) {
    printf("%d: %d: ", k, s->set_n[k]);
    for (i = 0; i < s->set_n[k]; i ++) {
      printf("(0x%lx,%d) ", s->s[k][i], (int)s->count[k][i]);
    }
    printf("\n");
  }
}

static int ohist_64_find_insertion_index(const ohist_64_t *s,
					 uint64_t idx,
					 int k,
					 int start,
					 int end)
{
  int c;

  if (start > end) {
    return start;
  }

  if (start == end) {
    if (idx == s->s[k][start]) {
      return start;
    } else if (idx < s->s[k][start]) {
      return start;
    } else {
      return start + 1;
    }
  }

  c = (start + end)/2;

  if (s->s[k][c] < idx) {
    return ohist_64_find_insertion_index(s, idx, k, c + 1, end);
  } else if (s->s[k][c] > idx) {
    return ohist_64_find_insertion_index(s, idx, k, start, c - 1);
  } else if (s->s[k][c] == idx) {
    return c;
  } else {
    return -1;
  }
}
