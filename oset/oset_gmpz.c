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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "oset_gmpz.h"

#include "slog.h"

static const int INCREMENT = 1024;

struct _oset_gmpz {
  int size;

  int n;
  mpz_t *s;
  int *depth;
};

static int oset_gmpz_expand(oset_gmpz_t *s);

oset_gmpz_t *
oset_gmpz_create(void)
{
  oset_gmpz_t *t;
  int i;

  t = (oset_gmpz_t *)malloc(sizeof(oset_gmpz_t));
  if (t == NULL) {
    return NULL;
  }

  t->size = INCREMENT;
  t->n = 0;
  t->s = (mpz_t *)malloc(sizeof(mpz_t) * INCREMENT);
  if (t->s == NULL) {
    return NULL;
  }
  t->depth = (int *)malloc(sizeof(int) * INCREMENT);
  if (t->depth == NULL) {
    return NULL;
  }
  
  for (i = 0; i < t->size; i ++) {
    mpz_init(t->s[i]);
  }
  memset(t->depth, 0, sizeof(int) * INCREMENT);

  return t;
}

void
oset_gmpz_delete(oset_gmpz_t *s)
{
  int i;

  if (s != NULL) {
    free(s->depth);
    
    for (i = 0; i < s->size; i ++) {
      mpz_clear(s->s[i]); 
    }
    free(s->s);

    free(s);
  }
}

int
oset_gmpz_insert(oset_gmpz_t *s, mpz_t i, int depth)
{
  int j;
  int ii;
  int c;

  if (mpz_sgn(i) < 0) {
    gmp_fprintf(stderr, "oset_gmpz_insert: invalid index %Zd\n", i);
    return -1;
  }

  if (depth <= 0) {
    ERROR("invalid depth %d", depth);
    return -1;
  }

  if (s->n == s->size) {
    if (oset_gmpz_expand(s) < 0) {
      return -1;
    }
  }

  if (s->n == 0 || mpz_cmp(i, s->s[s->n - 1]) > 0) {

    /* End add */
    mpz_set(s->s[s->n], i);
    s->depth[s->n] = depth;
    s->n ++;

    return 1;

  } else if (mpz_cmp(i, s->s[0]) < 0) {

    /* Start add */
    for (j = s->n; j > 0; j --) {
      mpz_set(s->s[j], s->s[j - 1]);
      s->depth[j] = s->depth[j - 1];
    }
    mpz_set(s->s[0], i);
    s->depth[0] = depth;
    s->n ++;

    return 1;
    
  } else {

    ii = -1;
    for (j = 0; j < s->n; j ++) {
      
      c = mpz_cmp(s->s[j], i);
      if (c == 0) {
	return 0;
      } else if (c > 0) {
	ii = j;
	break;
      }
    }

    if (ii >= 0) {
      for (j = s->n; j > ii; j --) {
	mpz_set(s->s[j], s->s[j - 1]);
	s->depth[j] = s->depth[j - 1];
      }
      mpz_set(s->s[ii], i);
      s->depth[ii] = depth;
      s->n ++;

      return 1;
    }

  }


  ERROR("error inserting");
  return 0;
  
}

int 
oset_gmpz_remove(oset_gmpz_t *s, mpz_t i)
{
  int j;
  int di;
  int c;

  di = -1;
  for (j = 0; j < s->n; j ++) {
    c = mpz_cmp(s->s[j], i);
    if (c == 0) {
      di = j;
      break;
    } 

    if (c > 0) {
      ERROR("warning: asked to remove non-existing element");
      return 0;
    }
  }

  if (di >= 0) { 
    for (j = di; j < (s->n - 1); j ++) {
      mpz_set(s->s[j], s->s[j + 1]);
      s->depth[j] = s->depth[j + 1];
    }
    s->n --;
    return 1;
  }

  gmp_fprintf(stderr, "oset_gmpz_remove: warning: element to remove not found %Zd\n", i);
  oset_gmpz_dump(s);
  return 0;
}

int
oset_gmpz_count(const oset_gmpz_t *s)
{
  return s->n;
}

int
oset_gmpz_nth_element(const oset_gmpz_t *s, int n, mpz_t e, int *depth)
{
  if (n < 0 || n >= s->n) {
    ERROR("index out of range");
    return -1;
  }

  mpz_set(e, s->s[n]);
  if (depth != NULL) {
    *depth = s->depth[n];
  }
  return 0;
}

int 
oset_gmpz_is_element(const oset_gmpz_t *s, mpz_t i)
{
  int j;
  int c;

  for (j = 0; j < s->n; j ++) {
    c = mpz_cmp(i, s->s[j]);
    if (c == 0) {
      return -1;
    }

    if (c < 0) {
      return 0;
    }
  }

  return 0;
}

int
oset_gmpz_clone(oset_gmpz_t *dst, const oset_gmpz_t *src)
{
  int i;

  dst->n = 0;
  while (src->n > dst->size) {
    /*
     * If necessary, increase the destination size to match the source. 
     */
    if (oset_gmpz_expand(dst) < 0) {
      ERROR("failed to expand destination to fit source");
      return -1;
    }
  }

  dst->n = src->n;
  for (i = 0; i < src->n; i ++) {
    mpz_set(dst->s[i], src->s[i]);
    dst->depth[i] = src->depth[i];
  }

  return 0;
}

int 
oset_gmpz_intersection(oset_gmpz_t *dst, const oset_gmpz_t *a, const oset_gmpz_t *b)
{
  int ai;
  int bi;
  int c;
  
  dst->n = 0;
  ai = 0;
  bi = 0;

  while (ai < a->n && bi < b->n) {

    c = mpz_cmp(a->s[ai], b->s[bi]);

    if (c == 0) {

      if (dst->n == dst->size) {
	if (oset_gmpz_expand(dst) < 0) {
	  return -1;
	}
      }

      mpz_set(dst->s[dst->n], a->s[ai]);
      dst->depth[dst->n] = a->depth[ai];
      dst->n ++;

      ai ++;
      bi ++;
    } else if (c < 0) {
      ai ++;
    } else {
      bi ++;
    }
  }

  return dst->n;
}

static int oset_gmpz_expand(oset_gmpz_t *s)
{
  int new_size;
  mpz_t *new_s;
  int *new_depth;
  int i;

  new_size = s->size + INCREMENT;
  new_s = (mpz_t *)malloc(sizeof(mpz_t) * new_size);
  if (new_s == NULL) {
    return -1;
  }
  for (i = 0; i < new_size; i ++) {
    mpz_init(new_s[i]);
  }

  new_depth = (int *)malloc(sizeof(int) * new_size);
  if (new_depth == NULL) {
    return -1;
  }

  for (i = 0; i < s->n; i ++) {
    mpz_set(new_s[i], s->s[i]);
    new_depth[i] = s->depth[i];
  }

  free(s->depth);
  
  for (i = 0; i < new_size; i ++) {
    mpz_clear(s->s[i]);
  }
  free(s->s);

  s->s = new_s;
  s->depth = new_depth;
  s->size = new_size;

  return 0;
}

int
oset_gmpz_inorder(const oset_gmpz_t *s) 
{
  int i;
  int c;

  for (i = 1; i < s->n; i ++) {
    c = mpz_cmp(s->s[i - 1], s->s[i]);
    if (c >= 0) {
      return 0;
    }
  }
  return -1;
}

int 
oset_gmpz_weighted_choose(const oset_gmpz_t *s, double alpha, double u, double *sum_weights, double *weight)
{
  int sum;
  int i;
  double t;
  double w;

  sum = oset_gmpz_weighted_sum(s, alpha);

  t = 0.0;
  for (i = 0; i < s->n; i ++) {
    w = pow((double)s->depth[i], alpha);
    t += w/sum;
    if (t > u) {
      *weight = w;
      *sum_weights = sum;
      return i;
    }
  }

  return -1;
}

int 
oset_gmpz_inverse_weighted_choose(const oset_gmpz_t *s, double alpha, double u, double *sum_weights, double *weight)
{
  double sum;
  int i;
  double t;
  double w;

  sum = oset_gmpz_inverse_weighted_sum(s, alpha);

  t = 0.0;
  for (i = 0; i < s->n; i ++) {
    w = pow((double)s->depth[i], -alpha);
    t += w/sum;
    if (t > u) {
      *weight = w;
      *sum_weights = sum;
      return i;
    }
  }

  return -1;
}

double
oset_gmpz_weighted_sum(const oset_gmpz_t *s, double alpha)
{
  int i;
  double sum;

  sum = 0.0;
  for (i = 0; i < s->n; i ++) {
    sum += pow((double)s->depth[i], alpha);
  }
  
  return sum;
}

double
oset_gmpz_weight(const oset_gmpz_t *s, double alpha, mpz_t i)
{
  int j;

  for (j = 0; j < s->n; j ++) {
    if (mpz_cmp(s->s[j], i) == 0) {
      return pow((double)s->depth[j], alpha);
    }
  }

  return 0.0;
}

double
oset_gmpz_inverse_weighted_sum(const oset_gmpz_t *s, double alpha)
{
  int i;
  double sum;

  sum = 0.0;
  for (i = 0; i < s->n; i ++) {
    sum += pow((double)s->depth[i], -alpha);
  }
  
  return sum;
}

double
oset_gmpz_inverse_weight(const oset_gmpz_t *s, double alpha, mpz_t i)
{
  int j;

  for (j = 0; j < s->n; j ++) {
    if (mpz_cmp(s->s[j], i) == 0) {
      return pow((double)s->depth[j], -alpha);
    }
  }

  return 0.0;
}

void
oset_gmpz_dump(const oset_gmpz_t *s)
{
  int i;

  for (i = 0; i < s->n; i ++) {
    gmp_printf("(%Zd,%d) ", s->s[i], s->depth[i]);
  }
  printf("\n");
}
