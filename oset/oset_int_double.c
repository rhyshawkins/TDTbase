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

#include "oset_int_double.h"

#include "slog.h"

#define USE_BSEARCH

static const int INCREMENT = 1024;

struct _oset_int_double {
  int size;

  int n;
  int *s;
  double *v;
};

static int oset_int_double_expand(oset_int_double_t *s);

static int oset_int_double_find_exact(const oset_int_double_t *s, 
			       int target,
			       int start,
			       int end);

static int oset_int_double_find_insertion_index(const oset_int_double_t *s,
					 int idx,
					 int start,
					 int end);

oset_int_double_t *
oset_int_double_create(void)
{
  oset_int_double_t *t;

  t = (oset_int_double_t *)malloc(sizeof(oset_int_double_t));
  if (t == NULL) {
    return NULL;
  }

  t->size = INCREMENT;
  t->n = 0;
  t->s = (int *)malloc(sizeof(int) * INCREMENT);
  if (t->s == NULL) {
    return NULL;
  }
  t->v = (double *)malloc(sizeof(double) * INCREMENT);
  if (t->v == NULL) {
    return NULL;
  }
  
  memset(t->s, 0, sizeof(int) * INCREMENT);
  memset(t->v, 0, sizeof(double) * INCREMENT);

  return t;
}

void
oset_int_double_destroy(oset_int_double_t *s)
{
  if (s != NULL) {
    free(s->v);
    free(s->s);
    free(s);
  }
}

int
oset_int_double_clear(oset_int_double_t *s)
{
  if (s == NULL) {
    return -1;
  }

  s->n = 0;
  memset(s->s, 0, sizeof(int) * s->size);
  memset(s->v, 0, sizeof(double) * s->size);

  return 0;
}

int
oset_int_double_insert(oset_int_double_t *s,
		       int i,
		       double value,
		       oset_int_double_multiple_insertion_action_t insert_action)
{
  int j;
  int ii;

  if (s->n == s->size) {
    if (oset_int_double_expand(s) < 0) {
      return -1;
    }
  }

  ii = oset_int_double_find_insertion_index(s,
				     i,
				     0,
				     s->n - 1);

  if (ii < 0) {

    /* 
     * Element already present
     */
    switch (insert_action) {
    case OSET_INT_DOUBLE_OVERWRITE:
    case OSET_INT_DOUBLE_SUM:
      ii = oset_int_double_find_exact(s, 
				      i,
				      0,
				      s->n - 1);
      if (ii < 0) {
	ERROR("Failed to find element");
	return -1;
      }

      if (insert_action == OSET_INT_DOUBLE_OVERWRITE) {
	s->v[ii] = value;
      } else {
	s->v[ii] += value;
      }
      break;

    default:
      break;
    }

    return 0;
  }

  for (j = s->n; j > ii; j --) {
    s->s[j] = s->s[j - 1];
    s->v[j] = s->v[j - 1];
  }
  s->s[ii] = i;
  s->v[ii] = value;
  s->n ++;
  return 1;

}

int 
oset_int_double_remove(oset_int_double_t *s, int i)
{
  int j;
  int di;

  di = oset_int_double_find_exact(s, 
			   i,
			   0,
			   s->n - 1);

  if (di >= 0) { 
    for (j = di; j < (s->n - 1); j ++) {
      s->s[j] = s->s[j + 1];
      s->v[j] = s->v[j + 1];
    }
    s->n --;
    if (s->n < 0) {
      ERROR("error removing: %d, %d", i, s->n);
      return -1;
    }
    return 1;
  }

  return 0;
}

int
oset_int_double_count(const oset_int_double_t *s)
{
  return s->n;
}

int 
oset_int_double_nth_element(const oset_int_double_t *s, int n, int *idx, double *value)
{
  if (n < 0 || n >= s->n) {
    return -1;
  }

  if (idx != NULL) {
    *idx = s->s[n];
  }
  if (value != NULL) {
    *value = s->v[n];
  }

  return s->s[n];
}

int 
oset_int_double_is_element(const oset_int_double_t *s, int i)
{

  return (oset_int_double_find_exact(s, 
			      i,
			      0, 
			      s->n - 1) >= 0);

}

void
oset_int_double_dump(const oset_int_double_t *s)
{
  int i;

  for (i = 0; i < s->n; i ++) {
    printf("(%d,%f) ", s->s[i], s->v[i]);
  }
  printf("\n");
}

static int
oset_int_double_expand(oset_int_double_t *s)
{
  int new_size;
  int *new_s;
  double *new_v;
  int i;

  new_size = s->size + INCREMENT;
  new_s = (int *)malloc(sizeof(int) * new_size);
  if (new_s == NULL) {
    return -1;
  }

  new_v = (double *)malloc(sizeof(double) * new_size);
  if (new_v == NULL) {
    return -1;
  }

  for (i = 0; i < s->n; i ++) {
    new_s[i] = s->s[i];
    new_v[i] = s->v[i];
  }

  free(s->v);
  free(s->s);
  s->s = new_s;
  s->v = new_v;
  s->size = new_size;

  return 0;
}

static int
oset_int_double_find_exact(const oset_int_double_t *s, 
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
    return oset_int_double_find_exact(s, target, c + 1, end);
  } else if (s->s[c] > target) {
    return oset_int_double_find_exact(s, target, start, c - 1);
  } else {
    return c;
  }
}

static int
oset_int_double_find_insertion_index(const oset_int_double_t *s,
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
    return oset_int_double_find_insertion_index(s, idx, c + 1, end);
  } else if (s->s[c] > idx) {
    return oset_int_double_find_insertion_index(s, idx, start, c - 1);
  } else {
    return -1;
  }
}
