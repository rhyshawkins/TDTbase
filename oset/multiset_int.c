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

#include "multiset_int.h"

#include "slog.h"

static const int DEPTH_INCREMENT = 16;
static const int SET_INCREMENT = 1024;

struct _multiset_int {
  int depth_size;

  int *set_size;
  int *set_n;
  int **s;
};

static int multiset_int_expand_depth(multiset_int_t *s, int depth);
static int multiset_int_expand_set(multiset_int_t *s, int depth);

static int multiset_int_find_exact(int *s, 
				   int target,
				   int start,
				   int end);

static int multiset_int_find_insertion_index(int *s,
					     int idx,
					     int start,
					     int end);

multiset_int_t *
multiset_int_create(void)
{
  multiset_int_t *s;
  int i;

  s = malloc(sizeof(multiset_int_t));
  if (s == NULL) {
    ERROR("failed to allocate memory");
    return NULL;
  }
  
  s->depth_size = DEPTH_INCREMENT;

  s->set_size = malloc(sizeof(int) * DEPTH_INCREMENT);
  if (s->set_size == NULL) {
    ERROR("failed to allocate set size");
    return NULL;
  }
  for (i = 0; i < s->depth_size; i ++) {
    s->set_size[i] = SET_INCREMENT;
  }

  s->set_n = malloc(sizeof(int) * DEPTH_INCREMENT);
  if (s->set_n == NULL) {
    ERROR("failed to allocate set n");
    return NULL;
  }
  for (i = 0; i < s->depth_size; i ++) {
    s->set_n[i] = 0;
  }

  s->s = malloc(sizeof(int *) * DEPTH_INCREMENT);
  if (s->s == NULL) {
    ERROR("failed to allocate set array");
    return NULL;
  }

  for (i = 0; i < s->depth_size; i ++) {
    s->s[i] = malloc(sizeof(int) * SET_INCREMENT);
    if (s->s[i] == NULL) {
      ERROR("failed to allocate set");
      return NULL;
    }

    memset(s->s[i], 0, sizeof(int) * SET_INCREMENT);
  }

  return s;
}


void
multiset_int_destroy(multiset_int_t *s)
{
  int i;

  if (s != NULL) {

    for (i = 0; i < s->depth_size; i ++) {
      free(s->s[i]);
    }
    free(s->s);

    free(s->set_n);
    free(s->set_size);
    free(s);
  }
}

int
multiset_int_clear(multiset_int_t *s)
{
  int i;
  
  if (s == NULL) {
    return -1;
  }

  for (i = 0; i < s->depth_size; i ++) {
    s->set_n[i] = 0;
    memset(s->s[i], 0, sizeof(int) * s->set_size[i]);
  }

  return 0;
}

int
multiset_int_insert(multiset_int_t *s, int index, int depth)
{
  int ii;
  int j;

  if (depth >= s->depth_size) {
    if (multiset_int_expand_depth(s, depth) < 0) {
      return -1;
    }
  }

  if (s->set_n[depth] == s->set_size[depth]) {
    if (multiset_int_expand_set(s, depth) < 0) {
      return -1;
    }
  }

  ii = multiset_int_find_insertion_index(s->s[depth], 
					 index,
					 0,
					 s->set_n[depth] - 1);
  /* Negative means index was already present */
  if (ii < 0) {
    return 0;
  }

  /* Shift tail indices back 1 to make room */
  for (j = s->set_n[depth]; j > ii; j --) {
    s->s[depth][j] = s->s[depth][j - 1];
  }

  s->s[depth][ii] = index;
  s->set_n[depth] ++;

  return 1;
}

int
multiset_int_remove(multiset_int_t *s, int index, int depth)
{
  int di;
  int j;

  if (depth < 0 || depth >= s->depth_size) {
    return -1;
  }

  di = multiset_int_find_exact(s->s[depth], 
			       index,
			       0,
			       s->set_n[depth] - 1);

  if (di >= 0) { 
    for (j = di; j < (s->set_n[depth] - 1); j ++) {
      s->s[depth][j] = s->s[depth][j + 1];
    }
    s->set_n[depth] --;
    return 1;
  }

  return 0;
}

int
multiset_int_total_count(multiset_int_t *s)
{
  int i;
  int c;

  if (s == NULL) {
    return -1;
  }
  
  c = 0;
  for (i = 0; i < s->depth_size; i ++) {
    c += s->set_n[i];
  }

  return c;
}

int 
multiset_int_restricted_total_count(multiset_int_t *s, int maxdepth)
{
  int i;
  int c;
  int depthlimit;

  c = 0;
  depthlimit = s->depth_size - 1;
  if (maxdepth >= 0 && maxdepth < depthlimit) {
    depthlimit = maxdepth;
  }
  for (i = 0; i <= depthlimit; i ++) {
    c += s->set_n[i];
  }

  return c;
}

int 
multiset_int_depth_count(multiset_int_t *s, int depth)
{
  if (s == NULL ||
      depth < 0 ||
      depth >= s->depth_size) {
    return -1;
  }

  return s->set_n[depth];
}

int
multiset_int_nonempty_count(multiset_int_t *s, int maxdepth)
{
  int i;
  int c;
  int depthlimit;

  c = 0;
  depthlimit = s->depth_size - 1;
  if (maxdepth >= 0 && maxdepth < depthlimit) {
    depthlimit = maxdepth;
  }
  for (i = 0; i <= depthlimit; i ++) {
    if (s->set_n[i] > 0) {
      c ++;
    }
  }

  return c;
}

int 
multiset_int_is_element(multiset_int_t *s, int index, int depth)
{
  if (s == NULL ||
      depth < 0 ||
      depth >= s->depth_size) {
    return 0;
  }

  if (s->set_n[depth] == 0) {
    return 0;
  }

  return (multiset_int_find_exact(s->s[depth], 
				  index,
				  0, 
				  s->set_n[depth] - 1) >= 0);
}

int
multiset_int_nth_element(multiset_int_t *s, int depth, int i, int *index)
{
  if (depth < 0 || depth >= s->depth_size) {
    return -1;
  }
  if (i < 0 || i >= s->set_n[depth]) {
    return -1;
  }

  *index = s->s[depth][i];
  return 0;
}

int 
multiset_int_choose_depth(multiset_int_t *s,
			  double u,
			  int maxdepth,
			  int *depth,
			  int *ndepths)
{
  int i;
  int j;
  int cdepths;
  int depthlimit;

  cdepths = 0;

  depthlimit = s->depth_size - 1;
  if (maxdepth >= 0 && maxdepth < depthlimit) {
    depthlimit = maxdepth;
  }

  for (i = 0; i <= depthlimit; i ++) {
    if (s->set_n[i] > 0) {
      cdepths ++;
    }
  }
  if (cdepths == 0) {
    ERROR("no non-empty sets");
    return -1;
  }

  j = (int)(u * (double)cdepths);
  for (i = 0; i < s->depth_size; i ++) {
    if (s->set_n[i] > 0) {
      if (j == 0) {
	*depth = i;
	*ndepths = cdepths;
	return 0;
      } else { 
	j --;
      }
    }
  }

  ERROR("failed to choose");
  return -1;
}

int
multiset_int_choose_index(multiset_int_t *s,
			  int depth,
			  double u,
			  int *index,
			  int *nindices)
{
  int j;

  if (s == NULL ||
      depth < 0 ||
      depth >= s->depth_size) {
    return -1;
  }  

  if (s->set_n[depth] == 0) {
    return -1;
  }
  
  j = (int)(u * (double)s->set_n[depth]);
  *index = s->s[depth][j];
  *nindices = s->set_n[depth];

  return 0;
}

int
multiset_int_choose_index_globally(multiset_int_t *s,
				   double u,
				   int maxdepth,
				   int *index,
				   int *depth,
				   int *nindices)
{
  int depthlimit;
  int cindices;
  int i;
  int j;

  cindices = 0;

  depthlimit = s->depth_size - 1;
  if (maxdepth >= 0 && maxdepth < depthlimit) {
    depthlimit = maxdepth;
  }

  for (i = 0; i <= depthlimit; i ++) {
    cindices += s->set_n[i];
  }

  if (cindices == 0) {
    return -1;
  }
  
  j = (int)(u * (double)cindices);

  for (i = 0; i <= depthlimit; i ++) {

    if (s->set_n[i] > j) {
      *index = s->s[i][j];
      *depth = i;
      *nindices = cindices;
      return 0;
    } else {
      j -= s->set_n[i];
    }
  }

  ERROR("failed to recover index %d %f", cindices, u);
  return -1;
}

int 
multiset_int_choose_index_weighted(multiset_int_t *s,
				   double u,
				   int maxdepth,
				   double depthweight,
				   int *index,
				   int *depth,
				   double *prob)
{
  int depthlimit;
  int i;
  int j;
  double sum;
  double v;
  double dv;

  sum = 0.0;

  depthlimit = s->depth_size - 1;
  if (maxdepth >= 0 && maxdepth < depthlimit) {
    depthlimit = maxdepth;
  }

  for (i = 0; i <= depthlimit; i ++) {
    sum += (double)s->set_n[i] * pow((double)(i + 1), depthweight);
  }

  v = sum * u;

  for (i = 0; i <= depthlimit; i ++) {
    dv = (double)s->set_n[i] * pow((double)(i + 1), depthweight);
    if (v < dv) {
      /* This depth */
      j = (int)(v/dv * (double)s->set_n[i]);
      *index = s->s[i][j];
      *depth = i;
      *prob = pow((double)(i + 1), depthweight)/sum;

      return 0;
    } else {
      v -= dv;
    }
  }

  return -1;
}

int
multiset_int_reverse_choose_index_weighted(multiset_int_t *s,
					   int maxdepth,
					   double depthweight,
					   int index,
					   int depth,
					   double *prob)
{
  int depthlimit;
  int i;
  double sum;


  if (multiset_int_is_element(s, index, depth)) {

    sum = 0.0;
    depthlimit = s->depth_size - 1;
    if (maxdepth >= 0 && maxdepth < depthlimit) {
      depthlimit = maxdepth;
    }
    
    for (i = 0; i <= depthlimit; i ++) {
      sum += (double)s->set_n[i] * pow((double)(i + 1), depthweight);
    }
    
    *prob = pow((double)(depth + 1), depthweight)/sum;

    return 0;
  }

  return -1;
}

void
multiset_int_dump(const multiset_int_t *s)
{
  int d;
  int i;

  for (d = 0; d < s->depth_size; d ++) {
    if (s->set_n[d] > 0) {
      printf("Depth %d:\n  ", d);
      for (i = 0; i < s->set_n[d]; i ++) {
	printf("%d ", s->s[d][i]);
      }
      printf("\n");
    }
  }
}

int
multiset_int_write(const multiset_int_t *s, FILE *fp)
{
  int d;
  int i;

  fprintf(fp, "%d\n", s->depth_size);

  for (d = 0; d < s->depth_size; d ++) {
   
    fprintf(fp, "%d %d\n", d, s->set_n[d]);

    for (i = 0; i < s->set_n[d]; i ++) {
      fprintf(fp, "%d\n", s->s[d][i]);
    }
  }

  return 0;
}

int 
multiset_int_read(multiset_int_t *s, FILE *fp)
{
  int maxd;
  int d;
  int j;

  int di;
  int dc;

  multiset_int_clear(s);

  if (fscanf(fp, "%d\n", &maxd) != 1) {
    ERROR("failed to read no. depths");
    return -1;
  }

  if (multiset_int_expand_depth(s, maxd - 1) < 0) {
    ERROR("failed to expand depth");
    return -1;
  }

  for (d = 0; d < maxd; d ++) {

    if (fscanf(fp, "%d %d", &di, &dc) != 2) {
      ERROR("failed to read depth header");
      return -1;
    }

    if (di != d) {
      ERROR("depth mismatch");
      return -1;
    }

    s->set_n[d] = 0;
    for (j = 0; j < dc; j ++, s->set_n[d] ++) {
      if (j == s->set_size[d]) {
	if (multiset_int_expand_set(s, d) < 0) {
	  ERROR("failed to expand set");
	  return -1;
	}
      }

      if (fscanf(fp, "%d\n", &(s->s[d][j])) != 1) {
	ERROR("failed to read set element");
	return -1;
      }
    }
  }

  return 0;
    
}

/*
 * Internal functions
 */

static int multiset_int_expand_depth(multiset_int_t *s, int depth)
{
  /*
   * Increase no. of depths to at least depth + 1
   */
  int newdepth_size;
  int *newset_size;
  int *newset_n;
  int **new_s;

  int i;

  newdepth_size = s->depth_size;
  while (newdepth_size <= depth) {
    newdepth_size += DEPTH_INCREMENT;
  }

  if (newdepth_size == s->depth_size) {
    /* No expansion necessary */
    return 0;
  }

  /* Temporarily disable depth expansion */
  printf("error: newdepth size = %d, current = %d, requested %d\n", newdepth_size, s->depth_size, depth);
  return -1;

  newset_size = malloc(sizeof(int) * newdepth_size);
  if (newset_size == NULL) {
    ERROR("failed to allocate new set size");
    return -1;
  }

  newset_n = malloc(sizeof(int) * newdepth_size);
  if (newset_n == NULL) {
    ERROR("failed to allocate new set n");
    return -1;
  }

  new_s = malloc(sizeof(int*) * newdepth_size);
  if (new_s == NULL) {
    ERROR("failed to allocate new set array");
    return -1;
  }

  for (i = 0; i < s->depth_size; i ++) {
    newset_size[i] = s->set_size[i];
    newset_n[i] = s->set_n[i];
    new_s[i] = s->s[i];
  }

  free(s->set_size);
  free(s->set_n);
  free(s->s);

  s->set_size = newset_size;
  s->set_n = newset_n;
  s->s = new_s;
  s->depth_size = newdepth_size;

  return 0;
}

static int multiset_int_expand_set(multiset_int_t *s, int depth)
{
  int new_size;
  int *new_s;
  int i;

  new_size = s->set_size[depth] + SET_INCREMENT;
  new_s = (int *)malloc(sizeof(int) * new_size);
  if (new_s == NULL) {
    ERROR("failed to allocate new set");
    return -1;
  }
  memset(new_s, 0, sizeof(int) * new_size);

  for (i = 0; i < s->set_n[depth]; i ++) {
    new_s[i] = s->s[depth][i];
  }

  free(s->s[depth]);
  s->s[depth] = new_s;
  s->set_size[depth] = new_size;

  return 0;

}

static int multiset_int_find_exact(int *s, 
				   int target,
				   int start,
				   int end)
{
  /* Return index of element in s->s[] array that matches target */
  int c;

  if (start > end) {
    return -1;
  } else if (start == end) {
    
    if (s[start] == target) {
      return start;
    } else {
      return -1;
    }
  }

  c = (start + end)/2;

  if (s[c] < target) {
    return multiset_int_find_exact(s, target, c + 1, end);
  } else if (s[c] > target) {
    return multiset_int_find_exact(s, target, start, c - 1);
  } else {
    return c;
  }
}

static int multiset_int_find_insertion_index(int *s,
					     int idx,
					     int start,
					     int end)
{
  int c;

  if (start > end) {
    return start;
  }

  if (start == end) {
    if (idx == s[start]) {
      return -1;
    } else if (idx < s[start]) {
      return start;
    } else {
      return start + 1;
    }
  }

  c = (start + end)/2;

  if (s[c] < idx) {
    return multiset_int_find_insertion_index(s, idx, c + 1, end);
  } else if (s[c] > idx) {
    return multiset_int_find_insertion_index(s, idx, start, c - 1);
  } else {
    return -1;
  }
}

