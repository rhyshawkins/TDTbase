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

#include "multiset_int_double.h"

#include "slog.h"

static const int DEPTH_INCREMENT = 16;
static const int SET_INCREMENT = 1024;

struct _multiset_int_double {
  int depth_size;

  int *set_size;
  int *set_n;
  int **s;
  double **v;
};

static int multiset_int_double_expand_set(multiset_int_double_t *s, int depth, int minsize);

static int multiset_int_double_find_exact(int *s, 
					  int target,
					  int start,
					  int end);

static int multiset_int_double_find_insertion_index(int *s,
						    int idx,
						    int start,
						    int end);

multiset_int_double_t *
multiset_int_double_create(void)
{
  multiset_int_double_t *s;
  int i;

  s = malloc(sizeof(multiset_int_double_t));
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

  s->v = malloc(sizeof(double *) * DEPTH_INCREMENT);
  if (s->v == NULL) {
    ERROR("failed to allocate value array");
    return NULL;
  }

  for (i = 0; i < s->depth_size; i ++) {
    s->v[i] = malloc(sizeof(double) * SET_INCREMENT);
    if (s->v[i] == NULL) {
      ERROR("failed to allocate set");
      return NULL;
    }

    memset(s->v[i], 0, sizeof(double) * SET_INCREMENT);
  }


  return s;
}


void
multiset_int_double_destroy(multiset_int_double_t *s)
{
  int i;

  if (s != NULL) {

    for (i = 0; i < s->depth_size; i ++) {
      free(s->v[i]);
    }
    free(s->v);

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
multiset_int_double_clear(multiset_int_double_t *s)
{
  int i;
  
  if (s == NULL) {
    return -1;
  }

  for (i = 0; i < s->depth_size; i ++) {
    s->set_n[i] = 0;
    memset(s->s[i], 0, sizeof(int) * s->set_size[i]);
    memset(s->v[i], 0, sizeof(double) * s->set_size[i]);
  }

  return 0;
}

int
multiset_int_double_clone(multiset_int_double_t *dest,
			  const multiset_int_double_t *src)
{
  int d;
  int j;

  if (dest->depth_size != src->depth_size) {
    ERROR("depth size mismatch");
    return -1;
  }

  for (d = 0; d < dest->depth_size; d ++) {

    if (dest->set_size[d] < src->set_size[d]) {
      if (multiset_int_double_expand_set(dest, d, src->set_size[d]) < 0) {
	ERROR("failed to expand dest set");
	return -1;
      }
    }

    for (j = 0; j < src->set_n[d]; j ++) {
      dest->s[d][j] = src->s[d][j];
      dest->v[d][j] = src->v[d][j];
    }
    dest->set_n[d] = src->set_n[d];
  }

  return 0;
}

int
multiset_int_double_insert(multiset_int_double_t *s, int index, int depth, double value)
{
  int ii;
  int j;

  if (depth >= s->depth_size) {
    return -1;
  }

  if (s->set_n[depth] == s->set_size[depth]) {
    if (multiset_int_double_expand_set(s, depth, s->set_n[depth] + 1) < 0) {
      return -1;
    }
  }

  ii = multiset_int_double_find_insertion_index(s->s[depth], 
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
    s->v[depth][j] = s->v[depth][j - 1];
  }

  s->s[depth][ii] = index;
  s->v[depth][ii] = value;
  s->set_n[depth] ++;

  return 1;
}

int
multiset_int_double_get(const multiset_int_double_t *s, int index, int depth, double *value)
{
  int di;

  di = multiset_int_double_find_exact(s->s[depth], 
				      index,
				      0,
				      s->set_n[depth] - 1);
  
  if (di >= 0) {
    *value = s->v[depth][di];
    return 0;
  }

  return -1;
}

int
multiset_int_double_set(multiset_int_double_t *s, int index, int depth, double value)
{
  int di;

  di = multiset_int_double_find_exact(s->s[depth], 
			       index,
			       0,
			       s->set_n[depth] - 1);

  if (di >= 0) {
    s->v[depth][di] = value;
    return 0;
  }

  return -1;
}

int
multiset_int_double_remove(multiset_int_double_t *s, int index, int depth)
{
  int di;
  int j;

  if (depth < 0 || depth >= s->depth_size) {
    return -1;
  }

  di = multiset_int_double_find_exact(s->s[depth], 
			       index,
			       0,
			       s->set_n[depth] - 1);

  if (di >= 0) { 
    for (j = di; j < (s->set_n[depth] - 1); j ++) {
      s->s[depth][j] = s->s[depth][j + 1];
      s->v[depth][j] = s->v[depth][j + 1];
    }
    s->set_n[depth] --;
    return 1;
  }

  return 0;
}

int
multiset_int_double_total_count(const multiset_int_double_t *s)
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
multiset_int_double_restricted_total_count(const multiset_int_double_t *s, int maxdepth)
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
multiset_int_double_depth_count(const multiset_int_double_t *s, int depth)
{
  if (s == NULL ||
      depth < 0 ||
      depth >= s->depth_size) {
    return -1;
  }

  return s->set_n[depth];
}

int
multiset_int_double_nonempty_count(const multiset_int_double_t *s, int maxdepth)
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
multiset_int_double_is_element(const multiset_int_double_t *s, int index, int depth)
{
  if (s == NULL ||
      depth < 0 ||
      depth >= s->depth_size) {
    return 0;
  }

  if (s->set_n[depth] == 0) {
    return 0;
  }

  return (multiset_int_double_find_exact(s->s[depth], 
				  index,
				  0, 
				  s->set_n[depth] - 1) >= 0);
}

int
multiset_int_double_nth_element(const multiset_int_double_t *s, int depth, int i, int *index, double *value)
{
  if (depth < 0 || depth >= s->depth_size) {
    return -1;
  }
  if (i < 0 || i >= s->set_n[depth]) {
    return -1;
  }

  *index = s->s[depth][i];
  *value = s->v[depth][i];
  return 0;
}

int 
multiset_int_double_choose_depth(const multiset_int_double_t *s,
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
multiset_int_double_choose_index(const multiset_int_double_t *s,
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
multiset_int_double_choose_index_globally(const multiset_int_double_t *s,
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
multiset_int_double_choose_index_weighted(const multiset_int_double_t *s,
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
multiset_int_double_reverse_choose_index_weighted(const multiset_int_double_t *s,
						  int maxdepth,
						  double depthweight,
						  int index,
						  int depth,
						  double *prob)
{
  int depthlimit;
  int i;
  double sum;


  if (multiset_int_double_is_element(s, index, depth)) {

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
multiset_int_double_dump(const multiset_int_double_t *s)
{
  int d;
  int i;

  for (d = 0; d < s->depth_size; d ++) {
    if (s->set_n[d] > 0) {
      printf("Depth %d:\n  ", d);
      for (i = 0; i < s->set_n[d]; i ++) {
	printf("(%d, %f) ", s->s[d][i], s->v[d][i]);
      }
      printf("\n");
    }
  }
}

int
multiset_int_double_write(const multiset_int_double_t *s, FILE *fp)
{
  int d;
  int i;

  fprintf(fp, "%d\n", s->depth_size);

  for (d = 0; d < s->depth_size; d ++) {
   
    fprintf(fp, "%d %d\n", d, s->set_n[d]);

    for (i = 0; i < s->set_n[d]; i ++) {
      fprintf(fp, "%d %.9g\n", s->s[d][i], s->v[d][i]);
    }
  }

  return 0;
}

int 
multiset_int_double_read(multiset_int_double_t *s, FILE *fp)
{
  int maxd;
  int d;
  int j;

  int di;
  int dc;

  multiset_int_double_clear(s);

  if (fscanf(fp, "%d\n", &maxd) != 1) {
    ERROR("failed to read no. depths");
    return -1;
  }

  if (maxd > s->depth_size) {
    ERROR("maximum depth exceeded %d > %d", maxd, s->depth_size);
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

    if (dc > s->set_size[d]) {
      if (multiset_int_double_expand_set(s, d, dc) < 0) {
	ERROR("failed to expand set");
	return -1;
      }
    }
    
    for (j = 0; j < dc; j ++) {

      if (fscanf(fp, "%d %lf\n", &(s->s[d][j]), &(s->v[d][j])) != 2) {
	ERROR("failed to read set element");
	return -1;
      }
    }

    s->set_n[d] = dc;
  }

  return 0;
    
}
int
multiset_int_double_write_binary(multiset_int_double_t *s,
				 multiset_int_double_write_t write_function,
				 void *fp)
{
  int d;
  int i;

  if (write_function(&(s->depth_size), sizeof(int), 1, fp) != 1) {
    ERROR("failed to write header");
    return -1;
  }

  for (d = 0; d < s->depth_size; d ++) {

    if (write_function(&(d), sizeof(int), 1, fp) != 1) {
      ERROR("failed to write depth");
      return -1;
    }
    if (write_function(&(s->set_n[d]), sizeof(int), 1, fp) != 1) {
      ERROR("failed to write depth count");
      return -1;
    }

    for (i = 0; i < s->set_n[d]; i ++) {

      if (write_function(&(s->s[d][i]), sizeof(int), 1, fp) != 1) {
	ERROR("failed to write index");
	return -1;
      }
      
      if (write_function(&(s->v[d][i]), sizeof(double), 1, fp) != 1) {
	ERROR("failed to write value");
	return -1;
      }
    }
  }

  return 0;
}

int
multiset_int_double_read_binary(multiset_int_double_t *s,
				multiset_int_double_read_t read_function,
				void *fp)
{
  int depth;
  int d;
  int di;
  int i;
  
  multiset_int_double_clear(s);

  if (read_function(&depth, sizeof(int), 1, fp) != 1) {
    ERROR("failed to read header");
    return -1;
  }

  if (depth > s->depth_size) {
    ERROR("depth mismatch: %d > %d", depth, s->depth_size);
    return -1;
  }

  for (d = 0; d < depth; d ++) {

    if (read_function(&(di), sizeof(int), 1, fp) != 1) {
      ERROR("failed to read depth");
      return -1;
    }

    if (di != d) {
      ERROR("depth index mismatch (%d != %d)", d, di);
      return -1;
    }

    if (read_function(&(s->set_n[d]), sizeof(int), 1, fp) != 1) {
      ERROR("failed to read depth count");
      return -1;
    }

    if (s->set_n[d] > s->set_size[d]) {
      if (multiset_int_double_expand_set(s, d, s->set_n[d]) < 0) {
	ERROR("failed to expand set");
	return -1;
      }
    }
								  
    for (i = 0; i < s->set_n[d]; i ++) {

      if (read_function(&(s->s[d][i]), sizeof(int), 1, fp) != 1) {
	ERROR("failed to read index");
	return -1;
      }
      
      if (read_function(&(s->v[d][i]), sizeof(double), 1, fp) != 1) {
	ERROR("failed to read value");
	return -1;
      }
    }
  }

  return 0;
}



/*
 * Internal functions
 */

static int multiset_int_double_expand_set(multiset_int_double_t *s, int depth, int minsize)
{
  int new_size;
  int *new_s;
  double *new_v;
  int i;

  if (s->set_size[depth] >= minsize) {
    return 0;
  }

  new_size = s->set_size[depth];
  while (new_size < minsize) {
    new_size += SET_INCREMENT;
  }

  
  new_s = (int *)malloc(sizeof(int) * new_size);
  if (new_s == NULL) {
    ERROR("failed to allocate new set");
    return -1;
  }

  new_v = (double *)malloc(sizeof(double) * new_size);
  if (new_v == NULL) {
    ERROR("failed to allocate new value set");
    return -1;
  }

  for (i = 0; i < s->set_n[depth]; i ++) {
    new_s[i] = s->s[depth][i];
    new_v[i] = s->v[depth][i];
  }

  free(s->s[depth]);
  free(s->v[depth]);

  s->s[depth] = new_s;
  s->v[depth] = new_v;

  s->set_size[depth] = new_size;

  return 0;

}

static int multiset_int_double_find_exact(int *s, 
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
    return multiset_int_double_find_exact(s, target, c + 1, end);
  } else if (s[c] > target) {
    return multiset_int_double_find_exact(s, target, start, c - 1);
  } else {
    return c;
  }
}

static int multiset_int_double_find_insertion_index(int *s,
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
    return multiset_int_double_find_insertion_index(s, idx, c + 1, end);
  } else if (s[c] > idx) {
    return multiset_int_double_find_insertion_index(s, idx, start, c - 1);
  } else {
    return -1;
  }
}

