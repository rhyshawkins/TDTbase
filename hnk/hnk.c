//
//    HNK Library : A library for computing combinations of arrangements of
//    general trees for the Trans-dimensional Tree algorithm. See
//
//      R Hawkins and M Sambridge, "Geophysical imaging using trans-dimensional trees",
//      Geophysical Journal International, 2015, 203:2, 972 - 1000,
//      https://doi.org/10.1093/gji/ggv326
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
#include <stdarg.h>

#include "hnk.h"

#include "slog.h"

struct _hnk {
  int refcount;
  
  int maxh;
  int maxk;
  int maxsplit;
  int nsplits;
  
  int *maxk_at_h;
  int *maxk_at_h_storage;

  /*
   * Memoization of counts [maxh + 1][maxk + 1]
   */
  mpz_t **counts;

  /*
   * Memoization of power of 2 splits [maxh + 1][maxk + 1][log2 max split]
   */
  mpz_t ***split_counts;

  int nsubtree;
  hnk_t **subtree;

  hnk_compute_maxk_at_h_t cmaxk;
  hnk_compute_hnk_t chnk;

  /*
   * Memoization of ratios [maxh + 1][maxk + 1]
   */
  double **ratios;

  mpf_t r;
  mpf_t n;
  mpf_t d;

  mpz_t a;
  mpz_t b;
};

static int 
memoize_count(hnk_t *t,
	      int h,
	      int k);

static int
memoize_ratio(hnk_t *t,
	      int h,
	      int k);

static int 
nsplits(int maxsplit)
{
  int i = 2;
  int j = 1;
  while (i < maxsplit) {
    i *= 2;
    j ++;
  }
  return j;
}

static int splitindex(int split)
{
  return nsplits(split) - 1;
}

int
hnk_naggregatesplits(int nsubtree)
{
  if (nsubtree <= 0) {
    return -1;
  }
  
  /* Next power of 2 */
  nsubtree --;
  nsubtree |= nsubtree >> 1;
  nsubtree |= nsubtree >> 2;
  nsubtree |= nsubtree >> 4;
  nsubtree |= nsubtree >> 8;
  nsubtree |= nsubtree >> 16;
  nsubtree ++;
  
  return nsubtree;
}

static int
memoize_split_count(hnk_t *t,
		    hnk_t *subtree,
		    int h,
		    int k,
		    int nsplit);

static int 
memoize_aggregate_split_count(hnk_t *t,
			      hnk_t **subtree,
			      int nsubtree,
			      int index,
			      int h,
			      int k,
			      int nsplit);

hnk_t *
hnk_create(int maxh,
	   int maxk,
	   int maxsplit,
	   hnk_compute_maxk_at_h_t cmaxk,
	   hnk_compute_hnk_t chnk,
	   hnk_t *subtree)
{
  hnk_t *t;
  int i;

  t = (hnk_t*)malloc(sizeof(hnk_t));
  if (t == NULL) {
    return NULL;
  }

  t->refcount = 0;
  t->maxh = maxh;
  t->maxk = maxk;
  t->maxsplit = maxsplit;
  t->nsplits = nsplits(maxsplit);

  t->maxk_at_h = (int*)malloc(sizeof(int) * (maxh + 1));
  if (t->maxk_at_h == NULL) {
    return NULL;
  }

  t->maxk_at_h_storage = (int*)malloc(sizeof(int) * (maxh + 1));
  if (t->maxk_at_h_storage == NULL) {
    return NULL;
  }

  for (i = 0; i <= maxh; i ++) {
    t->maxk_at_h[i] = -1;
    t->maxk_at_h_storage[i] = -1;
  }
  
  t->counts = (mpz_t**)malloc(sizeof(mpz_t*) * (maxh + 1));
  if (t->counts == NULL) {
    return NULL;
  }

  for (i = 0; i <= maxh; i ++) {
    t->counts[i] = NULL;
  }

  t->split_counts = (mpz_t***)malloc(sizeof(mpz_t**) * (maxh + 1));
  if (t->split_counts == NULL) {
    return NULL;
  }

  for (i = 0; i <= maxh; i ++) {
    t->split_counts[i] = NULL;
  }

  t->ratios = (double**)malloc(sizeof(double*) * (maxh + 1));
  for (i = 0; i <= maxh; i ++) {
    t->ratios[i] = NULL;
  }
  
  t->cmaxk = cmaxk;
  t->chnk = chnk;

  t->subtree = NULL;
  t->nsubtree = 0;

  if (subtree != NULL) {
    t->subtree = malloc(sizeof(hnk_t *));
    if (t->subtree == NULL) {
      ERROR("failed to allocate subtree array");
      return NULL;
    }

    t->subtree[0] = subtree;
    t->nsubtree = 1;
  }

  mpf_init(t->r);
  mpf_init(t->n);
  mpf_init(t->d);

  mpz_init(t->a);
  mpz_init(t->b);

  return t;
}

hnk_t *
hnk_create_aggregate(int maxh,
		     int maxk,
		     int maxsplit,
		     int nsubtrees,
		     ...)
{
  hnk_t *t;
  int i;

  t = hnk_create_aggregate_empty(maxh, maxk, maxsplit, nsubtrees);
  if (t == NULL) {
    return NULL;
  }
  
  va_list ap;
  va_start(ap, nsubtrees); 
  for (i = 0; i < nsubtrees; i ++) {
    t->subtree[i] = va_arg(ap, hnk_t*);
    t->subtree[i]->refcount ++;
  }
  va_end(ap);

  return t;
}

hnk_t *
hnk_create_aggregate_empty(int maxh,
			   int maxk,
			   int maxsplit,
			   int nsubtrees)
{
  hnk_t *t;
  int i;

  t = (hnk_t*)malloc(sizeof(hnk_t));
  if (t == NULL) {
    return NULL;
  }

  t->refcount = 0;
  t->maxh = maxh;
  t->maxk = maxk;
  t->maxsplit = maxsplit;
  t->nsplits = hnk_naggregatesplits(maxsplit);

  t->maxk_at_h = (int*)malloc(sizeof(int) * (maxh + 1));
  if (t->maxk_at_h == NULL) {
    return NULL;
  }

  t->maxk_at_h_storage = (int*)malloc(sizeof(int) * (maxh + 1));
  if (t->maxk_at_h_storage == NULL) {
    return NULL;
  }

  for (i = 0; i <= maxh; i ++) {
    t->maxk_at_h[i] = -1;
    t->maxk_at_h_storage[i] = -1;
  }

  t->counts = (mpz_t**)malloc(sizeof(mpz_t*) * (maxh + 1));
  if (t->counts == NULL) {
    return NULL;
  }

  for (i = 0; i <= maxh; i ++) {
    t->counts[i] = NULL;
  }

  t->split_counts = (mpz_t***)malloc(sizeof(mpz_t**) * (maxh + 1));
  if (t->split_counts == NULL) {
    return NULL;
  }

  for (i = 0; i <= maxh; i ++) {
    t->split_counts[i] = NULL;
  }

  t->ratios = (double**)malloc(sizeof(double*) * (maxh + 1));
  for (i = 0; i <= maxh; i ++) {
    t->ratios[i] = NULL;
  }
  
  t->cmaxk = hnk_aggregate_maxk_at_h;
  t->chnk = hnk_aggregate_hnk;

  t->nsubtree = nsubtrees;
  t->subtree = malloc(sizeof(hnk_t *) * nsubtrees);

  if (t->subtree == NULL) {
    ERROR("failed to allocate subtree array");
    return NULL;
  }

  mpf_init(t->r);
  mpf_init(t->n);
  mpf_init(t->d);

  mpz_init(t->a);
  mpz_init(t->b);

  return t;
}

int
hnk_aggregate_set_subtree(hnk_t *t,
			  int i,
			  hnk_t *subtree)
{
  if (t == NULL ||
      i < 0 ||
      i >= t->nsubtree) {
    return -1;
  }
  
  t->subtree[i] = subtree;
  t->subtree[i]->refcount ++;

  return 0;
}

int hnk_aggregate_maxk_at_h(hnk_t *t,
			    int h,
			    int nroot)
{
  int i;
  int s;
  int mh;

  if (t == NULL) {
    return -1;
  }

  if (h == 0) {
    return nroot * 1;
  }

  s = nroot * 1;
  for (i = 0; i < t->nsubtree; i ++) {
    if (t->subtree[i]->maxh < (h - 1)) {
      mh = hnk_get_maxk_at_h(t->subtree[i], t->subtree[i]->maxh);
    } else {
      mh = hnk_get_maxk_at_h(t->subtree[i], h - 1);
    }

    if (mh < 0) {
      ERROR("failed to compute subtree maxk");
      return -1;
    }

    s += mh;
  }
  return s;
}

int hnk_aggregate_hnk(hnk_t *t,
		      int h,
		      int k,
		      mpz_t hnk)
{
  return hnk_aggregate_split(t, 
			     t->subtree,
			     t->nsubtree,
			     0,
			     h,
			     k,
			     t->nsubtree,
			     hnk);
}

static int hnk_save_stream(hnk_t *t,
			   FILE *fp)
{
  int i;
  int j;
  int k;
  int n;

  int maxk_storage;

  if (t == NULL) {
    ERROR("NULL hnk passed");
    return -1;
  }
  
  fwrite(&(t->maxh), sizeof(int), 1, fp);
  fwrite(&(t->maxk), sizeof(int), 1, fp);
  fwrite(&(t->maxsplit), sizeof(int), 1, fp);

  for (i = 0; i <= t->maxh; i ++) {
    if (t->maxk_at_h[i] < 0) {
      t->maxk_at_h[i] = (t->cmaxk)(t, i, 1);
    }
  }    
  fwrite(t->maxk_at_h, sizeof(int), t->maxh + 1, fp);

  /*
   * Counts
   */
  for (i = 0; i <= t->maxh; i ++) {
    if (t->counts[i] == NULL) {
      j = 0;
      fwrite(&j, sizeof(int), 1, fp);
    } else {
      maxk_storage = hnk_get_maxk_at_h_storage(t, i);
      fwrite(&maxk_storage, sizeof(int), 1, fp);

      for (j = 0; j <= maxk_storage; j ++) {
	mpz_out_raw(fp, t->counts[i][j]);
      }
    }
  }

  /*
   * Split counts
   */
  n = t->nsplits;
  for (i = 0; i <= t->maxh; i ++) {
    if (t->split_counts[i] == NULL) {
      j = 0;
      fwrite(&j, sizeof(int), 1, fp);
    } else {
      maxk_storage = hnk_get_maxk_at_h_storage(t, i);
      fwrite(&maxk_storage, sizeof(int), 1, fp);
      
      for (j = 0; j <= maxk_storage; j ++) {
	if (t->split_counts[i][j] == NULL) {
	  k = 0;
	  fwrite(&k, sizeof(int), 1, fp);
	} else {
	  fwrite(&n, sizeof(int), 1, fp);
	  
	  for (k = 0; k < n; k ++) {
	    mpz_out_raw(fp, t->split_counts[i][j][k]);
	  }
	}
      }
    }
  }

  /*
   * Ratios
   */
  for (i = 0; i <= t->maxh; i ++) {
    if (t->ratios[i] == NULL) {
      j = 0;
      fwrite(&j, sizeof(int), 1, fp);
    } else {
      maxk_storage = hnk_get_maxk_at_h_storage(t, i);
      fwrite(&maxk_storage, sizeof(int), 1, fp);
      
      for (j = 0; j <= maxk_storage; j ++) {
	
	fwrite(&t->ratios[i][j], sizeof(double), 1, fp);
	
      }
    }
  }
  

  /*
   * Subtrees
   */
  for (i = 0; i < t->nsubtree; i ++) {

    if (t->subtree[i] == t) {
      /*
       * Skip recursive subtrees
       */
      j = 0;
      fwrite(&j, sizeof(int), 1, fp);
      continue;
      
    } else {

      j = 1;
      fwrite(&j, sizeof(int), 1, fp);
      
      if (hnk_save_stream(t->subtree[i], fp) < 0) {
	ERROR("failed to write subtree");
	return -1;
      }
    }
  }

  return 0;
}

int hnk_save(hnk_t *t,
	     const char *filename)
{
  FILE *fp;
  
  fp = fopen(filename, "w");
  if (fp == NULL) {
    ERROR("failed to create file");
    return -1;
  }

  if (hnk_save_stream(t, fp) < 0) {
    ERROR("failed to save to file");
    return -1;
  }

  fclose(fp);
  return 0;

  fclose(fp);
  
  return 0;
}

static int hnk_restore_stream(hnk_t *t,
			      FILE *fp)
{
  int i;
  int j;
  int k;
  int n;

  int maxh;
  int maxk;
  int maxsplit;
  int maxk_storage;
  
  if (fread(&maxh, sizeof(int), 1, fp) != 1) {
    ERROR("failed to read maxh");
    return -1;
  }
  if (maxh != t->maxh) {
    ERROR("maxh mismatch: %d != %d", maxh, t->maxh);
    return -1;
  }
  
  if (fread(&maxk, sizeof(int), 1, fp) != 1) {
    ERROR("failed to read maxk");
    return -1;
  }
  if (maxk != t->maxk) {
    ERROR("maxk mismatch: %d != %d", maxk, t->maxk);
    return -1;
  }

  if (fread(&maxsplit, sizeof(int), 1, fp) != 1) {
    ERROR("failed to read maxsplit");
    return -1;
  }
  if (maxsplit != t->maxsplit) {
    ERROR("maxsplit mismatch: %d != %d", maxsplit, t->maxsplit);
    return -1;
  }

  if (fread(t->maxk_at_h, sizeof(int), t->maxh + 1, fp) != (t->maxh + 1)) {
    ERROR("failed to read maxk at h");
    return -1;
  }
  for (i = 0; i <= t->maxh; i ++) {
    if (t->maxk_at_h[i] < 0) {
      ERROR("maxk at h not set correctly");
      return -1;
    }
  }

  /*
   * Counts
   */
  for (i = 0; i <= t->maxh; i ++) {

    if (fread(&n, sizeof(int), 1, fp) != 1) {
      ERROR("failed to read n counts");
      return -1;
    }

    if (n > 0) {

      maxk_storage = hnk_get_maxk_at_h_storage(t, i);
      if (n != maxk_storage) {
	ERROR("mismatch in maxk at h: %d != %d", n, t->maxk_at_h[i]);
	return -1;
      }

      if (t->counts[i] == NULL) {
	/*
	 * Allocate counts
	 */
	t->counts[i] = (mpz_t*)malloc(sizeof(mpz_t) * (maxk_storage + 1));
	if (t->counts[i] == NULL) {
	  ERROR("failed to allocate table entry");
	  return -1;
	}
    
	for (j = 0; j <= maxk_storage; j ++) {
	  mpz_init_set_si(t->counts[i][j], -1);
	}
      }

      for (j = 0; j <= maxk_storage; j ++) {
	mpz_inp_raw(t->counts[i][j], fp);
      }
    }
  }
  
    
  /*
   * Split counts
   */
  n = t->nsplits;
  for (i = 0; i <= t->maxh; i ++) {
    if (fread(&j, sizeof(int), 1, fp) != 1) {
      ERROR("failed to read number of k split counts");
      return -1;
    }

    if (j > 0) {

      maxk_storage = hnk_get_maxk_at_h_storage(t, i);
      if (j != maxk_storage) {
	ERROR("splitcounts maxk mismatch %d != %d", j, t->maxk);
	return -1;
      }
      
      /*
       * Allocate if missing
       */
      if (t->split_counts[i] == NULL) {
	t->split_counts[i] = (mpz_t**)malloc(sizeof(mpz_t*) * (maxk_storage + 1));
	if (t->split_counts[i] == NULL) {
	  ERROR("failed to allocate split counts");
	  return -1;
	}

	for (j = 0; j <= maxk_storage; j ++) {
	  t->split_counts[i][j] = NULL;
	}
      }

      for (j = 0; j <= maxk_storage; j ++) {
	if (fread(&k, sizeof(int), 1, fp) != 1) {
	  ERROR("failed to read number of n split counts");
	  return -1;
	}

	if (k > 0) {

	  if (k != n) {
	    ERROR("splitcounts nsplits mismatch %d != %d", k, n);
	    return -1;
	  }

	  /*
	   * Allocate if missing
	   */
	  if (t->split_counts[i][j] == NULL) {
	    t->split_counts[i][j] = (mpz_t*)malloc(sizeof(mpz_t) * n);
	    if (t->split_counts[i][j] == NULL) {
	      ERROR("failed to allocate split array");
	      return -1;
	    }

	    for (k = 0; k < n; k ++) {
	      mpz_init_set_si(t->split_counts[i][j][k], -1);
	    }
	  }
	  
	  for (k = 0; k < n; k ++) {
	    mpz_inp_raw(t->split_counts[i][j][k], fp);
	  }
	}
      }
    }
  }

  /*
   * Ratios
   */
  for (i = 0; i <= t->maxh; i ++) {

    if (fread(&j, sizeof(int), 1, fp) != 1) {
      ERROR("failed to read ratio count");
      return -1;
    }
    
    if (j > 0) {

      maxk_storage = hnk_get_maxk_at_h_storage(t, i);
      if (j != maxk_storage) {
	ERROR("ratio size mismatch %d != %d", j, t->maxk_at_h[i]);
	return -1;
      }
      
      if (t->ratios[i] == NULL) {
	t->ratios[i] = (double*)malloc(sizeof(double) * (maxk_storage + 1));
	if (t->ratios[i] == NULL) {
	  ERROR("failed to allocate table entry");
	  return -1;
	}
	
	for (j = 0; j <= maxk_storage; j ++) {
	  t->ratios[i][j] = -1.0;
	}
      }
	
      for (j = 0; j <= maxk_storage; j ++) {
	
	if (fread(&t->ratios[i][j], sizeof(double), 1, fp) != 1) {
	  ERROR("failed to read ratio entry");
	  return -1;
	}
	
      }
    }
  }
  

  /*
   * Subtrees
   */
  for (i = 0; i < t->nsubtree; i ++) {

    if (fread(&j, sizeof(int), 1, fp) != 1) {
      ERROR("failed to read subtree mask");
      return -1;
    }
    
    if (j == 1) {
      
      if (hnk_restore_stream(t->subtree[i], fp) < 0) {
	ERROR("failed to read subtree");
	return -1;
      }
      
    } else {
      if (j != 0) {
	ERROR("flag incorrect: %d", j);
	return -1;
      }
    }
  }

  return 0;
}

int hnk_restore(hnk_t *t,
		const char *filename)
{
  FILE *fp;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    ERROR("failed to open file");
    return -1;
  }

  if (hnk_restore_stream(t, fp) < 0) {
    ERROR("failed to restore from file");
    return -1;
  }

  fclose(fp);
  return 0;
}


void
hnk_destroy(hnk_t *t)
{
  int i;
  int j;
  int k;
  int ns;
  int maxk_storage;

  if (t != NULL) {

    t->refcount --;
    if (t->refcount <= 0) {

      for (i = 0; i <= t->maxh; i ++) {

	if (t->counts[i] != NULL) {
	  maxk_storage = hnk_get_maxk_at_h_storage(t, i);
	  for (j = 0; j <= maxk_storage; j ++) {
	    mpz_clear(t->counts[i][j]);
	  }
	  free(t->counts[i]);
	}
      }
      
      free(t->counts);
      
      ns = t->nsplits;
      for (i = 0; i <= t->maxh; i ++) {
	if (t->split_counts[i] != NULL) {
	  maxk_storage = hnk_get_maxk_at_h_storage(t, i);
	  for (j = 0; j <= maxk_storage; j ++) {
	    if (t->split_counts[i][j] != NULL) {
	      for (k = 0; k < ns; k ++) {
		mpz_clear(t->split_counts[i][j][k]);
	      }
	      free(t->split_counts[i][j]);
	    }
	  }
	  free(t->split_counts[i]);
	}
      }
      
      free(t->split_counts);
      
      for (i = 0; i <= t->maxh; i ++) {
	if (t->ratios[i] != NULL) {
	  free(t->ratios[i]);
	}
      }
      free(t->ratios);
      
      free(t->maxk_at_h_storage);
      free(t->maxk_at_h);
      
      mpf_clear(t->r);
      mpf_clear(t->n);
      mpf_clear(t->d);
      
      mpz_clear(t->a);
      mpz_clear(t->b);
      
      if (t->nsubtree > 0) {
	for (i = 0; i < t->nsubtree; i ++) {
	  hnk_destroy(t->subtree[i]);
	}
	free(t->subtree);
      }
      
      free(t);
    }
  }
}

int
hnk_get_maxk_at_h(hnk_t *t, 
		  int h)
{
  if (h < 0 || h > t->maxh) {
    /* ERROR("hnk_get_maxk_at_h: invalid height %d\n", h); */
    return -1;
  }

  if (t->maxk_at_h[h] < 0) {

    t->maxk_at_h[h] = (t->cmaxk)(t, h, 1);

  }

  return t->maxk_at_h[h];
}

int
hnk_get_maxk_at_h_storage(hnk_t *t, 
			  int h)
{
  if (h < 0 || h > t->maxh) {
    return -1;
  }

  if (t->maxk_at_h_storage[h] < 0) {

    t->maxk_at_h_storage[h] = (t->cmaxk)(t, h, 1);
    if (t->maxk_at_h_storage[h] > t->maxk) {
      t->maxk_at_h_storage[h] = t->maxk;
    }

  }

  return t->maxk_at_h_storage[h];
}

int
hnk_get_hnk(hnk_t *t, 
	    int h, 
	    int k, 
	    mpz_t hnk)
{
  int maxk;

  if (h < 0 || h > t->maxh) {
    ERROR("invalid height (%d > %d)", h, t->maxh);
    return -1;
  }

  maxk = hnk_get_maxk_at_h(t, h);
  if (maxk < 0) {
    ERROR("failed to get maxk");
    return -1;
  }
  
  if (k < 0 || k > maxk) {
    mpz_set_ui(hnk, 0);
    return 0;
  }

  if (memoize_count(t, h, k) < 0) {
    ERROR("failed to memoize");
    return -1;
  }

  mpz_set(hnk, t->counts[h][k]);

  return 0;
}

int
hnk_is_hnk_memoized(hnk_t *t,
		    int h,
		    int k)
{
  if (h < 0 || h > t->maxh) {
    return 0;
  }
  
  if (t->counts[h] == NULL) {
    return 0;
  }

  if (k < 0 || k > t->maxk_at_h[h]) {
    return 0;
  }

  if (mpz_sgn(t->counts[h][k]) < 0) {
    return 0;
  }

  return 1;
}

int
hnk_highest_memoized_k(hnk_t *t,
		       int h)
{
  int i;
  int maxk_storage;
  
  if (h < 0 || h > t->maxh) {
    return 0;
  }
  
  if (t->counts[h] == NULL) {
    return 0;
  }

  maxk_storage = hnk_get_maxk_at_h_storage(t, h);
  if (maxk_storage < 0) {
    return -1;
  }

  for (i = maxk_storage; i >= 0; i --) {

    if (mpz_sgn(t->counts[h][i]) > 0) {
      return i;
    }
  }

  return 0;
}


/*
 * Returns ratio of hnk(k+1)/hnk(k) for given k
 */
int
hnk_get_kplus1_ratio(hnk_t *t, 
		     int h, 
		     int k, 
		     double *ratio)
{
  int maxk;

  if (h < 0 || h > t->maxh) {
    ERROR("h out of range (%d)", h);
    return -1;
  }

  maxk = hnk_get_maxk_at_h(t, h);
  if (maxk < 0) {
    ERROR("failed to get maxk");
    return -1;
  }
  if (k < 0 || k > maxk) {
    *ratio = 0.0;
    return 0;
  }

  if (memoize_ratio(t,
		    h,
		    k) < 0) {
    ERROR("failed to memoize");
    return -1;
  }

  *ratio = t->ratios[h][k];

  return 0;
}

int
hnk_get_subtree_maxk_at_h(hnk_t *t, 
			  int h)
{
  if (t == NULL ||
      t->nsubtree == 0 ||
      t->subtree[0] == NULL) {
    return -1;
  }

  return hnk_get_maxk_at_h(t->subtree[0], h);
}

int
hnk_get_subtree_hnk(hnk_t *t, 
		    int h, 
		    int k, 
		    mpz_t hnk)
{
  if (t == NULL ||
      t->nsubtree == 0 ||
      t->subtree[0] == NULL) {
    return -1;
  }

  return hnk_get_hnk(t->subtree[0], h, k, hnk);
}

hnk_t *
hnk_get_subtree(hnk_t *t)
{
  if (t->nsubtree == 0) {
    return NULL;
  }
  return t->subtree[0];
}

int
hnk_save_text(hnk_t *t, 
	      const char *filename)
{
  return -1;
}
	   
hnk_t *
hnk_load_text(const char *filename, 
	      hnk_compute_maxk_at_h_t cmaxk,
	      hnk_compute_hnk_t chnk,
	      hnk_t *subtree)
{
  return NULL;
}

static int 
memoize_count(hnk_t *t,
	      int h,
	      int k)
{
  int maxk;
  int i;
  int maxk_storage;

  if (h < 0 || h > t->maxh) {
    ERROR("invalid h");
    return -1;
  }

  maxk = hnk_get_maxk_at_h(t, h);
  if (maxk < 0) {
    ERROR("failed to get maxk");
    return -1;
  }
  
  if (k < 0 || k > maxk) {
    return -1;
  }

  if (t->counts[h] == NULL) {
    
    /*
     * Initialise table now for this height
     */
    maxk_storage = hnk_get_maxk_at_h_storage(t, h);
    if (maxk_storage < 0) {
      ERROR("failed to get maxk for storage");
      return -1;
    }
    
    t->counts[h] = (mpz_t*)malloc(sizeof(mpz_t) * (maxk_storage + 1));
    if (t->counts[h] == NULL) {
      ERROR("failed to allocate table entry");
      return -1;
    }
    
    for (i = 0; i <= maxk_storage; i ++) {
      mpz_init_set_si(t->counts[h][i], -1);
    }
    
  }

  if (mpz_sgn(t->counts[h][k]) < 0) {
    /* 
     * Not been set yet
     */
    if ((t->chnk)(t, h, k, t->counts[h][k]) < 0) {
      ERROR("failed to compute hnk for %d, %d",
	    h, k);
      return -1;
    }
  }

  return 0;
}

static int
memoize_ratio(hnk_t *t,
	      int h,
	      int k)
{
  int maxk;
  int i;
  int maxk_storage;

  if (h < 0 || h > t->maxh) {
    ERROR("h out of range (%d)", h);
    return -1;
  }

  maxk = hnk_get_maxk_at_h(t, h);
  if (maxk < 0) {
    ERROR("failed to get maxk");
    return -1;
  }
  
  if (k < 0 || k > maxk) {
    ERROR("k out of range (%d)", k);
    return -1;
  }

  if (t->ratios[h] == NULL) {
    
    /*
     * Initialise table now for this height
     */
    maxk_storage = hnk_get_maxk_at_h_storage(t, h);
    if (maxk_storage < 0) {
      ERROR("failed to get maxk for storage");
      return -1;
    }

    t->ratios[h] = (double*)malloc(sizeof(double) * (maxk_storage + 1));
    if (t->ratios[h] == NULL) {
      ERROR("failed to allocate table entry");
      return -1;
    }
    
    for (i = 0; i <= maxk_storage; i ++) {
      t->ratios[h][i] = -1.0;
    }
    
  }

  if (t->ratios[h][k] < 0.0) {
    /* 
     * Not been set yet
     */
    if (memoize_count(t, h, k) < 0) { 
      ERROR("failed to memoize den: %d", k);
      return -1;
    }
    mpf_set_z(t->d, t->counts[h][k]);

    if (k == maxk) {
      mpf_set_ui(t->n, 0);
    } else if (memoize_count(t, h, k + 1) < 0) {
      ERROR("failed to memoize num: %d", k);
      return -1;
    } else {
      mpf_set_z(t->n, t->counts[h][k + 1]);
    }

    mpf_div(t->r, t->n, t->d);

    t->ratios[h][k] = mpf_get_d(t->r);
  }

  return 0;
}

static int
memoize_split_count(hnk_t *t,
		    hnk_t *subtree,
		    int h,
		    int k,
		    int nsplit)
{
  int j;
  int j0;
  int jmax;
  int maxk;
  mpz_t a, b, c, s;

  int si;
  int ns;

  int maxk_storage;

  if (h > t->maxh ||
      k > t->maxk) {
    ERROR("h %d > maxh %d || k %d > maxk %d", h, t->maxh, k, t->maxk);
    return -1;
  }

  if ((nsplit % 2) != 0) {
    ERROR("not a power of 2 (%d)", nsplit);
    return -1;
  }

  if (t->split_counts[h] == NULL) {

    maxk_storage = hnk_get_maxk_at_h_storage(t, h);
    if (maxk_storage < 0) {
      ERROR("failed to get maxk for storage");
      return -1;
    }
    
    t->split_counts[h] = (mpz_t **)malloc(sizeof(mpz_t*) * (maxk_storage + 1));
    if (t->split_counts[h] == NULL) {
      ERROR("failed to allocate k array");
      return -1;
    }
    
    for (j = 0; j <= maxk_storage; j ++) {
      t->split_counts[h][j] = NULL;
    }
  }

  if (t->split_counts[h][k] == NULL) {
    ns = t->nsplits;

    t->split_counts[h][k] = (mpz_t*)malloc(sizeof(mpz_t) * ns);
    if (t->split_counts[h][k] == NULL) {
      ERROR("failed to allocate split array");
      return -1;
    }

    for (j = 0; j < ns; j ++) {
      mpz_init_set_si(t->split_counts[h][k][j], -1);
    }
  }

  si = splitindex(nsplit);

  if (mpz_sgn(t->split_counts[h][k][si]) < 0) {

    /*
     * Need to compute
     */

    if (nsplit == 2) {

      maxk = hnk_get_maxk_at_h(subtree, h - 1);
      if (maxk < 0) {
	ERROR("failed to get maxk (nsplit = 2)");
	return -1;
      }

      mpz_init_set_ui(s, 0);
      mpz_init(a);
      mpz_init(b);
      mpz_init(c);
      
      /*
       * Small optimisation: skip first few arrangements that we know will fail.
       */
      j0 = 0;
      if ((k - 1) > maxk) {
	j0 = k - 1 - maxk;
      }
      
      for (j = j0; j <= k - 1; j ++) {
	
	if (hnk_get_hnk(subtree, h - 1, j, a) < 0) {
	  ERROR("failed to get a (%d %d)", h - 1, j);
	  return -1;
	}
	
	if (hnk_get_hnk(subtree, h - 1, k - 1 - j, b) < 0) {
	  ERROR("failed to get b");
	  return -1;
	}
	
	mpz_mul(c, a, b);
	mpz_add(a, s, c);
	mpz_set(s, a);
      }
      
      mpz_set(t->split_counts[h][k][si], s);
      
      mpz_clear(a);
      mpz_clear(b);
      mpz_clear(c);
      mpz_clear(s);

    } else {
      
      maxk = hnk_get_maxk_at_h(subtree, h - 1);
      if (maxk < 0) {
	ERROR("failed to get maxk");
	return -1;
      }

      /* printf("Splitting %d: %d %d:\n", nsplit, h, k); */
      mpz_init_set_ui(s, 0);
      mpz_init(a);
      mpz_init(b);
      mpz_init(c);
      
      /*
       * Small optimisation: skip first few arrangements that we know will fail.
       */
      j0 = 0;
      if ((k - 1) > (maxk*nsplit/2)) {
      	j0 = k - 1 - (maxk*nsplit/2);
      }
      
      /*
       * Large optimisation: skip last few arrangements that wont fit in first subtree.
       */
      jmax = k - 1;
      if (jmax > (maxk*nsplit/2)) {
      	jmax = (maxk*nsplit/2);
      }

      for (j = j0; j <= jmax; j ++) {
	
	if (hnk_general_split(t,
			      subtree,
			      h,
			      j + 1, 
			      nsplit/2,
			      a) < 0) {
	  ERROR("failed to get sub tree");
	  return -1;
	}
	/* if (hnk_get_hnk(subtree, h - 1, j, a) < 0) { */
	/*   ERROR("hnk_general_split: failed to get a\n"); */
	/*   return -1; */
	/* } */
	
	if (hnk_general_split(t,
			      subtree,
			      h,
			      k - j, /* Note the lack of -1 here as we've already accounted for the used k */
			      nsplit/2,
			      b) < 0) {
	  ERROR("failed to get sub tree");
	  return -1;
	}
	//gmp_printf("    (%d) %d %d = %Zd x %Zd\n", h, j, k - j, a, b);
	
	mpz_mul(c, a, b);
	mpz_add(a, c, s);
	mpz_set(s, a);
      }
      
      mpz_set(t->split_counts[h][k][si], s);
      
      mpz_clear(a);
      mpz_clear(b);
      mpz_clear(c);
      mpz_clear(s);
    }
  }

  return 0;
}

int
hnk_general_split(hnk_t *t,
		  hnk_t *subtree,
		  int h,
		  int k,
		  int nsplit,
		  mpz_t count)
{
  int j;
  int j0;
  int jmax;
  int maxk;
  mpz_t a, b, c, s;

  if (k == 0 ||
      k == 1) {
    mpz_set_ui(count, 1);
    return 0;
  }

  maxk = hnk_get_maxk_at_h(subtree, h - 1);
  if (maxk < 0) {
    ERROR("failed to get maxk");
    return -1;
  }
  
  if (k > ((nsplit * maxk) + 1)) {
    mpz_set_ui(count, 0);
    return 0;
  } else if (k == ((nsplit * maxk) + 1)) {
    mpz_set_ui(count, 1);
    return 0;
  }

  switch (nsplit) {
    
  case 0:
  case 1:
    ERROR("invalid nsplit: %d", nsplit);
    return -1;

  case 2:

    if (memoize_split_count(t, 
			    subtree,
			    h,
			    k,
			    nsplit) < 0) {
      return -1;
    }

    mpz_set(count, t->split_counts[h][k][splitindex(nsplit)]);
    
    return 0;

  default:

    if ((nsplit % 2) == 1) {

      mpz_init_set_ui(s, 0);
      mpz_init(a);
      mpz_init(b);
      mpz_init(c);
      
      /*
       * Small optimisation: skip first few arrangements that we know will fail.
       */
      j0 = 0;
      if ((k - 1) > (nsplit - 1)*maxk) {
	j0 = k - 1 - (nsplit - 1)*maxk;
      }
      
      /*
       * Large optimisation: skip last few arrangements that wont fit in first subtree.
       */
      jmax = k - 1;
      if (jmax > maxk) {
    	jmax = maxk;
      }
      
      for (j = j0; j <= jmax; j ++) {
	
	if (hnk_get_hnk(subtree, h - 1, j, a) < 0) {
	  ERROR("failed to get a (%d %d)", h - 1, j);
	  return -1;
	}
	
	if (hnk_general_split(t,
			      subtree,
			      h,
			      k - j, /* Note the lack of -1 here as we've already accounted for the used k */
			      nsplit - 1,
			      b) < 0) {
	  ERROR("failed to get sub tree");
	  return -1;
	}
	
	mpz_mul(c, a, b);
	mpz_add(a, c, s);
	mpz_set(s, a);
      }
      
      mpz_set(count, s);
      
      mpz_clear(a);
      mpz_clear(b);
      mpz_clear(c);
      mpz_clear(s);
      return 0;

    } else {
      
      if (memoize_split_count(t, 
			      subtree,
			      h,
			      k,
			      nsplit) < 0) {
	return -1;
      }
      
      mpz_set(count, t->split_counts[h][k][splitindex(nsplit)]);
      
      return 0;

    }
  }
  
  ERROR("unreachable");
  return -1;
}

int
hnk_aggregate_split(hnk_t *t,
		    hnk_t **subtree,
		    int nsubtree,
		    int index,
		    int h,
		    int k,
		    int nsplit,
		    mpz_t count)
{
  int j;
  int j0;
  int jmax;
  int maxk;
  int maxkleft;
  int maxkright;
  mpz_t a, b, c, s;
  int i;
  int ileft;
  int sh;

  if (nsubtree < 2 ||
      nsplit < 2) {
    ERROR("invalid number of splits %d %d\n",
	  nsubtree, nsplit);
    return -1;
  }
      
  if (k == 0 ||
      k == 1) {
    /* printf("  hnk_aggregate_split: early exit 1\n"); */
    mpz_set_ui(count, 1);
    return 0;
  }

  if (h == 0 && k > nsplit) {
    mpz_set_ui(count, 0);
    return 0;
  }

  /*
   * Determine where to split and the maximum k allowed each side
   */
  maxkleft = 0;
  maxkright = 0;

  if (nsubtree % 2 == 0) {
    ileft = nsubtree/2;
  } else {
    ileft = nsubtree/2 + 1;
  }

  for (i = 0; i < nsubtree; i ++) {
    sh = h - 1;
    if (sh > subtree[i]->maxh) {
      sh = subtree[i]->maxh;
    }

    maxk = hnk_get_maxk_at_h(subtree[i], sh);
    if (maxk < 0) {
      ERROR("failed to get maxk");
      return -1;
    }

    if (i < ileft) {
      maxkleft += maxk;
    } else {
      maxkright += maxk;
    }
  }

  /*
   * Edge cases with known answers
   */
  if (k > (maxkleft + maxkright + 1)) {
    /* printf("  hnk_aggregate_split: early exit 0 (maxk)\n"); */
    mpz_set_ui(count, 0);
    return 0;
  } else if (k == (maxkleft + maxkright + 1)) {
    /* printf("  hnk_aggregate_split: early exit 0 (maxk)\n"); */
    mpz_set_ui(count, 1);
    return 0;
  }


  if (nsplit == 3) {

    /*
     * We need to compute the first 2 as a split, and the last member of the 
     * aggregate as an individual tree.
     *
     * We tree the left hand first 2 subtrees + the root as a single tree and
     * the right hand subtree as a distinct tree. Hence there is always 1
     * node int the left (ie j0 = 1) and we specify the height of the left
     * tree as h and right right as h - 1.
     */

    mpz_init_set_ui(s, 0);
    mpz_init(a);
    mpz_init(b);
    mpz_init(c);
      
    /*
     * Small optimisation: skip first few arrangements that we know will fail.
     */
    j0 = 1;
    if ((k - 1) > maxkright) {
      j0 = k - 1 - maxkright;
    }
    
    /*
     * Large optimisation: skip last few arrangements that wont fit in first subtree.
     */
    jmax = k;
    if (jmax > (maxkleft + 1)) {
      jmax = maxkleft + 1;
    }

    for (j = j0; j <= jmax; j ++) {

      if (hnk_aggregate_split(t,
			      subtree,
			      2,
			      2*index + 1,
			      h,     /* Must be h not h - 1 here as we're using this sub-trees root */
			      j, 
			      2,
			      a) < 0) {
	ERROR("failed to get sub tree (%d %d)", h - 1, index);
	return -1;
      }

      sh = h - 1;
      if (sh > subtree[2]->maxh) {
	sh = subtree[2]->maxh;
      }
      if (hnk_get_hnk(subtree[2], sh, k - j, b) < 0) {
	ERROR("failed to get right");
	return -1;
      }
      

      mpz_mul(c, a, b);
      mpz_add(a, c, s);
      mpz_set(s, a);
    }

    mpz_set(count, s);
    
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(c);
    mpz_clear(s);
    return 0;

  } else {


    if (memoize_aggregate_split_count(t, 
				      subtree,
				      nsubtree,
				      index,
				      h,
				      k,
				      nsplit) < 0) {
      return -1;
    }

    /* gmp_printf("  memoize: %d %d (nsub = %d, index = %d) %Zd\n", h, k, nsubtree, index, t->split_counts[h][k][index]); */
    mpz_set(count, t->split_counts[h][k][index]);
    
    return 0;
  }
  
  ERROR("unreachable");
  return -1;
}

static int 
memoize_aggregate_split_count(hnk_t *t,
			      hnk_t **subtree,
			      int nsubtree,
			      int index,
			      int h,
			      int k,
			      int nsplit)
{
  int j;
  int j0;
  int jmax;
  mpz_t a, b, c, s;

  int si;
  int ns;

  int i;
  int ileft;
  int maxk;
  int maxkleft;
  int maxkright;
  int sh;

  if (h > t->maxh ||
      k > t->maxk) {
    ERROR("h %d > maxh %d || k %d > maxk %d", h, t->maxh, k, t->maxk);
    return -1;
  }

  if (t->split_counts[h] == NULL) {
    t->split_counts[h] = (mpz_t **)malloc(sizeof(mpz_t*) * (t->maxk + 1));
    if (t->split_counts[h] == NULL) {
      ERROR("failed to allocate k array");
      return -1;
    }
    for (j = 0; j <= t->maxk; j ++) {
      t->split_counts[h][j] = NULL;
    }
  }

  if (t->split_counts[h][k] == NULL) {
    ns = t->nsplits;

    t->split_counts[h][k] = (mpz_t*)malloc(sizeof(mpz_t) * ns);
    if (t->split_counts[h][k] == NULL) {
      ERROR("failed to allocate split array");
      return -1;
    }

    for (j = 0; j < ns; j ++) {
      mpz_init_set_si(t->split_counts[h][k][j], -1);
    }
  }

  si = index;
  if (si < 0 || si >= t->nsplits) {
    ERROR("index out of range %d (%d)",
	  si, t->nsplits);
    return -1;
  }

  if (mpz_sgn(t->split_counts[h][k][si]) < 0) {

    /*
     * Need to compute, first compute the split 
     */

    maxkleft = 0;
    maxkright = 0;
    
    if (nsubtree % 2 == 0) {
      ileft = nsubtree/2;
    } else {
      ileft = nsubtree/2 + 1;
    }
    
    for (i = 0; i < nsubtree; i ++) {
      sh = h - 1;
      if (sh > subtree[i]->maxh) {
	sh = subtree[i]->maxh;
      }
      maxk = hnk_get_maxk_at_h(subtree[i], sh);
      if (i < ileft) {
	maxkleft += maxk;
      } else {
	maxkright += maxk;
      }
    }

    if (nsplit == 2) {

      mpz_init_set_ui(s, 0);
      mpz_init(a);
      mpz_init(b);
      mpz_init(c);
      
      /*
       * Small optimisation: skip first few arrangements that we know will fail.
       */
      j0 = 0;
      if ((k - 1) > maxkright) {
	j0 = k - 1 - maxkright;
      }
      
      for (j = j0; j <= k - 1; j ++) {

	sh = h - 1;
	if (sh > subtree[0]->maxh) {
	  sh = subtree[0]->maxh;
	}
	if (hnk_get_hnk(subtree[0], sh, j, a) < 0) {
	  ERROR("failed to get a (%p, %d %d)", subtree[0], sh, j);
	  return -1;
	}
	
	sh = h - 1;
	if (sh > subtree[1]->maxh) {
	  sh = subtree[1]->maxh;
	}
	if (hnk_get_hnk(subtree[1], sh, k - 1 - j, b) < 0) {
	  ERROR("failed to get b");
	  return -1;
	}

	mpz_mul(c, a, b);
	mpz_add(a, s, c);
	mpz_set(s, a);
      }
      
      mpz_set(t->split_counts[h][k][si], s);
      
      mpz_clear(a);
      mpz_clear(b);
      mpz_clear(c);
      mpz_clear(s);

    } else {
      
      /* printf("Splitting %d: %d %d:\n", nsplit, h, k); */
      mpz_init_set_ui(s, 0);
      mpz_init(a);
      mpz_init(b);
      mpz_init(c);
      
      /*
       * Small optimisation: skip first few arrangements that we know will fail.
       */
      j0 = 0;
      if ((k - 1) > maxkright) {
      	j0 = k - 1 - maxkright;
      }
      
      /*
       * Large optimisation: skip last few arrangements that wont fit in first subtree.
       */
      jmax = k - 1;
      if (jmax > maxkleft) {
      	jmax = maxkleft;
      }

      for (j = j0; j <= jmax; j ++) {
	
	if (hnk_aggregate_split(t,
				subtree,
				ileft,
				2*index + 1,
				h,
				j + 1, 
				ileft,
				a) < 0) {
	  ERROR("failed to get sub tree");
	  return -1;
	}
	/* if (hnk_get_hnk(subtree, h - 1, j, a) < 0) { */
	/*   ERROR("hnk_general_split: failed to get a\n"); */
	/*   return -1; */
	/* } */

	if (hnk_aggregate_split(t,
				subtree + ileft,
				nsubtree - ileft,
				2*index + 2,
				h,
				k - j, /* Note the lack of -1 here as we've already accounted for the used k */
				nsubtree - ileft,
				b) < 0) {
	  ERROR("failed to get sub tree");
	  return -1;
	}
	
	mpz_mul(c, a, b);
	mpz_add(a, c, s);
	mpz_set(s, a);
      }
      
      mpz_set(t->split_counts[h][k][si], s);
      
      mpz_clear(a);
      mpz_clear(b);
      mpz_clear(c);
      mpz_clear(s);
    }
  }

  return 0;
}

