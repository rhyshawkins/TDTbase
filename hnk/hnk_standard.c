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

#include "hnk_standard.h"

static int n_ary_tree_maxk_at_h(int n, hnk_t *t, int h, int nroot);
static int n_ary_tree_hnk(int n, hnk_t *t, int h, int k, mpz_t hnk);

/*
 * Unary Tree Functions
 */
int 
unary_tree_maxk_at_h(hnk_t *t, int h, int nroot)
{
  return n_ary_tree_maxk_at_h(1, t, h, nroot);
}

int
unary_tree_hnk(hnk_t *t, int h, int k, mpz_t hnk)
{
  return n_ary_tree_hnk(1, t, h, k, hnk);
}

/*
 * Binary Tree Functions
 */
int
binary_tree_maxk_at_h(hnk_t *t, int h, int nroot)
{
  return n_ary_tree_maxk_at_h(2, t, h, nroot);
}

int
binary_tree_hnk(hnk_t *t, int h, int k, mpz_t hnk)
{
  return n_ary_tree_hnk(2, t, h, k, hnk);
}

hnk_t *
hnk_create_binary_tree(int maxh,
		       int maxk)
{
  return hnk_create(maxh, maxk, 2, binary_tree_maxk_at_h, binary_tree_hnk, NULL);
}

/*
 * Ternary Tree Functions
 */
int
ternary_tree_maxk_at_h(hnk_t *t, int h, int nroot)
{
  return n_ary_tree_maxk_at_h(3, t, h, nroot);
}

int
ternary_tree_hnk(hnk_t *t, int h, int k, mpz_t hnk)
{
  return n_ary_tree_hnk(3, t, h, k, hnk);
}

hnk_t *
hnk_create_ternary_tree(int maxh,
			int maxk)
{
  return hnk_create(maxh, maxk, 3, ternary_tree_maxk_at_h, ternary_tree_hnk, NULL);
}

/*
 * Quaternary Tree Functions
 */
int
quaternary_tree_maxk_at_h(hnk_t *t, int h, int nroot)
{
  return n_ary_tree_maxk_at_h(4, t, h, nroot);
}

int
quaternary_tree_hnk(hnk_t *t, int h, int k, mpz_t hnk)
{
  return n_ary_tree_hnk(4, t, h, k, hnk);
}

hnk_t *
hnk_create_quaternary_tree(int maxh,
			   int maxk)
{
  return hnk_create(maxh, maxk, 4, quaternary_tree_maxk_at_h, quaternary_tree_hnk, NULL);
}

int
senary_tree_maxk_at_h(hnk_t *t, int h, int nroot)
{
  return n_ary_tree_maxk_at_h(6, t, h, nroot);
}

int
senary_tree_hnk(hnk_t *t, int h, int k, mpz_t hnk)
{
  return n_ary_tree_hnk(6, t, h, k, hnk);
}

/*
 * Septenary Tree Functions
 */
int
septenary_tree_maxk_at_h(hnk_t *t, int h, int nroot)
{
  return n_ary_tree_maxk_at_h(7, t, h, nroot);
}

int
septenary_tree_hnk(hnk_t *t, int h, int k, mpz_t hnk)
{
  return n_ary_tree_hnk(7, t, h, k, hnk);
}

/*
 * Octary Tree functions
 */
int
octary_tree_maxk_at_h(hnk_t *t, int h, int nroot)
{
  return n_ary_tree_maxk_at_h(8, t, h, nroot);
}

int
octary_tree_hnk(hnk_t *t, int h, int k, mpz_t hnk)
{
  return n_ary_tree_hnk(8, t, h, k, hnk);
}

hnk_t *
hnk_create_octary_tree(int maxh,
		       int maxk)
{
  return hnk_create(maxh, maxk, 8, octary_tree_maxk_at_h, octary_tree_hnk, NULL);
}

/*
 * Nonary Tree functions
 */

int 
nonary_tree_maxk_at_h(hnk_t *t, int h, int nroot)
{
  return n_ary_tree_maxk_at_h(9, t, h, nroot);
}

int
nonary_tree_hnk(hnk_t *t, int h, int k, mpz_t hnk)
{
  return n_ary_tree_hnk(9, t, h, k, hnk);
}


/*
 * The local functions that do the work
 */
static int n_ary_tree_maxk_at_h(int n, hnk_t *t, int h, int nroot)
{
  hnk_t *sub;

  if (n < 1) {
    return -1;
  }

  if (h > 31) {
    return -1;
  }

  if (h < 0) {
    return 0;
  }

  if (h == 0) {
    return nroot * 1;
  } 
    
  sub = hnk_get_subtree(t);
  if (sub == NULL) {
    sub = t;
  }
  
  return nroot * (1 + n * hnk_get_maxk_at_h(sub, h - 1));
}

static int n_ary_tree_hnk(int n, hnk_t *t, int h, int k, mpz_t hnk)
{
  hnk_t *sub;

  if (n < 1) {
    return -1;
  }

  sub = hnk_get_subtree(t);
  if (sub == NULL) {
    sub = t;
  }

  if (n == 1) {
    if (h == 0 || k == 1 || k == 0) {

      if (k < 0 || k > 1) {
	mpz_set_ui(hnk, 0);
      } else {
	mpz_set_ui(hnk, 1);
      }
      return 0;

    } else {

      return hnk_get_hnk(sub, h - 1, k - 1, hnk);

    }
  } else {

    return hnk_general_split(t, sub, h, k, n, hnk);

  }
}
