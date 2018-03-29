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
#ifndef hnk_standard_h
#define hnk_standard_h

#include <gmp.h>

#include "hnk.h"

/*
 * Unary Tree Functions
 */
int 
unary_tree_maxk_at_h(hnk_t *t, int h, int nroot);

int
unary_tree_hnk(hnk_t *t, int h, int k, mpz_t hnk);

/*
 * Binary Tree Functions
 */
int
binary_tree_maxk_at_h(hnk_t *t, int h, int nroot);

int
binary_tree_hnk(hnk_t *t, int h, int k, mpz_t hnk);

hnk_t *
hnk_create_binary_tree(int maxh,
		       int maxk);

/*
 * Ternary Tree Functions
 */
int
ternary_tree_maxk_at_h(hnk_t *t, int h, int nroot);

int
ternary_tree_hnk(hnk_t *t, int h, int k, mpz_t hnk);

hnk_t *
hnk_create_ternary_tree(int maxh,
			int maxk);

/*
 * Quaternary Tree Functions
 */
int
quaternary_tree_maxk_at_h(hnk_t *t, int h, int nroot);

int
quaternary_tree_hnk(hnk_t *t, int h, int k, mpz_t hnk);

hnk_t *
hnk_create_quaternary_tree(int maxh,
			   int maxk);

/*
 * Senary Tree Functions
 */
int
senary_tree_maxk_at_h(hnk_t *t, int h, int nroot);

int
senary_tree_hnk(hnk_t *t, int h, int k, mpz_t hnk);

/*
 * Septenary Tree Functions
 */
int
septenary_tree_maxk_at_h(hnk_t *t, int h, int nroot);

int
septenary_tree_hnk(hnk_t *t, int h, int k, mpz_t hnk);

/*
 * Octary Tree functions
 */
int
octary_tree_maxk_at_h(hnk_t *t, int h, int nroot);

int
octary_tree_hnk(hnk_t *t, int h, int k, mpz_t hnk);

hnk_t *
hnk_create_octary_tree(int maxh,
		       int maxk);

/*
 * Nonary Tree functions
 */
int 
nonary_tree_maxk_at_h(hnk_t *t, int h, int nroot);

int
nonary_tree_hnk(hnk_t *t, int h, int k, mpz_t hnk);


  /*
Nullary means 0-ary.
Unary means 1-ary.
Binary means 2-ary.
Ternary means 3-ary.
Quaternary means 4-ary.
Quinary means 5-ary.
Senary means 6-ary.
Septenary means 7-ary.
Octary means 8-ary.
Nonary means 9-ary.
  */

#endif /* hnk_standard_h */
