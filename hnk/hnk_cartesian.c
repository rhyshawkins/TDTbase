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

#include "hnk.h"
#include "hnk_cartesian.h"
#include "hnk_standard.h"

#include "slog.h"

/*
 * 1,2 Tree for 1D
 */

hnk_t *
hnk_cartesian_12_create(int maxh,
			int maxk)
{
  hnk_t *sub;

  sub = hnk_create(maxh - 1, 
		   maxk - 1, 
		   2,
		   binary_tree_maxk_at_h, 
		   binary_tree_hnk,
		   NULL);
  if (sub == NULL) {
    ERROR("failed to create binary sub tree");
    return NULL;
  }

  return hnk_create(maxh,
		    maxk,
		    1,
		    unary_tree_maxk_at_h,
		    unary_tree_hnk,
		    sub);
}

hnk_t *
hnk_cartesian_13_create(int maxh,
			int maxk)
{
  hnk_t *sub;

  sub = hnk_create(maxh - 1, 
		   maxk - 1, 
		   3,
		   ternary_tree_maxk_at_h, 
		   ternary_tree_hnk,
		   NULL);
  if (sub == NULL) {
    ERROR("failed to create ternary sub tree");
    return NULL;
  }

  return hnk_create(maxh,
		    maxk,
		    1,
		    unary_tree_maxk_at_h,
		    unary_tree_hnk,
		    sub);
}

/*
 * 3,4 Tree for 2D 
 */

hnk_t *
hnk_cartesian_34_create(int maxh,
			int maxk)
{
  hnk_t *sub;

  sub = hnk_create(maxh - 1, 
		   maxk - 1, 
		   4,
		   quaternary_tree_maxk_at_h, 
		   quaternary_tree_hnk,
		   NULL);
  if (sub == NULL) {
    ERROR("failed to create quad sub tree");
    return NULL;
  }

  return hnk_create(maxh,
		    maxk,
		    3,
		    ternary_tree_maxk_at_h,
		    ternary_tree_hnk,
		    sub);
}

/*
 * 78 Tree for 3D Grids
 */

hnk_t *
hnk_cartesian_78_create(int maxh,
			int maxk)
{
  hnk_t *sub;

  sub = hnk_create(maxh - 1, 
		   maxk - 1, 
		   8,
		   octary_tree_maxk_at_h, 
		   octary_tree_hnk,
		   NULL);
  if (sub == NULL) {
    ERROR("failed to create oct sub tree");
    return NULL;
  }

  return hnk_create(maxh,
		    maxk,
		    7,
		    septenary_tree_maxk_at_h,
		    septenary_tree_hnk,
		    sub);

}
  


hnk_t *
hnk_cartesian_79_create(int maxh,
			int maxk)
{
  hnk_t *sub;

  sub = hnk_create(maxh - 1, 
		   maxk - 1, 
		   9,
		   nonary_tree_maxk_at_h, 
		   nonary_tree_hnk,
		   NULL);
  if (sub == NULL) {
    ERROR("failed to create nonary sub tree");
    return NULL;
  }

  return hnk_create(maxh,
		    maxk,
		    7,
		    septenary_tree_maxk_at_h,
		    septenary_tree_hnk,
		    sub);

}

