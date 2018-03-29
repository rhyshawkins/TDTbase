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

#include "hnk_s2.h"
#include "hnk_standard.h"

#include "slog.h"

hnk_t *
hnk_s2_icosahedron_create(int maxh,
			  int maxk)
{
  /* S2 face with icosahedron subdivision:
   * 
   * 20 Triangles from root
   * Each Triangle has 4 children
   * 
   */

  hnk_t *h4;
  hnk_t *root;

  h4 = hnk_create_quaternary_tree(maxh - 1, maxk - 1);
  if (h4 == NULL) {
    ERROR("failed to create quaternary subtree");
    return NULL;
  }

  root = hnk_create_aggregate(maxh, maxk,
			      20, 20,
			      h4, h4, h4, h4, h4,
			      h4, h4, h4, h4, h4,
			      h4, h4, h4, h4, h4,
			      h4, h4, h4, h4, h4);
  if (root == NULL) {
    ERROR("failed to create aggregate tree");
    return NULL;
  }

  return root;
}

