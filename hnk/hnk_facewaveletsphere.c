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

#include "hnk_facewaveletsphere.h"
#include "hnk_standard.h"

#include "slog.h"

hnk_t *
hnk_facewaveletsphere_icosahedron_create(int maxh,
					 int maxk)
{
  /* S2 face wavelet with icosahedron:
   * 
   * 20 Triangles from root
   * Each Triangle has 3 children
   * 
   */

  hnk_t *h3;
  hnk_t *root;

  h3 = hnk_create_ternary_tree(maxh - 1, maxk - 1);
  if (h3 == NULL) {
    ERROR("failed to create ternary subtree");
    return NULL;
  }

  root = hnk_create_aggregate(maxh, maxk,
			      20, 20,
			      h3, h3, h3, h3, h3,
			      h3, h3, h3, h3, h3,
			      h3, h3, h3, h3, h3,
			      h3, h3, h3, h3, h3);
  if (root == NULL) {
    ERROR("failed to create aggregate tree");
    return NULL;
  }

  return root;
}

hnk_t *
hnk_facewaveletsphereshell_icosahedron_create(int maxh,
					      int maxk,
					      int radial_degree)
{
  hnk_t *senary;
  hnk_t *h26;
  hnk_t *h16;
  hnk_t *root;

  int i;
  
  if (maxh != (radial_degree + 1)) {
    /* Not supported for now. */
    ERROR("varying radial degree not supported, require radial_degree = maxh - 1");
    return NULL;
  }

  senary = hnk_create(maxh - 1, maxk - 1, 6, senary_tree_maxk_at_h, senary_tree_hnk, NULL);
  if (senary == NULL) {
    ERROR("failed to create senary tree");
    return NULL;
  }

  h26 = hnk_create_aggregate_empty(maxh - 1, maxk - 1, 8, 8);
  if (h26 == NULL) {
    ERROR("failed to create h26 tree");
    return NULL;
  }

  if (hnk_aggregate_set_subtree(h26, 0, h26) < 0 ||
      hnk_aggregate_set_subtree(h26, 1, h26) < 0) {
    ERROR("failed to set h26 subtrees");
    return NULL;
  }

  for (i = 2; i < 8; i ++) {
    if (hnk_aggregate_set_subtree(h26, i, senary) < 0) {
      ERROR("failed to set h26 subtrees %d", i);
      return NULL;
    }
  }

  h16 = hnk_create_aggregate_empty(maxh - 1, maxk - 1, 7, 7);
  if (h16 == NULL) {
    ERROR("failed to create h16 tree");
    return NULL;
  }

  if (hnk_aggregate_set_subtree(h16, 0, h26) < 0) {
    ERROR("failed to set h16 subtree");
    return NULL;
  }

  for (i = 1; i < 7; i ++) {
    if (hnk_aggregate_set_subtree(h16, i, senary) < 0) {
      ERROR("failed to set h16 subtree %d", i);
      return NULL;
    }
  }

  root = hnk_create_aggregate(maxh, maxk, 20, 20,
			      h16, h16, h16, h16, h16,
			      h16, h16, h16, h16, h16,
			      h16, h16, h16, h16, h16,
			      h16, h16, h16, h16, h16);
  if (root == NULL) {
    ERROR("failed to create top level tree");
    return NULL;
  }

  return root;
}
