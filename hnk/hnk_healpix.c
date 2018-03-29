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
#include <stdlib.h>

#include "hnk_healpix.h"
#include "hnk_standard.h"

#include "slog.h"

hnk_t *
hnk_healpix2d_create(int maxh,
		     int maxk)
{
  hnk_t *h4;
  hnk_t *root;

  h4 = hnk_create_quaternary_tree(maxh - 1, maxk - 1);
  if (h4 == NULL) {
    ERROR("Failed to create quaternary subtree");
    return NULL;
  }

  root = hnk_create_aggregate(maxh, maxk,
			      12, 12,
			      h4, h4, h4, h4, h4, h4,
			      h4, h4, h4, h4, h4, h4);
  if (root == NULL) {
    ERROR("Failed to create aggregate tree");
    return NULL;
  }

  return root;
}

hnk_t *
hnk_healpix3d_create(int maxh,
		     int maxk)
{
  hnk_t *h8;
  hnk_t *root;

  h8 = hnk_create_quaternary_tree(maxh - 1, maxk - 1);
  if (h8 == NULL) {
    ERROR("Failed to create quaternary subtree");
    return NULL;
  }

  root = hnk_create_aggregate(maxh, maxk,
			      12, 12,
			      h8, h8, h8, h8, h8, h8,
			      h8, h8, h8, h8, h8, h8);
  if (root == NULL) {
    ERROR("Failed to create aggregate tree");
    return NULL;
  }

  return root;
}
