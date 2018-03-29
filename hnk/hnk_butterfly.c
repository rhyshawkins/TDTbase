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

#include "hnk_butterfly.h"

#include "hnk_standard.h"
#include "hnk_cartesian.h"

#include "slog.h"

/*
 * For S2
 */
hnk_t *
hnk_butterfly_tetrahedron_create(int maxh,
				 int maxk)
{
  hnk_t *pole;
  hnk_t *cart34;
  hnk_t *b;
  /*
   * Tetrahedron tree is of the form:
   *
   *              pole - 0 - pole
   *                    / \
   *                 3,4   3,4
   */

  pole = hnk_create(0, 1, 1, unary_tree_maxk_at_h, unary_tree_hnk, NULL);
  if (pole == NULL) {
    ERROR("failed to create pole");
    return NULL;
  }

  cart34 = hnk_cartesian_34_create(maxh - 1, maxk - 1);
  if (cart34 == NULL) {
    ERROR("failed to create 3,4 sub tree");
    return NULL;
  }

  b = hnk_create_aggregate(maxh, maxk, 4, 4, pole, cart34, cart34, pole);
  if (b == NULL) {
    ERROR("failed to create aggregate tree");
    return NULL;
  }

  return b;
}


hnk_t *
hnk_butterfly_octahedron_create(int maxh,
				int maxk)
{
  hnk_t *pole;
  hnk_t *cart34;
  hnk_t *b;
  /*
   * Octahedron tree is of the form:
   * 
   *                pole   pole
   *                    \ /
   *               3,4 - 0 - 3,4
   *                    / \
   *                 3,4   3,4
   */

  pole = hnk_create(0, 1, 1, unary_tree_maxk_at_h, unary_tree_hnk, NULL);
  if (pole == NULL) {
    ERROR("failed to create pole");
    return NULL;
  }

  cart34 = hnk_cartesian_34_create(maxh - 1, maxk - 1);
  if (cart34 == NULL) {
    ERROR("failed to create 3,4 sub tree");
    return NULL;
  }

  /* b = hnk_create_aggregate(maxh, maxk, 6, 6, pole, cart34, cart34, cart34, cart34, pole); */
  b = hnk_create_aggregate(maxh, maxk, 6, 6, cart34, cart34, pole, cart34, cart34, pole);
  if (b == NULL) {
    ERROR("failed to create aggregate tree");
    return NULL;
  }

  return b;
}

hnk_t *
hnk_butterfly_icosahedron_create(int maxh,
				 int maxk)
{
  hnk_t *pole;
  hnk_t *cart34;
  hnk_t *b;

  /*
   * Icosahedron tree is of the form:
   * 
   *                pole   pole
   *                    \ /
   *                     0      
   *                    /  
   *                 3,4 x 10
   */

  pole = hnk_create(0, 1, 1, unary_tree_maxk_at_h, unary_tree_hnk, NULL);
  if (pole == NULL) {
    ERROR("failed to create pole");
    return NULL;
  }

  cart34 = hnk_cartesian_34_create(maxh - 1, maxk - 1);
  if (cart34 == NULL) {
    ERROR("failed to create 3,4 sub tree");
    return NULL;
  }

  b = hnk_create_aggregate(maxh, maxk, 12, 12,
			   pole,
			   cart34, cart34, cart34, cart34, cart34, 
			   cart34, cart34, cart34, cart34, cart34, 
			   pole);
  if (b == NULL) {
    ERROR("failed to create aggregate tree");
    return NULL;
  }

  return b;
}

/*
 * For thick shells
 */
hnk_t *
hnk_butterfly_tetrahedron_shell_create(int maxh,
				       int maxk,
				       int radial_degree)
{
  hnk_t *pole;
  hnk_t *cart78;
  hnk_t *b;
  
  if (maxh == (radial_degree + 1)) {

    /*
     * Simplest case
     */

    pole = hnk_cartesian_12_create(maxh - 1, maxk - 1);
    if (pole == NULL) {
      ERROR("failed to create pole");
      return NULL;
    }

    cart78 = hnk_cartesian_78_create(maxh - 1, maxk - 1);
    if (cart78 == NULL) {
      ERROR("failed to create 7,8 sub tree");
      return NULL;
    }
    
    b = hnk_create_aggregate(maxh, maxk, 4, 4,
			     pole,
			     cart78, cart78,
			     pole);
    if (b == NULL) {
      ERROR("failed to create aggregate tree");
      return NULL;
    }

    return b;

  } else {

    ERROR("unimplemented");
    return NULL;

  }
  
}

hnk_t *
hnk_butterfly_octahedron_shell_create(int maxh,
				      int maxk,
				      int radial_degree)
{

  hnk_t *pole;
  hnk_t *cart78;
  hnk_t *b;
  
  if (maxh == (radial_degree + 1)) {

    /*
     * Simplest case
     */

    pole = hnk_cartesian_12_create(maxh - 1, maxk - 1);
    if (pole == NULL) {
      ERROR("failed to create pole");
      return NULL;
    }

    cart78 = hnk_cartesian_78_create(maxh - 1, maxk - 1);
    if (cart78 == NULL) {
      ERROR("failed to create 7,8 sub tree");
      return NULL;
    }
    
    b = hnk_create_aggregate(maxh, maxk, 6, 6,
			     pole,
			     cart78, cart78, cart78, cart78,
			     pole);
    if (b == NULL) {
      ERROR("failed to create aggregate tree");
      return NULL;
    }

    return b;

  } else {

    ERROR("unimplemented");
    return NULL;

  }

}


hnk_t *
hnk_butterfly_icosahedron_shell_create(int maxh,
				       int maxk,
				       int radial_degree)
{
  hnk_t *pole;
  hnk_t *cart78;
  hnk_t *b;
  
  if (maxh == (radial_degree + 1)) {

    /*
     * Simplest case
     */

    pole = hnk_cartesian_12_create(maxh - 1, maxk - 1);
    if (pole == NULL) {
      ERROR("failed to create pole");
      return NULL;
    }

    cart78 = hnk_cartesian_78_create(maxh - 1, maxk - 1);
    if (cart78 == NULL) {
      ERROR("failed to create 7,8 sub tree");
      return NULL;
    }
    
    b = hnk_create_aggregate(maxh, maxk, 12, 12,
			     pole,
			     cart78, cart78, cart78, cart78, cart78,
			     cart78, cart78, cart78, cart78, cart78,
			     pole);
    if (b == NULL) {
      ERROR("failed to create aggregate tree");
      return NULL;
    }

    return b;

  } else {

    ERROR("unimplemented");
    return NULL;

  }
}
