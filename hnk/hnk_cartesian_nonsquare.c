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

#include "hnk_cartesian_nonsquare.h"

#include "hnk.h"
#include "hnk_cartesian.h"
#include "hnk_standard.h"

#include "slog.h"

static hnk_t *
hnk_cartesian_nonsquare_2D_create_subquadtree(int maxh_width,
					      int maxh_height,
					      int maxk);

hnk_t *
hnk_cartesian_nonsquare_2D_create(int maxh_width,
				  int maxh_height,
				  int maxk)
{
  hnk_t *sh = NULL;
  hnk_t *sv = NULL;
  hnk_t *sd = NULL;
  int maxh;
  
  if (maxh_width < 1 ||
      maxh_height < 1 ||
      maxh_width > 16 ||
      maxh_height > 16) {
    ERROR("height/width out of range (%d %d)\n",
	  maxh_width, maxh_height);
    return NULL;
  }


  /*
   * Just return a standard 3,4 tree if it's actually square
   */
  if (maxh_width == maxh_height) {

    return hnk_cartesian_34_create(maxh_width,
				   maxk);

  }

  /*
   * Reorder so that width > height, as this simplifies the code and for counting arrangements
   * this doesn't change anything.
   */
  if (maxh_width < maxh_height) {
    maxh = maxh_height;
    maxh_height = maxh_width;
    maxh_width = maxh;
  } else {
    maxh = maxh_width;
  }

  if (maxh_width < maxh_height) {
    ERROR("width < height\n");
    return NULL;
  }
  

  /* 
   * Horizontal just a quad
   */
  sh = hnk_cartesian_nonsquare_2D_create_subquadtree(maxh_width - 1,
						     maxh_height,
						     maxk - 1);
  
  /*
   * Diagonal just a quad
   */
  sd = hnk_create_quaternary_tree(maxh_height - 1,
				  maxk - 1);
  
  /*
   * Vertical a staggered quad
   */
  sv = hnk_create_quaternary_tree(maxh_height - 1,
				  maxk - 1);

  if (sh == NULL ||
      sd == NULL ||
      sv == NULL) {
    ERROR("failed to construct subtrees (%p %p %p)\n",
	  sh, sd, sv);
    return NULL;
  }
    
  return hnk_create_aggregate(maxh,
			      maxk,
			      3,
			      3,
			      sh,
			      sd,
			      sv);
}

static hnk_t *
hnk_cartesian_nonsquare_2D_create_subquadtree(int maxh_width,
					      int maxh_height,
					      int maxk)
{
  hnk_t *sh;
  hnk_t *sv;

  if (maxh_width == maxh_height) {

    return hnk_create_quaternary_tree(maxh_width,
				      maxk);

  }

  if (maxh_width == 0) {
    return hnk_create_binary_tree(maxh_height, maxk);
  }

  if (maxh_height == 0) {
    return hnk_create_binary_tree(maxh_width, maxk);
  }

  if (maxh_width < maxh_height) {
    ERROR("width < height (%d < %d)",
	  maxh_width, maxh_height);
    return NULL;
  }
  

  sh = hnk_cartesian_nonsquare_2D_create_subquadtree(maxh_width - 1,
						     maxh_height,
						     maxk - 1);
  if (sh == NULL) {
    return NULL;
  }
  
  sv = hnk_create_quaternary_tree(maxh_height - 1,
				  maxk - 1);
  
  return hnk_create_aggregate(maxh_width,
			      maxk,
			      4,
			      4,
			      sh,
			      sh,
			      sv,
			      sv);

}


hnk_t *
hnk_cartesian_nonsquare_2D_create_sub(int maxh_width,
				      int maxh_height,
				      int maxk)
{
  int maxh_min;
  int subcount;

  hnk_t *h34;
  hnk_t *h;

  int i;

  maxh_min = maxh_width;
  if (maxh_height < maxh_min) {
    maxh_min = maxh_height;
  }

  subcount =
    (1 << (maxh_width - maxh_min)) *
    (1 << (maxh_height - maxh_min));

  if (subcount == 1) {
    return hnk_cartesian_34_create(maxh_width, maxk);
  }

  h34 = hnk_cartesian_34_create(maxh_min, maxk - 1);
  if (h34 == NULL) {
    return NULL;
  }

  h = hnk_create_aggregate_empty(maxh_min + 1,
				 maxk,
				 subcount,
				 subcount);
  if (h == NULL) {
    return NULL;
  }

  for (i = 0; i < subcount; i ++) {
    if (hnk_aggregate_set_subtree(h, i, h34) < 0) {
      return NULL;
    }
  }

  return h;
}


/*
 * 3D Non-square
 */

static hnk_t *
hnk_cartesian_nonsquare_3D_create_suboctree(int maxh_width,
					    int maxh_height,
					    int maxh_depth,
					    int maxk);

static void
reorder(int *a, int *b)
{
  int t;
  
  if ((*a) < (*b)) {
    t = (*a);
    *a = (*b);
    *b = t;
  }
}

hnk_t *
hnk_cartesian_nonsquare_3D_create(int maxh_width,
				  int maxh_height,
				  int maxh_depth,
				  int maxk)
{
  /*
   * The 7 top-level branches
   */
  hnk_t *sht = NULL;
  hnk_t *svt = NULL;
  hnk_t *sdt = NULL;
  hnk_t *sbb = NULL;
  hnk_t *shb = NULL;
  hnk_t *svb = NULL;
  hnk_t *sdb = NULL;
  
  int maxh;
  
  if (maxh_width < 1 ||
      maxh_height < 1 ||
      maxh_depth < 1 ||
      maxh_width > 16 ||
      maxh_height > 16 ||
      maxh_depth > 16) {
    ERROR("height/width/depth out of range (%d %d %d)\n",
	  maxh_width, maxh_height, maxh_depth);
    return NULL;
  }


  /*
   * Just return a standard 7,8 tree if it's actually square
   */
  if (maxh_width == maxh_height &&
      maxh_width == maxh_depth) {

    return hnk_cartesian_78_create(maxh_width,
				   maxk);

  }

  /*
   * Reorder components so that depth is smallest, width is largest. The tree structure is simply
   * rearranged if the order is flipped so the counting is unaffected. We do this to reduce the
   * number of edge cases we need to check. From here we know that width >= height >= depth and
   * that they are not all equal.
   */
  reorder(&maxh_width, &maxh_height);
  reorder(&maxh_height, &maxh_depth);
  reorder(&maxh_width, &maxh_height);

  /*
   * Double check
   */
  if (maxh_width < maxh_height ||
      maxh_height < maxh_depth) {
    ERROR("bug: failed to reorder heights\n");
    return NULL;
  }
  

  if (maxh_width == maxh_height) {

    /*
     * (width == height) > depth
     *
     * ??t subtrees are staggered, ??b trees are octrees
     */

    sbb = hnk_create_octary_tree(maxh_depth - 1,
				 maxk - 1);
    shb = sbb;
    svb = sbb;
    sdb = sbb;

    /* All ??t trees the same */
    sht = hnk_cartesian_nonsquare_3D_create_suboctree(maxh_width - 1,
						      maxh_height - 1,
						      maxh_depth,
						      maxk - 1);
    svt = sht;
    sdt = sht;

  } else {

    if (maxh_height == maxh_depth) {

      /* 
       * width > (height == depth)
       *
       * sht subtree staggered, rest are octrees
       */

      sht = hnk_cartesian_nonsquare_3D_create_suboctree(maxh_width - 1,
							maxh_height,
							maxh_depth,
							maxk - 1);

      sdt = hnk_create_octary_tree(maxh_height - 1,
				   maxk - 1);

      svt = sdt;
      sbb = sdt;
      shb = sdt;
      sdb = sdt;
      svb = sdt;

    } else {

      /*
       * width > height > depth
       *
       * s?b are octrees, s?t are staggered
       */

      sht = hnk_cartesian_nonsquare_3D_create_suboctree(maxh_width - 1,
							maxh_height,
							maxh_depth,
							maxk - 1);

      /* Note the use of height in width, this is because the primary constraint in this
       * direction is the height so we want to use that constraint.
       */
      sdt = hnk_cartesian_nonsquare_3D_create_suboctree(maxh_height - 1,
							maxh_height - 1,
							maxh_depth,
							maxk - 1);
      svt = sdt;

      sbb = hnk_create_octary_tree(maxh_depth - 1,
				   maxk - 1);
      shb = sbb;
      sdb = sbb;
      svb = sbb;
    }
  }

  maxh = maxh_width;

  if (sht == NULL ||
      sdt == NULL ||
      svt == NULL ||
      sbb == NULL ||
      shb == NULL ||
      sdb == NULL ||
      svb == NULL) {
    ERROR("failed to create one or more subtrees (%p %p %p %p %p %p %p)\n",
	  sht, sdt, svt, sbb, shb, sdb, svb);
    return NULL;
  }
  
  return hnk_create_aggregate(maxh,
			      maxk,
			      7,
			      7,
			      sht,
			      sdt,
			      svt,
			      sbb,
			      shb,
			      sdb,
			      svb);
}

static hnk_t *
hnk_cartesian_nonsquare_3D_create_suboctree(int maxh_width,
					    int maxh_height,
					    int maxh_depth,
					    int maxk)
{
  hnk_t *sh;
  hnk_t *sv;
  hnk_t *sd;

  /* 
   * Straight octree now
   */
  if (maxh_width == maxh_height &&
      maxh_width == maxh_depth) {

    return hnk_create_octary_tree(maxh_width,
				  maxk);

  }

  /*
   * Double checker ordering of heights is still correct
   */
  if (maxh_width < maxh_height ||
      maxh_height < maxh_depth) {
    ERROR("width,height,depth out of order (%d %d %d)",
	  maxh_width, maxh_height, maxh_depth);
    return NULL;
  }

  if (maxh_width == 0) {
    if (maxh_height == 0) {
      return hnk_create_binary_tree(maxh_depth, maxk);
    } else if (maxh_depth == 0) {
      return hnk_create_binary_tree(maxh_height, maxk);
    } else {
      return hnk_cartesian_nonsquare_2D_create_subquadtree(maxh_height, maxh_depth, maxk);
    }
  }

  if (maxh_height == 0) {
    if (maxh_width == 0) {
      return hnk_create_binary_tree(maxh_depth, maxk);
    } else if (maxh_depth == 0) {
      return hnk_create_binary_tree(maxh_width, maxk);
    } else {
      return hnk_cartesian_nonsquare_2D_create_subquadtree(maxh_width, maxh_depth, maxk);
    }
  }

  if (maxh_depth == 0) {
    if (maxh_width == 0) {
      return hnk_create_binary_tree(maxh_height, maxk);
    } else if (maxh_height == 0) {
      return hnk_create_binary_tree(maxh_width, maxk);
    } else {
      return hnk_cartesian_nonsquare_2D_create_subquadtree(maxh_width, maxh_height, maxk);
    }
  }

  if (maxh_width > maxh_height) {

    if (maxh_height > maxh_depth) {

      /* width > height > depth */
      sd = hnk_create_octary_tree(maxh_depth - 1, maxk);

      sv = hnk_cartesian_nonsquare_3D_create_suboctree(maxh_height - 1,
						       maxh_height - 1,
						       maxh_depth,
						       maxk);

      sh = hnk_cartesian_nonsquare_3D_create_suboctree(maxh_width - 1,
						       maxh_height,
						       maxh_depth,
						       maxk);

      return hnk_create_aggregate(maxh_width,
				  maxk,
				  8,
				  8,
				  sh, sh, sv, sv, sd, sd, sd, sd);

      


    } else {

      /* width > (height == depth) */

      sd = hnk_create_octary_tree(maxh_depth - 1, maxk);
      sv = sd;

      sh = hnk_cartesian_nonsquare_3D_create_suboctree(maxh_width - 1,
						       maxh_height,
						       maxh_depth,
						       maxk);
      
      return hnk_create_aggregate(maxh_width,
				  maxk,
				  8,
				  8,
				  sh, sh, sv, sv, sd, sd, sd, sd);

    }

  } else {

    /* (width == height) > depth */

    sd = hnk_create_octary_tree(maxh_depth - 1, maxk);

    sh = hnk_cartesian_nonsquare_3D_create_suboctree(maxh_width - 1,
						     maxh_height - 1,
						     maxh_depth,
						     maxk);
    sv = sh;
    
    return hnk_create_aggregate(maxh_width,
				maxk,
				8,
				8,
				sh, sh, sv, sv, sd, sd, sd, sd);

  }
      
  ERROR("unreachable");
  return NULL;
}

hnk_t *
hnk_cartesian_nonsquare_3D_create_sub(int maxh_width,
				      int maxh_height,
				      int maxh_depth,
				      int maxk)
{
  int maxh_min;
  int subcount;

  hnk_t *h78;
  hnk_t *h;

  int i;

  maxh_min = maxh_width;
  if (maxh_height < maxh_min) {
    maxh_min = maxh_height;
  }
  if (maxh_depth < maxh_min) {
    maxh_min = maxh_depth;
  }

  subcount =
    (1 << (maxh_width - maxh_min)) *
    (1 << (maxh_height - maxh_min)) *
    (1 << (maxh_depth - maxh_min));

  if (subcount == 1) {
    return hnk_cartesian_78_create(maxh_width, maxk);
  }

  h78 = hnk_cartesian_78_create(maxh_min, maxk - 1);
  if (h78 == NULL) {
    return NULL;
  }

  h = hnk_create_aggregate_empty(maxh_min + 1,
				 maxk,
				 subcount,
				 subcount);
  if (h == NULL) {
    return NULL;
  }

  for (i = 0; i < subcount; i ++) {
    if (hnk_aggregate_set_subtree(h, i, h78) < 0) {
      return NULL;
    }
  }

  return h;
}

hnk_t *
hnk_cartesian_nonsquare_3D_create_spectral(int nbase,
					   int maxh,
					   int maxk)
{
  hnk_t *base;
  hnk_t *oct;
  int i;

  oct = hnk_create_octary_tree(maxh - 1, maxk - 1);
  if (oct == NULL) {
    return NULL;
  }
  
  base = hnk_create_aggregate_empty(maxh, maxk, nbase, nbase);
  if (base == NULL) {
    return NULL;
  }

  for (i = 0; i < nbase; i ++) {
    if (hnk_aggregate_set_subtree(base, i, oct) < 0) {
      return NULL;
    }
  }

  return base;
}
