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
#ifndef hnk_cartesian_h
#define hnk_cartesian_h

#include "hnk.h"

/*
 * 1,2 Tree for 1D
 */

hnk_t *
hnk_cartesian_12_create(int maxh,
			int maxk);

/*
 * 1,3 Tree for 3D Triangular subdivision (poles)
 */

hnk_t *
hnk_cartesian_13_create(int maxh,
			int maxk);

/*
 * 3,4 Tree for 2D
 */

hnk_t *
hnk_cartesian_34_create(int maxh,
			int maxk);
  
/*
 * 7,8 Tree for 3D
 */

hnk_t *
hnk_cartesian_78_create(int maxh,
			int maxk);
  

/*
 * 7,9 Tree for 3D Triangular subdivision
 */

hnk_t *
hnk_cartesian_79_create(int maxh,
			int maxk);

#endif /* hnk_cartesian_h */
