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
#ifndef hnk_cartesian_nonsquare_h
#define hnk_cartesian_nonsquare_h

#include "hnk.h"

/*
 * 2D Non-square
 */
hnk_t *
hnk_cartesian_nonsquare_2D_create(int maxh_width,
				  int maxh_height,
				  int maxk);

hnk_t *
hnk_cartesian_nonsquare_2D_create_sub(int maxh_width,
				      int maxh_height,
				      int maxk);


/*
 * 3D Non-square
 */
hnk_t *
hnk_cartesian_nonsquare_3D_create(int maxh_width,
				  int maxh_height,
				  int maxh_depth,
				  int maxk);

hnk_t *
hnk_cartesian_nonsquare_3D_create_sub(int maxh_width,
				      int maxh_height,
				      int maxh_depth,
				      int maxk);

/*
 * 3D Spectral Element Cell meshes
 */
hnk_t *
hnk_cartesian_nonsquare_3D_create_spectral(int nbase,
					   int maxh,
					   int maxk);

#endif /* hnk_cartesian_nonsquare_h */
