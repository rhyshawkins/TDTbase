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
#ifndef hnk_healpix_h
#define hnk_healpix_h

#include "hnk.h"

/*
 * These are subdivision trees that start with 12 top level children
 * followed by 4 (2d) or 8 (3d) thereafter.
 */
hnk_t *
hnk_healpix2d_create(int maxh,
		     int maxk);

hnk_t *
hnk_healpix3d_create(int maxh,
		     int maxk);


#endif /* hnk_healpix_h */
