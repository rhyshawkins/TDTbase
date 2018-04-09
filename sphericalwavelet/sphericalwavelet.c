//
//    Spherical Subdivision/Wavelet library
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
#include <string.h>

#include "sphericalwavelet.h"

sphericalwavelet_t *
sphericalwavelet_create(int degree)
{
  return NULL;
}

void
sphericalwavelet_destroy(sphericalwavelet_t *s)
{
}

/*
 * Coefficient count for tetrahedron subdivision
 */
int sw_tetra_ncoeff(int depth)
{
  // Progression is 4 6 24
  return -1;
}

int sw_tetra_totalcoeff(int depth)
{
  return -1;
}

/*
 * Indexing for tetrahedron subdivision
 */
int sw_tetra_depth_of_index(int index)
{
  return -1;
}
