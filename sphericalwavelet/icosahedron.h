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

#ifndef icosahedron_h
#define icosahedron_h

#include "manifold.h"

manifold_t *
icosahedron_create(int degree);

/*
 * Counting functions
 */
int 
icosahedron_nvertices(int depth);

int 
icosahedron_nedges(int depth);

int 
icosahedron_ntriangles(int depth);

/*
 * Edge subtended angle -> resolution at depth
 */
double
icosahedron_angle(int depth);



#endif /* icosahedron_h */
