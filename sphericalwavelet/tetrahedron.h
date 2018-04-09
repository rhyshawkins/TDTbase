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

#ifndef tetrahedron_h
#define tetrahedron_h

#include "manifold.h"

manifold_t *
tetrahedron_create(int degree);

/*
 * Counting functions
 */
int 
tetrahedron_nvertices(int depth);

int 
tetrahedron_nedges(int depth);

int 
tetrahedron_ntriangles(int depth);



#endif /* tetrahedron_h */
