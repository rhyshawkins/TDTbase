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

#ifndef octahedron_h
#define octahedron_h

#include "manifold.h"

manifold_t *
octahedron_create(int degree);

void
octahedron_destroy(manifold_t *o);

int
octahedron_parent(manifold_t *o, int index);

int
octahedron_children(manifold_t *o, int index, int *child_indices);

/*
 * Printing functions
 */
void
octahedron_print_triangles(manifold_t *o, int depth);

/*
 * Counting functions
 */
int 
octahedron_nvertices(int depth);

int 
octahedron_nedges(int depth);

int 
octahedron_ntriangles(int depth);



#endif /* octahedron_h */
