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

#ifndef sphericalwavelet_h
#define sphericalwavelet_h

typedef struct _sphericalwavelet sphericalwavelet_t;

typedef enum {
  SPHERICALWAVELET_TETRAHEDRON,
  SPHERICALWAVELET_OCTAHEDRON,
  SPHERICALWAVELET_ICOSAHEDRON
} sphericalwavelet_base_t;

sphericalwavelet_t *
sphericalwavelet_create(int degree);

void
sphericalwavelet_destroy(sphericalwavelet_t *s);

/*
 * Coefficient count for tetrahedron subdivision
 */
int sw_tetra_ncoeff(int depth);
int sw_tetra_totalcoeff(int depth);

/*
 * Indexing for tetrahedron subdivision
 */
int sw_tetra_depth_of_index(int index);


#endif /* sphericalwavelet_h */
