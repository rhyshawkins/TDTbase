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

#ifndef face_subdivision_h
#define face_subdivision_h

#include "manifold.h"

/*
 * This ia a non-wavelet subdivision scheme. For a given manifold we define the forward transform
 * as the mean of the 4 sub-triangles is the value of the parent triangle, and the 4 differences
 * from the mean are set to the 4 sub-triangles. This is similar to a Haar wavelet.
 */

/*
 * Full transforms
 */
int
face_subdivision_forward(manifold_t *m,
			 double *coeff);

int
face_subdivision_inverse(manifold_t *m,
			 double *coeff);

/*
 * Individual steps
 */
int
face_subdivision_forward_step(manifold_t *m,
			      double *coeff,
			      int depth);

int
face_subdivision_inverse_step(manifold_t *m,
			      double *coeff,
			      int depth);

#endif /* face_subdivision_h */
