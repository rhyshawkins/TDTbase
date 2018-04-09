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

#ifndef face_wavelet_h
#define face_wavelet_h

#include "manifold.h"

/*
 * Full transforms
 */
int
face_wavelet_biohaar_forward(manifold_t *m,
			     double *coeff);

int
face_wavelet_biohaar_inverse(manifold_t *m,
			     double *coeff);

/*
 * Individual steps
 */
int
face_wavelet_biohaar_forward_step(manifold_t *m,
			      double *coeff,
			      int depth);

int
face_wavelet_biohaar_inverse_step(manifold_t *m,
				  double *coeff,
				  int depth);

/*
 * Shell transforms : radial degree == manifold degree
 */

typedef int (*face_wavelet_radial_forward_step_t)(double *s,
						  int width,
						  int stride,
						  double *work);
typedef int (*face_wavelet_radial_inverse_step_t)(double *s,
						  int width,
						  int stride,
						  double *work);

  
int
face_wavelet_biohaar_shell_forward(manifold_t *m,
				   double *coeff,
				   int ncoeff,
				   double *workspace,
				   face_wavelet_radial_forward_step_t radial_forward_step);

int
face_wavelet_biohaar_shell_inverse(manifold_t *m,
				   double *coeff,
				   int ncoeff, 
				   double *workspace,
				   face_wavelet_radial_inverse_step_t radial_inverse_step);


#endif /* face_wavelet_h */
