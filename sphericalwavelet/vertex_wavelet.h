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

#ifndef vertex_wavelet_h
#define vertex_wavelet_h

#include "manifold.h"

/*
 * Full transforms
 */
int
vertex_wavelet_butterfly_forward(manifold_t *m,
				 double *coeff);

int
vertex_wavelet_butterfly_inverse(manifold_t *m,
				 double *coeff);

int
vertex_wavelet_butterfly_forward_lifted(manifold_t *m,
					double *coeff);

int
vertex_wavelet_butterfly_inverse_lifted(manifold_t *m,
					double *coeff);

/*
 * Individual steps
 */
int
vertex_wavelet_butterfly_forward_step(manifold_t *m,
				      double *coeff,
				      int depth);

int
vertex_wavelet_butterfly_forward_lifted_step(manifold_t *m,
					     double *coeff,
					     int depth);

int
vertex_wavelet_butterfly_inverse_step(manifold_t *m,
				      double *coeff,
				      int depth);

int
vertex_wavelet_butterfly_inverse_lifted_step(manifold_t *m,
					     double *coeff,
					     int depth);

/*
 * Shell transforms : radial degree == manifold degree
 */
int
vertex_wavelet_butterfly_shell_forward(manifold_t *m,
				       double *coeff,
				       int ncoeff,
				       double *workspace);

int
vertex_wavelet_butterfly_shell_inverse(manifold_t *m,
				       double *coeff,
				       int ncoeff, 
				       double *workspace);

int
vertex_wavelet_butterfly_shell_forward_lifted(manifold_t *m,
					      double *coeff,
					      int ncoeff,
					      double *workspace);

int
vertex_wavelet_butterfly_shell_inverse_lifted(manifold_t *m,
					      double *coeff,
					      int ncoeff,
					      double *workspace);



#endif /* vertex_wavelet_h */
