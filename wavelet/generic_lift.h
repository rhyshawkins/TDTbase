//
//    Wavelet transform library
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

#ifndef generic_lift_h
#define generic_lift_h


typedef int (*generic_lift_forward1d_step_t)(double *s,
					     int width,
					     int stride,
					     double *work);

typedef int (*generic_lift_inverse1d_step_t)(double *s,
					     int width,
					     int stride,
					     double *work);

/*
 * 1D Full Transform
 */
int
generic_lift_forward1d(double *s,
		       int width,
		       int stride,
		       double *work,
		       generic_lift_forward1d_step_t transform);

int
generic_lift_inverse1d(double *s,
		       int width,
		       int stride,
		       double *work,
		       generic_lift_inverse1d_step_t transform);

/*
 * 2D Step
 */
int
generic_lift_inverse2d_step(double *s,
			    int width,
			    int height,
			    int stride,
			    double *work,
			    generic_lift_forward1d_step_t row_transform,
			    generic_lift_forward1d_step_t col_transform);

/*
 * 2D Full Transform
 */
int
generic_lift_forward2d(double *s,
		       int width,
		       int height,
		       int stride,
		       double *work,
		       generic_lift_forward1d_step_t row_transform,
		       generic_lift_forward1d_step_t col_transform,
		       int subtile);

int
generic_lift_inverse2d(double *s,
		       int width,
		       int height,
		       int stride,
		       double *work,
		       generic_lift_inverse1d_step_t row_transform,
		       generic_lift_inverse1d_step_t col_transform,
		       int subtile);

/*
 * 3D Full Transform
 */
int
generic_lift_forward3d(double *s,
		       int width,
		       int height,
		       int depth,
		       int stride,
		       int slicestride,
		       double *work,
		       generic_lift_forward1d_step_t row_transform,
		       generic_lift_forward1d_step_t col_transform,
		       generic_lift_forward1d_step_t dep_transform,
		       int subtile);


int
generic_lift_inverse3d(double *s,
		       int width,
		       int height,
		       int depth,
		       int stride,
		       int slicestride,
		       double *work,
		       generic_lift_inverse1d_step_t row_transform,
		       generic_lift_inverse1d_step_t col_transform,
		       generic_lift_inverse1d_step_t dep_transform,
		       int subtile);


#endif /* generic_lift_h */
