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

#ifndef haar_lift_h
#define haar_lift_h

/*
 * 1D Forward
 */
int
haar_lift_forward1d_haar(double *s,
			 int width,
			 int stride,
			 double *work);

int
haar_lift_forward1d_haar_step(double *s,
			      int width,
			      int stride,
			      double *work);

/*
 * 1D Inverse
 */
int
haar_lift_inverse1d_haar(double *s,
			 int width,
			 int stride,
			 double *work);

int
haar_lift_inverse1d_haar_step(double *s,
			      int width,
			      int stride,
			      double *work);

/*
 * 2D Forward
 */
int 
haar_lift_forward2d_haar(double *s, 
			 int width,
			 int height,
			 int stride,
			 double *work,
			 int subtile);

int 
haar_lift_forward2d_haar_step(double *s, 
			      int width,
			      int height,
			      int stride,
			      double *work);

/*
 * 2D Inverse
 */
int
haar_lift_inverse2d_haar(double *s,
			 int width,
			 int height,
			 int stride,
			 double *work,
			 int subtile);

int 
haar_lift_inverse2d_haar_step(double *s, 
			      int width,
			      int height,
			      int stride,
			      double *work);

/*
 * 3D Forward
 */
int 
haar_lift_forward3d_haar(double *s,
			 int width,
			 int height,
			 int depth,
			 int rowstride,
			 int slicestride,
			 double *work,
			 int subtile);

int 
haar_lift_forward3d_haar_step(double *s,
			      int width,
			      int height,
			      int depth,
			      int rowstride,
			      int slicestride,
			      double *work);

int 
haar_lift_forward3d_haar_2dstep(double *s,
				int width,
				int height,
				int stride,
				int rowstride,
				double *work);


/*
 * 3D Inverse
 */
int 
haar_lift_inverse3d_haar(double *s,
			 int width,
			 int height,
			 int depth,
			 int rowstride,
			 int slicestride,
			 double *work,
			 int subtile);

int 
haar_lift_inverse3d_haar_step(double *s,
			      int width,
			      int height,
			      int depth,
			      int rowstride,
			      int slicestride,
			      double *work);

int 
haar_lift_inverse3d_haar_2dstep(double *s,
				int width,
				int height,
				int stride,
				int rowstride,
				double *work);

#endif /* haar_lift_h */
