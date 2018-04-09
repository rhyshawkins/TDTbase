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

#ifndef cdf97_lift_periodic_h
#define cdf97_lift_periodic_h

/*
 * 1D Forward
 */
int
cdf97_lift_periodic_forward1d_cdf97(double *s,
			   int width,
			   int stride,
			   double *work);

int
cdf97_lift_periodic_forward1d_cdf97_step(double *s,
				int width,
				int stride,
				double *work);

/*
 * 1D Inverse
 */
int
cdf97_lift_periodic_inverse1d_cdf97(double *s,
			   int width,
			   int stride,
			   double *work);

int
cdf97_lift_periodic_inverse1d_cdf97_step(double *s,
				int width,
				int stride,
				double *work);

/*
 * 2D Forward
 */
int 
cdf97_lift_periodic_forward2d_cdf97(double *s, 
				    int width,
				    int height,
				    int stride,
				    double *work,
				    int subtile);

int 
cdf97_lift_periodic_forward2d_cdf97_step(double *s, 
					 int width,
					 int height,
					 int stride,
					 double *work);

/*
 * 2D Inverse
 */
int
cdf97_lift_periodic_inverse2d_cdf97(double *s,
				    int width,
				    int height,
				    int stride,
				    double *work,
				    int subtile);

int 
cdf97_lift_periodic_inverse2d_cdf97_step(double *s, 
					 int width,
					 int height,
					 int stride,
					 double *work);

/*
 * 3D Forward
 */
int 
cdf97_lift_periodic_forward3d_cdf97(double *s,
				    int width,
				    int height,
				    int depth,
				    int rowstride,
				    int slicestride,
				    double *work,
				    int subtile);

int 
cdf97_lift_periodic_forward3d_cdf97_step(double *s,
					 int width,
					 int height,
					 int depth,
					 int rowstride,
					 int slicestride,
					 double *work);

int 
cdf97_lift_periodic_forward3d_cdf97_2dstep(double *s,
					   int width,
					   int height,
					   int stride,
					   int rowstride,
					   double *work);


/*
 * 3D Inverse
 */
int 
cdf97_lift_periodic_inverse3d_cdf97(double *s,
				    int width,
				    int height,
				    int depth,
				    int rowstride,
				    int slicestride,
				    double *work,
				    int subtile);

int 
cdf97_lift_periodic_inverse3d_cdf97_step(double *s,
				int width,
				int height,
				int depth,
				int rowstride,
				int slicestride,
				double *work);

int 
cdf97_lift_periodic_inverse3d_cdf97_2dstep(double *s,
				  int width,
				  int height,
				  int stride,
				  int rowstride,
				  double *work);

#endif /* cdf97_lift_periodic_h */
