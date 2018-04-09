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

#ifndef radialbasisfunction_h
#define radialbasisfunction_h

typedef double (*radialbasisfunction2D_t)(double cx, double cy, double r,
					  double x, double y);

typedef double (*radialbasisfunction3D_t)(double cx, double cy, double cz, double rxy, double rz,
					  double x, double y, double z);

/*
 * Linear
 */
double rbf_linear1D(double x, double r);

double rbf_linear2D(double cx, double cy, double r,
		    double x, double y);

double rbf_linear2DS2(double clon, double clat, double r,
		      double lon, double lat);

double rbf_linear3D(double cx, double cy, double dz, double rxy, double rz,
		    double x, double y, double z);

/*
 * Quadratic
 */
double rbf_quadratic2D(double cx, double cy, double r,
		       double x, double y);

double rbf_quadratic3D(double cx, double cy, double dz, double rxy, double rz,
		       double x, double y, double z);

/*
 * Piecewise Cubic smooth, 0 gradients at radius = 0, r
 */
double rbf_cubic2D(double cx, double cy, double r,
		   double x, double y);

double rbf_cubic3D(double cx, double cy, double dz, double rxy, double rz,
		   double x, double y, double z);

/*
 * Cosine: Peak of 1.0 at radius = 0, 
 */
double rbf_cosine2D(double cx, double cy, double r,
		    double x, double y);

double rbf_cosine3D(double cx, double cy, double dz, double rxy, double rz,
		    double x, double y, double z);

/*
 * Gaussian: Un-normalized, sigma = r
 */
double rbf_gaussian2D(double cx, double cy, double r,
		      double x, double y);

double rbf_gaussian3D(double cx, double cy, double dz, double rxy, double rz,
		      double x, double y, double z);

/*
 * Lanczos: 
 */
double rbf_lanczos2D(double cx, double cy, double r,
		     double x, double y);

double rbf_lanczos3D(double cx, double cy, double dz, double rxy, double rz,
		     double x, double y, double z);

#endif /* radialbasisfunction_h */
