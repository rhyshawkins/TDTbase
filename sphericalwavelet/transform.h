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

#ifndef transform_h
#define transform_h

void
transform_from_axis_angle(double x, double y, double z,
			  double theta,
			  double *transform);

void
transform_compute_axis_angle(double x1, double y1, double z1,
			     double x2, double y2, double z2,
			     double *ax, double *ay, double *az,
			     double *theta);
			     
#endif /* transform_h */
