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

#ifndef daubechies_h
#define daubechies_h

int
daubechies2d_forward_d4(double *s,
		      int width,
		      int height,
		      int stride,
		      double *work);

int
daubechies2d_inverse_d4(double *s,
		      int width,
		      int height,
		      int stride,
		      double *work);

int
daubechies1d_forward_d4_step(double *s,
			     int width,
			     int step,
			     double *work);

int
daubechies1d_inverse_d4_step(double *s,
			     int width,
			     int step,
			     double *work);


#endif /* daubechies_h */
