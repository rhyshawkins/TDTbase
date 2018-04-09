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

#ifndef cdf97_lift_impulse_h
#define cdf97_lift_impulse_h

/*
 * Size computations
 */
int cdf97_lift_impulse_1dsize(int degree, int impulse_degree);

/*
 * 1D Impulse Functions
 */
int cdf97_lift_impulse_1d(int degree, 
			  int impulse_degree, 
			  double impulse, 
			  double *v, 
			  int size,
			  int offset);


#endif /* cdf97_lift_impulse_h */
