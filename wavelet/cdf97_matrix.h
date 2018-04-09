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

#ifndef cdf97_matrix_h
#define cdf97_matrix_h

#include <gsl/gsl_matrix.h>

int
cdf97_matrix_forward2d_create_col_step(int j_max, int j, gsl_matrix **_m);

int
cdf97_matrix_inverse2d_create_col_step(int j_max, int j, gsl_matrix **_m);

int
cdf97_matrix_forward2d_create_row_step(int j_max, int j, gsl_matrix **_m);

int
cdf97_matrix_inverse2d_create_row_step(int j_max, int j, gsl_matrix **_m);

int
cdf97_matrix_forward2d_create_step(int j_max, int j, gsl_matrix **_m);

int
cdf97_matrix_inverse2d_create_step(int j_max, int j, gsl_matrix **_m);

int
cdf97_matrix_forward2d_create(int j_max, gsl_matrix **m);

int
cdf97_matrix_inverse2d_create(int j_max, gsl_matrix **_m);


#endif /* cdf97_matrix_h */
