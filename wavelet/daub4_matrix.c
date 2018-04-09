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

#include <stdio.h>
#include <math.h>

#include <gsl/gsl_blas.h>

#include "daub4_matrix.h"
#include "generic_matrix.h"
#include "boundary.h"

/* #define H0 0.48296291314453414337487159986 */
/* #define H1 0.83651630373780790557529378092 */
/* #define H2 0.22414386804201338102597276224 */
/* #define H3 -0.12940952255126038117444941881 */

/*
 * These coefficients are scaled by 1/sqrt(2) from the true orthongal terms to match
 * the normalization of the lifting version.
 */
#define H0 0.34150635094610965
#define H1 0.5915063509461096
#define H2 0.15849364905389032
#define H3 -9.150635094610965e-2

#define G0 H3
#define G1 -H2
#define G2 H1
#define G3 -H0

static const int NFORWARD_SCALING = 4;
static const double FORWARD_SCALING[4] = {
  H0, H1, H2, H3};
static const int FORWARD_SCALING_OFFSET[4] = {
  0, 1, 2, 3};

static const int NFORWARD_WAVELET = 4;
static const double FORWARD_WAVELET[4] = {
  G0, G1, G2, G3};
static const int FORWARD_WAVELET_OFFSET[4] = {
  -1, 0, 1, 2};

#define HB0 0.6830127018922194
#define HB1 1.1830127018922194
#define HB2 0.3169872981077807
#define HB3 -0.18301270189221933

#define GB0 HB3
#define GB1 -HB2
#define GB2 HB1
#define GB3 -HB0

static const int NINVERSE_SCALING = 4;
static const double INVERSE_SCALING[4] = {
  HB2, HB1, HB0, HB3};
static const int INVERSE_SCALING_OFFSET[4] = {
  -2, -1, 0, 1};

static const int NINVERSE_WAVELET = 4;
static const double INVERSE_WAVELET[4] = {
  GB0, GB3, GB2, GB1};
static const int INVERSE_WAVELET_OFFSET[4] = {
  -3, -2, -1, 0};

int
daub4_matrix_forward2d_create_col_step(int j_max, int j, gsl_matrix **_m)
{
  return generic_matrix_create_forward2d_col_substep(j_max,
						     j,
						     NFORWARD_SCALING,
						     FORWARD_SCALING,
						     FORWARD_SCALING_OFFSET,
						     NFORWARD_WAVELET,
						     FORWARD_WAVELET,
						     FORWARD_WAVELET_OFFSET,
						     _m,
						     wavelet_boundary_periodic);
}

int
daub4_matrix_inverse2d_create_col_step(int j_max, int j, gsl_matrix **_m)
{
  return generic_matrix_create_inverse2d_col_substep(j_max,
						     j,
						     NINVERSE_SCALING,
						     INVERSE_SCALING,
						     INVERSE_SCALING_OFFSET,
						     NINVERSE_WAVELET,
						     INVERSE_WAVELET,
						     INVERSE_WAVELET_OFFSET,
						     _m,
						     wavelet_boundary_periodic);
}

int
daub4_matrix_forward2d_create_row_step(int j_max, int j, gsl_matrix **_m)
{
  return generic_matrix_create_forward2d_row_substep(j_max,
						     j,
						     NFORWARD_SCALING,
						     FORWARD_SCALING,
						     FORWARD_SCALING_OFFSET,
						     NFORWARD_WAVELET,
						     FORWARD_WAVELET,
						     FORWARD_WAVELET_OFFSET,
						     _m,
						     wavelet_boundary_periodic);
}

int
daub4_matrix_inverse2d_create_row_step(int j_max, int j, gsl_matrix **_m)
{
  return generic_matrix_create_inverse2d_row_substep(j_max,
						     j,
						     NINVERSE_SCALING,
						     INVERSE_SCALING,
						     INVERSE_SCALING_OFFSET,
						     NINVERSE_WAVELET,
						     INVERSE_WAVELET,
						     INVERSE_WAVELET_OFFSET,
						     _m,
						     wavelet_boundary_periodic);
}

int
daub4_matrix_forward2d_create_step(int j_max, int j, gsl_matrix **_m)
{
  return generic_matrix_create_forward2d_substep(j_max,
						 j,
						 NFORWARD_SCALING,
						 FORWARD_SCALING,
						 FORWARD_SCALING_OFFSET,
						 NFORWARD_WAVELET,
						 FORWARD_WAVELET,
						 FORWARD_WAVELET_OFFSET,
						 _m,
						 wavelet_boundary_periodic);
}

int
daub4_matrix_inverse2d_create_step(int j_max, int j, gsl_matrix **_m)
{
  return generic_matrix_create_inverse2d_substep(j_max,
						 j,
						 NINVERSE_SCALING,
						 INVERSE_SCALING,
						 INVERSE_SCALING_OFFSET,
						 NINVERSE_WAVELET,
						 INVERSE_WAVELET,
						 INVERSE_WAVELET_OFFSET,
						 _m,
						 wavelet_boundary_periodic);
}

int
daub4_matrix_forward2d_create(int j_max, gsl_matrix **_m)
{
  return generic_matrix_create_forward2d(j_max,
					 NFORWARD_SCALING,
					 FORWARD_SCALING,
					 FORWARD_SCALING_OFFSET,
					 NFORWARD_WAVELET,
					 FORWARD_WAVELET,
					 FORWARD_WAVELET_OFFSET,
					 _m,
					 wavelet_boundary_periodic);
}

int
daub4_matrix_inverse2d_create(int j_max, gsl_matrix **_m)
{
  return generic_matrix_create_inverse2d(j_max,
					 NINVERSE_SCALING,
					 INVERSE_SCALING,
					 INVERSE_SCALING_OFFSET,
					 NINVERSE_WAVELET,
					 INVERSE_WAVELET,
					 INVERSE_WAVELET_OFFSET,
					 _m,
					 wavelet_boundary_periodic);

}


