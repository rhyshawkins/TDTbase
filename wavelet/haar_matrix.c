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

#include "haar_matrix.h"
#include "generic_matrix.h"
#include "boundary.h"

/*
 * This is a Haar wavelet normalized on the forward transform only which means
 * that the scaling terms are "averages" of the next higher resolution level.
 * This is different from the normalization using 1/sqrt(2) which then has
 * the H2 matrix [H0 H1; G0 G1] = [HB0 HB1; GB0 GB1].
 */
#define H0 0.5
#define H1 0.5

#define G0 0.5
#define G1 -0.5

static const int NFORWARD_SCALING = 2;
static const double FORWARD_SCALING[2] = {
  H0, H1};
static const int FORWARD_SCALING_OFFSET[2] = {
  0, 1};

static const int NFORWARD_WAVELET = 2;
static const double FORWARD_WAVELET[2] = {
  G0, G1};
static const int FORWARD_WAVELET_OFFSET[2] = {
  -1, 0};

#define HB0 1.0
#define HB1 1.0

#define GB0 1.0
#define GB1 -1.0

static const int NINVERSE_SCALING = 2;
static const double INVERSE_SCALING[2] = {
  HB0, HB1};
static const int INVERSE_SCALING_OFFSET[2] = {
  0, 1};

static const int NINVERSE_WAVELET = 2;
static const double INVERSE_WAVELET[2] = {
  GB0, GB1};
static const int INVERSE_WAVELET_OFFSET[2] = {
  -1, 0};

int
haar_matrix_forward2d_create_col_step(int j_max, int j, gsl_matrix **_m)
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
						     wavelet_boundary_reflect);
}

int
haar_matrix_inverse2d_create_col_step(int j_max, int j, gsl_matrix **_m)
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
						     wavelet_boundary_reflect);
}

int
haar_matrix_forward2d_create_row_step(int j_max, int j, gsl_matrix **_m)
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
						     wavelet_boundary_reflect);
}

int
haar_matrix_inverse2d_create_row_step(int j_max, int j, gsl_matrix **_m)
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
						     wavelet_boundary_reflect);
}

int
haar_matrix_forward2d_create_step(int j_max, int j, gsl_matrix **_m)
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
						 wavelet_boundary_reflect);
}

int
haar_matrix_inverse2d_create_step(int j_max, int j, gsl_matrix **_m)
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
						 wavelet_boundary_reflect);
}

int
haar_matrix_forward2d_create(int j_max, gsl_matrix **_m)
{
  return generic_matrix_create_forward2d(j_max,
					 NFORWARD_SCALING,
					 FORWARD_SCALING,
					 FORWARD_SCALING_OFFSET,
					 NFORWARD_WAVELET,
					 FORWARD_WAVELET,
					 FORWARD_WAVELET_OFFSET,
					 _m,
					 wavelet_boundary_reflect);
}

int
haar_matrix_inverse2d_create(int j_max, gsl_matrix **_m)
{
  return generic_matrix_create_inverse2d(j_max,
					 NINVERSE_SCALING,
					 INVERSE_SCALING,
					 INVERSE_SCALING_OFFSET,
					 NINVERSE_WAVELET,
					 INVERSE_WAVELET,
					 INVERSE_WAVELET_OFFSET,
					 _m,
					 wavelet_boundary_reflect);

}


