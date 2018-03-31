
#include <stdio.h>

#include <gsl/gsl_blas.h>

#include "cdf97_matrix.h"
#include "generic_matrix.h"
#include "boundary.h"

//static const double H[5] = {0.602949018236,  0.266864118443, -0.078223266529, -0.016864118443, 0.026748757411};
//static const double G[4] = {0.557543526229, -0.295635881557, -0.028771763114,  0.045635881557};

#define H0 0.602949018236
#define H1 0.266864118443
#define H2 -0.078223266529
#define H3 -0.016864118443
#define H4 0.026748757411

#define G0 0.557543526229
#define G1 -0.295635881557
#define G2 -0.028771763114
#define G3 0.045635881557

static const int NFORWARD_SCALING = 9;
static const double FORWARD_SCALING[9] = {
  H4, H3, H2, H1, H0, H1, H2, H3, H4};
static const int FORWARD_SCALING_OFFSET[9] = {
  -4, -3, -2, -1, 0, 1, 2, 3, 4};

static const int NFORWARD_WAVELET = 7;
static const double FORWARD_WAVELET[7] = {
  G3, G2, G1, G0, G1, G2, G3};
static const int FORWARD_WAVELET_OFFSET[7] = {
  -3, -2, -1, 0, 1, 2, 3};
  //  -4, -3, -2, -1, 0, 1, 2};

//static const double HB[4] = {1.115087052458, 0.591271763114, -5.7543526228e-2, -9.1271763114e-2};
//static const double GB[5] = {1.205898036472, -0.533728236886, -0.156446533058, 3.3728236886e-2, 5.3497514822e-2};

#define HB0 1.115087052458
#define HB1 0.591271763114
#define HB2 -5.7543526228e-2
#define HB3 -9.1271763114e-2

#define GB0 1.205898036472
#define GB1 -0.533728236886
#define GB2 -0.156446533058
#define GB3 3.3728236886e-2
#define GB4 5.3497514822e-2

static const int NINVERSE_SCALING = 7;
static const double INVERSE_SCALING[7] = {
  GB3, HB2, GB1, HB0, GB1, HB2, GB3};
static const int INVERSE_SCALING_OFFSET[7] = {
  -3, -2, -1, 0, 1, 2, 3};

static const int NINVERSE_WAVELET = 9;
static const double INVERSE_WAVELET[9] = {
  GB4, HB3, GB2, HB1, GB0, HB1, GB2, HB3, GB4};
static const int INVERSE_WAVELET_OFFSET[9] = {
  -4, -3, -2, -1, 0, 1, 2, 3, 4};

int
cdf97_matrix_forward2d_create_col_step(int j_max, int j, gsl_matrix **_m)
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
cdf97_matrix_inverse2d_create_col_step(int j_max, int j, gsl_matrix **_m)
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
cdf97_matrix_forward2d_create_row_step(int j_max, int j, gsl_matrix **_m)
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
cdf97_matrix_inverse2d_create_row_step(int j_max, int j, gsl_matrix **_m)
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
cdf97_matrix_forward2d_create_step(int j_max, int j, gsl_matrix **_m)
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
cdf97_matrix_inverse2d_create_step(int j_max, int j, gsl_matrix **_m)
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
cdf97_matrix_forward2d_create(int j_max, gsl_matrix **_m)
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
cdf97_matrix_inverse2d_create(int j_max, gsl_matrix **_m)
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


