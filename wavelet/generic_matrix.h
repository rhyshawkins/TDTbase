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

#ifndef generic_matrix_h
#define generic_matrix_h

#include <gsl/gsl_matrix.h>

/*
 * Helper functions
 */

int
generic_matrix_interlace_index(int i, int width);

int
generic_matrix_deinterlace_index(int i, int width);

/*
 * Forward Transforms
 */

typedef int (*generic_matrix_edge_func_t)(int i, int width);

int
generic_matrix_fill_forward2d_col_A(int j_max,
				    int si,
				    int nscaling,
				    const double *scaling,
				    const int *soffset,
				    gsl_matrix *m,
				    int row_offset,
				    int col_offset,
				    generic_matrix_edge_func_t edge);

int
generic_matrix_fill_forward2d_col_B(int j_max,
				    int wi,
				    int nwavelet,
				    const double *wavelet,
				    const int *woffset,
				    gsl_matrix *m,
				    int row_offset,
				    int col_offset,
				    generic_matrix_edge_func_t edge);

int
generic_matrix_fill_forward2d_col_step(int j_max,
				       int nscaling,
				       const double *scaling,
				       const int *soffset,
				       int nwavelet,
				       const double *wavelet,
				       const int *woffset,
				       gsl_matrix *m,
				       generic_matrix_edge_func_t edge);

int
generic_matrix_create_forward2d_col_step(int j_max,
					 int nscaling,
					 const double *scaling,
					 const int *soffset,
					 int nwavelet,
					 const double *wavelet,
					 const int *woffset,
					 gsl_matrix **_m,
					 generic_matrix_edge_func_t edge);

int
generic_matrix_create_forward2d_col_substep(int j_max,
					    int j0,
					    int nscaling,
					    const double *scaling,
					    const int *soffset,
					    int nwavelet,
					    const double *wavelet,
					    const int *woffset,
					    gsl_matrix **_m,
					    generic_matrix_edge_func_t edge);

int
generic_matrix_fill_forward2d_row_A(int j_max,
				    int nscaling,
				    const double *scaling,
				    const int *soffset,
				    gsl_matrix *m,
				    int row_offset,
				    int col_offset,
				    generic_matrix_edge_func_t edge);

int
generic_matrix_fill_forward2d_row_B(int j_max,
				    int nwavelet,
				    const double *wavelet,
				    const int *woffset,
				    gsl_matrix *m,
				    int row_offset,
				    int col_offset,
				    generic_matrix_edge_func_t edge);

int
generic_matrix_fill_forward2d_row_step(int j_max,
				       int nscaling,
				       const double *scaling,
				       const int *soffset,
				       int nwavelet,
				       const double *wavelet,
				       const int *woffset,
				       gsl_matrix *m,
				       generic_matrix_edge_func_t edge);

int
generic_matrix_create_forward2d_row_step(int j_max,
					 int nscaling,
					 const double *scaling,
					 const int *soffset,
					 int nwavelet,
					 const double *wavelet,
					 const int *woffset,
					 gsl_matrix **_m,
					 generic_matrix_edge_func_t edge);

int
generic_matrix_create_forward2d_row_substep(int j_max,
					    int j0,
					    int nscaling,
					    const double *scaling,
					    const int *soffset,
					    int nwavelet,
					    const double *wavelet,
					    const int *woffset,
					    gsl_matrix **_m,
					    generic_matrix_edge_func_t edge);

int
generic_matrix_fill_forward2d_step(int jmax,
				   int nscaling,
				   const double *scaling,
				   const int *soffset,
				   int nwavelet,
				   const double *wavelet,
				   const int *woffset,
				   gsl_matrix *m,
				   generic_matrix_edge_func_t edge);

int
generic_matrix_create_forward2d_step(int jmax,
				     int nscaling,
				     const double *scaling,
				     const int *soffset,
				     int nwavelet,
				     const double *wavelet,
				     const int *woffset,
				     gsl_matrix **m,
				     generic_matrix_edge_func_t edge);

int
generic_matrix_create_forward2d_substep(int jmax,
					int j0,
					int nscaling,
					const double *scaling,
					const int *soffset,
					int nwavelet,
					const double *wavelet,
					const int *woffset,
					gsl_matrix **m,
					generic_matrix_edge_func_t edge);


int
generic_matrix_create_forward2d(int j_max,
				int nscaling,
				const double *scaling,
				const int *soffset,
				int nwavelet,
				const double *wavelet,
				const int *woffset,
				gsl_matrix **m,
				generic_matrix_edge_func_t edge);


/*
 * Inverse Transforms
 */

int
generic_matrix_fill_inverse2d_col_A(int j_max,
				    int si,
				    int nscaling,
				    const double *scaling,
				    const int *soffset,
				    gsl_matrix *m,
				    int row_offset,
				    int col_offset,
				    generic_matrix_edge_func_t edge);

int
generic_matrix_fill_inverse2d_col_B(int j_max,
				    int wi,
				    int nwavelet,
				    const double *wavelet,
				    const int *woffset,
				    gsl_matrix *m,
				    int row_offset,
				    int col_offset,
				    generic_matrix_edge_func_t edge);

int
generic_matrix_fill_inverse2d_col_step(int j_max,
				       int nscaling,
				       const double *scaling,
				       const int *soffset,
				       int nwavelet,
				       const double *wavelet,
				       const int *woffset,
				       gsl_matrix *m,
				       generic_matrix_edge_func_t edge);

int
generic_matrix_create_inverse2d_col_step(int j_max,
					 int nscaling,
					 const double *scaling,
					 const int *soffset,
					 int nwavelet,
					 const double *wavelet,
					 const int *woffset,
					 gsl_matrix **_m,
					 generic_matrix_edge_func_t edge);

int
generic_matrix_create_inverse2d_col_substep(int j_max,
					    int j0,
					    int nscaling,
					    const double *scaling,
					    const int *soffset,
					    int nwavelet,
					    const double *wavelet,
					    const int *woffset,
					    gsl_matrix **_m,
					    generic_matrix_edge_func_t edge);

int
generic_matrix_fill_inverse2d_row_A(int j_max,
				    int nscaling,
				    const double *scaling,
				    const int *soffset,
				    gsl_matrix *m,
				    int row_offset,
				    int col_offset,
				    generic_matrix_edge_func_t edge);

int
generic_matrix_fill_inverse2d_row_B(int j_max,
				    int nwavelet,
				    const double *wavelet,
				    const int *woffset,
				    gsl_matrix *m,
				    int row_offset,
				    int col_offset,
				    generic_matrix_edge_func_t edge);

int
generic_matrix_fill_inverse2d_row_step(int j_max,
				       int nscaling,
				       const double *scaling,
				       const int *soffset,
				       int nwavelet,
				       const double *wavelet,
				       const int *woffset,
				       gsl_matrix *m,
				       generic_matrix_edge_func_t edge);

int
generic_matrix_create_inverse2d_row_step(int j_max,
					 int nscaling,
					 const double *scaling,
					 const int *soffset,
					 int nwavelet,
					 const double *wavelet,
					 const int *woffset,
					 gsl_matrix **_m,
					 generic_matrix_edge_func_t edge);

int
generic_matrix_create_inverse2d_row_substep(int j_max,
					    int j0,
					    int nscaling,
					    const double *scaling,
					    const int *soffset,
					    int nwavelet,
					    const double *wavelet,
					    const int *woffset,
					    gsl_matrix **_m,
					    generic_matrix_edge_func_t edge);

int
generic_matrix_fill_inverse2d_step(int jmax,
				   int nscaling,
				   const double *scaling,
				   const int *soffset,
				   int nwavelet,
				   const double *wavelet,
				   const int *woffset,
				   gsl_matrix *m,
				   generic_matrix_edge_func_t edge);

int
generic_matrix_create_inverse2d_step(int jmax,
				     int nscaling,
				     const double *scaling,
				     const int *soffset,
				     int nwavelet,
				     const double *wavelet,
				     const int *woffset,
				     gsl_matrix **m,
				     generic_matrix_edge_func_t edge);

int
generic_matrix_create_inverse2d_substep(int jmax,
					int j0,
					int nscaling,
					const double *scaling,
					const int *soffset,
					int nwavelet,
					const double *wavelet,
					const int *woffset,
					gsl_matrix **m,
					generic_matrix_edge_func_t edge);


int
generic_matrix_create_inverse2d(int j_max,
				int nscaling,
				const double *scaling,
				const int *soffset,
				int nwavelet,
				const double *wavelet,
				const int *woffset,
				gsl_matrix **m,
				generic_matrix_edge_func_t edge);

#endif /* generic_matrix_h */
