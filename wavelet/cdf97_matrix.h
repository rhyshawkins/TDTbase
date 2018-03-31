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
