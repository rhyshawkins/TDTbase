#ifndef daub4_matrix_h
#define daub4_matrix_h

#include <gsl/gsl_matrix.h>

int
daub4_matrix_forward2d_create_col_step(int j_max, int j, gsl_matrix **_m);

int
daub4_matrix_inverse2d_create_col_step(int j_max, int j, gsl_matrix **_m);

int
daub4_matrix_forward2d_create_row_step(int j_max, int j, gsl_matrix **_m);

int
daub4_matrix_inverse2d_create_row_step(int j_max, int j, gsl_matrix **_m);

int
daub4_matrix_forward2d_create_step(int j_max, int j, gsl_matrix **_m);

int
daub4_matrix_inverse2d_create_step(int j_max, int j, gsl_matrix **_m);

int
daub4_matrix_forward2d_create(int j_max, gsl_matrix **m);

int
daub4_matrix_inverse2d_create(int j_max, gsl_matrix **_m);


#endif /* daub4_matrix_h */
