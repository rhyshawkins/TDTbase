#ifndef face_wavelet_h
#define face_wavelet_h

#include "manifold.h"

/*
 * Full transforms
 */
int
face_wavelet_biohaar_forward(manifold_t *m,
			     double *coeff);

int
face_wavelet_biohaar_inverse(manifold_t *m,
			     double *coeff);

/*
 * Individual steps
 */
int
face_wavelet_biohaar_forward_step(manifold_t *m,
			      double *coeff,
			      int depth);

int
face_wavelet_biohaar_inverse_step(manifold_t *m,
				  double *coeff,
				  int depth);

/*
 * Shell transforms : radial degree == manifold degree
 */

typedef int (*face_wavelet_radial_forward_step_t)(double *s,
						  int width,
						  int stride,
						  double *work);
typedef int (*face_wavelet_radial_inverse_step_t)(double *s,
						  int width,
						  int stride,
						  double *work);

  
int
face_wavelet_biohaar_shell_forward(manifold_t *m,
				   double *coeff,
				   int ncoeff,
				   double *workspace,
				   face_wavelet_radial_forward_step_t radial_forward_step);

int
face_wavelet_biohaar_shell_inverse(manifold_t *m,
				   double *coeff,
				   int ncoeff, 
				   double *workspace,
				   face_wavelet_radial_inverse_step_t radial_inverse_step);


#endif /* face_wavelet_h */
