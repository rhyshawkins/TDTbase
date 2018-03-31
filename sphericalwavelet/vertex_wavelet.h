#ifndef vertex_wavelet_h
#define vertex_wavelet_h

#include "manifold.h"

/*
 * Full transforms
 */
int
vertex_wavelet_butterfly_forward(manifold_t *m,
				 double *coeff);

int
vertex_wavelet_butterfly_inverse(manifold_t *m,
				 double *coeff);

int
vertex_wavelet_butterfly_forward_lifted(manifold_t *m,
					double *coeff);

int
vertex_wavelet_butterfly_inverse_lifted(manifold_t *m,
					double *coeff);

/*
 * Individual steps
 */
int
vertex_wavelet_butterfly_forward_step(manifold_t *m,
				      double *coeff,
				      int depth);

int
vertex_wavelet_butterfly_forward_lifted_step(manifold_t *m,
					     double *coeff,
					     int depth);

int
vertex_wavelet_butterfly_inverse_step(manifold_t *m,
				      double *coeff,
				      int depth);

int
vertex_wavelet_butterfly_inverse_lifted_step(manifold_t *m,
					     double *coeff,
					     int depth);

/*
 * Shell transforms : radial degree == manifold degree
 */
int
vertex_wavelet_butterfly_shell_forward(manifold_t *m,
				       double *coeff,
				       int ncoeff,
				       double *workspace);

int
vertex_wavelet_butterfly_shell_inverse(manifold_t *m,
				       double *coeff,
				       int ncoeff, 
				       double *workspace);

int
vertex_wavelet_butterfly_shell_forward_lifted(manifold_t *m,
					      double *coeff,
					      int ncoeff,
					      double *workspace);

int
vertex_wavelet_butterfly_shell_inverse_lifted(manifold_t *m,
					      double *coeff,
					      int ncoeff,
					      double *workspace);



#endif /* vertex_wavelet_h */
