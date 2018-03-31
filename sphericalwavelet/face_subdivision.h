#ifndef face_subdivision_h
#define face_subdivision_h

#include "manifold.h"

/*
 * This ia a non-wavelet subdivision scheme. For a given manifold we define the forward transform
 * as the mean of the 4 sub-triangles is the value of the parent triangle, and the 4 differences
 * from the mean are set to the 4 sub-triangles. This is similar to a Haar wavelet.
 */

/*
 * Full transforms
 */
int
face_subdivision_forward(manifold_t *m,
			 double *coeff);

int
face_subdivision_inverse(manifold_t *m,
			 double *coeff);

/*
 * Individual steps
 */
int
face_subdivision_forward_step(manifold_t *m,
			      double *coeff,
			      int depth);

int
face_subdivision_inverse_step(manifold_t *m,
			      double *coeff,
			      int depth);

#endif /* face_subdivision_h */
