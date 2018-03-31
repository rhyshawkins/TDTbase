#ifndef daubechies_h
#define daubechies_h

int
daubechies2d_forward_d4(double *s,
		      int width,
		      int height,
		      int stride,
		      double *work);

int
daubechies2d_inverse_d4(double *s,
		      int width,
		      int height,
		      int stride,
		      double *work);

int
daubechies1d_forward_d4_step(double *s,
			     int width,
			     int step,
			     double *work);

int
daubechies1d_inverse_d4_step(double *s,
			     int width,
			     int step,
			     double *work);


#endif /* daubechies_h */
