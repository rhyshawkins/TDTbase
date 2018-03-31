#ifndef cdf97_lift_impulse_h
#define cdf97_lift_impulse_h

/*
 * Size computations
 */
int cdf97_lift_impulse_1dsize(int degree, int impulse_degree);

/*
 * 1D Impulse Functions
 */
int cdf97_lift_impulse_1d(int degree, 
			  int impulse_degree, 
			  double impulse, 
			  double *v, 
			  int size,
			  int offset);


#endif /* cdf97_lift_impulse_h */
