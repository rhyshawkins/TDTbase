#ifndef generic_lift_h
#define generic_lift_h


typedef int (*generic_lift_forward1d_step_t)(double *s,
					     int width,
					     int stride,
					     double *work);

typedef int (*generic_lift_inverse1d_step_t)(double *s,
					     int width,
					     int stride,
					     double *work);

/*
 * 1D Full Transform
 */
int
generic_lift_forward1d(double *s,
		       int width,
		       int stride,
		       double *work,
		       generic_lift_forward1d_step_t transform);

int
generic_lift_inverse1d(double *s,
		       int width,
		       int stride,
		       double *work,
		       generic_lift_inverse1d_step_t transform);

/*
 * 2D Step
 */
int
generic_lift_inverse2d_step(double *s,
			    int width,
			    int height,
			    int stride,
			    double *work,
			    generic_lift_forward1d_step_t row_transform,
			    generic_lift_forward1d_step_t col_transform);

/*
 * 2D Full Transform
 */
int
generic_lift_forward2d(double *s,
		       int width,
		       int height,
		       int stride,
		       double *work,
		       generic_lift_forward1d_step_t row_transform,
		       generic_lift_forward1d_step_t col_transform,
		       int subtile);

int
generic_lift_inverse2d(double *s,
		       int width,
		       int height,
		       int stride,
		       double *work,
		       generic_lift_inverse1d_step_t row_transform,
		       generic_lift_inverse1d_step_t col_transform,
		       int subtile);

/*
 * 3D Full Transform
 */
int
generic_lift_forward3d(double *s,
		       int width,
		       int height,
		       int depth,
		       int stride,
		       int slicestride,
		       double *work,
		       generic_lift_forward1d_step_t row_transform,
		       generic_lift_forward1d_step_t col_transform,
		       generic_lift_forward1d_step_t dep_transform,
		       int subtile);


int
generic_lift_inverse3d(double *s,
		       int width,
		       int height,
		       int depth,
		       int stride,
		       int slicestride,
		       double *work,
		       generic_lift_inverse1d_step_t row_transform,
		       generic_lift_inverse1d_step_t col_transform,
		       generic_lift_inverse1d_step_t dep_transform,
		       int subtile);


#endif /* generic_lift_h */
