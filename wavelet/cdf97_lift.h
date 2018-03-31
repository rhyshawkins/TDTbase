#ifndef cdf97_lift_h
#define cdf97_lift_h

/*
 * 1D Forward
 */
int
cdf97_lift_forward1d_cdf97(double *s,
			   int width,
			   int stride,
			   double *work);

int
cdf97_lift_forward1d_cdf97_step(double *s,
				int width,
				int stride,
				double *work);

/*
 * 1D Inverse
 */
int
cdf97_lift_inverse1d_cdf97(double *s,
			   int width,
			   int stride,
			   double *work);

int
cdf97_lift_inverse1d_cdf97_step(double *s,
				int width,
				int stride,
				double *work);

/*
 * 2D Forward
 */
int 
cdf97_lift_forward2d_cdf97(double *s, 
			   int width,
			   int height,
			   int stride,
			   double *work,
			   int subtile);

int 
cdf97_lift_forward2d_cdf97_step(double *s, 
				int width,
				int height,
				int stride,
				double *work);

/*
 * 2D Inverse
 */
int
cdf97_lift_inverse2d_cdf97(double *s,
			   int width,
			   int height,
			   int stride,
			   double *work,
			   int subtile);

int 
cdf97_lift_inverse2d_cdf97_step(double *s, 
				int width,
				int height,
				int stride,
				double *work);

/*
 * 3D Forward
 */
int 
cdf97_lift_forward3d_cdf97(double *s,
			   int width,
			   int height,
			   int depth,
			   int rowstride,
			   int slicestride,
			   double *work,
			   int subtile);


int 
cdf97_lift_forward3d_cdf97_step(double *s,
				int width,
				int height,
				int depth,
				int rowstride,
				int slicestride,
				double *work);

int 
cdf97_lift_forward3d_cdf97_2dstep(double *s,
				  int width,
				  int height,
				  int stride,
				  int rowstride,
				  double *work);


/*
 * 3D Inverse
 */
int 
cdf97_lift_inverse3d_cdf97(double *s,
			   int width,
			   int height,
			   int depth,
			   int rowstride,
			   int slicestride,
			   double *work,
			   int subtile);

int 
cdf97_lift_inverse3d_cdf97_step(double *s,
				int width,
				int height,
				int depth,
				int rowstride,
				int slicestride,
				double *work);

int 
cdf97_lift_inverse3d_cdf97_2dstep(double *s,
				  int width,
				  int height,
				  int stride,
				  int rowstride,
				  double *work);

#endif /* cdf97_lift_h */
