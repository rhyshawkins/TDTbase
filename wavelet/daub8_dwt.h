#ifndef daub8_h
#define daub8_h

/*
 * 1D Forward
 */
int
daub8_dwt_forward1d_daub8(double *s,
			  int width,
			  int stride,
			  double *work);

int
daub8_dwt_forward1d_daub8_step(double *s,
			      int width,
			      int stride,
			      double *work);

/*
 * 1D Inverse
 */
int
daub8_dwt_inverse1d_daub8(double *s,
			 int width,
			 int stride,
			 double *work);

int
daub8_dwt_inverse1d_daub8_step(double *s,
			      int width,
			      int stride,
			      double *work);

/*
 * 2D Forward
 */
int 
daub8_dwt_forward2d_daub8(double *s, 
			  int width,
			  int height,
			  int stride,
			  double *work,
			  int subtile);

int 
daub8_dwt_forward2d_daub8_step(double *s, 
			      int width,
			      int height,
			      int stride,
			      double *work);

/*
 * 2D Inverse
 */
int
daub8_dwt_inverse2d_daub8(double *s,
			  int width,
			  int height,
			  int stride,
			  double *work,
			  int subtile);

int 
daub8_dwt_inverse2d_daub8_step(double *s, 
			      int width,
			      int height,
			      int stride,
			      double *work);

/*
 * 3D Forward
 */
int 
daub8_dwt_forward3d_daub8(double *s,
			  int width,
			  int height,
			  int depth,
			  int rowstride,
			  int slicestride,
			  double *work,
			  int subtile);

int 
daub8_dwt_forward3d_daub8_step(double *s,
			       int width,
			       int height,
			       int depth,
			       int rowstride,
			       int slicestride,
			       double *work);

int 
daub8_dwt_forward3d_daub8_2dstep(double *s,
				 int width,
				 int height,
				 int stride,
				 int rowstride,
				 double *work);


/*
 * 3D Inverse
 */
int 
daub8_dwt_inverse3d_daub8(double *s,
			  int width,
			  int height,
			  int depth,
			  int rowstride,
			  int slicestride,
			  double *work,
			  int subtile);

int 
daub8_dwt_inverse3d_daub8_step(double *s,
			      int width,
			      int height,
			      int depth,
			      int rowstride,
			      int slicestride,
			      double *work);

int 
daub8_dwt_inverse3d_daub8_2dstep(double *s,
				int width,
				int height,
				int stride,
				int rowstride,
				double *work);

#endif /* daub8_h */
