
#include "daub6_dwt.h"

#include "boundary.h"

/*
 * These coefficients are scaled by 1/sqrt(2) from the true orthongal terms to match
 * the normalization of the lifting version. This then requires a factor of 2 in the 
 * inverse transform.
 */
#define H0  0.235233603892081840 //  0.33267055295008261599851158914
#define H1  0.570558457915721800 //  0.80689150931109257649449360409
#define H2  0.325182500263116240 //  0.45987750211849157009515194215
#define H3 -0.095467207784163680 // -0.13501102001025458869638990670
#define H4 -0.060416104155198096 // -0.08544127388202666169281916918
#define H5  0.024908749868441864 //  0.03522629188570953660274066472

#define G0 H5
#define G1 -H4
#define G2 H3
#define G3 -H2
#define G4 H1
#define G5 -H0

#define HB0 (2.0*H0)
#define HB1 (2.0*H1)
#define HB2 (2.0*H2)
#define HB3 (2.0*H3)
#define HB4 (2.0*H4)
#define HB5 (2.0*H5)

#define GB0 (2.0*G0)
#define GB1 (2.0*G1)
#define GB2 (2.0*G2)
#define GB3 (2.0*G3)
#define GB4 (2.0*G4)
#define GB5 (2.0*G5)

/*
 * 1D Forward
 */
int
daub6_dwt_forward1d_daub6(double *s,
			  int width,
			  int stride,
			  double *work)
{
  int w;

  w = width;

  while (w > 1) {

    if (daub6_dwt_forward1d_daub6_step(s, w, stride, work) < 0) {
      return -1;
    }

    w >>= 1;
  }

  return 0;
}

int
daub6_dwt_forward1d_daub6_step(double *s,
			       int width,
			       int stride,
			       double *work)
{
  int i;
  
  for (i = 0; i < width; i ++) {
    work[i] = s[i*stride];
  }

  for (i = 0; i < width/2; i ++) {
    s[i*stride] =
      H0 * work[2*i] + 
      H1 * work[2*i + 1] + 
      H2 * work[wavelet_boundary_periodic(2*i + 2, width)] + 
      H3 * work[wavelet_boundary_periodic(2*i + 3, width)] +
      H4 * work[wavelet_boundary_periodic(2*i + 4, width)] + 
      H5 * work[wavelet_boundary_periodic(2*i + 5, width)];
				
    s[(width/2 + i)*stride] =
      G0 * work[2*i] + 
      G1 * work[2*i + 1] + 
      G2 * work[wavelet_boundary_periodic(2*i + 2, width)] + 
      G3 * work[wavelet_boundary_periodic(2*i + 3, width)] +
      G4 * work[wavelet_boundary_periodic(2*i + 4, width)] + 
      G5 * work[wavelet_boundary_periodic(2*i + 5, width)];
  }

  return 0;
}

/*
 * 1D Inverse
 */
int
daub6_dwt_inverse1d_daub6(double *s,
			  int width,
			  int stride,
			  double *work)
{
  int w;
  int levels;
  int i;

  w = width;
  levels = 0;

  while (w > 2) {
    levels ++;
    w >>= 1;
  }

  for (i = 0; i <= levels; i ++) {

    if (daub6_dwt_inverse1d_daub6_step(s, w, stride, work) < 0) {
      return -1;
    }

    w <<= 1;
  }

  return 0;
}

int
daub6_dwt_inverse1d_daub6_step(double *s,
			       int width,
			       int stride,
			       double *work)
{
  int i;

  for (i = 0; i < width/2; i ++) {
    work[2*i] = s[i*stride];
    work[2*i + 1] = s[(i + width/2)*stride];
  }
  
  for (i = 0; i < width/2; i ++) {
    s[(2*i) * stride] =
      HB4 * work[wavelet_boundary_periodic(2*i - 4, width)] + 
      HB1 * work[wavelet_boundary_periodic(2*i - 3, width)] + 
      HB2 * work[wavelet_boundary_periodic(2*i - 2, width)] + 
      HB3 * work[wavelet_boundary_periodic(2*i - 1, width)] + 
      HB0 * work[2*i + 0] + 
      HB5 * work[2*i + 1];
				
    s[(2*i + 1)*stride] =
      GB0 * work[wavelet_boundary_periodic(2*i - 4, width)] + 
      GB5 * work[wavelet_boundary_periodic(2*i - 3, width)] + 
      GB2 * work[wavelet_boundary_periodic(2*i - 2, width)] + 
      GB3 * work[wavelet_boundary_periodic(2*i - 1, width)] + 
      GB4 * work[2*i + 0] + 
      GB1 * work[2*i + 1];
  }

  return 0;
}

/*
 * 2D Forward
 */
int 
daub6_dwt_forward2d_daub6(double *s, 
			  int width,
			  int height,
			  int stride,
			  double *work,
			  int subtile)
{
  int w;
  int h;

  w = width;
  h = height;

  while (w > 1 && h > 1) {

    if (daub6_dwt_forward2d_daub6_step(s, w, h, stride, work) < 0) {
      return -1;
    }

    w >>= 1;
    h >>= 1;
  }

  if (!subtile) {

    while (w > 1) {
      /*
       * Left with a single row
       */
      if (daub6_dwt_forward1d_daub6_step(s, w, 1, work) < 0) {
	return -1;
      }
      
      w >>= 1;
    }
    
    while (h > 1) {
      /*
       * Left with a single column
       */
      if (daub6_dwt_forward1d_daub6_step(s, h, stride, work) < 0) {
	return -1;
      }
      
      h >>= 1;
    }
  }

  return 0;
}

int 
daub6_dwt_forward2d_daub6_step(double *s, 
			       int width,
			       int height,
			       int stride,
			       double *work)
{
  int i;

  /*
   * 1D Transform on Columns
   */
  for (i = 0; i < width; i ++) {
    
    daub6_dwt_forward1d_daub6_step(s + i,
				  height,
				  stride,
				  work);
  }

  /*
   * 1D Transform on Rows
   */
  for (i = 0; i < height; i ++) {

    daub6_dwt_forward1d_daub6_step(s + stride*i,
				  width,
				  1,
				  work);
  }
  
  return 0;
}

/*
 * 2D Inverse
 */
int
daub6_dwt_inverse2d_daub6(double *s,
			  int width,
			  int height,
			  int stride,
			  double *work,
			  int subtile)
{
  int w;
  int h;
  int levels;

  int wlevels;
  int hlevels;
  int i;

  w = width;
  h = height;
  levels = 0;
  wlevels = 0;
  hlevels = 0;

  if (!subtile) {
    while (w > 2 && h > 2) {
      levels ++;
      w >>= 1;
      h >>= 1;
    }
    
    while (w > 2) {
      wlevels ++;
      w >>= 1;
    }
    
    while (h > 2) {
      hlevels ++;
      h >>= 1;
    }
  }
  
  for (i = 0; i < wlevels; i ++) {
    
    if (daub6_dwt_inverse1d_daub6_step(s, w, 1, work) < 0) {
      return -1;
    }

    w <<= 1;
  }

  for (i = 0; i < hlevels; i ++) {

    if (daub6_dwt_inverse1d_daub6_step(s, h, stride, work) < 0) {
      return -1;
    }

    h <<= 1;
  }

  for (i = 0; i <= levels; i ++) {

    if (daub6_dwt_inverse2d_daub6_step(s, w, h, stride, work) < 0) {
      return -1;
    }

    w <<= 1;
    h <<= 1;
  }

  return 0;
}

int 
daub6_dwt_inverse2d_daub6_step(double *s, 
			      int width,
			      int height,
			      int stride,
			      double *work)
{
  int i;

  /*
   * 1D Inverse Transform on Rows
   */

  for (i = 0; i < height; i ++) {

    daub6_dwt_inverse1d_daub6_step(s + i*stride,
				  width,
				  1,
				  work);
    
  }

  /*
   * 1D Inverse Transform on Columns
   */

  for (i = 0; i < width; i ++) {

    daub6_dwt_inverse1d_daub6_step(s + i,
				  height,
				  stride,
				  work);
  }
      
  return 0;
}

/*
 * 3D Forward
 */
int 
daub6_dwt_forward3d_daub6(double *s,
			  int width,
			  int height,
			  int depth,
			  int rowstride,
			  int slicestride,
			  double *work,
			  int subtile)
{
  int w;
  int h;
  int d;

  w = width;
  h = height;
  d = depth;

  while (w > 1 && h > 1 && d > 1) {

    if (daub6_dwt_forward3d_daub6_step(s, w, h, d, rowstride, slicestride, work) < 0) {
      return -1;
    }

    w >>= 1;
    h >>= 1;
    d >>= 1;
  }

  if (!subtile) {
    if (d == 1) {
      while (w > 1 && h > 1) {

	if (daub6_dwt_forward3d_daub6_2dstep(s, w, h, 1, rowstride, work) < 0) {
	  return -1;
	}

	w >>= 1;
	h >>= 1;
      }

      while (w > 1) {
      
	if (daub6_dwt_forward1d_daub6_step(s, w, 1, work) < 0) {
	  return -1;
	}
      
	w >>= 1;
      }

      while (h > 1) {

	if (daub6_dwt_forward1d_daub6_step(s, h, rowstride, work) < 0) {
	  return -1;
	}

	h >>= 1;
      }

    } else if (h == 1) {

      while (w > 1 && d > 1) {

	if (daub6_dwt_forward3d_daub6_2dstep(s, w, d, 1, slicestride, work) < 0) {
	  return -1;
	}

	w >>= 1;
	d >>= 1;

      }

      while (w > 1) {

	if (daub6_dwt_forward1d_daub6_step(s, w, 1, work) < 0) {
	  return -1;
	}

	w >>= 1;
      }

      while (d > 1) {

	if (daub6_dwt_forward1d_daub6_step(s, d, slicestride, work) < 0) {
	  return -1;
	}

	d >>= 1;
      }

    } else if (w == 1) {

      while (h > 1 && d > 1) {

	if (daub6_dwt_forward3d_daub6_2dstep(s, h, d, rowstride, slicestride, work) < 0) {
	  return -1;
	}

	h >>= 1;
	d >>= 1;

      }

      while (h > 1) {

	if (daub6_dwt_forward1d_daub6_step(s, h, rowstride, work) < 0) {
	  return -1;
	}

	h >>= 1;
      }

      while (d > 1) {

	if (daub6_dwt_forward1d_daub6_step(s, d, slicestride, work) < 0) {
	  return -1;
	}

	d >>= 1;
      }

    }
  }
  
  return 0;
}
   

int 
daub6_dwt_forward3d_daub6_step(double *s,
			      int width,
			      int height,
			      int depth,
			      int rowstride,
			      int slicestride,
			      double *work)
{
  int i;
  int j;
  int o;

  /*
   * 1D Transform on Rows
   */
  for (i = 0; i < height; i ++) {
    for (j = 0; j < depth; j ++) {

      o = j*slicestride + i*rowstride;

      daub6_dwt_forward1d_daub6_step(s + o,
				    width,
				    1,
				    work);

    }
  }

  /*
   * 1D Transform on Columns
   */
  for (i = 0; i < width; i ++) {
    for (j = 0; j < depth; j ++) {
    
      o = j*slicestride + i;

      daub6_dwt_forward1d_daub6_step(s + o,
				    height,
				    rowstride,
				    work);
      
    }
  }

  /*
   * 1D Transforms on Slices
   */
  for (i = 0; i < width; i ++) {
    for (j = 0; j < height; j ++) {

      o = j*rowstride + i;

      daub6_dwt_forward1d_daub6_step(s + o,
				    depth,
				    slicestride,
				    work);

    }
  }

  return 0;
}

int 
daub6_dwt_forward3d_daub6_2dstep(double *s,
				int width,
				int height,
				int stride,
				int rowstride,
				double *work)
{
  int i;

  /*
   * 1D Transform on Rows
   */
  for (i = 0; i < height; i ++) {
    daub6_dwt_forward1d_daub6_step(s + i*rowstride,
				  width,
				  stride,
				  work);

  }

  /*
   * 1D Transform on Columns
   */
  for (i = 0; i < width; i ++) {
    daub6_dwt_forward1d_daub6_step(s + i*stride,
				  height,
				  rowstride,
				  work);

  }

  return 0;
}


/*
 * 3D Inverse
 */
int 
daub6_dwt_inverse3d_daub6(double *s,
			  int width,
			  int height,
			  int depth,
			  int rowstride,
			  int slicestride,
			  double *work,
			  int subtile)
{
  int w;
  int h;
  int d;
  int levels;

  int whlevels;
  int wdlevels;
  int hdlevels;

  int wlevels;
  int hlevels;
  int dlevels;

  int i;

  w = width;
  h = height;
  d = depth;
  levels = 0;

  whlevels = 0;
  wdlevels = 0;
  hdlevels = 0;

  wlevels = 0;
  hlevels = 0;
  dlevels = 0;

  while (w > 2 && h > 2 && d > 2) {
    levels ++;
    w >>= 1;
    h >>= 1;
    d >>= 1;
  }

  if (!subtile) {
    while (w > 2 && h > 2) {
      whlevels ++;
      w >>= 1;
      h >>= 1;
    }

    while (w > 2 && d > 2) {
      wdlevels ++;
      w >>= 1;
      d >>= 1;
    }

    while (h > 2 && d > 2) {
      hdlevels ++;
      h >>= 1;
      d >>= 1;
    }
  
    while (w > 2) {
      wlevels ++;
      w >>= 1;
    }

    while (h > 2) {
      hlevels ++;
      h >>= 1;
    }

    while (d > 2) {
      dlevels ++;
      d >>= 1;
    }

    /*
     * 1D expansion for non-square
     */
    for (i = 0; i < wlevels; i ++) {
      if (daub6_dwt_inverse1d_daub6_step(s, w, 1, work) < 0) {
	return -1;
      }

      w <<= 1;
    }

    for (i = 0; i < hlevels; i ++) {
      if (daub6_dwt_inverse1d_daub6_step(s, h, rowstride, work) < 0) {
	return -1;
      }

      h <<= 1;
    }

    for (i = 0; i < dlevels; i ++) {
      if (daub6_dwt_inverse1d_daub6_step(s, d, slicestride, work) < 0) {
	return -1;
      }

      d <<= 1;
    }

    /*
     * 2D expansion for non-square
     */
    for (i = 0; i < whlevels; i ++) {
      if (daub6_dwt_inverse3d_daub6_2dstep(s, w, h, 1, rowstride, work) < 0) {
	return -1;
      }

      w <<= 1;
      h <<= 1;
    }

    for (i = 0; i < wdlevels; i ++) {
      if (daub6_dwt_inverse3d_daub6_2dstep(s, w, d, 1, slicestride, work) < 0) {
	return -1;
      }

      w <<= 1;
      d <<= 1;
    }
  
    for (i = 0; i < hdlevels; i ++) {
      if (daub6_dwt_inverse3d_daub6_2dstep(s, h, d, rowstride, slicestride, work) < 0) {
	return -1;
      }

      h <<= 1;
      d <<= 1;
    }
  }
  
  /*
   * 3D expansion
   */
  for (i = 0; i <= levels; i ++) {

    if (daub6_dwt_inverse3d_daub6_step(s, w, h, d, rowstride, slicestride, work) < 0) {
      return -1;
    }

    w <<= 1;
    h <<= 1;
    d <<= 1;
  }

  return 0;
}

int 
daub6_dwt_inverse3d_daub6_step(double *s,
			      int width,
			      int height,
			      int depth,
			      int rowstride,
			      int slicestride,
			      double *work)
{
  int i;
  int j;
  int o;

  /*
   * 1D Transforms on Slices
   */
  for (i = 0; i < width; i ++) {
    for (j = 0; j < height; j ++) {

      o = j*rowstride + i;

      daub6_dwt_inverse1d_daub6_step(s + o,
				    depth,
				    slicestride,
				    work);

    }
  }

  /*
   * 1D Inverse Transform on Columns
   */
  for (i = 0; i < width; i ++) {
    for (j = 0; j < depth; j ++) {
    
      o = j*slicestride + i;

      daub6_dwt_inverse1d_daub6_step(s + o,
				    height,
				    rowstride,
				    work);
    }
  }

  /*
   * 1D Transform on Rows
   */
  for (i = 0; i < height; i ++) {
    for (j = 0; j < depth; j ++) {

      o = j*slicestride + i*rowstride;

      daub6_dwt_inverse1d_daub6_step(s + o,
				    width,
				    1,
				    work);

    }
  }

  return 0;
}

int 
daub6_dwt_inverse3d_daub6_2dstep(double *s,
				int width,
				int height,
				int stride,
				int rowstride,
				double *work)
{
  int i;

  /*
   * 1D Transform on Columns
   */
  for (i = 0; i < width; i ++) {
    daub6_dwt_inverse1d_daub6_step(s + i*stride,
				  height,
				  rowstride,
				  work);
    
  }

  /*
   * 1D Transform on Rows
   */
  for (i = 0; i < height; i ++) {
    daub6_dwt_inverse1d_daub6_step(s + i*rowstride,
				  width,
				  stride,
				  work);

  }

  return 0;
}

