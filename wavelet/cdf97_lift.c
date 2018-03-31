
#include <stdio.h>

#include "cdf97_lift.h"

#include "generic_lift.h"

static const double a1 = -1.58613434200;
static const double a2 = -0.05298011854;
static const double a3 =  0.88291107620;
static const double a4 =  0.44350685220;

static const double k1 = 0.81289306611596146; // 1.14960439886/sqrt(2)
static const double k2 = 0.61508705245700002; // (1/1.14960439886)/sqrt(2)

static const double ik1 = 1.230174104914;
static const double ik2 = 1.6257861322319229;

/*
 * 1D
 */

int
cdf97_lift_forward1d_cdf97(double *s,
			   int width,
			   int stride,
			   double *work)
{
  return generic_lift_forward1d(s, width, stride, work,
				cdf97_lift_forward1d_cdf97_step);
}

int
cdf97_lift_forward1d_cdf97_step(double *s,
				int width,
				int stride,
				double *work)
{
  int i;

  /*
   * Copy to workspace
   */
  for (i = 0; i < width; i ++) {
    work[i] = s[stride*i];
  }

  /*
   * Lifting steps
   */
  for (i = 1; i < (width - 1); i += 2) {
    work[i] += a1 * (work[i - 1] + work[i + 1]);
  }
  work[width - 1] += 2.0 * a1 * work[width - 2];
  
  for (i = 2; i < width; i += 2) {
    work[i] += a2 * (work[i - 1] + work[i + 1]);
  }
  work[0] += 2.0 * a2 * work[1];
  
  for (i = 1; i < (width - 1); i += 2) {
    work[i] += a3 * (work[i - 1] + work[i + 1]);
  }
  work[width - 1] += 2.0 * a3 * work[width - 2];
  
  for (i = 2; i < width; i += 2) {
    work[i] += a4 * (work[i - 1] + work[i + 1]);
  }
  work[0] += 2.0 * a4 * work[1];

  /*
   * Copy back and de-interleave
   */
  for (i = 0; i < width/2; i ++) {
    s[stride*i] = k1 * work[2*i];
    s[stride*(width/2 + i)] = k2 * work[2*i + 1];
  }

  return 0;
}

int
cdf97_lift_inverse1d_cdf97(double *s,
			   int width,
			   int stride,
			   double *work)
{
  return generic_lift_inverse1d(s, width, stride, work,
				cdf97_lift_inverse1d_cdf97_step);
}

int
cdf97_lift_inverse1d_cdf97_step(double *s,
				int width,
				int stride,
				double *work)
{
  int i;

  /*
   * Copy to workspace and interleave
   */
  for (i = 0; i < width/2; i ++) {
    work[2*i] = ik1 * s[i*stride];
    work[2*i + 1] = ik2 * s[(width/2 + i)*stride];
  }

  /*
   * Inverse lifting steps
   */
  for (i = 2; i < width; i += 2) {
    work[i] -= a4 * (work[i - 1] + work[i + 1]);
  }
  work[0] -= 2.0 * a4 * work[1];
  
  for (i = 1; i < (width - 1); i += 2) {
    work[i] -= a3 * (work[i - 1] + work[i + 1]);
  }
  work[width - 1] -= 2.0 * a3 * work[width - 2];
  
  for (i = 2; i < width; i += 2) {
    work[i] -= a2 * (work[i - 1] + work[i + 1]);
  }
  work[0] -= 2.0 * a2 * work[1];
  
  for (i = 1; i < (width - 1); i += 2) {
    work[i] -= a1 * (work[i - 1] + work[i + 1]);
  }
  work[width - 1] -= 2.0 * a1 * work[width - 2];
  
  /*
   * Copy back
   */
  for (i = 0; i < width; i ++) {
    s[i*stride] = work[i];
  }

  return 0;
}

/*
 * 2D
 */

int 
cdf97_lift_forward2d_cdf97(double *s, 
			   int width,
			   int height,
			   int stride,
			   double *work,
			   int subtile)
{
  return generic_lift_forward2d(s,
				width,
				height,
				stride,
				work,
				cdf97_lift_forward1d_cdf97_step,
				cdf97_lift_forward1d_cdf97_step,
				subtile);
}

int 
cdf97_lift_forward2d_cdf97_step(double *s,
				int width,
				int height,
				int stride,
				double *work)
{
  int i;

  /*
   * 1D Forward Transform on Rows
   */

  for (i = 0; i < height; i ++) {

    cdf97_lift_forward1d_cdf97_step(s + i*stride,
				    width,
				    1,
				    work);

  }

  /*
   * 1D Forward Transform on Columns
   */

  for (i = 0; i < width; i ++) {

    cdf97_lift_forward1d_cdf97_step(s + i,
				    height,
				    stride,
				    work);
  }
      
  return 0;
}



int
cdf97_lift_inverse2d_cdf97(double *s,
			   int width,
			   int height,
			   int stride,
			   double *work,
			   int subtile)
{
  return generic_lift_inverse2d(s,
				width,
				height,
				stride,
				work,
				cdf97_lift_inverse1d_cdf97_step,
				cdf97_lift_inverse1d_cdf97_step,
				subtile);
}

int 
cdf97_lift_inverse2d_cdf97_step(double *s,
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

    cdf97_lift_inverse1d_cdf97_step(s + i*stride,
				    width,
				    1,
				    work);

  }

  /*
   * 1D Inverse Transform on Columns
   */

  for (i = 0; i < width; i ++) {

    cdf97_lift_inverse1d_cdf97_step(s + i,
				    height,
				    stride,
				    work);
  }
      
  return 0;
}

int 
cdf97_lift_forward3d_cdf97(double *s,
			   int width,
			   int height,
			   int depth,
			   int rowstride,
			   int slicestride,
			   double *work,
			   int subtile)
{
  return generic_lift_forward3d(s,
				width,
				height,
				depth,
				rowstride,
				slicestride,
				work,
				cdf97_lift_forward1d_cdf97_step,
				cdf97_lift_forward1d_cdf97_step,
				cdf97_lift_forward1d_cdf97_step,
				subtile);
}

int 
cdf97_lift_forward3d_cdf97_step(double *s,
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

      cdf97_lift_forward1d_cdf97_step(s + o,
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

      cdf97_lift_forward1d_cdf97_step(s + o,
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

      cdf97_lift_forward1d_cdf97_step(s + o,
				      depth,
				      slicestride,
				      work);

    }
  }

  return 0;
}

int 
cdf97_lift_forward3d_cdf97_2dstep(double *s,
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
    cdf97_lift_forward1d_cdf97_step(s + i*rowstride,
				    width,
				    stride,
				    work);

  }

  /*
   * 1D Transform on Columns
   */
  for (i = 0; i < width; i ++) {
    cdf97_lift_forward1d_cdf97_step(s + i*stride,
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
cdf97_lift_inverse3d_cdf97(double *s,
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
      if (cdf97_lift_inverse1d_cdf97_step(s, w, 1, work) < 0) {
	return -1;
      }
      
      w <<= 1;
    }
    
    for (i = 0; i < hlevels; i ++) {
      if (cdf97_lift_inverse1d_cdf97_step(s, h, rowstride, work) < 0) {
	return -1;
      }
      
      h <<= 1;
    }
    
    for (i = 0; i < dlevels; i ++) {
      if (cdf97_lift_inverse1d_cdf97_step(s, d, slicestride, work) < 0) {
	return -1;
      }
      
      d <<= 1;
    }
    
    /*
     * 2D expansion for non-square
     */
    for (i = 0; i < whlevels; i ++) {
      if (cdf97_lift_inverse3d_cdf97_2dstep(s, w, h, 1, rowstride, work) < 0) {
	return -1;
      }
      
      w <<= 1;
      h <<= 1;
    }
    
    for (i = 0; i < wdlevels; i ++) {
      if (cdf97_lift_inverse3d_cdf97_2dstep(s, w, d, 1, slicestride, work) < 0) {
	return -1;
      }
      
      w <<= 1;
      d <<= 1;
    }
    
    for (i = 0; i < hdlevels; i ++) {
      if (cdf97_lift_inverse3d_cdf97_2dstep(s, h, d, rowstride, slicestride, work) < 0) {
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

    if (cdf97_lift_inverse3d_cdf97_step(s, w, h, d, rowstride, slicestride, work) < 0) {
      return -1;
    }

    w <<= 1;
    h <<= 1;
    d <<= 1;
  }

  return 0;
}

int 
cdf97_lift_inverse3d_cdf97_step(double *s,
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

      cdf97_lift_inverse1d_cdf97_step(s + o,
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

      cdf97_lift_inverse1d_cdf97_step(s + o,
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

      cdf97_lift_inverse1d_cdf97_step(s + o,
				      width,
				      1,
				      work);

    }
  }

  return 0;
}

int 
cdf97_lift_inverse3d_cdf97_2dstep(double *s,
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
    cdf97_lift_inverse1d_cdf97_step(s + i*stride,
				    height,
				    rowstride,
				    work);

  }

  /*
   * 1D Transform on Rows
   */
  for (i = 0; i < height; i ++) {
    cdf97_lift_inverse1d_cdf97_step(s + i*rowstride,
				    width,
				    stride,
				    work);

  }

  return 0;
}

