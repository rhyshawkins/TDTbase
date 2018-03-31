
#include <stdio.h>

#include "generic_lift.h"

/*
 * 1D Full Transform
 */
int
generic_lift_forward1d(double *s,
		       int width,
		       int stride,
		       double *work,
		       generic_lift_forward1d_step_t transform)
{
  int w;

  w = width;

  while (w > 1) {

    if (transform(s, w, stride, work) < 0) {
      return -1;
    }

    w >>= 1;
  }

  return 0;
}

int
generic_lift_inverse1d(double *s,
		       int width,
		       int stride,
		       double *work,
		       generic_lift_inverse1d_step_t transform)
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

    if (transform(s, w, stride, work) < 0) {
      return -1;
    }

    w <<= 1;
  }

  return 0;
}

/*
 * 2D Step transform
 */

int
generic_lift_inverse2d_step(double *s,
			    int width,
			    int height,
			    int stride,
			    double *work,
			    generic_lift_forward1d_step_t row_transform,
			    generic_lift_forward1d_step_t col_transform)
{
  int j;
  
  /*
   * 1D Transform on Rows
   */
  for (j = 0; j < height; j ++) {
    if (row_transform(s + j*stride,
		      width,
		      1,
		      work) < 0) {
      return -1;
    }
  }

  /*
   * 1D Transform on Columns
   */
  for (j = 0; j < width; j ++) {
    if (col_transform(s + j,
		      height,
		      stride,
		      work) < 0) {
      return -1;
    }
  }
  
  return 0;
}

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
		       int subtile)
{
  int w;
  int h;
  int i;

  w = width;
  h = height;

  while (w > 1 && h > 1) {

    /*
     * 1D Transform on Columns
     */
    for (i = 0; i < w; i ++) {
      if (col_transform(s + i,
			h,
			stride,
			work) < 0) {
	return -1;
      }
    }

    /*
     * 1D Transform on Rows
     */
    for (i = 0; i < h; i ++) {
      if (row_transform(s + i*stride,
			w,
			1,
			work) < 0) {
	return -1;
      }
    }

    w >>= 1;
    h >>= 1;
  }

  if (!subtile) {
    while (w > 1) {
      /*
       * Left with a single row
       */
      if (row_transform(s, w, 1, work) < 0) {
	return -1;
      }
      
      w >>= 1;
    }
    
    while (h > 1) {
      /*
       * Left with a single column
       */
      if (col_transform(s, h, stride, work) < 0) {
	return -1;
      }
      
      h >>= 1;
    }
  }
  
  return 0;
}

int
generic_lift_inverse2d(double *s,
		       int width,
		       int height,
		       int stride,
		       double *work,
		       generic_lift_inverse1d_step_t row_transform,
		       generic_lift_inverse1d_step_t col_transform,
		       int subtile)
{
  int w;
  int h;
  int levels;

  int wlevels;
  int hlevels;
  int i;
  int j;

  w = width;
  h = height;
  levels = 0;
  wlevels = 0;
  hlevels = 0;

  while (w > 2 && h > 2) {
    levels ++;
    w >>= 1;
    h >>= 1;
  }

  if (!subtile) {
    while (w > 2) {
      wlevels ++;
      w >>= 1;
    }
    
    while (h > 2) {
      hlevels ++;
      h >>= 1;
    }

    for (i = 0; i < wlevels; i ++) {
      
      if (row_transform(s, w, 1, work) < 0) {
	return -1;
      }
      
      w <<= 1;
    }
    
    for (i = 0; i < hlevels; i ++) {
      
      if (col_transform(s, h, stride, work) < 0) {
	return -1;
      }
      
      h <<= 1;
    }
  }
  
  for (i = 0; i <= levels; i ++) {

    /*
     * 1D Transform on Rows
     */
    for (j = 0; j < h; j ++) {
      if (row_transform(s + j*stride,
			w,
			1,
			work) < 0) {
	return -1;
      }
    }

    /*
     * 1D Transform on Columns
     */
    for (j = 0; j < w; j ++) {
      if (col_transform(s + j,
			h,
			stride,
			work) < 0) {
	return -1;
      }
    }


    w <<= 1;
    h <<= 1;
  }

  return 0;
}

/*
 * 3D Full Transform
 */
int
generic_lift_forward3d(double *s,
		       int width,
		       int height,
		       int depth,
		       int rowstride,
		       int slicestride,
		       double *work,
		       generic_lift_forward1d_step_t row_transform,
		       generic_lift_forward1d_step_t col_transform,
		       generic_lift_forward1d_step_t dep_transform,
		       int subtile)
{
  int w;
  int h;
  int d;

  int i;
  int j;
  int o;
  
  w = width;
  h = height;
  d = depth;

  /*
   * Full 3D Transform steps
   */
  while (w > 1 && h > 1 && d > 1) {

    /*
     * 1D Transform on Rows
     */
    for (i = 0; i < h; i ++) {
      for (j = 0; j < d; j ++) {

	o = j*slicestride + i*rowstride;
	if (row_transform(s + o,
			  w,
			  1,
			  work) < 0) {
	  return -1;
	}
      }
    }

    /*
     * 1D Transform on Columns
     */
    for (i = 0; i < w; i ++) {
      for (j = 0; j < d; j ++) {

	o = j*slicestride + i;

	if (col_transform(s + o,
			  h,
			  rowstride,
			  work) < 0) {
	  return -1;
	}
      }
    }
			
    /*
     * 1D Transform on Slices
     */
    for (i = 0; i < w; i ++) {
      for (j = 0; j < h; j ++) {

	o = j*rowstride + i;

	if (dep_transform(s + o,
			  d,
			  slicestride,
			  work) < 0) {
	  return -1;
	}
      }
    }
	
    w >>= 1;
    h >>= 1;
    d >>= 1;
  }

  if (!subtile) {
    /*
     * Remaining 2D Transform steps
     */
    if (d == 1) {
      while (w > 1 && h > 1) {
	
	/*
	 * 1D Transform on Rows
	 */
	for (i = 0; i < h; i ++) {
	  o = i*rowstride;
	  if (row_transform(s + o,
			    w,
			    1,
			    work) < 0) {
	    return -1;
	  }
	}
	
	/*
	 * 1D Transform on Columns
	 */
	for (i = 0; i < w; i ++) {
	  o = i;
	  
	  if (col_transform(s + o,
			    h,
			    rowstride,
			    work) < 0) {
	    return -1;
	  }
	}
	
	w >>= 1;
	h >>= 1;
      }
      
    } else if (h == 1) {
      
      while (w > 1 && d > 1) {
	
	/*
	 * 1D Transform on Rows
	 */
	for (j = 0; j < d; j ++) {
	  
	  o = j*slicestride;
	  if (row_transform(s + o,
			    w,
			    1,
			    work) < 0) {
	    return -1;
	  }
	}
	
	
	/*
	 * 1D Transform on Slices
	 */
	for (i = 0; i < w; i ++) {
	  
	  o = i;
	  
	  if (dep_transform(s + o,
			    d,
			    slicestride,
			    work) < 0) {
	    return -1;
	  }
	}
	
	w >>= 1;
	d >>= 1;
	
      }
      
    } else if (w == 1) {
      
      while (h > 1 && d > 1) {
	
	
	/*
	 * 1D Transform on Columns
	 */
	for (j = 0; j < d; j ++) {
	  
	  o = j*slicestride;
	  
	  if (col_transform(s + o,
			    h,
			    rowstride,
			    work) < 0) {
	    return -1;
	  }
	}
	
	/*
	 * 1D Transform on Slices
	 */
	for (j = 0; j < h; j ++) {
	  
	  o = j*rowstride;
	  
	  if (dep_transform(s + o,
			    d,
			    slicestride,
			    work) < 0) {
	    return -1;
	  }
	}
	
	h >>= 1;
	d >>= 1;
	
      }
      
    }
    
    /*
     * Remaining 1D Transform steps
     */
    while (w > 1) {
      
      if (row_transform(s, w, 1, work) < 0) {
	return -1;
      }
      
      w >>= 1;
    }
    
    while (h > 1) {
      
      if (col_transform(s, h, rowstride, work) < 0) {
	return -1;
      }
      
      h >>= 1;
    }
    
    while (d > 1) {
      
      if (dep_transform(s, d, slicestride, work) < 0) {
	return -1;
      }
      
      d >>= 1;
    }
  }
  
  return 0;
}


int
generic_lift_inverse3d(double *s,
		       int width,
		       int height,
		       int depth,
		       int rowstride,
		       int slicestride,
		       double *work,
		       generic_lift_inverse1d_step_t row_transform,
		       generic_lift_inverse1d_step_t col_transform,
		       generic_lift_inverse1d_step_t dep_transform,
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
  int j;
  int k;

  int o;

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
      if (row_transform(s, w, 1, work) < 0) {
	return -1;
      }
      
      w <<= 1;
    }
    
    for (i = 0; i < hlevels; i ++) {
      if (col_transform(s, h, rowstride, work) < 0) {
	return -1;
      }
      
      h <<= 1;
    }
    
    for (i = 0; i < dlevels; i ++) {
      if (dep_transform(s, d, slicestride, work) < 0) {
	return -1;
      }
      
      d <<= 1;
    }
    
    /*
     * 2D expansion for non-square
     */
    for (i = 0; i < whlevels; i ++) {

      for (j = 0; j < w; j ++) {
	if (col_transform(s + j, h, rowstride, work) < 0) {
	  return -1;
	}
      }

      for (j = 0; j < h; j ++) {
	if (row_transform(s + j*rowstride, w, 1, work) < 0) {
	  return -1;
	}
      }

      w <<= 1;
      h <<= 1;
    }
    
    for (i = 0; i < wdlevels; i ++) {

      for (j = 0; j < w; j ++) {
	if (dep_transform(s + j, d, slicestride, work) < 0) {
	  return -1;
	}
      }

      for (j = 0; j < d; j ++) {
	if (row_transform(s + j*slicestride, w, 1, work) < 0) {
	  return -1;
	}
      }

      w <<= 1;
      d <<= 1;
    }
    
    for (i = 0; i < hdlevels; i ++) {

      for (j = 0; j < h; j ++) {
	if (dep_transform(s + j * rowstride, d, slicestride, work) < 0) {
	  return -1;
	}
      }

      for (j = 0; j < d; j ++) {
	if (col_transform(s + j*slicestride, h, rowstride, work) < 0) {
	  return -1;
	}
      }

      h <<= 1;
      d <<= 1;
    }
  }
  
  /*
   * 3D expansion
   */
  for (i = 0; i <= levels; i ++) {

    /*
     * 1D Transform on Slices
     */
    for (j = 0; j < w; j ++) {
      for (k = 0; k < h; k ++) {
	o = k * rowstride + j;

	if (dep_transform(s + o, d, slicestride, work) < 0) {
	  return -1;
	}
      }
    }
	
    /*
     * 1D Transform on Columns
     */
    for (j = 0; j < w; j ++) {
      for (k = 0; k < d; k ++) {
	o = k * slicestride + j;

	if (col_transform(s + o, h, rowstride, work) < 0) {
	  return -1;
	}
      }
    }

    /*
     * 1D Transform on Rows
     */
    for (j = 0; j < h; j ++) {
      for (k = 0; k < d; k ++) {
	o = k * slicestride + j * rowstride;

	if (col_transform(s + o, w, 1, work) < 0) {
	  return -1;
	}
      }
    }

    w <<= 1;
    h <<= 1;
    d <<= 1;
  }

  return 0;
}

