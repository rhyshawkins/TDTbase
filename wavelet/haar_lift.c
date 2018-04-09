//
//    Wavelet transform library
//    
//    Copyright (C) 2014 - 2018 Rhys Hawkins
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#include "haar_lift.h"

/*
 * 1D Forward
 */
int
haar_lift_forward1d_haar(double *s,
			 int width,
			 int stride,
			 double *work)
{
  int w;

  w = width;

  while (w > 1) {

    if (haar_lift_forward1d_haar_step(s, w, stride, work) < 0) {
      return -1;
    }

    w >>= 1;
  }

  return 0;
}

int
haar_lift_forward1d_haar_step(double *s,
			      int width,
			      int stride,
			      double *work)
{
  int i;

  /*
   * Copy to workspace 
   */
  for (i = 0; i < width; i ++) {
    work[i] = s[stride * i];
  }

  /*
   * Lifting steps
   */
  for (i = 1; i < width; i += 2) {
    work[i] -= work[i - 1];
  }

  for (i = 0; i < width; i += 2) {
    work[i] += 0.5 * work[i + 1];
  }

  /*
   * Copy back and de-interleave
   */
  for (i = 0; i < width/2; i ++) {
    s[stride*i] = work[2*i];
    s[stride*(width/2 + i)] = -0.5*work[2*i + 1];
  }

  return 0;
}

/*
 * 1D Inverse
 */
int
haar_lift_inverse1d_haar(double *s,
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

    if (haar_lift_inverse1d_haar_step(s, w, stride, work) < 0) {
      return -1;
    }

    w <<= 1;
  }

  return 0;
}
  

int
haar_lift_inverse1d_haar_step(double *s,
			      int width,
			      int stride,
			      double *work)
{
  int i;

  /*
   * Copy to workspace and interleave
   */
  for (i = 0; i < width/2; i ++) {
    work[2*i] = s[i*stride];
    work[2*i + 1] = -2.0 * s[(width/2 + i)*stride];
  }

  /*
   * Lifting steps
   */
  for (i = 0; i < width; i += 2) {
    work[i] -= 0.5 * work[i + 1];
  }
  

  for (i = 1; i < width; i += 2) {
    work[i] += work[i - 1];
  }

  /*
   * Copy back
   */
  for (i = 0; i < width; i ++) {
    s[i*stride] = work[i];
  }

  return 0;
}

/*
 * 2D Forward
 */
int 
haar_lift_forward2d_haar(double *s, 
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

    if (haar_lift_forward2d_haar_step(s, w, h, stride, work) < 0) {
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
      if (haar_lift_forward1d_haar_step(s, w, 1, work) < 0) {
	return -1;
      }
      
      w >>= 1;
    }
    
    while (h > 1) {
      /*
       * Left with a single column
       */
      if (haar_lift_forward1d_haar_step(s, h, stride, work) < 0) {
	return -1;
      }
      
      h >>= 1;
    }
  }

  return 0;
}


int 
haar_lift_forward2d_haar_step(double *s, 
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
    
    haar_lift_forward1d_haar_step(s + i,
				  height,
				  stride,
				  work);
  }

  /*
   * 1D Transform on Rows
   */
  for (i = 0; i < height; i ++) {

    haar_lift_forward1d_haar_step(s + stride*i,
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
haar_lift_inverse2d_haar(double *s,
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
      
      if (haar_lift_inverse1d_haar_step(s, w, 1, work) < 0) {
	return -1;
      }
      
      w <<= 1;
    }
    
    for (i = 0; i < hlevels; i ++) {
      
      if (haar_lift_inverse1d_haar_step(s, h, stride, work) < 0) {
	return -1;
      }

      h <<= 1;
    }
  }

  for (i = 0; i <= levels; i ++) {

    if (haar_lift_inverse2d_haar_step(s, w, h, stride, work) < 0) {
      return -1;
    }

    w <<= 1;
    h <<= 1;
  }

  return 0;
}

int 
haar_lift_inverse2d_haar_step(double *s, 
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

    haar_lift_inverse1d_haar_step(s + i*stride,
				  width,
				  1,
				  work);
    
  }

  /*
   * 1D Inverse Transform on Columns
   */

  for (i = 0; i < width; i ++) {

    haar_lift_inverse1d_haar_step(s + i,
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
haar_lift_forward3d_haar(double *s,
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

    if (haar_lift_forward3d_haar_step(s, w, h, d, rowstride, slicestride, work) < 0) {
      return -1;
    }

    w >>= 1;
    h >>= 1;
    d >>= 1;
  }

  if (!subtile) {

    if (d == 1) {
      while (w > 1 && h > 1) {
	
	if (haar_lift_forward3d_haar_2dstep(s, w, h, 1, rowstride, work) < 0) {
	  return -1;
	}
	
	w >>= 1;
	h >>= 1;
      }
      
      while (w > 1) {
	
	if (haar_lift_forward1d_haar_step(s, w, 1, work) < 0) {
	  return -1;
	}
	
	w >>= 1;
      }
      
      while (h > 1) {
	
	if (haar_lift_forward1d_haar_step(s, h, rowstride, work) < 0) {
	  return -1;
	}
	
	h >>= 1;
      }
      
    } else if (h == 1) {
      
      while (w > 1 && d > 1) {
	
	if (haar_lift_forward3d_haar_2dstep(s, w, d, 1, slicestride, work) < 0) {
	  return -1;
	}
	
	w >>= 1;
	d >>= 1;
	
      }

      while (w > 1) {
	
	if (haar_lift_forward1d_haar_step(s, w, 1, work) < 0) {
	  return -1;
	}
	
	w >>= 1;
      }
      
      while (d > 1) {
	
	if (haar_lift_forward1d_haar_step(s, d, slicestride, work) < 0) {
	  return -1;
	}
	
	d >>= 1;
      }
      
    } else if (w == 1) {
      
      while (h > 1 && d > 1) {
	
	if (haar_lift_forward3d_haar_2dstep(s, h, d, rowstride, slicestride, work) < 0) {
	  return -1;
	}
	
	h >>= 1;
	d >>= 1;
	
      }
      
      while (h > 1) {
	
	if (haar_lift_forward1d_haar_step(s, h, rowstride, work) < 0) {
	  return -1;
	}
	
	h >>= 1;
      }
      
      while (d > 1) {
	
	if (haar_lift_forward1d_haar_step(s, d, slicestride, work) < 0) {
	  return -1;
	}
	
	d >>= 1;
      }
    }
  }

  return 0;
}

int 
haar_lift_forward3d_haar_step(double *s,
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

      haar_lift_forward1d_haar_step(s + o,
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

      haar_lift_forward1d_haar_step(s + o,
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

      haar_lift_forward1d_haar_step(s + o,
				    depth,
				    slicestride,
				    work);

    }
  }

  return 0;
}

int 
haar_lift_forward3d_haar_2dstep(double *s,
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
    haar_lift_forward1d_haar_step(s + i*rowstride,
				  width,
				  stride,
				  work);

  }

  /*
   * 1D Transform on Columns
   */
  for (i = 0; i < width; i ++) {
    haar_lift_forward1d_haar_step(s + i*stride,
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
haar_lift_inverse3d_haar(double *s,
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
      if (haar_lift_inverse1d_haar_step(s, w, 1, work) < 0) {
	return -1;
      }

      w <<= 1;
    }

    for (i = 0; i < hlevels; i ++) {
      if (haar_lift_inverse1d_haar_step(s, h, rowstride, work) < 0) {
	return -1;
      }

      h <<= 1;
    }

    for (i = 0; i < dlevels; i ++) {
      if (haar_lift_inverse1d_haar_step(s, d, slicestride, work) < 0) {
	return -1;
      }

      d <<= 1;
    }

    /*
     * 2D expansion for non-square
     */
    for (i = 0; i < whlevels; i ++) {
      if (haar_lift_inverse3d_haar_2dstep(s, w, h, 1, rowstride, work) < 0) {
	return -1;
      }

      w <<= 1;
      h <<= 1;
    }

    for (i = 0; i < wdlevels; i ++) {
      if (haar_lift_inverse3d_haar_2dstep(s, w, d, 1, slicestride, work) < 0) {
	return -1;
      }

      w <<= 1;
      d <<= 1;
    }
  
    for (i = 0; i < hdlevels; i ++) {
      if (haar_lift_inverse3d_haar_2dstep(s, h, d, rowstride, slicestride, work) < 0) {
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

    if (haar_lift_inverse3d_haar_step(s, w, h, d, rowstride, slicestride, work) < 0) {
      return -1;
    }

    w <<= 1;
    h <<= 1;
    d <<= 1;
  }

  return 0;
}

int 
haar_lift_inverse3d_haar_step(double *s,
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

      haar_lift_inverse1d_haar_step(s + o,
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

      haar_lift_inverse1d_haar_step(s + o,
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

      haar_lift_inverse1d_haar_step(s + o,
				    width,
				    1,
				    work);

    }
  }

  return 0;
}

int 
haar_lift_inverse3d_haar_2dstep(double *s,
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
    haar_lift_inverse1d_haar_step(s + i*stride,
				  height,
				  rowstride,
				  work);
    
  }

  /*
   * 1D Transform on Rows
   */
  for (i = 0; i < height; i ++) {
    haar_lift_inverse1d_haar_step(s + i*rowstride,
				  width,
				  stride,
				  work);

  }

  return 0;
}

