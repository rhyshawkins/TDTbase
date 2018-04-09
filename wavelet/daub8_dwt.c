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

#include "daub8_dwt.h"

#include "boundary.h"

/*
 * These coefficients are scaled by 1/sqrt(2) from the true orthongal terms to match
 * the normalization of the lifting version. This then requires a factor of 2 in the 
 * inverse transform.
 */
#define H0  0.162901714025649180  //  0.23037781330889650086329118304,
#define H1  0.505472857545914400  //  0.71484657055291564708992195527,
#define H2  0.446100069123379800  //  0.63088076792985890788171633830,
#define H3 -0.019787513117822320  // -0.02798376941685985421141374718,
#define H4 -0.132253583684519870  // -0.18703481171909308407957067279,
#define H5  0.021808150237088625  //  0.03084138183556076362721936253,
#define H6  0.023251800535490877  //  0.03288301166688519973540751355,
#define H7 -0.007493494665180735  // -0.01059740178506903210488320852

#define G0 H7
#define G1 -H6
#define G2 H5
#define G3 -H4
#define G4 H3
#define G5 -H2
#define G6 H1
#define G7 -H0

#define HB0 (2.0*H0)
#define HB1 (2.0*H1)
#define HB2 (2.0*H2)
#define HB3 (2.0*H3)
#define HB4 (2.0*H4)
#define HB5 (2.0*H5)
#define HB6 (2.0*H6)
#define HB7 (2.0*H7)

#define GB0 (2.0*G0)
#define GB1 (2.0*G1)
#define GB2 (2.0*G2)
#define GB3 (2.0*G3)
#define GB4 (2.0*G4)
#define GB5 (2.0*G5)
#define GB6 (2.0*G6)
#define GB7 (2.0*G7)

/*
 * 1D Forward
 */
int
daub8_dwt_forward1d_daub8(double *s,
			  int width,
			  int stride,
			  double *work)
{
  int w;

  w = width;

  while (w > 1) {

    if (daub8_dwt_forward1d_daub8_step(s, w, stride, work) < 0) {
      return -1;
    }

    w >>= 1;
  }

  return 0;
}

int
daub8_dwt_forward1d_daub8_step(double *s,
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
      H5 * work[wavelet_boundary_periodic(2*i + 5, width)] +
      H6 * work[wavelet_boundary_periodic(2*i + 6, width)] +
      H7 * work[wavelet_boundary_periodic(2*i + 7, width)];
				
    s[(width/2 + i)*stride] =
      G0 * work[2*i] + 
      G1 * work[2*i + 1] + 
      G2 * work[wavelet_boundary_periodic(2*i + 2, width)] + 
      G3 * work[wavelet_boundary_periodic(2*i + 3, width)] +
      G4 * work[wavelet_boundary_periodic(2*i + 4, width)] + 
      G5 * work[wavelet_boundary_periodic(2*i + 5, width)] +
      G6 * work[wavelet_boundary_periodic(2*i + 6, width)] +
      G7 * work[wavelet_boundary_periodic(2*i + 7, width)];
      
  }

  return 0;
}

/*
 * 1D Inverse
 */
int
daub8_dwt_inverse1d_daub8(double *s,
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

    if (daub8_dwt_inverse1d_daub8_step(s, w, stride, work) < 0) {
      return -1;
    }

    w <<= 1;
  }

  return 0;
}

int
daub8_dwt_inverse1d_daub8_step(double *s,
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
      HB6 * work[wavelet_boundary_periodic(2*i - 6, width)] +
      GB6 * work[wavelet_boundary_periodic(2*i - 5, width)] +
      HB4 * work[wavelet_boundary_periodic(2*i - 4, width)] + 
      GB4 * work[wavelet_boundary_periodic(2*i - 3, width)] + 
      HB2 * work[wavelet_boundary_periodic(2*i - 2, width)] + 
      GB2 * work[wavelet_boundary_periodic(2*i - 1, width)] + 
      HB0 * work[2*i + 0] + 
      GB0 * work[2*i + 1];
				
    s[(2*i + 1)*stride] =
      HB7 * work[wavelet_boundary_periodic(2*i - 6, width)] + 
      GB7 * work[wavelet_boundary_periodic(2*i - 5, width)] + 
      HB5 * work[wavelet_boundary_periodic(2*i - 4, width)] + 
      GB5 * work[wavelet_boundary_periodic(2*i - 3, width)] + 
      HB3 * work[wavelet_boundary_periodic(2*i - 2, width)] + 
      GB3 * work[wavelet_boundary_periodic(2*i - 1, width)] + 
      HB1 * work[2*i + 0] + 
      GB1 * work[2*i + 1];
  }

  return 0;
}

/*
 * 2D Forward
 */
int 
daub8_dwt_forward2d_daub8(double *s, 
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

    if (daub8_dwt_forward2d_daub8_step(s, w, h, stride, work) < 0) {
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
      if (daub8_dwt_forward1d_daub8_step(s, w, 1, work) < 0) {
	return -1;
      }
      
      w >>= 1;
    }
    
    while (h > 1) {
      /*
       * Left with a single column
       */
      if (daub8_dwt_forward1d_daub8_step(s, h, stride, work) < 0) {
	return -1;
      }
      
      h >>= 1;
    }
  }

  return 0;
}

int 
daub8_dwt_forward2d_daub8_step(double *s, 
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
    
    daub8_dwt_forward1d_daub8_step(s + i,
				  height,
				  stride,
				  work);
  }

  /*
   * 1D Transform on Rows
   */
  for (i = 0; i < height; i ++) {

    daub8_dwt_forward1d_daub8_step(s + stride*i,
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
daub8_dwt_inverse2d_daub8(double *s,
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
      
      if (daub8_dwt_inverse1d_daub8_step(s, w, 1, work) < 0) {
	return -1;
      }
      
      w <<= 1;
    }
    
    for (i = 0; i < hlevels; i ++) {
      
      if (daub8_dwt_inverse1d_daub8_step(s, h, stride, work) < 0) {
	return -1;
      }
      
      h <<= 1;
    }
  }

  for (i = 0; i <= levels; i ++) {

    if (daub8_dwt_inverse2d_daub8_step(s, w, h, stride, work) < 0) {
      return -1;
    }

    w <<= 1;
    h <<= 1;
  }

  return 0;
}

int 
daub8_dwt_inverse2d_daub8_step(double *s, 
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

    daub8_dwt_inverse1d_daub8_step(s + i*stride,
				  width,
				  1,
				  work);
    
  }

  /*
   * 1D Inverse Transform on Columns
   */

  for (i = 0; i < width; i ++) {

    daub8_dwt_inverse1d_daub8_step(s + i,
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
daub8_dwt_forward3d_daub8(double *s,
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

    if (daub8_dwt_forward3d_daub8_step(s, w, h, d, rowstride, slicestride, work) < 0) {
      return -1;
    }

    w >>= 1;
    h >>= 1;
    d >>= 1;
  }

  if (!subtile) {
    if (d == 1) {
      while (w > 1 && h > 1) {

	if (daub8_dwt_forward3d_daub8_2dstep(s, w, h, 1, rowstride, work) < 0) {
	  return -1;
	}

	w >>= 1;
	h >>= 1;
      }

      while (w > 1) {
      
	if (daub8_dwt_forward1d_daub8_step(s, w, 1, work) < 0) {
	  return -1;
	}
      
	w >>= 1;
      }

      while (h > 1) {

	if (daub8_dwt_forward1d_daub8_step(s, h, rowstride, work) < 0) {
	  return -1;
	}

	h >>= 1;
      }

    } else if (h == 1) {

      while (w > 1 && d > 1) {

	if (daub8_dwt_forward3d_daub8_2dstep(s, w, d, 1, slicestride, work) < 0) {
	  return -1;
	}

	w >>= 1;
	d >>= 1;

      }

      while (w > 1) {

	if (daub8_dwt_forward1d_daub8_step(s, w, 1, work) < 0) {
	  return -1;
	}

	w >>= 1;
      }

      while (d > 1) {

	if (daub8_dwt_forward1d_daub8_step(s, d, slicestride, work) < 0) {
	  return -1;
	}

	d >>= 1;
      }

    } else if (w == 1) {

      while (h > 1 && d > 1) {

	if (daub8_dwt_forward3d_daub8_2dstep(s, h, d, rowstride, slicestride, work) < 0) {
	  return -1;
	}

	h >>= 1;
	d >>= 1;

      }

      while (h > 1) {

	if (daub8_dwt_forward1d_daub8_step(s, h, rowstride, work) < 0) {
	  return -1;
	}

	h >>= 1;
      }

      while (d > 1) {

	if (daub8_dwt_forward1d_daub8_step(s, d, slicestride, work) < 0) {
	  return -1;
	}

	d >>= 1;
      }

    }
  }
  
  return 0;
}
   

int 
daub8_dwt_forward3d_daub8_step(double *s,
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

      daub8_dwt_forward1d_daub8_step(s + o,
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

      daub8_dwt_forward1d_daub8_step(s + o,
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

      daub8_dwt_forward1d_daub8_step(s + o,
				    depth,
				    slicestride,
				    work);

    }
  }

  return 0;
}

int 
daub8_dwt_forward3d_daub8_2dstep(double *s,
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
    daub8_dwt_forward1d_daub8_step(s + i*rowstride,
				  width,
				  stride,
				  work);

  }

  /*
   * 1D Transform on Columns
   */
  for (i = 0; i < width; i ++) {
    daub8_dwt_forward1d_daub8_step(s + i*stride,
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
daub8_dwt_inverse3d_daub8(double *s,
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
      if (daub8_dwt_inverse1d_daub8_step(s, w, 1, work) < 0) {
	return -1;
      }

      w <<= 1;
    }

    for (i = 0; i < hlevels; i ++) {
      if (daub8_dwt_inverse1d_daub8_step(s, h, rowstride, work) < 0) {
	return -1;
      }

      h <<= 1;
    }

    for (i = 0; i < dlevels; i ++) {
      if (daub8_dwt_inverse1d_daub8_step(s, d, slicestride, work) < 0) {
	return -1;
      }

      d <<= 1;
    }

    /*
     * 2D expansion for non-square
     */
    for (i = 0; i < whlevels; i ++) {
      if (daub8_dwt_inverse3d_daub8_2dstep(s, w, h, 1, rowstride, work) < 0) {
	return -1;
      }

      w <<= 1;
      h <<= 1;
    }

    for (i = 0; i < wdlevels; i ++) {
      if (daub8_dwt_inverse3d_daub8_2dstep(s, w, d, 1, slicestride, work) < 0) {
	return -1;
      }

      w <<= 1;
      d <<= 1;
    }
  
    for (i = 0; i < hdlevels; i ++) {
      if (daub8_dwt_inverse3d_daub8_2dstep(s, h, d, rowstride, slicestride, work) < 0) {
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

    if (daub8_dwt_inverse3d_daub8_step(s, w, h, d, rowstride, slicestride, work) < 0) {
      return -1;
    }

    w <<= 1;
    h <<= 1;
    d <<= 1;
  }

  return 0;
}

int 
daub8_dwt_inverse3d_daub8_step(double *s,
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

      daub8_dwt_inverse1d_daub8_step(s + o,
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

      daub8_dwt_inverse1d_daub8_step(s + o,
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

      daub8_dwt_inverse1d_daub8_step(s + o,
				    width,
				    1,
				    work);

    }
  }

  return 0;
}

int 
daub8_dwt_inverse3d_daub8_2dstep(double *s,
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
    daub8_dwt_inverse1d_daub8_step(s + i*stride,
				  height,
				  rowstride,
				  work);
    
  }

  /*
   * 1D Transform on Rows
   */
  for (i = 0; i < height; i ++) {
    daub8_dwt_inverse1d_daub8_step(s + i*rowstride,
				  width,
				  stride,
				  work);

  }

  return 0;
}

