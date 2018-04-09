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

#include "daubechies.h"

static const double H1 = 0.48296291314453414337487159986;
static const double H2 = 0.83651630373780790557529378092;
static const double H3 = 0.22414386804201338102597276224;
static const double H4 = -0.12940952255126038117444941881;

int
daubechies2d_forward_d4(double *s,
			int width,
			int height,
			int stride,
			double *work)
{
  return -1;
}

int
daubechies2d_inverse_d4(double *s,
			int width,
			int height,
			int stride,
			double *work)
{
  return -1;
}

int
daubechies1d_forward_d4_step(double *s,
			     int width,
			     int step,
			     double *work)
{
  int i;
  int j;
  int mod;
  int hw;

  mod = width * step;

  hw = width/2;

  for (i = 0, j = 0; i < hw; i ++, j += (2 * step)) {
    work[i] = 
      H1 * s[j] +
      H2 * s[j + step] +
      H3 * s[(j + 2*step) % mod] +
      H4 * s[(j + 3*step) % mod];

    work[hw + i] = 
      H4 * s[j] -
      H3 * s[j + step] +
      H2 * s[(j + 2*step) % mod] -
      H1 * s[(j + 3*step) % mod];
  }

  for (i = 0, j = 0; i < width; i ++, j += step) {
    s[j] = work[i];
  }

  return 0;
}

int
daubechies1d_inverse_d4_step(double *s,
			     int width,
			     int step,
			     double *work)
{
  /* for (i = 0, j = 0; i < hw; i ++, j += (2 * step)) { */

  /*   /\* Even *\/ */
  /*   work[2*i] =  */
  /*     H[2] * s[(i + width*step - 1) % fw] +  */
  /*     H1 * s[ */
  

  /* return -1; */
  
}
