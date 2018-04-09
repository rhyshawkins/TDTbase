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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cdf97_lift_impulse.h"

static const double a1 = -1.58613434200;
static const double a2 = -0.05298011854;
static const double a3 =  0.88291107620;
static const double a4 =  0.44350685220;

static const double k1 = 0.81289306611596146; // 1/1.230174104914
static const double k2 = 0.61508705245700002; // 1.230174104914/2

static const double ik1 = 1.230174104914;
static const double ik2 = 1.6257861322319229;

static int cdf97_lift_inverse1d_interlaced(double *v,
					   int width,
					   int stride);

int cdf97_lift_impulse_1dsize(int degree, int impulse_depth)
{
  if (impulse_depth < 0 || impulse_depth > degree) {
    return -1;
  }

  if (impulse_depth == degree) {
    return 9; /* Width of cdf97 filter */
  } else {
    return cdf97_lift_impulse_1dsize(degree, impulse_depth + 1) + (1 << (degree - impulse_depth)) * 7;
  }
}

int cdf97_lift_impulse_1d(int degree, 
			  int impulse_degree, 
			  double impulse, 
			  double *v, 
			  int size,
			  int offset)
{
  memset(v, 0, sizeof(double) * size);
  v[offset] = impulse;
  
  return cdf97_lift_inverse1d_interlaced(v, size, 1);
}

static int cdf97_lift_inverse1d_interlaced(double *v,
					   int width,
					   int stride)
{
  int i;

  for (i = 0; i < width; i ++) {
    if (i % 2 == 0) {
      v[i] *= ik1;
    } else {
      v[i] *= ik2;
    }
  }

  for (i = 2; i < width; i += 2) {
    v[stride*i] -= a4 * (v[stride*(i - 1)] + v[stride*(i + 1)]);
  }
  v[0] -= 2.0 * a4 * v[stride];
  
  for (i = 1; i < (width - 1); i += 2) {
    v[stride*i] -= a3 * (v[stride*(i - 1)] + v[stride*(i + 1)]);
  }
  v[stride*(width - 1)] -= 2.0 * a3 * v[stride*(width - 2)];
  
  for (i = 2; i < width; i += 2) {
    v[stride*i] -= a2 * (v[stride*(i - 1)] + v[stride*(i + 1)]);
  }
  v[0] -= 2.0 * a2 * v[stride];
  
  for (i = 1; i < (width - 1); i += 2) {
    v[stride*i] -= a1 * (v[stride*(i - 1)] + v[stride*(i + 1)]);
  }
  v[stride*(width - 1)] -= 2.0 * a1 * v[stride*(width - 2)];

  return 0;
}
