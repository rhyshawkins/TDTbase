//
//    Spherical Subdivision/Wavelet library
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

#include <math.h>

#include "transform.h"
#include "vertex3.h"

void
transform_from_axis_angle(double x, double y, double z,
			  double theta,
			  double *transform)
{
  double c, s, C;
  
  c = cos(theta);
  s = sin(theta);
  C = 1 - c;

  /* xxC+c & xyC-zs & xzC+ys */
  /* yxC+zs & yyC+c & yzC-xs */
  /* zxC-ys & zyC+xs & zzC+c */
    
  transform[0] = x*x*C + c;
  transform[1] = x*y*C - z*s;
  transform[2] = x*z*C + y*s;
  transform[3] = 0.0;

  transform[4] = y*x*C + z*s;
  transform[5] = y*y*C + c;
  transform[6] = y*z*C - x*s;
  transform[7] = 0.0;

  transform[8] = z*x*C - y*s;
  transform[9] = z*y*C + x*s;
  transform[10] = z*z*C + c;
  transform[11] = 0.0;

  transform[12] = 0.0;
  transform[13] = 0.0;
  transform[14] = 0.0;
  transform[15] = 1.0;
}

void
transform_compute_axis_angle(double x1, double y1, double z1,
			     double x2, double y2, double z2,
			     double *ax, double *ay, double *az,
			     double *theta)
{
  double nx, ny, nz;
  double nl;
  double al;
  double bl;

  vertex3_cross(x1, y1, z1,
		x2, y2, z2,
		&nx, &ny, &nz);

  nl = sqrt(nx*nx + ny*ny + nz*nz);

  if (nl > 0.0) {

    al = sqrt(x1*x1 + y1*y1 + z1*z1);
    bl = sqrt(x2*x2 + y2*y2 + z2*z2);
    
    *ax = nx/nl;
    *ay = ny/nl;
    *az = nz/nl;

    *theta = asin(nl/(al * bl));

  } else {
    /* Return 0 rotation */
    *ax = 0.0;
    *ay = 0.0;
    *az = 1.0;
    *theta = 0.0;
  }
}
