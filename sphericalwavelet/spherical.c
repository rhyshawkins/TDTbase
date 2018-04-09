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

#include "spherical.h"

double
spherical_gcdist(double lon1, double lat1, double lon2, double lat2, double r)
{
  double rlon1, rlat1;
  double rlon2, rlat2;
  double dlon, dlat, dsigma;
  double sdlon2;
  double sdlat2;
  
  rlon1 = lon1 * M_PI/180.0;
  rlat1 = lat1 * M_PI/180.0;
    
  rlon2 = lon2 * M_PI/180.0;
  rlat2 = lat2 * M_PI/180.0;

  dlat = fabs(rlat2 - rlat1);
  dlon = fabs(rlon2 - rlon1);

  sdlat2 = sin(dlat/2.0);
  sdlat2 *= sdlat2;

  sdlon2 = sin(dlon/2.0);
  sdlon2 *= sdlon2;
  
  /*
   * Haversine formula for better accuracy for small angles
   */
  dsigma = 2.0 * asin(sqrt(sdlat2 + cos(rlat1)*cos(rlat2)*sdlon2));

  /*
   * Standard cosine rule computation
   */
  return r * dsigma;
}
