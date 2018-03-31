
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
