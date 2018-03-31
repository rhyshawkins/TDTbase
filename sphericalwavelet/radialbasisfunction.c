
#include <math.h>

#include "radialbasisfunction.h"

#include "spherical.h"

static double
radius2D(double cx, double cy, double x, double y);

static double
sphericalradius2D(double clon, double clat, double lon, double lat);

/*
 * Linear
 */
double rbf_linear1D(double x, double r)
{
  if (x < 0.0 || x > r) {
    return 0.0;
  } else {
    return 1.0 - x/r;
  }
}
double rbf_linear2D(double cx, double cy, double r,
		    double x, double y)
{
  return rbf_linear1D(radius2D(cx, cy, x, y), r);
}

double rbf_linear2DS2(double clon, double clat, double r,
		      double lon, double lat)
{
  return rbf_linear1D(sphericalradius2D(clon, clat, lon, lat), r);
}

double rbf_linear3D(double cx, double cy, double cz, double rxy, double rz,
		    double x, double y, double z)
{
  return
    rbf_linear1D(fabs(cz - z), rz) *
    rbf_linear1D(radius2D(cx, cy, x, y), rxy);
}

/*
 * Quadratic
 */
double rbf_quadratic2D(double cx, double cy, double r,
		       double x, double y)
{
  double R;

  R = radius2D(cx, cy, x, y);

  if (R > r) {
    return 0.0;
  } else {
    return -(R - r)*(R + r)/(r*r);
  }
}

double rbf_quadratic3D(double cx, double cy, double cz, double rxy, double rz,
		       double x, double y, double z)
{
  double R;
  double Rz;

  Rz = fabs(cz - z);
  if (Rz > rz) {
    return 0.0;
  } else {
    R = radius2D(cx, cy, x, y);
    if (R > rxy) {
      return 0.0;
    } else {
      return
	(R - rxy)*(R + rxy)/(rxy*rxy) *
	(Rz - rz)*(R + rz)/(rz*rz);
    }
  }
}

/*
 * Piecewise Cubic smooth, 0 gradients at radius = 0, r
 */
double rbf_cubic2D(double cx, double cy, double r,
		   double x, double y)
{
  double a;
  double b;
  double R;

  R = radius2D(cx, cy, x, y);
  if (R > r) {
    return 0.0;
  } else {
    a = 2.0/(r * r * r);
    b = -3.0/(r * r);

    return a*R*R*R + b*R*R + 1.0;
  }
}
  
double rbf_cubic3D(double cx, double cy, double cz, double rxy, double rz,
		   double x, double y, double z)
{
  double az;
  double bz;
  double Rz;
  
  double a;
  double b;
  double R;

  Rz = fabs(cz - z);
  if (Rz > rz) {
    return 0.0;
  } else {
    
    R = radius2D(cx, cy, x, y);
    if (R > rxy) {
      return 0.0;
    } else {

      az = 2.0/(rz * rz * rz);
      bz = -3.0/(rz * rz);
      
      a = 2.0/(rxy * rxy * rxy);
      b = -3.0/(rxy * rxy);
      
      return
	(az*Rz*Rz*Rz + bz*Rz*Rz + 1.0) * 
	(a*R*R*R + b*R*R + 1.0);
    }
  }
}

/*
 * Cosine: Peak of 1.0 at radius = 0, 
 */
double rbf_cosine2D(double cx, double cy, double r,
		    double x, double y)
{
  double R;
  
  R = radius2D(cx, cy, x, y);
  if (R > r) {
    return 0.0;
  } else {
    return 0.5*cos(M_PI * R/r) + 0.5;
  }
}

double rbf_cosine3D(double cx, double cy, double cz, double rxy, double rz,
		    double x, double y, double z)
{
  double R;
  double Rz;

  Rz = fabs(cz - z);
  if (Rz > rz) {
    return 0.0;
  } else {
    
    R = radius2D(cx, cy, x, y);
    if (R > rxy) {
      return 0.0;
    } else {
      return
	(0.5*cos(M_PI * Rz/rz) + 0.5) * 
	(0.5*cos(M_PI * R/rxy) + 0.5);
    }
  }
}

/*
 * Gaussian: Un-normalized, sigma = r
 */
double rbf_gaussian2D(double cx, double cy, double r,
		      double x, double y)
{
  double R;

  R = radius2D(cx, cy, x, y);
  if (R > r) {
    return 0.0;
  } else {
    return exp(-0.5 * 9.0 * (R*R)/(r*r));
  }
}

double rbf_gaussian3D(double cx, double cy, double cz, double rxy, double rz,
		      double x, double y, double z)
{
  double R;
  double Rz;

  Rz = fabs(cz - z);
  if (Rz > rz) {
    return 0.0;
  } else {
    
    R = radius2D(cx, cy, x, y);
    if (R > rxy) {
      return 0.0;
    } else {
      return
	exp(-0.5 * 9.0 * (Rz*Rz)/(rz*rz)) * 
	exp(-0.5 * 9.0 * (R*R)/(rxy*rxy));
    }
  }
}

/*
 * Lanczos: 
 */

static double lanczos(double x, double a)
{
  if (x <= 0.0) {
    return 1.0;
  } else if (x < a) {
    return a * sin(M_PI * x) * sin(M_PI * x/a)/(M_PI * M_PI * x * x);
  } else {
    return 0.0;
  }
}

double rbf_lanczos2D(double cx, double cy, double r,
		     double x, double y)
{
  double R;

  R = radius2D(cx, cy, x, y);
  if (R > r) {
    return 0.0;
  } else {
    return lanczos(R, r);
  }
}

double rbf_lanczos3D(double cx, double cy, double cz, double rxy, double rz,
		     double x, double y, double z)
{
  double R;
  double Rz;

  Rz = fabs(cz - z);
  if (Rz > rz) {
    return 0.0;
  } else {
    
    R = radius2D(cx, cy, x, y);
    if (R > rxy) {
      return 0.0;
    } else {
      return
	lanczos(Rz, rz) *
	lanczos(R, rxy);
    }
  }
}

static double
radius2D(double cx, double cy, double x, double y)
{
  double dx;
  double dy;

  dx = cx - x;
  dy = cy - y;

  return sqrt(dx*dx + dy*dy);
}

static double
sphericalradius2D(double clon, double clat, double lon, double lat)
{
  return spherical_gcdist(clon, clat, lon, lat, 1.0);
}

