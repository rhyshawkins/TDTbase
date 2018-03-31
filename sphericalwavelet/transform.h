#ifndef transform_h
#define transform_h

void
transform_from_axis_angle(double x, double y, double z,
			  double theta,
			  double *transform);

void
transform_compute_axis_angle(double x1, double y1, double z1,
			     double x2, double y2, double z2,
			     double *ax, double *ay, double *az,
			     double *theta);
			     
#endif /* transform_h */
