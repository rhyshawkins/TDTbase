
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "icosahedron.h"

int main(int argc, char *argv[])
{
  int i;
  double d;
  int nt;

  FILE *fp;
  
  fp = fopen("icosahedron_resolution.txt", "w");
  if (fp == NULL) {
    fprintf(stderr, "error: failed to create file\n");
    return -1;
  }
  
  fprintf(fp, "Degree : %15s %15s %15s %15s %15s\n", "Triangles", "Radian", "Degree", "Earth km", "SH Deg.");
  
  for (i = 1; i <= 10; i ++) {
    d = icosahedron_angle(i);
    if (d < 0.0) {
      break;
    }
    nt = icosahedron_ntriangles(i);
      
    fprintf(fp, "    %02d : %15d %15.9f %15.9f %15.9f %15d\n", i, nt, d, d * 180.0/M_PI, d * 6371.0, (int)ceil(2.0*M_PI/d));

  }

  fclose(fp);

  return 0;
}
    
  
