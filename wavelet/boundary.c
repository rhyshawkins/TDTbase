
#include "boundary.h"

int
wavelet_boundary_periodic(int i, int width)
{
  if (i < 0) {
    return wavelet_boundary_periodic(i + width, width);
  } else {
    return i % width;
  }
}

int
wavelet_boundary_reflect(int i, int width)
{
  if (i < 0) {
    return wavelet_boundary_reflect(-i, width);
  } else if (i >= width) {
    return wavelet_boundary_reflect(2*width - 2 - i, width);
  } else {
    return i;
  }
}



