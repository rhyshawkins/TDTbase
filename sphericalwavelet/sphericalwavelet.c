
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sphericalwavelet.h"

sphericalwavelet_t *
sphericalwavelet_create(int degree)
{
  return NULL;
}

void
sphericalwavelet_destroy(sphericalwavelet_t *s)
{
}

/*
 * Coefficient count for tetrahedron subdivision
 */
int sw_tetra_ncoeff(int depth)
{
  // Progression is 4 6 24
  return -1;
}

int sw_tetra_totalcoeff(int depth)
{
  return -1;
}

/*
 * Indexing for tetrahedron subdivision
 */
int sw_tetra_depth_of_index(int index)
{
  return -1;
}
