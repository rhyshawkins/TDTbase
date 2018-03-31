
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wavelet2d_tile.h"

#include "slog.h"

wavelet2d_tile_t *
wavelet2d_tile_create(int size, int overlap)
{
  wavelet2d_tile_t *w;
  int i;

  if (size <= 0 || 
      (size & (size - 1))) {
    ERROR("invalid size");
    return NULL;
  }
  
  if (overlap <= 0 || overlap >= size) {
    ERROR("invalid overlap");
    return NULL;
  }

  w = malloc(sizeof(wavelet2d_tile_t));
  if (w == NULL) {
    ERROR("failed to allocate structure");
    return NULL;
  }

  w->size = size;
  w->overlap = overlap;
  w->rowstride = size + 2*overlap;

  w->coeff = malloc(sizeof(double) * w->rowstride * w->rowstride);
  if (w->coeff == NULL) {
    ERROR("failed to allocate coefficients");
    return NULL;
  }
  memset(w->coeff, 0, sizeof(double) * w->rowstride * w->rowstride);
  
  for (i = 0; i < NEIGHBOUR_COUNT; i ++) {
    w->neighbours[i] = NULL;
    w->edges[i] = EDGE_REFLECT;
  }

  return w;
}

void
wavelet2d_tile_destroy(wavelet2d_tile_t *w)
{
  if (w != NULL) {
    free(w->coeff);
    free(w);
  }
}


