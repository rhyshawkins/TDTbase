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


