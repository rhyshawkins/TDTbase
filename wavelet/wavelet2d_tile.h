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

#ifndef wavelet2d_tile_h
#define wavelet2d_tile_h

typedef struct wavelet2d_tile_ wavelet2d_tile_t;

typedef enum {
  NEIGHBOUR_LEFT = 0,
  NEIGHBOUR_TOP,
  NEIGHBOUR_RIGHT,
  NEIGHBOUR_BOTTOM,
  NEIGHBOUR_TOPLEFT,
  NEIGHBOUR_TOPRIGHT,
  NEIGHBOUR_BOTTOMRIGHT,
  NEIGHBOUR_BOTTOMLEFT,
  NEIGHBOUR_COUNT
} wavelet2d_neighbour_t;
  
typedef enum {
  EDGE_REFLECT = 0,
  EDGE_NEIGHBOUR_LEFT,
  EDGE_NEIGHBOUR_LEFT_REVERSE,
  EDGE_NEIGHBOUR_TOP,
  EDGE_NEIGHBOUR_TOP_REVERSE,
  EDGE_NEIGHBOUR_RIGHT,
  EDGE_NEIGHBOUR_RIGHT_REVERSE,
  EDGE_NEIGHBOUR_BOTTOM,
  EDGE_NEIGHBOUR_BOTTOM_REVERSE,
  
  EDGE_COUNT
} wavelet2d_edge_t;

struct wavelet2d_tile_ {

  int size;
  int overlap;
  int rowstride;

  wavelet2d_tile_t *neighbours[NEIGHBOUR_COUNT];
  wavelet2d_edge_t edges[NEIGHBOUR_COUNT];
  
  double *coeff;
};

wavelet2d_tile_t *
wavelet2d_tile_create(int size, int overlap);

void
wavelet2d_tile_destroy(wavelet2d_tile_t *w);

#endif /* wavelet2d_tile_h */
