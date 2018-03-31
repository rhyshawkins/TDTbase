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
