#ifndef sphericalwavelet_h
#define sphericalwavelet_h

typedef struct _sphericalwavelet sphericalwavelet_t;

typedef enum {
  SPHERICALWAVELET_TETRAHEDRON,
  SPHERICALWAVELET_OCTAHEDRON,
  SPHERICALWAVELET_ICOSAHEDRON
} sphericalwavelet_base_t;

sphericalwavelet_t *
sphericalwavelet_create(int degree);

void
sphericalwavelet_destroy(sphericalwavelet_t *s);

/*
 * Coefficient count for tetrahedron subdivision
 */
int sw_tetra_ncoeff(int depth);
int sw_tetra_totalcoeff(int depth);

/*
 * Indexing for tetrahedron subdivision
 */
int sw_tetra_depth_of_index(int index);


#endif /* sphericalwavelet_h */
