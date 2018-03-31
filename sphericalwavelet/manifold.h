#ifndef manifold_h
#define manifold_h

#include "vertex3.h"
#include "edge.h"
#include "triangle.h"

typedef int (*vertex_count_t)(int);
typedef int (*edge_count_t)(int);
typedef int (*triangle_count_t)(int);

typedef enum {
  MANIFOLD_NORMALIZE = -1,
  MANIFOLD_ASIS = 0
} manifold_normalization_t;

typedef struct _manifold manifold_t;
struct _manifold {

  int degree;

  /* Functions for computing the number of vertices/edges/triangles at each subdivision level */
  vertex_count_t nverticesatdepth;
  edge_count_t nedgesatdepth;
  triangle_count_t ntrianglesatdepth;

  /* Vertices are index from 0 .. total number of vertices across all depths */
  int nvertices;
  vertex3_t *vertices;

  /* Edges and triangles are index by (depth, 0..number of edges/triangles per depth - 1) */
  int *nedges;
  edge_t **edges;

  int ntotaltriangles;
  int *ntriangles;
  triangle_t **triangles;
};

manifold_t *
manifold_create(int degree,
		vertex_count_t nverticesatdepth,
		edge_count_t nedgesatdepth,
		triangle_count_t ntrianglesatdepth);

void
manifold_destroy(manifold_t *m);

int
manifold_subdivide(manifold_t *m,
		   int depth);

int
manifold_build_neighbors(manifold_t *m);

int
manifold_compute_areas(manifold_t *m,
		       int depth);

int
manifold_save_geo(manifold_t *m,
		  const char *filename);

int
manifold_valid(manifold_t *m);

int
manifold_get_children_vertices(manifold_t *m, int pi, int *ci);

int
manifold_get_parent(manifold_t *m, int ci);

/*
 * Lookup functions
 */

int
manifold_find_enclosing_triangle(manifold_t *m,
				 triangle_workspace_t *workspace,
				 double lon, double lat,
				 int *ti,
				 double *ba, double *bb, double *bc);

int
manifold_find_nearest_vertex(manifold_t *m,
			     triangle_workspace_t *workspace,
			     double lon, double lat,
			     int *vi);

int
manifold_compute_barycentre_coordinates(manifold_t *m,
					triangle_workspace_t *workspace,
					double lon,
					double lat,
					int *va, int *vb, int *vc,
					double *ba, double *bb, double *bc);

/*
 * Helper functions for construction
 */
int
manifold_set_vertex(manifold_t *m,
		    int vi,
		    double x, double y, double z,
		    int depth,
		    int normalize);

int
manifold_set_edge(manifold_t *m,
		  int depth,
		  int ei,
		  int va, int vb);

int
manifold_set_triangle(manifold_t *m,
		      int depth,
		      int ti,
		      int va, int vb, int vc,
		      int eab, int ebc, int eca);


int
manifold_get_triangle_other_vertex(manifold_t *m,
				   int depth,
				   int ei,
				   int ti,
				   int *v);

int
manifold_get_edge_f_vertices(manifold_t *m,
			     int depth,
			     int ei,
			     int *f);
			    
int
manifold_get_edge_e_vertices(manifold_t *m,
			     int depth,
			     int ei,
			     int *e);

#endif /* manifold_h */
  
