#ifndef edge_h
#define edge_h

struct _edge {
  int a;
  int b;

  int parent;
  int child_edges[2];

  int triangles[2];
};
typedef struct _edge edge_t;

void
edge_init(edge_t *e);

int
edge_add_triangle(edge_t *e, int ti);

#endif /* edge_h */
