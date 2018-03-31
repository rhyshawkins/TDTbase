
#include <stdio.h>

#include "manifold.h"

#include "slog.h"

static int
manifold_all_vertices_set(manifold_t *m);

static int
manifold_all_edges_set(manifold_t *m);

static int
manifold_all_parent_child_valid(manifold_t *m);

int
manifold_valid(manifold_t *m)
{
  int ec;

  ec = 0;
  
  ec += manifold_all_vertices_set(m);

  ec += manifold_all_edges_set(m);

  ec += manifold_all_parent_child_valid(m);

  return (ec == 0);
}

static int
manifold_all_vertices_set(manifold_t *m)
{
  int i;
  int j;
  int ec;
  int expected_children;

  ec = 0;
  for (i = 0; i < m->nvertices; i ++) {

    
    if (m->vertices[i].depth < 0) {
      ERROR("vertex %d depth not set", i);
      ec ++;

    } else {

      if (m->vertices[i].depth > 0) {
	if (m->vertices[i].depth == m->degree) {
	  expected_children = 0;
	} else {
	  expected_children = 4;
	}
	
	if (m->vertices[i].parent < 0) {
	  ERROR("vertex %d parent not set", i);
	  ec ++;
	}

	if (m->vertices[i].v[0] < 0 ||
	    m->vertices[i].v[1] < 0) {
	  ERROR("vertex %d v vertices not set (%d %d)",
		i,
		m->vertices[i].v[0],
		m->vertices[i].v[1]);
	  ec ++;
	}
	  
	if (m->vertices[i].f[0] < 0 ||
	    m->vertices[i].f[1] < 0) {
	  ERROR("vertex %d f vertices not set (%d %d)",
		i,
		m->vertices[i].f[0],
		m->vertices[i].f[1]);
	  ec ++;
	}

	if (m->vertices[i].e[0] < 0 ||
	    m->vertices[i].e[1] < 0 ||
	    m->vertices[i].e[2] < 0 ||
	    m->vertices[i].e[3] < 0) {
	  ERROR("vertex %d e vertices not set (%d %d %d %d)",
		i,
		m->vertices[i].e[0],
		m->vertices[i].e[1],
		m->vertices[i].e[2],
		m->vertices[i].e[3]);
	  ec ++;
	}

      } else {
	if (i == 0 || i == 1) {
	  expected_children = 0;
	} else {
	  if (m->degree > 0) {
	    expected_children = 3;
	  } else {
	    expected_children = 0;
	  }
	}
      }

      for (j = 0; j < expected_children; j ++) {
	if (m->vertices[i].children[j] < 0) {
	  ERROR("vertex %d missing child %d", i, j);
	  ec ++;
	}
      }

      for (j = expected_children; j < 4; j ++) {
	if (m->vertices[i].children[j] >= 0) {
	  ERROR("vertex %d has extra child %d", i, j);
	  ec ++;
	}
      }
    }
  }

  return ec;
}

static int
manifold_all_edges_set(manifold_t *m)
{
  int d;
  int i;
  int ec;

  ec = 0;

  for (d = 0; d <= m->degree; d ++) {

    for (i = 0; i < m->nedges[d]; i ++) {

      if (m->edges[d][i].a < 0 ||
	  m->edges[d][i].b < 0) {
	ERROR("invalid edge %d (%d %d)",
	      i,
	      m->edges[d][i].a,
	      m->edges[d][i].b);
	ec ++;

      } else {

	if (m->edges[d][i].triangles[0] < 0 ||
	    m->edges[d][i].triangles[1] < 0) {
	  ERROR("edge %d missing triangles (%d %d)",
		i,
		m->edges[d][i].triangles[0],
		m->edges[d][i].triangles[1]);
	  ec ++;
	}

	if (d < m->degree) {
	  if (m->edges[d][i].child_edges[0] < 0 ||
	      m->edges[d][i].child_edges[1] < 0) {
	    ERROR("edge %d missing child_edges (%d %d)",
		  i,
		  m->edges[d][i].child_edges[0],
		  m->edges[d][i].child_edges[1]);
	    ec ++;
	  }
	}
	  
	
      }
    }
  }

  return ec;
}

static int
manifold_all_parent_child_valid(manifold_t *m)
{
  int ec;
  int i;
  int j;

  ec = 0;
  for (i = 2; i < m->nvertices; i ++) {

    for (j = 0; j < 4; j ++) {

      if (m->vertices[i].children[j] >= 0) {

	if (m->vertices[m->vertices[i].children[j]].parent != i) {
	  ERROR("parent %d isn't referenced by child %d",
		i, m->vertices[i].children[j]);
	  ec ++;
	}
      }
    }
  }

  return ec;
}

