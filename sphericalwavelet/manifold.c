
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>

#include "manifold.h"

#include "slog.h"


manifold_t *
manifold_create(int degree,
		vertex_count_t nverticesatdepth,
		edge_count_t nedgesatdepth,
		triangle_count_t ntrianglesatdepth)
{
  manifold_t *m;
  int i;

  m = malloc(sizeof(manifold_t));
  if (m == NULL) {
    return NULL;
  }

  m->degree = degree;
  m->nverticesatdepth = nverticesatdepth;
  m->nedgesatdepth = nedgesatdepth;
  m->ntrianglesatdepth = ntrianglesatdepth;

  m->nvertices = m->nverticesatdepth(m->degree);
  m->vertices = malloc(sizeof(vertex3_t) * m->nvertices);
  if (m->vertices == NULL) {
    ERROR("failed to allocate vertices");
    return NULL;
  }

  m->nedges = malloc(sizeof(int) * (degree + 1));
  if (m->nedges == NULL) {
    ERROR("failed to allocate nedges");
    return NULL;
  }

  m->ntriangles = malloc(sizeof(int) * (degree + 1));
  if (m->ntriangles == NULL) {
    ERROR("failed to allocate ntriangles");
    return NULL;
  }

  m->edges = malloc(sizeof(edge_t*) * (degree + 1));
  if (m->edges == NULL) {
    ERROR("failed to allocate edges");
    return NULL;
  }

  m->triangles = malloc(sizeof(edge_t*) * (degree + 1));
  if (m->triangles == NULL) {
    ERROR("failed to allocate triangles");
    return NULL;
  }

  m->ntotaltriangles = 0;
  for (i = 0; i <= degree; i ++) {
    m->nedges[i] = m->nedgesatdepth(i);
    m->ntriangles[i] = m->ntrianglesatdepth(i);
    m->ntotaltriangles += m->ntriangles[i];
    
    m->edges[i] = malloc(sizeof(edge_t) * m->nedges[i]);
    if (m->edges[i] == NULL) {
      ERROR("failed to create edges for depth %d", i);
      return NULL;
    }
    memset(m->edges[i], -1, sizeof(edge_t) * m->nedges[i]);

    m->triangles[i] = malloc(sizeof(triangle_t) * m->ntriangles[i]);
    if (m->triangles[i] == NULL) {
      ERROR("failed to create triangles for depth %d", i);
      return NULL;
    }
    memset(m->triangles[i], -1, sizeof(triangle_t) * m->ntriangles[i]);
  }

  return m;
}

void
manifold_destroy(manifold_t *m)
{
  int i;

  if (m != NULL) {
    for (i = 0; i <= m->degree; i ++) {
      free(m->triangles[i]);
      free(m->edges[i]);
    }

    free(m->triangles);
    free(m->ntriangles);
    free(m->edges);
    free(m->nedges);

    free(m->vertices);

    free(m);
  }
}


int
manifold_subdivide(manifold_t *m,
		   int depth)
{
  vertex3_t *va, *vb;
  vertex3_t *vm;

  int vi;
  int ei;
  int ti;
  int i;

  int via, vib, vic;
  int eiab, eibc, eica;

  edge_t *eab, *ebc, *eca;

  int viab, vibc, vica;
  int iabbc, ibcca, icaab;

  int ordered_ei[6];

  if (depth < 1 || depth > m->degree) {
    return -1;
  }

  /*
   * Split all the edges
   */
  vi = m->nverticesatdepth(depth - 1);
  ei = 0;

  for (i = 0; i < m->nedges[depth - 1]; i ++, vi ++) {

    /*
     * Get references to the end points and midpoint vertex
     */
    va = &(m->vertices[m->edges[depth - 1][i].a]);
    vb = &(m->vertices[m->edges[depth - 1][i].b]);

    vm = &(m->vertices[vi]);
    vertex3_initialize(vm);
    
    /*
     * Find and normalize the midpoint
     */
    vertex3_midpoint(va, vb, vm);
    vertex3_normalize(vm);
    vm->depth = depth;

    /*
     * Set the "v" vertices (easy)
     */
    vm->v[0] = m->edges[depth - 1][i].a;
    vm->v[1] = m->edges[depth - 1][i].b;

    /*
     * Set the parent child relation ships
     */
    if (va->depth == (depth - 1)) {
      /* vm is child of va */
      if (vertex3_add_child(va, vi) < 0) {
	ERROR("failed to add child to parent vertex");
	return -1;
      }
      vm->parent = vm->v[0];
    } else {
      /* vm is child of vb */
      if (vertex3_add_child(vb, vi) < 0) {
	ERROR("failed to add child to parent vertex");
	return -1;
      }
      vm->parent = vm->v[1];
    }

    /*
     * Set the f vertices
     */
    if (manifold_get_edge_f_vertices(m, depth - 1, i, vm->f) < 0) {
      ERROR("failed to generate f vertices");
      return -1;
    }

    /*
     * Set the e vertices
     */
    if (manifold_get_edge_e_vertices(m, depth - 1, i, vm->e) < 0) {
      ERROR("failed to generate e vertices");
      return -1;
    }

    /*
     * Add the 2 child edges
     */

    m->edges[depth - 1][i].child_edges[0] = ei;

    m->edges[depth][ei].a = m->edges[depth - 1][i].a;
    m->edges[depth][ei].b = vi;
    m->edges[depth][ei].parent = i;
    m->edges[depth][ei].child_edges[0] = -1;
    m->edges[depth][ei].child_edges[1] = -1;
    ei ++;
    
    m->edges[depth - 1][i].child_edges[1] = ei;

    m->edges[depth][ei].a = vi;
    m->edges[depth][ei].b = m->edges[depth - 1][i].b;
    m->edges[depth][ei].parent = i;
    m->edges[depth][ei].child_edges[0] = -1;
    m->edges[depth][ei].child_edges[1] = -1;
    ei ++;
  }

  /*
   * All new vertices created and all external edges created.
   * Now subdivide the triangles, and add the new internal edges
   */
  ti = 0;
  for (i = 0; i < m->ntriangles[depth - 1]; i ++) {

    via = m->triangles[depth - 1][i].a;
    vib = m->triangles[depth - 1][i].b;
    vic = m->triangles[depth - 1][i].c;

    eiab = m->triangles[depth - 1][i].ab;
    eibc = m->triangles[depth - 1][i].bc;
    eica = m->triangles[depth - 1][i].ca;

    /*
     * Get parent edges, these must have child edges set.
     */
    eab = &(m->edges[depth - 1][eiab]);
    ebc = &(m->edges[depth - 1][eibc]);
    eca = &(m->edges[depth - 1][eica]);

    if (eab->a == via) {
      ordered_ei[0] = eab->child_edges[0];
      ordered_ei[1] = eab->child_edges[1];
    } else if (eab->b == via) {
      ordered_ei[0] = eab->child_edges[1];
      ordered_ei[1] = eab->child_edges[0];
    } else {
      ERROR("edge/triangle mismatch a: %d %d: %d\n",
	      eab->a, eab->b, via);
      ERROR("%d %d (%d %d %d)\n", depth - 1, i, via, vib, vic);
      return -1;
    }

    if (ebc->a == vib) {
      ordered_ei[2] = ebc->child_edges[0];
      ordered_ei[3] = ebc->child_edges[1];
    } else if (ebc->b == vib) {
      ordered_ei[2] = ebc->child_edges[1];
      ordered_ei[3] = ebc->child_edges[0];
    } else {
      ERROR("edge/triangle mismatch b: %d %d: %d\n",
	      ebc->a, ebc->b, vib);
      ERROR("%d (%d %d %d)\n", i, via, vib, vic);
      return -1;
    }

    if (eca->a == vic) {
      ordered_ei[4] = eca->child_edges[0];
      ordered_ei[5] = eca->child_edges[1];
    } else if (eca->b == vic) {
      ordered_ei[4] = eca->child_edges[1];
      ordered_ei[5] = eca->child_edges[0];
    } else {
      ERROR("edge/triangle mismatch c: %d %d: %d\n",
	      eca->a, eca->b, vic);
      ERROR("%d (%d %d %d)\n", i, via, vib, vic);
      return -1;
    }

    viab = m->edges[depth][eab->child_edges[0]].b;
    vibc = m->edges[depth][ebc->child_edges[0]].b;
    vica = m->edges[depth][eca->child_edges[0]].b;

    /*
     * Construct the tree internal edges
     */
    iabbc = ei;
    if (manifold_set_edge(m, depth, ei,
			  viab, vibc) < 0) {
      ERROR("failed to set edge iabbc");
      return -1;
    }
    ei ++;

    ibcca = ei;
    if (manifold_set_edge(m, depth, ei,
			  vibc, vica) < 0) {
      ERROR("failed to set edge ibcca");
      return -1;
    }
    ei ++;

    icaab = ei;
    if (manifold_set_edge(m, depth, ei,
			  vica, viab) < 0) {
      ERROR("failed to set edge icaab");
      return -1;
    }
    ei ++;

    /*
     * Now construct the 4 sub triangles
     */
    if (manifold_set_triangle(m, depth, ti,
			      via, viab, vica,
			      ordered_ei[0], icaab, ordered_ei[5]) < 0) {
      ERROR("failed to set triangle 0");
      return -1;
    }
    m->triangles[depth][ti].parent = i;
    m->triangles[depth - 1][i].child_triangles[0] = ti;
    ti ++;

    if (manifold_set_triangle(m, depth, ti,
			      viab, vib, vibc,
			      ordered_ei[1], ordered_ei[2], iabbc) < 0) {
      ERROR("failed to set triangle 1");
      return -1;
    }
    m->triangles[depth][ti].parent = i;
    m->triangles[depth - 1][i].child_triangles[1] = ti;
    ti ++;

    if (manifold_set_triangle(m, depth, ti,
			      vica, vibc, vic,
			      ibcca, ordered_ei[3], ordered_ei[4]) < 0) {
      ERROR("failed to set triangle 2");
      return -1;
    }
    m->triangles[depth][ti].parent = i;
    m->triangles[depth - 1][i].child_triangles[2] = ti;
    ti ++;

    if (manifold_set_triangle(m, depth, ti,
			      vica, viab, vibc,
			      icaab, iabbc, ibcca) < 0) {
      ERROR("failed to set triangle 3");
      return -1;
    }
    m->triangles[depth][ti].parent = i;
    m->triangles[depth - 1][i].child_triangles[3] = ti;
    ti ++;
  }

  if (ei != m->nedgesatdepth(depth)) {
    ERROR("incorrect number of edges generated %d != %d @ %d",
	    ei, m->nedgesatdepth(depth), depth);
    return -1;
  }
  
  if (ti != m->ntrianglesatdepth(depth)) {
    ERROR("incorrect number of triangles generated %d != %d @ %d",
	    ti, m->ntrianglesatdepth(depth), depth);
    return -1;
  }
    
  return 0;
}

int
manifold_build_neighbors(manifold_t *m)
{
  int j;
  int v0;
  int v1;
  
  for (j = 0; j < m->nedges[m->degree]; j ++) {

    v0 = m->edges[m->degree][j].a;
    v1 = m->edges[m->degree][j].b;
      
    if (vertex3_add_neighbor(&(m->vertices[v0]), v1) < 0) {
      return -1;
    }
      
    if (vertex3_add_neighbor(&(m->vertices[v1]), v0) < 0) {
      return -1;
    }

  }

  return 0;
}


int
manifold_compute_areas(manifold_t *m,
		       int depth)
{
  int i;
  double a;
  int vstart, vend;
  double total_area;
  
  if (depth == m->degree) {
    /*
     * Lowest level => approximate directly, it is assumed that the area has been initialized to 0.
     */

    total_area = 0.0;
    for (i = 0; i < m->ntriangles[depth]; i ++) {

      a = triangle_area(&(m->triangles[depth][i]),
			m->vertices);
      m->triangles[depth][i].area = a;
      
      total_area += a;

      m->vertices[m->triangles[depth][i].a].area += a/3.0;
      m->vertices[m->triangles[depth][i].b].area += a/3.0;
      m->vertices[m->triangles[depth][i].c].area += a/3.0;
      
    }

  } else {
    /* 
     * Higher levels accumulate
     */
    vstart = m->nverticesatdepth(depth);
    vend = m->nverticesatdepth(depth + 1);

    for (i = vstart; i < vend; i ++) {
      a = m->vertices[i].area;

      if (m->vertices[i].v[0] < 0 ||
	  m->vertices[i].v[1] < 0) {
	ERROR("v vertices unset");
	return -1;
      }
      m->vertices[m->vertices[i].v[0]].area += a/2.0;
      m->vertices[m->vertices[i].v[1]].area += a/2.0;
      
      if (m->vertices[i].f[0] < 0 ||
	  m->vertices[i].f[1] < 0) {
	ERROR("f vertices unset");
	return -1;
      }
      m->vertices[m->vertices[i].f[0]].area += a/4.0;
      m->vertices[m->vertices[i].f[1]].area += a/4.0;

      if (m->vertices[i].e[0] < 0 ||
	  m->vertices[i].e[1] < 0 ||
	  m->vertices[i].e[2] < 0 ||
	  m->vertices[i].e[3] < 0) {
	ERROR("e vertices unset");
	return -1;
      }
      m->vertices[m->vertices[i].e[0]].area -= a/16.0;
      m->vertices[m->vertices[i].e[1]].area -= a/16.0;
      m->vertices[m->vertices[i].e[2]].area -= a/16.0;
      m->vertices[m->vertices[i].e[3]].area -= a/16.0;
    }

    for (i = vstart; i < vend; i ++) {
      m->vertices[i].area =
	m->vertices[i].area/
	(m->vertices[m->vertices[i].v[0]].area + m->vertices[m->vertices[i].v[1]].area);
    }

    /*
     * Sum the child triangles into parent
     */
    vend = m->ntrianglesatdepth(depth);
    for (i = 0; i < vend; i ++) {
      m->triangles[depth][i].area =
	m->triangles[depth + 1][m->triangles[depth][i].child_triangles[0]].area + 
	m->triangles[depth + 1][m->triangles[depth][i].child_triangles[1]].area + 
	m->triangles[depth + 1][m->triangles[depth][i].child_triangles[2]].area + 
	m->triangles[depth + 1][m->triangles[depth][i].child_triangles[3]].area;
    }
  }
    
  return 0;
}


int
manifold_save_geo(manifold_t *m,
		  const char *filename)
{
  FILE *fp;
  int i;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    ERROR("failed to create file");
    return -1;
  }

  fprintf(fp, "PGEOMETRY V5\n");
  fprintf(fp, "NPoints %d NPrims %d\n", 
	  m->nvertices,
	  m->ntriangles[m->degree]);

  fprintf(fp, 
	  "NPointGroups 0 NPrimGroups 0\n"
	  "NPointAttrib 0 NVertexAttrib 0 NPrimAttrib 0 NAttrib 0\n");

  for (i = 0; i < m->nvertices; i ++) {

    fprintf(fp, "%f %f %f 1.0\n", 
	    m->vertices[i].x, m->vertices[i].y, m->vertices[i].z);
    
  }

  fprintf(fp, "Run %d Poly\n", m->ntriangles[m->degree]);

  for (i = 0; i < m->ntriangles[m->degree]; i ++) {
    
    fprintf(fp, " 3 < %d %d %d\n",
	    m->triangles[m->degree][i].a,
	    m->triangles[m->degree][i].b,
	    m->triangles[m->degree][i].c);

  }

  fprintf(fp, "beginExtra\n");
  fprintf(fp, "endExtra\n");

  fclose(fp);

  return 0;
}

int
manifold_get_child_vertices(manifold_t *m, int pi, int *ci)
{
  int i;
  
  if (m == NULL ||
      pi < 0 || pi >= m->nvertices) {
    ERROR("invalid parameters");
    return -1;
  }

  for (i = 0; i < 4; i ++) {
    if (m->vertices[pi].children[i] < 0) {
      break;
    }
    ci[i] = m->vertices[pi].children[i];
  }

  return i;
}

int
manifold_get_parent(manifold_t *m, int ci)
{
  if (m == NULL ||
      ci < 0 || ci >= m->nvertices) {
    ERROR("invalid parameters");
    return -1;
  }

  return m->vertices[ci].parent;
}

/*
 * Lookup functions
 */
int
manifold_find_enclosing_triangle(manifold_t *m,
				 triangle_workspace_t *workspace,
				 double lon, double lat,
				 int *ti,
				 double *ba, double *bb, double *bc)
{
  int t;
  int ct;
  double ta, tb, tc;
  double px, py, pz;

  int i;
  int depth;

  double epsilon;
  
  /*
   * Convert lon/lat to cartesian
   */
  vertex3_sphtocart(lon, lat,
		    &px, &py, &pz);

  /*
   * Linear search of the top level triangles
   */
  t = -1;
  for (i = 0; i < m->ntriangles[0]; i ++) {
    if (triangle_point_in_triangle(workspace,
				   &(m->triangles[0][i]),
				   m->vertices,
				   px, py, pz,
				   &ta, &tb, &tc,
				   DEFAULT_TRIANGLE_EPSILON)) {
      t = i;
      break;
    }
  }

  if (t < 0) {
    ERROR("failed to find enclosing triangle at top level");
    return -1;
  }

  /*
   * Now progress down the subdivision levels
   */
  depth = 0;
  while (depth < m->degree) {

    ct = -1;
    epsilon = DEFAULT_TRIANGLE_EPSILON;
    
    while (ct < 0) {
      for (i = 0; i < 4; i ++) {
	if (triangle_point_in_triangle(workspace,
				       &(m->triangles[depth + 1][m->triangles[depth][t].child_triangles[i]]),
				       m->vertices,
				       px, py, pz,
				       &ta, &tb, &tc,
				       epsilon)) {
	  ct = m->triangles[depth][t].child_triangles[i];
	  break;
	}
      }

      if (ct < 0) {
	epsilon *= 2.0;
      }
    }
    
    /* if (ct < 0) { */
    /*   ERROR("invalid child triangle (%d/%d)", depth, m->degree); */
    /*   for (i = 0; i < 4; i ++) { */
    /* 	ct = triangle_point_in_triangle(workspace, */
    /* 					&(m->triangles[depth + 1][m->triangles[depth][t].child_triangles[i]]), */
    /* 					m->vertices, */
    /* 					px, py, pz, */
    /* 					&ta, &tb, &tc); */
    /* 	ERROR("    %20.17f %20.17f %20.17f %d\n", ta, tb, tc, ct); */
    /*   } */
    /*   return -1; */
    /* } */

    t = ct;
    depth ++;
  }

  *ti = t;
  *ba = ta;
  *bb = tb;
  *bc = tc;

  return 0;
}

int
manifold_find_nearest_vertex(manifold_t *m,
			     triangle_workspace_t *workspace,
			     double lon, double lat,
			     int *vi)
{
  int ti;
  double ba, bb, bc;
  
  if (manifold_find_enclosing_triangle(m, workspace, lon, lat, &ti, &ba, &bb, &bc) < 0) {
    return -1;
  }

  if (ba > bb) {
    if (ba > bc) {
      *vi = m->triangles[m->degree][ti].a;
    } else {
      *vi = m->triangles[m->degree][ti].c;
    }
  } else {
    if (bb > bc) {
      *vi = m->triangles[m->degree][ti].b;
    } else {
      *vi = m->triangles[m->degree][ti].c;
    }
  }

  return 0;
}

int
manifold_compute_barycentre_coordinates(manifold_t *m,
					triangle_workspace_t *workspace,
					double lon,
					double lat,
					int *va, int *vb, int *vc,
					double *ba, double *bb, double *bc)
{
  int ti;
  
  if (manifold_find_enclosing_triangle(m, workspace, lon, lat, &ti, ba, bb, bc) < 0) {
    return -1;
  }

  *va = m->triangles[m->degree][ti].a;
  *vb = m->triangles[m->degree][ti].b;
  *vc = m->triangles[m->degree][ti].c;

  return 0;
}


int
manifold_set_vertex(manifold_t *m,
		    int vi,
		    double x, double y, double z,
		    int depth,
		    int normalize)
{
  vertex3_t *v;

  if (m == NULL ||
      vi < 0 || vi >= m->nvertices) {
    return -1;
  }

  v = m->vertices + vi;

  vertex3_initialize(v);

  v->x = x;
  v->y = y;
  v->z = z;

  v->depth = depth;

  if (normalize) {
    vertex3_normalize(v);
  }

  return 0;
}

int
manifold_set_edge(manifold_t *m,
		  int depth,
		  int ei,
		  int va, int vb)
{
  edge_t *e;

  if (m == NULL ||
      depth < 0 || depth > m->degree ||
      ei < 0 || ei >= m->nedges[depth]) {
    return -1;
  }

  e = m->edges[depth] + ei;

  edge_init(e);

  e->a = va;
  e->b = vb;

  return 0;
}

int
manifold_set_triangle(manifold_t *m,
		      int depth,
		      int ti,
		      int va, int vb, int vc,
		      int eab, int ebc, int eca)
{
  triangle_t *t;
  edge_t *e;
  
  if (m == NULL ||
      depth < 0 || depth > m->degree ||
      ti < 0 || ti >= m->ntriangles[depth]) {
    return -1;
  }

  t = m->triangles[depth] + ti;

  triangle_init(t);

  t->a = va;
  t->b = vb;
  t->c = vc;

  t->ab = eab;
  if (eab < 0 || eab >= m->nedges[depth]) {
    ERROR("edge ab out of range: %d (%d)", eab, m->nedges[depth]);
    return -1;
  }
  e = m->edges[depth] + eab;
  if (edge_add_triangle(e, ti) < 0) {
    ERROR("unable to add triangle to edge ab");
    return -1;
  }
  
  t->bc = ebc;
  if (ebc < 0 || ebc >= m->nedges[depth]) {
    ERROR("edge bc out of range: %d (%d)", ebc, m->nedges[depth]);
    return -1;
  }
  e = m->edges[depth] + ebc;
  if (edge_add_triangle(e, ti) < 0) {
    ERROR("unable to add triangle to edge bc");
    return -1;
  }

  t->ca = eca;
  if (eca < 0 || eca >= m->nedges[depth]) {
    ERROR("edge ca out of range: %d (%d)", eca, m->nedges[depth]);
    return -1;
  }
  e = m->edges[depth] + eca;
  if (edge_add_triangle(e, ti) < 0) {
    ERROR("unable to add triangle to edge ca");
    return -1;
  }

  return 0;
}

int
manifold_get_triangle_other_vertex(manifold_t *m,
				   int depth,
				   int ei,
				   int ti,
				   int *v)
{
  int va;
  int vb;
  triangle_t *t;
  
  va = m->edges[depth][ei].a;
  vb = m->edges[depth][ei].b;

  t = &(m->triangles[depth][ti]);
  if (t->a != va && t->a != vb) {
    *v = t->a;
    return 0;
  } else if (t->b != va && t->b != vb) {
    *v = t->b;
    return 0;
  } else if (t->c != va && t->c != vb) {
    *v = t->c;
    return 0;
  }

  ERROR("failed to find vertex");
  return -1;
}

int
manifold_get_edge_f_vertices(manifold_t *m,
			     int depth,
			     int ei,
			     int *f)
{
  int v;
  edge_t *e;

  e = &(m->edges[depth][ei]);
  
  if (e->triangles[0] < 0) {
    ERROR("triangle 0 missing");
    return -1;
  }
  if (manifold_get_triangle_other_vertex(m, depth, ei, e->triangles[0], &v) < 0) {
    ERROR("triangle 0 failed to get vertex");
    return -1;
  }
  f[0] = v;
  
  if (e->triangles[1] < 0) {
    ERROR("triangle 1 missing");
    return -1;
  }
  if (manifold_get_triangle_other_vertex(m, depth, ei, e->triangles[1], &v) < 0) {
    ERROR("triangle 1 failed to get vertex");
    return -1;
  }
  f[1] = v;

  /*
   * This can happen legitimately for the tetrahedron level 0
   */
  /* if (f[0] == f[1]) { */
  /*   ERROR("manifold_get_edge_f_vertices: duplicate vertices (triangles %d %d)\n", */
  /* 	    e->triangles[0], e->triangles[1]); */
  /*   return -1; */
  /* } */

  return 0;
}
			    
int
manifold_get_edge_e_vertices(manifold_t *m,
			     int depth,
			     int ei,
			     int *ev)
{
  int f[2];
  edge_t *e;
  triangle_t *t;

  int ea;
  int eb;

  e = &(m->edges[depth][ei]);
  
  /*
   * First triangle
   */
  if (e->triangles[0] < 0) {
    ERROR("triangle 0 missing");
    return -1;
  }

  t = &(m->triangles[depth][e->triangles[0]]);
  if (t->ab == ei) {
    ea = t->bc;
    eb = t->ca;
  } else if (t->bc == ei) {
    ea = t->ab;
    eb = t->ca;
  } else if (t->ca == ei) {
    ea = t->ab;
    eb = t->ca;
  } else {
    ERROR("triangle/edge inconsistency in triangle 0");
    return -1;
  }

  if (manifold_get_edge_f_vertices(m, depth, ea, f) < 0) {
    ERROR("failed to get f vertices 0, ea");
    return -1;
  }

  if (f[0] != e->a && f[0] != e->b) {
    ev[0] = f[0];
  } else if (f[1] != e->a && f[1] != e->b) {
    ev[0] = f[1];
  } else if (f[0] == f[1]) {
    /* This can happen for tetrahedron */
    ev[0] = f[0];
  } else {
    ERROR("failed to get value f vertices 0, ea");
    ERROR("edge %d %d: f %d %d: triangles %d %d",
	    e->a, e->b,
	    f[0], f[1],
	    e->triangles[0], e->triangles[1]);
    return -1;
  }

  if (manifold_get_edge_f_vertices(m, depth, eb, f) < 0) {
    ERROR("failed to get f vertices 0, eb");
    return -1;
  }

  if (f[0] != e->a && f[0] != e->b) {
    ev[1] = f[0];
  } else if (f[1] != e->a && f[1] != e->b) {
    ev[1] = f[1];
  } else if (f[0] == f[1]) {
    ev[1] = f[0];
  } else {
    ERROR("failed to get value f vertices 0, eb");
    ERROR("edge %d %d: f %d %d: triangles %d %d",
	    e->a, e->b,
	    f[0], f[1],
	    e->triangles[0], e->triangles[1]);
    return -1;
  }

  /*
   * 2nd Triangle 
   */
  if (e->triangles[1] < 0) {
    ERROR("triangle 1 missing");
    return -1;
  }

  t = &(m->triangles[depth][e->triangles[1]]);
  if (t->ab == ei) {
    ea = t->bc;
    eb = t->ca;
  } else if (t->bc == ei) {
    ea = t->ab;
    eb = t->ca;
  } else if (t->ca == ei) {
    ea = t->ab;
    eb = t->ca;
  } else {
    ERROR("triangle/edge inconsistency in triangle 1");
    return -1;
  }

  if (manifold_get_edge_f_vertices(m, depth, ea, f) < 0) {
    ERROR("failed to get f vertices 1, ea");
    return -1;
  }

  if (f[0] != e->a && f[0] != e->b) {
    ev[2] = f[0];
  } else if (f[1] != e->a && f[1] != e->b) {
    ev[2] = f[1];
  } else if (f[0] == f[1]) {
    ev[2] = f[0];
  } else {
    ERROR("failed to get value f vertices 1, ea");
    ERROR("edge %d %d: f %d %d: triangles %d %d",
	    e->a, e->b,
	    f[0], f[1],
	    e->triangles[0], e->triangles[1]);
    return -1;
  }

  if (manifold_get_edge_f_vertices(m, depth, eb, f) < 0) {
    ERROR("failed to get f vertices 1, eb");
    return -1;
  }

  if (f[0] != e->a && f[0] != e->b) {
    ev[3] = f[0];
  } else if (f[1] != e->a && f[1] != e->b) {
    ev[3] = f[1];
  } else if (f[0] == f[1]) {
    ev[3] = f[0];
  } else {
    ERROR("failed to get value f vertices 1, eb");
    ERROR("edge %d %d: f %d %d: triangles %d %d",
	  e->a, e->b,
	  f[0], f[1],
	  e->triangles[0], e->triangles[1]);
    return -1;
  }

  return 0;
}

double
manifold_area_for_depth(manifold_t *m,
			int depth)
{
  double a;
  int i;
  int vstart;
  int vend;
  
  if (depth == 0) {
    vstart = 0;
  } else {
    vstart = m->nverticesatdepth(depth - 1);
  }
  
  vend = m->nverticesatdepth(depth);

  a = 0.0;
  for (i = vstart; i < vend; i ++) {
    a += m->vertices[i].area;
  }

  return a;
}
