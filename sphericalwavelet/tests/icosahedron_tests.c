
#include <stdio.h>
#include <stdlib.h>
#include <check.h>

#include <math.h>

#include "icosahedron.h"

START_TEST(test_icosahedron_counting)
{
  ck_assert(icosahedron_nvertices(0) == 12);
  ck_assert(icosahedron_nedges(0) == 30);
  ck_assert(icosahedron_ntriangles(0) == 20);

  ck_assert(icosahedron_nvertices(1) == 42);
  ck_assert(icosahedron_nedges(1) == 120);
  ck_assert(icosahedron_ntriangles(1) == 80);

  ck_assert(icosahedron_nvertices(2) == 162);
  ck_assert(icosahedron_nedges(2) == 480);
  ck_assert(icosahedron_ntriangles(2) == 320);

  ck_assert(icosahedron_nvertices(3) == 642);
  ck_assert(icosahedron_nedges(3) == 1920);
  ck_assert(icosahedron_ntriangles(3) == 1280);
}
END_TEST

START_TEST(test_icosahedron_create)
{
  manifold_t *o;
  char filename[256];
  int i;

  for (i = 0; i < 5; i ++) {
    
  
    o = icosahedron_create(i);

    ck_assert(o != NULL);
    
    ck_assert(manifold_valid(o));

    sprintf(filename, "icosahedron_%d.geo", i);
    ck_assert(manifold_save_geo(o, filename) >= 0);

    manifold_destroy(o);
  }

}
END_TEST

START_TEST(test_icosahedron_edges)
{
  manifold_t *o;

  o = icosahedron_create(2);

  ck_assert(o != NULL);
  
  ck_assert(manifold_valid(o));

  /* printf("d1: edge0: %d %d\n", o->edges[1][0].a, o->edges[1][0].b); */
  /* printf("d1: edge1: %d %d\n", o->edges[1][1].a, o->edges[1][1].b); */
  /* printf("d1: edge2: %d %d\n", o->edges[1][2].a, o->edges[1][2].b); */
  /* printf("d1: edge3: %d %d\n", o->edges[1][3].a, o->edges[1][3].b); */
  /* printf("d1: edge4: %d %d\n", o->edges[1][4].a, o->edges[1][4].b); */
  /* printf("d1: edge5: %d %d\n", o->edges[1][5].a, o->edges[1][5].b); */


  manifold_destroy(o);
}
END_TEST

START_TEST(test_icosahedron_neighbors)
{
  manifold_t *o;
  int i;
  int nc;
  int vi;

  int nv5;
  
  o = icosahedron_create(3);

  ck_assert(o != NULL);
  
  ck_assert(manifold_valid(o));

  nv5 = 0;
  
  for (vi = 0; vi < o->nvertices; vi ++) {

    nc = 0;
    for (i = 0; i < 6; i ++) {
      if (o->vertices[vi].n[i] >= 0) {
	nc ++;
      } else {
	break;
      }
    }
    ck_assert(nc >= 5);

    if (nc == 5) {
      nv5 ++;
    }
  }

  /* Special case for icosahedron : there should be 12 vertices with only 5 neighbors */
  ck_assert(nv5 == 12);

  manifold_destroy(o);
}
END_TEST

Suite *
icosahedron_suite (void)
{
  Suite *s = suite_create ("Icosahedron");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");

  tcase_add_test (tc_core, test_icosahedron_counting);
  tcase_add_test (tc_core, test_icosahedron_create);

  tcase_add_test (tc_core, test_icosahedron_edges);
  tcase_add_test (tc_core, test_icosahedron_neighbors);
  
  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = icosahedron_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
