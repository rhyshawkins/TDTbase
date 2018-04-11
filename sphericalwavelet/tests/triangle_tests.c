

#include <stdio.h>
#include <stdlib.h>
#include <check.h>

#include "octahedron.h"
#include "triangle.h"

START_TEST(test_triangle_point_in_triangle)
{
  triangle_workspace_t *workspace;
  triangle_t triangle;
  vertex3_t vertices[3];

  double px, py, pz;
  double ba, bc, bb;
  
  workspace = triangle_point_in_triangle_create_workspace();
  ck_assert(workspace != NULL);

  vertex3_initialize(&(vertices[0]));
  vertex3_sphtocart(0.0, 0.0,
		    &(vertices[0].x), &(vertices[0].y), &(vertices[0].z));
  vertex3_initialize(&(vertices[1]));
  vertex3_sphtocart(0.0, 90.0,
		    &(vertices[1].x), &(vertices[1].y), &(vertices[1].z));
  vertex3_initialize(&(vertices[2]));
  vertex3_sphtocart(90.0, 0.0,
		    &(vertices[2].x), &(vertices[2].y), &(vertices[2].z));

  triangle_init(&triangle);
  triangle.a = 0;
  triangle.b = 1;
  triangle.c = 2;

  /*
   * Outside test
   */
  vertex3_sphtocart(-1.0, 45.0,
		    &px, &py, &pz);
  ck_assert(triangle_point_in_triangle(workspace,
				       &triangle,
				       vertices,
				       px, py, pz,
				       &ba, &bb, &bc) == 0);

  /*
   * Central test
   */
  vertex3_sphtocart(45.0, 45.0,
		    &px, &py, &pz);
  ck_assert(triangle_point_in_triangle(workspace,
				       &triangle,
				       vertices,
				       px, py, pz,
				       &ba, &bb, &bc));

  /*
   * Vertex test
   */
  vertex3_sphtocart(0.0, 0.0,
		    &px, &py, &pz);
  ck_assert(triangle_point_in_triangle(workspace,
				       &triangle,
				       vertices,
				       px, py, pz,
				       &ba, &bb, &bc));
  ck_assert(ba == 1.0);
  
  /*
   * Vertex test 1
   */
  vertex3_sphtocart(0.0, 90.0,
		    &px, &py, &pz);
  ck_assert(triangle_point_in_triangle(workspace,
				       &triangle,
				       vertices,
				       px, py, pz,
				       &ba, &bb, &bc));
  ck_assert(bb == 1.0);

  /*
   * Vertex test 2
   */
  vertex3_sphtocart(90.0, 0.0,
		    &px, &py, &pz);
  ck_assert(triangle_point_in_triangle(workspace,
				       &triangle,
				       vertices,
				       px, py, pz,
				       &ba, &bb, &bc));
  ck_assert(bc == 1.0);

  /*
   * Edge test 0
   */
  vertex3_sphtocart(0.0, 45.0,
		    &px, &py, &pz);
  ck_assert(triangle_point_in_triangle(workspace,
				       &triangle,
				       vertices,
				       px, py, pz,
				       &ba, &bb, &bc));
  ck_assert(bc == 0.0);
  
  /*
   * Edge test 1
   */
  vertex3_sphtocart(90.0, 45.0,
		    &px, &py, &pz);
  triangle_point_in_triangle(workspace,
					    &triangle,
					    vertices,
					    px, py, pz,
					    &ba, &bb, &bc);
  ck_assert(ba == 0.0);
  
  /*
   * Edge test 2
   */
  vertex3_sphtocart(45.0, 0.0,
		    &px, &py, &pz);
  ck_assert(triangle_point_in_triangle(workspace,
				       &triangle,
				       vertices,
				       px, py, pz,
				       &ba, &bb, &bc));
  ck_assert(bb == 0.0);

  triangle_point_in_triangle_free_workspace(workspace);
}
END_TEST

START_TEST(test_triangle_in_octahedron)
{
  manifold_t *o;
  triangle_workspace_t *workspace;
  int vi;
  
  o = octahedron_create(2);
  ck_assert(o != NULL);

  workspace = triangle_point_in_triangle_create_workspace();
  ck_assert(workspace != NULL);

  /* North pole */
  ck_assert(manifold_find_nearest_vertex(o, workspace,
					 45.0, 89.0,
					 &vi) == 0);
  ck_assert(vi == 0);

  ck_assert(manifold_find_nearest_vertex(o, workspace,
					 0.0, 90.0,
					 &vi) == 0);
  ck_assert(vi == 0);

  ck_assert(manifold_find_nearest_vertex(o, workspace,
					 -90.0, 89.0,
					 &vi) == 0);
  ck_assert(vi == 0);

  /* South pole */
  ck_assert(manifold_find_nearest_vertex(o, workspace,
					 45.0, -89.0,
					 &vi) == 0);
  ck_assert(vi == 1);

  ck_assert(manifold_find_nearest_vertex(o, workspace,
					 0.0, -90.0,
					 &vi) == 0);
  ck_assert(vi == 1);

  ck_assert(manifold_find_nearest_vertex(o, workspace,
					 -90.0, -89.0,
					 &vi) == 0);
  ck_assert(vi == 1);

  /* On/Near Equator Points */
  ck_assert(manifold_find_nearest_vertex(o, workspace,
					 0.0, 0.0,
					 &vi) == 0);
  ck_assert(vi == 2);
  
  ck_assert(manifold_find_nearest_vertex(o, workspace,
					 89.0, 1.0,
					 &vi) == 0);
  ck_assert(vi == 3);
  
  ck_assert(manifold_find_nearest_vertex(o, workspace,
					 179.0, -1.0,
					 &vi) == 0);
  ck_assert(vi == 4);
  
  ck_assert(manifold_find_nearest_vertex(o, workspace,
					 271.0, 1.0,
					 &vi) == 0);
  ck_assert(vi == 5);

  manifold_destroy(o);
}
END_TEST

START_TEST(test_triangle_in_icosahedron)
{
  manifold_t *o;
  triangle_workspace_t *workspace;
  int vi;
  int ti;
  double ba, bb, bc;
  
  o = octahedron_create(2);
  ck_assert(o != NULL);

  workspace = triangle_point_in_triangle_create_workspace();
  ck_assert(workspace != NULL);

  /* North pole */
  ck_assert(manifold_find_nearest_vertex(o, workspace,
					 45.0, 89.0,
					 &vi) == 0);
  ck_assert(vi == 0);

  ck_assert(manifold_find_nearest_vertex(o, workspace,
					 0.0, 90.0,
					 &vi) == 0);
  ck_assert(vi == 0);

  ck_assert(manifold_find_nearest_vertex(o, workspace,
					 -90.0, 89.0,
					 &vi) == 0);
  ck_assert(vi == 0);

  /* South pole */
  ck_assert(manifold_find_nearest_vertex(o, workspace,
					 45.0, -89.0,
					 &vi) == 0);
  ck_assert(vi == 1);

  ck_assert(manifold_find_nearest_vertex(o, workspace,
					 0.0, -90.0,
					 &vi) == 0);
  ck_assert(vi == 1);

  ck_assert(manifold_find_nearest_vertex(o, workspace,
					 -90.0, -89.0,
					 &vi) == 0);
  ck_assert(vi == 1);

  /* On/Near Equator Points */
  ck_assert(manifold_find_nearest_vertex(o, workspace,
					 0.0, 0.0,
					 &vi) == 0);
  ck_assert(vi == 2);
  
  ck_assert(manifold_find_nearest_vertex(o, workspace,
					 89.0, 1.0,
					 &vi) == 0);
  ck_assert(vi == 3);
  ck_assert(manifold_find_nearest_vertex(o, workspace,
					 90.0, 0.0,
					 &vi) == 0);
  ck_assert(vi == 3);
  ck_assert(manifold_find_enclosing_triangle(o, workspace,
					     90.0, 0.0,
					     &ti,
					     &ba, &bb, &bc) == 0);
	    
  
  
  ck_assert(manifold_find_nearest_vertex(o, workspace,
					 179.0, -1.0,
					 &vi) == 0);
  ck_assert(vi == 4);
  
  ck_assert(manifold_find_nearest_vertex(o, workspace,
					 271.0, 1.0,
					 &vi) == 0);
  ck_assert(vi == 5);

  manifold_destroy(o);
}
END_TEST

Suite *
triangle_suite (void)
{
  Suite *s = suite_create ("Triangle");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");

  tcase_add_test (tc_core, test_triangle_point_in_triangle);
  
  tcase_add_test (tc_core, test_triangle_in_octahedron);
  tcase_add_test (tc_core, test_triangle_in_icosahedron);

  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = triangle_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
