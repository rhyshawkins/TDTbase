
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <check.h>

#include "manifold.h"
#include "icosahedron.h"
#include "wavetreesphereface3d.h"

START_TEST (test_wavetreesphereface3d_create)
{
  static const int DEGREE = 1;
  wavetreesphereface3d_t *s;
  manifold_t *manifold;

  ck_assert(DEGREE > 0);
  
  manifold = icosahedron_create(DEGREE - 1);
  ck_assert(manifold != NULL);
  
  s = wavetreesphereface3d_create(manifold,
				  DEGREE, 
				  0.0);
  ck_assert(s != NULL);

  wavetreesphereface3d_destroy(s);
  manifold_destroy(manifold);
}
END_TEST

START_TEST (test_wavetreesphereface3d_sphereindices)
{
  static const int DEGREE = 3;
  wavetreesphereface3d_t *s;
  manifold_t *manifold;

  int t, r, d;
  int i;
  int rowstride;

  ck_assert(DEGREE > 0);
  
  manifold = icosahedron_create(DEGREE - 1);
  ck_assert(manifold != NULL);

  rowstride = manifold->ntotaltriangles;
  
  s = wavetreesphereface3d_create(manifold,
				  DEGREE, 
				  0.0);
  ck_assert(s != NULL);

  /*
   * Test root index
   */
  ck_assert(wavetreesphereface3d_index_to_sphereindices(s, 0, &t, &r, &d) == 0);
  ck_assert(t == 0);
  ck_assert(r == 0);
  ck_assert(d == 0);
  ck_assert(wavetreesphereface3d_sphereindices_to_index(s, t, r, d) == 0);

  /*
   * Test indices of degree 1 indices (1 .. 20)
   */
  for (i = 1; i <= 20; i ++) {
    ck_assert(wavetreesphereface3d_index_to_sphereindices(s, i, &t, &r, &d) == 0);

    /* printf("%2d: %2d %2d %2d (%2d)\n", */
    /* 	   i, t, r, d, */
    /* 	   wavetreesphereface3d_sphereindices_to_index(s, t, r, d)); */
    
    ck_assert(t == (i - 1));
    ck_assert(r == 0);
    ck_assert(d == 1);

    ck_assert(wavetreesphereface3d_sphereindices_to_index(s, t, r, d) == i);
  }

  /*
   * Test indices of degree 2 indices (21 .. 160)
   */

  for (i = 21; i <= 80; i ++) {
    /* Lateral Indices */
    ck_assert(wavetreesphereface3d_index_to_sphereindices(s, i, &t, &r, &d) == 0);

    /* printf("%2d: %2d %2d %2d (%2d)\n", */
    /* 	   i, t, r, d, */
    /* 	   wavetreesphereface3d_sphereindices_to_index(s, t, r, d)); */

    ck_assert(t == (i - 1));
    ck_assert(r == 0);
    ck_assert(d == 2);

    ck_assert(wavetreesphereface3d_sphereindices_to_index(s, t, r, d) == i);
  }

  for (i = 0; i < 80; i ++) {
    /* Radial/Lateral Indices */
    ck_assert(wavetreesphereface3d_index_to_sphereindices(s, i + 1 + rowstride, &t, &r, &d) == 0);

    /* printf("%2d: %2d %2d %2d (%2d)\n", */
    /* 	   i + 1 + rowstride, t, r, d, */
    /* 	   wavetreesphereface3d_sphereindices_to_index(s, t, r, d)); */

    ck_assert(t == i);
    ck_assert(r == 1);
    ck_assert(d == 2);

    ck_assert(wavetreesphereface3d_sphereindices_to_index(s, t, r, d) == i + 1 + rowstride);
  }
  
  wavetreesphereface3d_destroy(s);
  manifold_destroy(manifold);
}
END_TEST

START_TEST (test_wavetreesphereface3d_depth)
{
  static const int DEGREE = 3;

  manifold_t *manifold;
  wavetreesphereface3d_t *s;
  int i;
  int rowstride;

  ck_assert(DEGREE > 0);

  manifold = icosahedron_create(DEGREE - 1);
  ck_assert(manifold != NULL);
  
  s = wavetreesphereface3d_create(manifold,
				  DEGREE, 
				  0.0);
  ck_assert(s != NULL);

  ck_assert(wavetreesphereface3d_depthofindex(s, 0) == 0);
  
  for (i = 1; i < 21; i ++) {
    ck_assert(wavetreesphereface3d_depthofindex(s, i) == 1);
  }

  rowstride = wavetreesphereface3d_get_rowstride(s);
  ck_assert(rowstride > 0);

  /*
   * Depth 2 
   */
  for (i = 21; i < 101; i ++) {
    ck_assert(wavetreesphereface3d_depthofindex(s, i) == 2);
  }

  for (i = 1; i < 101; i ++) {
    ck_assert(wavetreesphereface3d_depthofindex(s, i + rowstride) == 2);
  }

  /*
   * Depth 3
   */
  for (i = 101; i < 421; i ++) {
    ck_assert(wavetreesphereface3d_depthofindex(s, i) == 3);
    ck_assert(wavetreesphereface3d_depthofindex(s, i + rowstride) == 3);
  }
    
  for (i = 1; i < 421; i ++) {
    ck_assert(wavetreesphereface3d_depthofindex(s, i + 2*rowstride) == 3);
    ck_assert(wavetreesphereface3d_depthofindex(s, i + 3*rowstride) == 3);
  }
    
  wavetreesphereface3d_destroy(s);
  manifold_destroy(manifold);
}
END_TEST

START_TEST (test_wavetreesphereface3d_childindices)
{
  static const int DEGREE = 3;

  manifold_t *manifold;
  wavetreesphereface3d_t *s;
  int indices[20];
  int i;
  int j;

  int rowstride;

  ck_assert(DEGREE > 0);

  manifold = icosahedron_create(DEGREE - 1);
  ck_assert(manifold != NULL);
  
  s = wavetreesphereface3d_create(manifold,
				  DEGREE, 
				  0.0);
  ck_assert(s != NULL);

  rowstride = wavetreesphereface3d_get_rowstride(s);
  ck_assert(rowstride > 0);

  /*
   * Root should give 20 children
   */
  ck_assert(wavetreesphereface3d_get_child_indices(s, 0, indices, 20) == 20);
  for (i = 0; i < 20; i ++) {
    ck_assert(indices[i] == (i + 1));
    ck_assert(wavetreesphereface3d_parent_index(s, indices[i]) == 0);
  }

  /*
   * Next level should give 7 children each
   */
  for (i = 0; i < 20; i ++) {
    ck_assert(wavetreesphereface3d_get_child_indices(s, i + 1, indices, 20) == 7);
    for (j = 0; j < 7; j ++) {
      /* printf("parent index: %d %d -> %d\n", */
      /* 	     j, */
      /* 	     indices[j], */
      /* 	     wavetreesphereface3d_parent_index(s, indices[j])); */
      ck_assert(wavetreesphereface3d_parent_index(s, indices[j]) == (i + 1));
    }
  }

  /*
   * Depth 2 but refinement of the radial top level triangles: 8 children
   */
  for (i = 0; i < 20; i ++) {
    ck_assert(wavetreesphereface3d_get_child_indices(s, i + 1 + rowstride, indices, 20) == 8);
    for (j = 0; j < 8; j ++) {
      /* printf("parent index: %d %d -> %d (%d)\n", */
      /* 	     j, */
      /* 	     indices[j], */
      /* 	     wavetreesphereface3d_parent_index(s, indices[j]), */
      /* 	     i + 1 + rowstride); */
      ck_assert(wavetreesphereface3d_parent_index(s, indices[j]) == (i + 1 + rowstride));
    }
  }

  /* Depth 2: 6 Children
   */
  for (i = 20; i < 100; i ++) {
    ck_assert(wavetreesphereface3d_get_child_indices(s, i + 1, indices, 20) == 6);
    for (j = 0; j < 6; j ++) {
      /* printf("parent index: %d %d -> %d (%d)\n", */
      /* 	     j, */
      /* 	     indices[j], */
      /* 	     wavetreesphereface3d_parent_index(s, indices[j]), */
      /* 	     i + 1); */
      ck_assert(wavetreesphereface3d_parent_index(s, indices[j]) == (i + 1));
    }
  }
  
  
  wavetreesphereface3d_destroy(s);
  manifold_destroy(manifold);
}
END_TEST

START_TEST(test_wavetreesphereface3d_birth)
{
  static const int DEGREE = 3;
  static const double DUMMY = 3.14;
  
  manifold_t *manifold;
  wavetreesphereface3d_t *s;

  double *a;
  int coeffsize;
  int size;
  int i;
  int indices[20];
  int nindices;
  int it, ir, id;
  
  ck_assert(DEGREE > 0);
  
  manifold = icosahedron_create(DEGREE - 1);
  ck_assert(manifold != NULL);
  
  s = wavetreesphereface3d_create(manifold,
				  DEGREE, 
				  0.0);
  ck_assert(s != NULL);
  
  size = wavetreesphereface3d_get_image_size(s);
  ck_assert(size > 0);
  
  coeffsize = wavetreesphereface3d_get_coeff_size(s);
  printf("coeffsize: %d\n", coeffsize);
  ck_assert(coeffsize == 649);
  
  a = malloc(sizeof(double) * size);
  ck_assert(a != NULL);
  
  /*
   * Initialize
   */
  ck_assert(wavetreesphereface3d_initialize(s, 1.0) >= 0);
  ck_assert(wavetreesphereface3d_coeff_count(s) == 1);
  
  ck_assert(wavetreesphereface3d_attachable_branches(s) == 12);
  ck_assert(wavetreesphereface3d_prunable_leaves(s) == 0);
  
  /* Set array to non-zero value */
  for (i = 0; i < size; i ++) {
    a[i] = DUMMY;
  }
  ck_assert(wavetreesphereface3d_map_to_array(s, a, size) >= 0);

  for (i = 0; i < 12; i ++) {
    ck_assert(a[i] == 1.0);
  }
  for (i = 12; i < size; i ++) {
    ck_assert(a[i] == DUMMY);
  }
  
  /*
   * Birth a node (north pole)
   */
  ck_assert(wavetreesphereface3d_propose_birth(s, 
					       1,
					       1,
					       2.0) >= 0);
  ck_assert(wavetreesphereface3d_commit(s) >= 0);
  ck_assert(wavetreesphereface3d_coeff_count(s) == 2);
  ck_assert(wavetreesphereface3d_attachable_branches(s) == 12);
  ck_assert(wavetreesphereface3d_prunable_leaves(s) == 1);
  

  /*
   * Birth a node (first non-polar point)
   */
  ck_assert(wavetreesphereface3d_propose_birth(s, 
					       3,
					       1,
					       3.0) >= 0);
  ck_assert(wavetreesphereface3d_commit(s) >= 0);
  ck_assert(wavetreesphereface3d_coeff_count(s) == 3);
  /* ck_assert(wavetreesphereface3d_attachable_branches(s) == 18); */
  ck_assert(wavetreesphereface3d_prunable_leaves(s) == 2);
  
  ck_assert(wavetreesphereface3d_index_to_sphereindices(s, 3, &it, &ir, &id) >= 0);
  printf("ii %d ik %d\n", it, ir);
  
  nindices = wavetreesphereface3d_get_child_indices(s,
						    3,
						    indices,
						    20);
  for (i = 0; i < nindices; i ++) {
    printf("Child 3: %d: %d\n", i, indices[i]);
  }
  
  
  /*
   * Birth a node (deeper first non-polar point)
   */
  ck_assert(wavetreesphereface3d_propose_birth(s, 
					       165,
					       2,
					       6.0) >= 0);
  ck_assert(wavetreesphereface3d_commit(s) >= 0);
  ck_assert(wavetreesphereface3d_coeff_count(s) == 4);
  /* ck_assert(wavetreesphereface3d_attachable_branches(s) == 18); */
  ck_assert(wavetreesphereface3d_prunable_leaves(s) == 2);
  
  nindices = wavetreesphereface3d_get_child_indices(s,
						    165,
						    indices,
						    15);
  for (i = 0; i < nindices; i ++) {
    printf("Child 165: %d: %d\n", i, indices[i]);
  }
  
  nindices = wavetreesphereface3d_get_child_indices(s,
						    327,
						    indices,
						    15);
  printf("327: %d children (parent %d)\n", nindices, wavetreesphereface3d_parent_index(s, 327));
  for (i = 0; i < nindices; i ++) {
    printf("Child 327: %d: %d\n", i, indices[i]);
  }
  
  nindices = wavetreesphereface3d_get_child_indices(s,
						    489,
						    indices,
						    15);
  
  wavetreesphereface3d_destroy(s);
  manifold_destroy(manifold);

  free(a);
}
END_TEST

#if 0
START_TEST(test_wavetreesphereface3d_death)
{
}
END_TEST

START_TEST(test_wavetreesphereface3d_value)
{
}
END_TEST

START_TEST(test_wavetreesphereface3d_image_mapping)
{
}
END_TEST

START_TEST(test_wavetreesphereface3d_saveload)
{
}
END_TEST
#endif /* 0 */

Suite *
wavetreesphereface3d_suite (void)
{
  Suite *s = suite_create ("Wave Tree Face Wavelet 3D");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");

  tcase_add_test (tc_core, test_wavetreesphereface3d_create);
  tcase_add_test (tc_core, test_wavetreesphereface3d_sphereindices);
  tcase_add_test (tc_core, test_wavetreesphereface3d_depth);
  tcase_add_test (tc_core, test_wavetreesphereface3d_childindices);

  tcase_add_test (tc_core, test_wavetreesphereface3d_birth);

  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = wavetreesphereface3d_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
