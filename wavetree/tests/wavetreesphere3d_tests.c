
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <check.h>

#include "manifold.h"
#include "icosahedron.h"
#include "wavetreesphere3d.h"

START_TEST (test_wavetreesphere3d_ncoefficients)
{
  static const int DEGREE = 1;
  wavetreesphere3d_t *s;
  manifold_t *manifold;

  ck_assert(DEGREE > 0);
  
  manifold = icosahedron_create(DEGREE - 1);
  ck_assert(manifold != NULL);
  
  s = wavetreesphere3d_create(manifold,
			      DEGREE, 
			      0.0);
  ck_assert(s != NULL);

  wavetreesphere3d_destroy(s);
  manifold_destroy(manifold);
}
END_TEST

START_TEST (test_wavetreesphere3d_childindices)
{
  static const int DEGREE = 3;

  manifold_t *manifold;
  wavetreesphere3d_t *s;
  int indices[12];
  int i;
  int j;
  int d;

  ck_assert(DEGREE > 0);

  manifold = icosahedron_create(DEGREE - 1);
  ck_assert(manifold != NULL);
  
  s = wavetreesphere3d_create(manifold,
			      DEGREE, 
			      0.0);
  ck_assert(s != NULL);

  ck_assert(wavetreesphere3d_get_child_indices(s, 0, indices, 12) == 12);
  ck_assert(indices[0] == 1);
  ck_assert(wavetreesphere3d_parent_index(s, indices[0]) == 0);
  
  ck_assert(indices[1] == 2);
  ck_assert(wavetreesphere3d_parent_index(s, indices[1]) == 0);

  ck_assert(indices[11] == 12);
  ck_assert(wavetreesphere3d_parent_index(s, indices[11]) == 0);

  /* Polar points -> 1 child */
  ck_assert(wavetreesphere3d_get_child_indices(s, 1, indices, 12) == 1);
  ck_assert(indices[0] == (manifold->nvertices + 1));
  ck_assert(wavetreesphere3d_parent_index(s, indices[0]) == 1);
  ck_assert(wavetreesphere3d_depthofindex(s, indices[0]) == 2);
  
  ck_assert(wavetreesphere3d_get_child_indices(s, 2, indices, 12) == 1);
  ck_assert(indices[0] == (manifold->nvertices + 2));
  ck_assert(wavetreesphere3d_parent_index(s, indices[0]) == 2);
  ck_assert(wavetreesphere3d_depthofindex(s, indices[0]) == 2);

  /* Non Polar (3 .. 12) -> 7 Children */
  for (j = 3; j < 13; j ++) {
    ck_assert(wavetreesphere3d_get_child_indices(s, j, indices, 12) == 7);
    ck_assert(wavetreesphere3d_depthofindex(s, j) == 1);
    for (i = 0; i < 7; i ++) {
      ck_assert(indices[i] > 0);
      ck_assert(wavetreesphere3d_parent_index(s, indices[i] == j));
      ck_assert(wavetreesphere3d_depthofindex(s, indices[i]) == 2);
    }
  }
  
  wavetreesphere3d_destroy(s);
  manifold_destroy(manifold);
}
END_TEST

START_TEST (test_wavetreesphere3d_2dindices)
{
  static const int DEGREE = 3;

  manifold_t *manifold;
  wavetreesphere3d_t *s;
  int i;
  int j;
  int k;
  int ii;
  int ij;
  int ik;

  int index;
  int width;

  ck_assert(DEGREE > 0);

  manifold = icosahedron_create(DEGREE - 1);
  ck_assert(manifold != NULL);
  
  s = wavetreesphere3d_create(manifold,
			      DEGREE, 
			      0.0);
  ck_assert(s != NULL);

  wavetreesphere3d_destroy(s);
  manifold_destroy(manifold);
}
END_TEST

START_TEST (test_wavetreesphere3d_parenting)
{
  static const int DEGREE = 3;

  manifold_t *manifold;
  wavetreesphere3d_t *s;
  int ncoeff;
  
  int i;
  int j;
  int k;
  int ii;
  int ij;
  int ik;
  int parent;
  int nchildren;
  int child_indices[15];
  int found;

  int index;
  int width;

  ck_assert(DEGREE > 0);

  manifold = icosahedron_create(DEGREE - 1);
  ck_assert(manifold != NULL);
  
  s = wavetreesphere3d_create(manifold,
			      DEGREE, 
			      0.0);
  ck_assert(s != NULL);

  ncoeff = wavetreesphere3d_get_coeff_size(s);
  ck_assert(ncoeff == 649);

  for (i = 1; i < ncoeff; i ++) {
    parent = wavetreesphere3d_parent_index(s, i);
    ck_assert(parent >= 0);

    nchildren = wavetreesphere3d_get_child_indices(s, parent, child_indices, 15);
    ck_assert(nchildren > 0);

    found = 0;
    for (j = 0; j < nchildren; j ++) {
      if (child_indices[j] == i) {
	found = -1;
	break;
      }
    }

    ck_assert(found);
  }

  wavetreesphere3d_destroy(s);
  manifold_destroy(manifold);
}
END_TEST

START_TEST (test_wavetreesphere3d_depth)
{
  static const int DEGREE = 3;

  manifold_t *manifold;
  wavetreesphere3d_t *s;
  int ncoeff;
  
  int i;
  int j;
  int k;
  int ii;
  int ij;
  int ik;
  int parent;
  int nchildren;
  int child_indices[15];
  int found;

  int index;
  int width;

  ck_assert(DEGREE > 0);

  manifold = icosahedron_create(DEGREE - 1);
  ck_assert(manifold != NULL);
  
  s = wavetreesphere3d_create(manifold,
			      DEGREE, 
			      0.0);
  ck_assert(s != NULL);

  ncoeff = wavetreesphere3d_get_coeff_size(s);
  ck_assert(ncoeff == 649);

  ck_assert(wavetreesphere3d_depthofindex(s, 0) == 0);

  for (i = 1; i < 13; i ++) {
    ck_assert(wavetreesphere3d_depthofindex(s, i) == 1);
  }

  for (i = 13; i < 43; i ++) {
    ck_assert(wavetreesphere3d_depthofindex(s, i) == 2);
  }
  for (i = 1; i < 43; i ++) {
    ck_assert(wavetreesphere3d_depthofindex(s, i + 162) == 2);
  }

  for (i = 43; i < 163; i ++) {
    ck_assert(wavetreesphere3d_depthofindex(s, i) == 3);
    ck_assert(wavetreesphere3d_depthofindex(s, i + 162) == 3);
  }
  
  for (i = 1; i < 163; i ++) {
    ck_assert(wavetreesphere3d_depthofindex(s, i + 2*162) == 3);
    ck_assert(wavetreesphere3d_depthofindex(s, i + 3*162) == 3);
  }

  wavetreesphere3d_destroy(s);
  manifold_destroy(manifold);
}
END_TEST

START_TEST(test_wavetreesphere3d_birth)
{
  static const int DEGREE = 3;
  static const double DUMMY = 3.14;
  
  manifold_t *manifold;
  wavetreesphere3d_t *s;

  double *a;
  int coeffsize;
  int size;
  int i;
  int indices[15];
  int nindices;
  int ii, ik;
  
  ck_assert(DEGREE > 0);
  
  manifold = icosahedron_create(DEGREE - 1);
  ck_assert(manifold != NULL);
  
  s = wavetreesphere3d_create(manifold,
			      DEGREE, 
			      0.0);
  ck_assert(s != NULL);

  size = wavetreesphere3d_get_image_size(s);
  ck_assert(size > 0);

  coeffsize = wavetreesphere3d_get_coeff_size(s);
  ck_assert(coeffsize == 649);

  a = malloc(sizeof(double) * size);
  ck_assert(a != NULL);

  /*
   * Initialize
   */
  ck_assert(wavetreesphere3d_initialize(s, 1.0) >= 0);
  ck_assert(wavetreesphere3d_coeff_count(s) == 1);

  ck_assert(wavetreesphere3d_attachable_branches(s) == 12);
  ck_assert(wavetreesphere3d_prunable_leaves(s) == 0);

  /* Set array to non-zero value */
  for (i = 0; i < size; i ++) {
    a[i] = DUMMY;
  }
  ck_assert(wavetreesphere3d_map_to_array(s, a, size) >= 0);

  for (i = 0; i < 12; i ++) {
    ck_assert(a[i] == 1.0);
  }
  for (i = 12; i < size; i ++) {
    ck_assert(a[i] == DUMMY);
  }

  /*
   * Birth a node (north pole)
   */
  ck_assert(wavetreesphere3d_propose_birth(s, 
					   1,
					   1,
					   2.0) >= 0);
  ck_assert(wavetreesphere3d_commit(s) >= 0);
  ck_assert(wavetreesphere3d_coeff_count(s) == 2);
  ck_assert(wavetreesphere3d_attachable_branches(s) == 12);
  ck_assert(wavetreesphere3d_prunable_leaves(s) == 1);
  

  /*
   * Birth a node (first non-polar point)
   */
  ck_assert(wavetreesphere3d_propose_birth(s, 
					   3,
					   1,
					   3.0) >= 0);
  ck_assert(wavetreesphere3d_commit(s) >= 0);
  ck_assert(wavetreesphere3d_coeff_count(s) == 3);
  /* ck_assert(wavetreesphere3d_attachable_branches(s) == 18); */
  ck_assert(wavetreesphere3d_prunable_leaves(s) == 2);

  ck_assert(wavetreesphere3d_sphere2dindices(s, 3, &ii, &ik) >= 0);
  printf("ii %d ik %d\n", ii, ik);

  nindices = wavetreesphere3d_get_child_indices(s,
						3,
						indices,
						15);
  for (i = 0; i < nindices; i ++) {
    printf("Child 3: %d: %d\n", i, indices[i]);
  }
  

  /*
   * Birth a node (deeper first non-polar point)
   */
  ck_assert(wavetreesphere3d_propose_birth(s, 
					   165,
					   2,
					   6.0) >= 0);
  ck_assert(wavetreesphere3d_commit(s) >= 0);
  ck_assert(wavetreesphere3d_coeff_count(s) == 4);
  /* ck_assert(wavetreesphere3d_attachable_branches(s) == 18); */
  ck_assert(wavetreesphere3d_prunable_leaves(s) == 2);
  
  nindices = wavetreesphere3d_get_child_indices(s,
						165,
						indices,
						15);
  for (i = 0; i < nindices; i ++) {
    printf("Child 165: %d: %d\n", i, indices[i]);
  }
  
  nindices = wavetreesphere3d_get_child_indices(s,
						327,
						indices,
						15);
  printf("327: %d children (parent %d)\n", nindices, wavetreesphere3d_parent_index(s, 327));
  for (i = 0; i < nindices; i ++) {
    printf("Child 327: %d: %d\n", i, indices[i]);
  }

  nindices = wavetreesphere3d_get_child_indices(s,
						489,
						indices,
						15);

  wavetreesphere3d_destroy(s);
  manifold_destroy(manifold);

  free(a);
}
END_TEST

START_TEST(test_wavetreesphere3d_death)
{
}
END_TEST

START_TEST(test_wavetreesphere3d_value)
{
}
END_TEST

START_TEST(test_wavetreesphere3d_image_mapping)
{
}
END_TEST

START_TEST(test_wavetreesphere3d_saveload)
{
}
END_TEST

Suite *
wavetreesphere3d_suite (void)
{
  Suite *s = suite_create ("Wave Tree 3D");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_wavetreesphere3d_ncoefficients); 
  tcase_add_test (tc_core, test_wavetreesphere3d_childindices); 
  tcase_add_test (tc_core, test_wavetreesphere3d_2dindices); 
  tcase_add_test (tc_core, test_wavetreesphere3d_parenting);
  tcase_add_test (tc_core, test_wavetreesphere3d_depth);
  
  /* tcase_add_test (tc_core, test_wavetreesphere3d_birth); */

  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = wavetreesphere3d_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_fork_status (sr, CK_NOFORK);

  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
