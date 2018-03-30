
#include <stdio.h>
#include <stdlib.h>
#include <check.h>

#include "multiset_int_double.h"

START_TEST (test_multiset_int_double_create)
{
  multiset_int_double_t *set;

  set = multiset_int_double_create();
  ck_assert_ptr_ne(set, NULL);

  multiset_int_double_destroy(set);
}
END_TEST

START_TEST(test_multiset_int_double_insert)
{
  int indices[] = {9, 4, 7, 8, 1, 3, 6};
  int depths[] = {3, 1, 2, 2, 1, 1, 2};
  double values[] = {0.5, 0.25, 0.33, 0.11, 0.78, 0.2, 0.6};
  int i;

  multiset_int_double_t *set;

  set = multiset_int_double_create();
  ck_assert_ptr_ne(set, NULL);

  for (i = 0; i < sizeof(indices)/sizeof(int); i ++) {

    ck_assert(multiset_int_double_insert(set, indices[i], depths[i], values[i]) > 0);

  }

  /* multiset_int_double_dump(set); */

  ck_assert(multiset_int_double_total_count(set) == 7);

  for (i = 0; i < sizeof(indices)/sizeof(int); i ++) {
    ck_assert(multiset_int_double_is_element(set, indices[i], depths[i]));
  }

  ck_assert(multiset_int_double_depth_count(set, 0) == 0);
  ck_assert(multiset_int_double_depth_count(set, 1) == 3);
  ck_assert(multiset_int_double_depth_count(set, 2) == 3);
  ck_assert(multiset_int_double_depth_count(set, 3) == 1);
  ck_assert(multiset_int_double_depth_count(set, 4) == 0);

  multiset_int_double_destroy(set);
}
END_TEST


START_TEST(test_multiset_int_double_remove)
{
  int indices[] = {9, 4, 7, 8, 1, 3, 6};
  int depths[] = {3, 1, 2, 2, 1, 1, 2};
  double values[] = {0.5, 0.25, 0.33, 0.11, 0.78, 0.2, 0.6};
  int i;

  multiset_int_double_t *set;

  set = multiset_int_double_create();
  ck_assert_ptr_ne(set, NULL);

  ck_assert(multiset_int_double_remove(set, 8, 1) == 0);

  ck_assert(multiset_int_double_total_count(set) == 0);

  for (i = 0; i < sizeof(indices)/sizeof(int); i ++) {

    ck_assert(multiset_int_double_insert(set, indices[i], depths[i], values[i]) > 0);

  }

  ck_assert(multiset_int_double_remove(set, 5, 1) == 0);

  ck_assert(multiset_int_double_total_count(set) == 7);

  for (i = 0; i < sizeof(indices)/sizeof(int); i ++) {
    ck_assert(multiset_int_double_is_element(set, indices[i], depths[i]));
  }

  /*
   * Remove middle
   */
  ck_assert(multiset_int_double_remove(set, 8, 2) == 1);

  ck_assert(multiset_int_double_total_count(set) == 6);
  
  ck_assert(multiset_int_double_is_element(set, 8, 2) == 0);

  for (i = 0; i < sizeof(indices)/sizeof(int); i ++) {
    if (indices[i] != 8) {
      ck_assert(multiset_int_double_is_element(set, indices[i], depths[i]));
    }
  }

  /* 2nd remove should return 0, ie nothing removed */
  ck_assert(multiset_int_double_remove(set, 8, 2) == 0); 

  /* Reinsert */
  ck_assert(multiset_int_double_insert(set, 8, 2, 1.0) > 0);
  ck_assert(multiset_int_double_total_count(set) == 7);


  /*
   * Remove first
   */
  ck_assert(multiset_int_double_remove(set, 1, 1) == 1);

  ck_assert(multiset_int_double_total_count(set) == 6);
  
  ck_assert(multiset_int_double_is_element(set, 1, 1) == 0);

  for (i = 0; i < sizeof(indices)/sizeof(int); i ++) {
    if (indices[i] != 1) {
      ck_assert(multiset_int_double_is_element(set, indices[i], depths[i]));
    }
  }

  /* 2nd remove should return 0, ie nothing removed */
  ck_assert(multiset_int_double_remove(set, 1, 1) == 0); 

  ck_assert(multiset_int_double_insert(set, 1, 1, 1.0) > 0);

  /*
   * Remove last
   */
  ck_assert(multiset_int_double_remove(set, 9, 3) == 1);

  ck_assert(multiset_int_double_total_count(set) == 6);
  
  ck_assert(multiset_int_double_is_element(set, 9, 3) == 0);

  for (i = 0; i < sizeof(indices)/sizeof(int); i ++) {
    if (indices[i] != 9) {
      ck_assert(multiset_int_double_is_element(set, indices[i], depths[i]));
    }
  }

  /* 2nd remove should return 0, ie nothing removed */
  ck_assert(multiset_int_double_remove(set, 9, 3) == 0); 

  multiset_int_double_destroy(set);
}
END_TEST

START_TEST(test_multiset_int_double_choice)
{
  int indices[] = {9, 4, 7, 8, 1, 3, 6};
  int depths[] = {3, 1, 2, 2, 1, 1, 2};
  double values[] = {0.5, 0.25, 0.33, 0.11, 0.78, 0.2, 0.6};
  int i;

  multiset_int_double_t *set;

  int depth;
  int ndepths;

  int index;
  int nindex;

  int maxdepth;

  set = multiset_int_double_create();
  ck_assert_ptr_ne(set, NULL);

  ck_assert(multiset_int_double_remove(set, 8, 2) == 0);

  ck_assert(multiset_int_double_total_count(set) == 0);

  for (i = 0; i < sizeof(indices)/sizeof(int); i ++) {

    ck_assert(multiset_int_double_insert(set, indices[i], depths[i], values[i]) > 0);

  }

  ck_assert(multiset_int_double_remove(set, 5, 1) == 0);
  ck_assert(multiset_int_double_total_count(set) == 7);

  for (i = 0; i < sizeof(indices)/sizeof(int); i ++) {
    ck_assert(multiset_int_double_is_element(set, indices[i], depths[i]));
  }

  /*
   * Choosing Depth 
   */
  maxdepth = 5;
  ck_assert(multiset_int_double_choose_depth(set, 0.0, maxdepth, &depth, &ndepths) == 0);
  ck_assert(depth == 1);
  ck_assert(ndepths == 3);

  ck_assert(multiset_int_double_choose_depth(set, 0.34, maxdepth, &depth, &ndepths) == 0);
  ck_assert(depth == 2);
  ck_assert(ndepths == 3);

  ck_assert(multiset_int_double_choose_depth(set, 0.67, maxdepth, &depth, &ndepths) == 0);
  ck_assert(depth == 3);
  ck_assert(ndepths == 3);

  maxdepth = 2;
  ck_assert(multiset_int_double_choose_depth(set, 0.0, maxdepth, &depth, &ndepths) == 0);
  ck_assert(depth == 1);
  ck_assert(ndepths == 2);

  ck_assert(multiset_int_double_choose_depth(set, 0.34, maxdepth, &depth, &ndepths) == 0);
  ck_assert(depth == 1);
  ck_assert(ndepths == 2);

  ck_assert(multiset_int_double_choose_depth(set, 0.67, maxdepth, &depth, &ndepths) == 0);
  ck_assert(depth == 2);
  ck_assert(ndepths == 2);

  /*
   * Choosing index at depth 1 
   */
  ck_assert(multiset_int_double_choose_index(set, 1, 0.0, &index, &nindex) == 0);
  ck_assert(index == 1);
  ck_assert(nindex == 3);

  ck_assert(multiset_int_double_choose_index(set, 1, 0.34, &index, &nindex) == 0);
  ck_assert(index == 3);
  ck_assert(nindex == 3);

  ck_assert(multiset_int_double_choose_index(set, 1, 0.67, &index, &nindex) == 0);
  ck_assert(index == 4);
  ck_assert(nindex == 3);

  /*
   * Choosing index at depth 2 
   */
  ck_assert(multiset_int_double_choose_index(set, 2, 0.0, &index, &nindex) == 0);
  ck_assert(index == 6);
  ck_assert(nindex == 3);

  ck_assert(multiset_int_double_choose_index(set, 2, 0.34, &index, &nindex) == 0);
  ck_assert(index == 7);
  ck_assert(nindex == 3);

  ck_assert(multiset_int_double_choose_index(set, 2, 0.67, &index, &nindex) == 0);
  ck_assert(index == 8);
  ck_assert(nindex == 3);

  /*
   * Choosing index at depth 3 
   */
  ck_assert(multiset_int_double_choose_index(set, 3, 0.0, &index, &nindex) == 0);
  ck_assert(index == 9);
  ck_assert(nindex == 1);

  ck_assert(multiset_int_double_choose_index(set, 3, 0.34, &index, &nindex) == 0);
  ck_assert(index == 9);
  ck_assert(nindex == 1);

  ck_assert(multiset_int_double_choose_index(set, 3, 0.67, &index, &nindex) == 0);
  ck_assert(index == 9);
  ck_assert(nindex == 1);

  multiset_int_double_destroy(set);
}
END_TEST

START_TEST(test_multiset_int_double_binary_readwrite)
{
  int indices[] = {9, 4, 7, 8, 1, 3, 6};
  int depths[] = {3, 1, 2, 2, 1, 1, 2};
  double values[] = {0.5, 0.25, 0.33, 0.11, 0.78, 0.2, 0.6};
  int i;
  
  multiset_int_double_t *out;
  multiset_int_double_t *in;

  FILE *fp;

  double value;

  out = multiset_int_double_create();
  ck_assert(out != NULL);

  in = multiset_int_double_create();
  ck_assert(in != NULL);
  
  for (i = 0; i < sizeof(indices)/sizeof(int); i ++) {

    ck_assert(multiset_int_double_insert(out, indices[i], depths[i], values[i]) >= 0);

  }

  fp = fopen("multiset_int_double_binary.dat", "w");
  ck_assert(fp != NULL);

  ck_assert(multiset_int_double_write_binary(out,
					     (multiset_int_double_write_t)fwrite,
					     (void*)fp) >= 0);

  fclose(fp);

  fp = fopen("multiset_int_double_binary.dat", "r");
  ck_assert(fp != NULL);

  ck_assert(multiset_int_double_read_binary(in,
					    (multiset_int_double_read_t)fread,
					    (void*)fp) >= 0);

  for (i = 0; i < sizeof(indices)/sizeof(int); i ++) {

    ck_assert(multiset_int_double_get(in, indices[i], depths[i], &value) >= 0);

    ck_assert(value == values[i]);
    
  }
  fclose(fp);
  
  multiset_int_double_destroy(out);
  multiset_int_double_destroy(in);
}
END_TEST

Suite *
multiset_int_double_suite (void)
{
  Suite *s = suite_create ("Multiset Int");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_multiset_int_double_create);
  tcase_add_test (tc_core, test_multiset_int_double_insert);
  tcase_add_test (tc_core, test_multiset_int_double_remove);
  tcase_add_test (tc_core, test_multiset_int_double_choice);

  tcase_add_test (tc_core, test_multiset_int_double_binary_readwrite);
  
  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = multiset_int_double_suite ();
  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
