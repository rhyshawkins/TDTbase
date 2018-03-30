
#include <stdio.h>
#include <stdlib.h>
#include <check.h>

#include "ttree_int.h"

START_TEST (test_ttree_int_create)
{
  ttree_int_t *set;

  set = ttree_int_create(8);
  ck_assert_ptr_ne(set, NULL);

  ttree_int_destroy(set);
}
END_TEST

START_TEST(test_ttree_int_insert)
{
  char sequence[][17] = { 
    "(.(..(....).).)",
    "((....).(....))",
    ""
  };

  int i;
  int k;
  int count;

  ttree_int_t *set;

  set = ttree_int_create(8);
  ck_assert_ptr_ne(set, NULL);

  k = 3;

  for (i = 0; strlen(sequence[i]) > 0; i ++) {

    ck_assert(ttree_int_insert(set, k, sequence[i], 1) >= 0);

  }

  for (i = 0; strlen(sequence[i]) > 0; i ++) {
    
    ck_assert(ttree_int_get(set, k, sequence[i], &count) >= 0);
    ck_assert(count == 1);

  }

  ttree_int_destroy(set);
}
END_TEST

struct iterate_test_data {
  char **sequences;
  int nseq;
  int count;
  int invalid_count;
};

static int iterate_test_cb(void *user, const char *string, int count)
{
  struct iterate_test_data *data = (struct iterate_test_data*)user;
  int i;
  int valid;

  valid = 0;
  for (i = 0; strlen(data->sequences[i]) > 0; i ++) {
    if (strcmp(data->sequences[i], string) == 0) {
      valid = -1;
      break;
    }
  }
  
  if (!valid) {
    fprintf(stderr, "iterate_test_cb: invalid string %s\n", string);
    data->invalid_count ++;
  }

  data->nseq ++;
  data->count += count;
  
  return 0;
}

START_TEST(test_ttree_int_iterate)
{
  char *sequence[] = { 
    "(.(..(....).).)",
    "((....).(....))",
    "(..(..(....).))",
    ""
  };

  int i;
  int k;

  ttree_int_t *set;

  struct iterate_test_data user;

  set = ttree_int_create(8);
  ck_assert_ptr_ne(set, NULL);

  k = 3;

  for (i = 0; strlen(sequence[i]) > 0; i ++) {

    ck_assert(ttree_int_insert(set, k, sequence[i], 1 + i) >= 0);

  }

  user.sequences = sequence;
  user.nseq = 0;
  user.count = 0;
  user.invalid_count = 0;
  ck_assert(ttree_int_iterate(set, k, iterate_test_cb, &user) >= 0);
  ck_assert(user.invalid_count == 0);
  ck_assert(user.nseq == 3);
  ck_assert(user.count == 6);
}
END_TEST

Suite *
ttree_int_suite (void)
{
  Suite *s = suite_create ("Ttree Int");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_ttree_int_create);
  tcase_add_test (tc_core, test_ttree_int_insert);
  tcase_add_test (tc_core, test_ttree_int_iterate);
  suite_add_tcase (s, tc_core);

  return s;
}

int main (void) 
{
  int number_failed;
  Suite *s = ttree_int_suite ();
  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
