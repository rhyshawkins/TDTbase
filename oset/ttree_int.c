
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ttree_int.h"

#include "slog.h"

typedef struct _ttree_node ttree_node_t;
struct _ttree_node {
  int count;
  char c;

  ttree_node_t *left;
  ttree_node_t *eq;
  ttree_node_t *right;
};
  
struct _ttree_int {

  int maxk;
  ttree_node_t **tree;

};

static void ttree_node_free(ttree_node_t *t);

static ttree_node_t *ttree_node_insert(ttree_node_t *t, const char *string, int incr);

static int ttree_node_get(ttree_node_t *t, const char *string, int *count);

static int ttree_node_iterate(ttree_node_t *t, ttree_int_iterate_t iterator, void *user, char *str, int len);

ttree_int_t *
ttree_int_create(int maxk)
{
  ttree_int_t *t;
  int i;

  t = malloc(sizeof(ttree_int_t));
  if (t == NULL) {
    return NULL;
  }

  t->maxk = maxk;
  t->tree = malloc(sizeof(ttree_node_t *) * (maxk + 1));
  if (t->tree == NULL) {
    return NULL;
  }

  for (i = 0; i < (maxk + 1); i ++) {
    t->tree[i] = NULL;
  }

  return t;
}

void
ttree_int_destroy(ttree_int_t *t)
{
  int i;

  if (t != NULL) {
    for (i = 0; i < (t->maxk + 1); i ++) {
      ttree_node_free(t->tree[i]);
    }
    free(t->tree);
    free(t);
  }
}

int
ttree_int_insert(ttree_int_t *t, int k, const char *string, int incr)
{
  ttree_node_t *p;

  if (k < 0 || k > t->maxk) {
    ERROR("invalid k %d", k);
    return -1;
  }
  
  if (string == NULL || string[0] == '\0') {
    ERROR("invalid string");
    return -1;
  }
  
  p = ttree_node_insert(t->tree[k], string, incr);
  if (p == NULL) {
    ERROR("failed to insert");
    return -1;
  }

  t->tree[k] = p;
  return 0;
}

int
ttree_int_get(ttree_int_t *t, int k, const char *string, int *count)
{
  return ttree_node_get(t->tree[k], string, count);
}

int
ttree_int_iterate(ttree_int_t *t, int k, ttree_int_iterate_t iterator, void *user)
{
  char buffer[1024];
  
  return ttree_node_iterate(t->tree[k], iterator, user, buffer, 0);
}

/*
 * Internal functions
 */
static void ttree_node_free(ttree_node_t *t)
{
  if (t != NULL) {
    ttree_node_free(t->left);
    ttree_node_free(t->eq);
    ttree_node_free(t->right);
    free(t);
  }
}

static ttree_node_t *ttree_node_insert(ttree_node_t *t, const char *string, int incr)
{
  ttree_node_t *p;

  if (t == NULL) {
    t = malloc(sizeof(ttree_node_t));
    if (t == NULL) {
      return NULL;
    }

    t->c = string[0];
    t->count = 0;
    t->left = NULL;
    t->eq = NULL;
    t->right = NULL;
  }

  if (t->c < string[0]) {

    p = ttree_node_insert(t->left, string, incr);
    if (p == NULL) {
      return NULL;
    }
    t->left = p;

  } else if (t->c > string[0]) {

    p = ttree_node_insert(t->right, string, incr);
    if (p == NULL) {
      return NULL;
    }
    t->right = p;

  } else {

    if (string[1] == '\0') {
      t->count += incr;
    } else {
      p = ttree_node_insert(t->eq, string + 1, incr);
      if (p == NULL) {
	return NULL;
      }
      t->eq = p;
    }
  }

  return t;
}

static int ttree_node_get(ttree_node_t *t, const char *string, int *count)
{
  if (t == NULL) {
    return -1;
  }

  if (t->c < string[0]) {
    return ttree_node_get(t->left, string, count);
  } else if (t->c > string[0]) {
    return ttree_node_get(t->right, string, count);
  } else {
    if (string[1] == '\0') {
      *count = t->count;
      return 0;
    } else {
      return ttree_node_get(t->eq, string + 1, count);
    }
  }
}

static int ttree_node_iterate(ttree_node_t *t, ttree_int_iterate_t iterator, void *user, char *str, int len)
{
  if (t == NULL) {
    return 0;
  }

  if (t->count == 0) {

    if (ttree_node_iterate(t->left, iterator, user, str, len) < 0) {
      return -1;
    }
    
    str[len] = t->c;
    if (ttree_node_iterate(t->eq, iterator, user, str, len + 1) < 0) {
      return -1;
    }

    if (ttree_node_iterate(t->right, iterator, user, str, len) < 0) {
      return -1;
    }

    return 0;

  } else {

    str[len] = t->c;
    str[len + 1] = '\0';
    return iterator(user, str, t->count);

  }
}
