
#ifndef tracking_h
#define tracking_h

#include <time.h>

struct tracking {
  int n;
  double mean;
  int started;
  struct timespec start;
  struct timespec end;
};
typedef struct tracking tracking_t;

void tracking_init(tracking_t *t);

int tracking_start(tracking_t *t);

int tracking_end(tracking_t *t);

void tracking_print(tracking_t *t, const char *label);

int tracking_samples(tracking_t *t);

double tracking_mean(tracking_t *t);

#endif /* tracking_h */
