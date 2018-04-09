//
//    Simple elapsed time tracking library
//    
//    Copyright (C) 2014 - 2018 Rhys Hawkins
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#include <stdio.h>

#include <time.h>

#include "tracking.h"

void tracking_init(tracking_t *t)
{
  t->n = 0;
  t->mean = 0.0;
  t->started = 0;
};

int tracking_start(tracking_t *t)
{
  if (t->started) {
    fprintf(stderr, "Tracking already started\n");
    return -1;
  }

  t->started = 1;
  clock_gettime(CLOCK_REALTIME, &t->start);

  return 0;
}

int tracking_end(tracking_t *t)
{
  long nsec;
  time_t dsec;
  double msec;
  double delta;

  if (t->started == 0) {
    fprintf(stderr, "Tracking not started\n");
    return -1;
  }

  clock_gettime(CLOCK_REALTIME, &t->end);

  dsec = t->end.tv_sec - t->start.tv_sec;
  nsec = t->end.tv_nsec - t->start.tv_nsec;

  msec = 
    (double)dsec * 1000000.0 + 
    (double)nsec / 1000.0;
    
  delta = msec - t->mean;
  t->n ++;
  t->mean += delta/(double)(t->n);

  t->started = 0;

  return 0;
}

void tracking_print(tracking_t *t, const char *label)
{
  fprintf(stderr, "%s mean time: %.3f us (%d samples)\n", label, t->mean, t->n);
}

int tracking_samples(tracking_t *t)
{
  return t->n;
}

double tracking_mean(tracking_t *t)
{
  return t->mean;
}

