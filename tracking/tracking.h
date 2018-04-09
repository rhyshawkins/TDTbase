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
