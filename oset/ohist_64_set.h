//
//    Ordered set library, used for maintaining arbitrary tree based models
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

#ifndef ohist_64_h
#define ohist_64_h

#include <stdint.h>

typedef struct _ohist_64_set ohist_64_set_t;

ohist_64_set_t *
ohist_64_set_create(int maxk);

void
ohist_64_set_destroy(ohist_64_set_t *s);

int
oset_clear(ohist_64_set_t *s);

int
ohist_64_set_insert(ohist_64_set_t *s, 
		    uint64_t idx, 
		    int k, 
		    const int *set,
		    int incr,
		    int *pset,
		    int *matching);

int
ohist_64_set_nelements(const ohist_64_set_t *s, int k);

int 
ohist_64_set_nth_element(const ohist_64_set_t *s, 
			 int k,
			 int n,
			 uint64_t *idx,
			 int *count);

void
ohist_64_set_dump(const ohist_64_set_t *s);

#endif /* ohist_64_set_h */

