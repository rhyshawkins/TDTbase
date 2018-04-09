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
#ifndef ttree_int_h
#define ttree_int_h

typedef struct _ttree_int ttree_int_t;

ttree_int_t *
ttree_int_create(int maxk);

void
ttree_int_destroy(ttree_int_t *t);

int
ttree_int_insert(ttree_int_t *t, int k, const char *string, int incr);

int
ttree_int_get(ttree_int_t *t, int k, const char *string, int *count);

typedef int (*ttree_int_iterate_t)(void *user, const char *string, int count);
int
ttree_int_iterate(ttree_int_t *t, int k, ttree_int_iterate_t iterator, void *user);

#endif /* ttree_int_h */
