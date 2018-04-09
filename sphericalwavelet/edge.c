//
//    Spherical Subdivision/Wavelet library
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
#include <stdlib.h>
#include <string.h>

#include "edge.h"

#include "slog.h"

void
edge_init(edge_t *e)
{
  e->a = -1;
  e->b = -1;

  e->parent = -1;
  e->child_edges[0] = -1;
  e->child_edges[1] = -1;

  e->triangles[0] = -1;
  e->triangles[1] = -1;
}

int
edge_add_triangle(edge_t *e, int ti)
{
  if (e == NULL || ti < 0) {
    ERROR("invalid parameters: %p %d",
	  e, ti);
    return -1;
  }
  
  if (e->triangles[0] < 0) {
    e->triangles[0] = ti;
    return 0;
  } else if (e->triangles[1] < 0) {
    e->triangles[1] = ti;
    return 0;
  }

  ERROR("no slots for triangle");
  return -1;
}
