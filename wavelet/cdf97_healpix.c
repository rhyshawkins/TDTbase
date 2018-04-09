//
//    Wavelet transform library
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

#include "cdf97_healpix.h"

#include "slog.h"

cdf97_healpix_t *
cdf97_healpix_create(int width)
{
  cdf97_healpix_t *r;
  int i;
  int size;

  r = malloc(sizeof(cdf97_healpix_t));
  if (r == NULL) {
    return NULL;
  }

  r->width = width;
  r->height = width;

  size = width * width;

  for (i = 0; i < 12; i ++) {
    r->tile[i] = malloc(sizeof(double) * size);
    if (r->tile[i] == NULL) {
      return NULL;
    }
    memset(r->tile[i], 0, sizeof(double) * size);
  }

  return r;
}

void
cdf97_healpix_destroy(cdf97_healpix_t *c)
{
  int i;

  if (c != NULL) {

    for (i = 0; i < 12; i ++) {
      free(c->tile[i]);
    }

    free(c);
  }
}

cdf97_healpix_workspace_t *
cdf97_healpix_workspace_create(int width)
{
  cdf97_healpix_workspace_t *r;
  int i;
  int size;

  r = malloc(sizeof(cdf97_healpix_workspace_t));
  if (r == NULL) {
    return NULL;
  }

  r->width = width;
  r->height = width;

  size = width * width;

  for (i = 0; i < 12; i ++) {
    r->tile[i] = malloc(sizeof(double) * size);
    if (r->tile[i] == NULL) {
      return NULL;
    }
    memset(r->tile[i], 0, sizeof(double) * size);
  }

  r->row = malloc(sizeof(double) * (width + 8));
  if (r->row == NULL) {
    return NULL;
  }
  memset(r->row, 0, sizeof(double) * (width + 8));

  return NULL;
}

void
cdf97_healpix_workspace_destroy(cdf97_healpix_workspace_t *c)
{
  int i;

  if (c != NULL) {

    free(c->row);

    for (i = 0; i < 12; i ++) {
      free(c->tile[i]);
    }

    free(c);
  }
}

int
cdf97_healpix_forward(cdf97_healpix_t *c, cdf97_healpix_workspace_t *workspace)
{
  int t;
  int w;

  int i;
  int j;

  int di;
  int dj;
  


  for (w = c->width; w >= 2; w /= 2) {

    /*
     * Do all the rows
     */
    cdf97_healpix_forward_rows(c, workspace, w);

    /*
     * Do all the cols
     */
    cdf97_healpix_forward_cols(c, workspace, w);

    /*
     * Results interleaved in the workspace object, de-interleave and copy back.
     */
    for (t = 0; t < 12; t ++) {

      for (j = 0; j < w; j ++) {

	if (j % 2 == 0) {
	  dj = j;
	} else {
	  dj = w/2 + j;
	}

	for (i = 0; i < w; i ++) {

	  if (i % 2 == 0) {
	    di = i;
	  } else {
	    di = w/2 + i;
	  }
	  
	  c->tile[t][dj * c->width + di] = workspace->tile[t][j * workspace->width + i];
	}
      }
    }
  }
    
  return 0;
}

int
cdf97_healpix_inverse(cdf97_healpix_t *c, cdf97_healpix_workspace_t *workspace)
{
  return -1;
}

int
cdf97_healpix_forward_rows(cdf97_healpix_t *c, cdf97_healpix_workspace_t *workspace, int width)
{
  int t;
  int r;
  int i;
  
  int left_t;
  int right_t;

  for (t = 0; t < 12; t ++) {

    switch (t) {
    case 0:
      left_t = 4;
      right_t = 1;
      break;
      
    case 1:
      left_t = 5;
      right_t = 2;
      break;

    case 2:
      left_t = 6;
      right_t = 3;
      break;

    case 3:
      left_t = 7;
      right_t = 0;
      break;

    default:
      ERROR("invalid tile %d", t);
      return -1;
    }

    for (r = 0; r < width; r ++) {

      /*
       * Copy left 4 points to workspace
       */
      for (i = 0; i < 4; i ++) {
	/* This is wrong when width is less than 4, need to write a recursive function */
	workspace->row[i] = c->tile[left_t][c->width * r + width - 5 + i];
      }
      
      /*
       * Copy row to workspace
       */
      for (i = 0; i < width; i ++) {
	workspace->row[i + 4] = c->tile[t][c->width * r + i];
      }

      /*
       * Copy right 4 points to workspace
       */

      /* 
       * Forward lifting transform
       */
      
      /*
       * Copy row to workspace
       */
      for (i = 0; i < width; i ++) {
	workspace->tile[t][c->width * r + i] = workspace->row[i + 4];
      }
    }
  }

  return -1;
}

int
cdf97_healpix_forward_cols(cdf97_healpix_t *c, cdf97_healpix_workspace_t *w, int width)
{
  return -1;
}

int cdf97_healpix_traverse_row(cdf97_healpix_t *c, 
			       int width,
			       int tile, 
			       int row, 
			       int dir, 
			       int *next_tile, 
			       int *type,
			       int *index,
			       int *next_dir)
{
  static const int LEFT_TILE[12] = {
    4, 5, 6, 7, 
    11, 8, 9, 10, 
    11, 8, 9, 10
  };

  static const int RIGHT_TILE[12] = {
    1, 2, 3, 0, 
    0, 1, 2, 3, 
    5, 6, 7, 4
  };

  static const int LEFT_TYPE[12] = {
    0, 0, 0, 0,
    0, 0, 0, 0,
    1, 1, 1, 1
  };

  static const int RIGHT_TYPE[12] = {
    1, 1, 1, 1,
    0, 0, 0, 0,
    0, 0, 0, 0
  };

  if (tile < 0 || tile >= 12) {
    ERROR("invalid tile");
    return -1;
  }

  if (dir < 0) {
    *next_tile = LEFT_TILE[tile];
    *type = LEFT_TYPE[tile];
    /* If we rotate from rows to cols we mirror the index */
    if (*type == 1) {
      *index = width - 1 - row;
    } else {
      *index = row;
    }
    *next_dir = dir;
  } else {
    *next_tile = RIGHT_TILE[tile];
    *type = RIGHT_TYPE[tile];
    /* If we rotate from rows to cols we mirror the index */
    if (*type == 1) {
      *index = width - 1 - row;
    } else {
      *index = row;
    }
    *next_dir = dir;
  }
    
  return 0;
}

int cdf97_healpix_traverse_col(cdf97_healpix_t *c, 
			       int width,
			       int tile, 
			       int col, 
			       int dir, 
			       int *next_tile, 
			       int *type,
			       int *index,
			       int *next_dir)
{
  static const int TOP_TILE[12] = {
    3, 0, 1, 2, 
    3, 0, 1, 2, 
    4, 5, 6, 7
  };

  static const int BOTTOM_TILE[12] = {
    5, 6, 7, 4,
    8, 9, 10, 11,
    9, 10, 11, 8
  };

  static const int TOP_TYPE[12] = {
    0, 0, 0, 0,
    1, 1, 1, 1,
    1, 1, 1, 1
  };

  static const int BOTTOM_TYPE[12] = {
    1, 1, 1, 1,
    1, 1, 1, 1,
    0, 0, 0, 0
  };

  if (tile < 0 || tile >= 12) {
    ERROR("invalid tile");
    return -1;
  }

  if (dir < 0) {
    *next_tile = TOP_TILE[tile];
    *type = TOP_TYPE[tile];
    /* If we rotate from cols to rows we mirror the index */
    if (*type == 0) {
      *index = width - 1 - col;
    } else {
      *index = col;
    }
    *next_dir = dir;
  } else {
    *next_tile = BOTTOM_TILE[tile];
    *type = BOTTOM_TYPE[tile];
    /* If we rotate from cols to rows we mirror the index */
    if (*type == 0) {
      *index = width - 1 - col;
    } else {
      *index = col;
    }
    *next_dir = dir;
  }
    
  return 0;
}

int cdf97_healpix_fill_rows(cdf97_healpix_t *c, 
			    int width,
			    int tile, 
			    int row, 
			    int dir, 
			    double *buffer, 
			    int n)
{
  int next_tile;
  int next_type;
  int next_index;
  int next_dir;

  int nc;
  int i;

  if (cdf97_healpix_traverse_row(c, 
				 width,
				 tile, 
				 row, 
				 dir, 
				 &next_tile, 
				 &next_type,
				 &next_index,
				 &next_dir) < 0) {
    return -1;
  }

  if (n > width) {
    nc = width;
  } else {
    nc = n;
  }

  for (i = 0; i < nc; i ++) {
    if (next_type == 0) {
      if (next_dir < 0) {
	buffer[i] = c->tile[next_tile][next_index * width + width - 1 - i];
      } else {
	buffer[i] = c->tile[next_tile][next_index * width + i];
      }
    } else {
      if (next_dir < 0) {
	buffer[i] = c->tile[next_tile][next_index + width*(width - 1 - i)];
      } else {
	buffer[i] = c->tile[next_tile][next_index + width*i];
      }
    }
  }

  if (nc < n) {
    if (next_type == 0) {
      return cdf97_healpix_fill_rows(c, width, next_tile, next_index, next_dir, buffer + nc, n - nc);
    } else {
      return cdf97_healpix_fill_cols(c, width, next_tile, next_index, next_dir, buffer + nc, n - nc);
    }
  }
}

int cdf97_healpix_fill_cols(cdf97_healpix_t *c, 
			    int width,
			    int tile, 
			    int col, 
			    int dir, 
			    double *buffer, 
			    int n)
{
  int next_tile;
  int next_type;
  int next_index;
  int next_dir;

  int nc;
  int i;

  if (cdf97_healpix_traverse_col(c, 
				 width,
				 tile, 
				 col, 
				 dir, 
				 &next_tile, 
				 &next_type,
				 &next_index,
				 &next_dir) < 0) {
    return -1;
  }

  if (n > width) {
    nc = width;
  } else {
    nc = n;
  }

  for (i = 0; i < nc; i ++) {
    if (next_type == 0) {
      if (next_dir < 0) {
	buffer[i] = c->tile[next_tile][next_index * width + width - 1 - i];
      } else {
	buffer[i] = c->tile[next_tile][next_index * width + i];
      }
    } else {
      if (next_dir < 0) {
	buffer[i] = c->tile[next_tile][next_index + width*(width - 1 - i)];
      } else {
	buffer[i] = c->tile[next_tile][next_index + width*i];
      }
    }
  }

  if (nc < n) {
    if (next_type == 0) {
      return cdf97_healpix_fill_rows(c, width, next_tile, next_index, next_dir, buffer + nc, n - nc);
    } else {
      return cdf97_healpix_fill_cols(c, width, next_tile, next_index, next_dir, buffer + nc, n - nc);
    }
  }
}

