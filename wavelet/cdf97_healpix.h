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

#ifndef cdf97_healpix_h
#define cdf97_healpix_h

typedef struct cdf97_healpix_ cdf97_healpix_t;
struct cdf97_healpix_ {

  int width;
  int height;

  double *tile[12];

};

typedef struct cdf97_healpix_workspace_ cdf97_healpix_workspace_t;
struct cdf97_healpix_workspace_ {
  int width;
  int height;

  double *row;
  double *tile[12];
};

cdf97_healpix_t *
cdf97_healpix_create(int width);

void
cdf97_healpix_destroy(cdf97_healpix_t *c);

cdf97_healpix_workspace_t *
cdf97_healpix_workspace_create(int width);

void
cdf97_healpix_workspace_destroy(cdf97_healpix_workspace_t *w);

int
cdf97_healpix_forward(cdf97_healpix_t *c, cdf97_healpix_workspace_t *w);

int
cdf97_healpix_inverse(cdf97_healpix_t *c, cdf97_healpix_workspace_t *w);

int
cdf97_healpix_forward_rows(cdf97_healpix_t *c, cdf97_healpix_workspace_t *w, int width);

int
cdf97_healpix_forward_cols(cdf97_healpix_t *c, cdf97_healpix_workspace_t *w, int width);

/*
 * Internal functions
 */

/*
 * These two functions compute where the specified tile/row/col/direction wraps to.
 * The return values are the next tile index, an integer type with 0 meaning that 
 * the transition is to a row in the next tile and 1 meaning a column. The index is
 * the index of the row or column and the next_dir represents the direction along
 * the row or column, +1 for positive starting at 0, -1 for negative starting at (size - 1)
 */
int cdf97_healpix_traverse_row(cdf97_healpix_t *c, 
			       int width,
			       int tile, 
			       int row, 
			       int dir, 
			       int *next_tile, 
			       int *type,
			       int *index,
			       int *next_dir);

int cdf97_healpix_traverse_col(cdf97_healpix_t *c, 
			       int width,
			       int tile, 
			       int col, 
			       int dir, 
			       int *next_tile, 
			       int *type,
			       int *index,
			       int *next_dir);

int cdf97_healpix_fill_rows(cdf97_healpix_t *c, 
			    int width,
			    int tile, 
			    int row, 
			    int dir, 
			    double *buffer, 
			    int n);

int cdf97_healpix_fill_cols(cdf97_healpix_t *c, 
			    int width, 
			    int tile, 
			    int col, 
			    int dir, 
			    double *buffer, 
			    int n);


#endif /* cdf97_healpix_h */
