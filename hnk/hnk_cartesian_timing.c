//
//    HNK Library : A library for computing combinations of arrangements of
//    general trees for the Trans-dimensional Tree algorithm. See
//
//      R Hawkins and M Sambridge, "Geophysical imaging using trans-dimensional trees",
//      Geophysical Journal International, 2015, 203:2, 972 - 1000,
//      https://doi.org/10.1093/gji/ggv326
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

#include "hnk.h"
#include "hnk_cartesian.h"

int main(int argc, char *argv[])
{
  int dim;
  int h;
  int kmax;

  hnk_t *t;
  mpz_t a;
  int i;

  dim = 3;
  h = 5;
  kmax = 1000;

  switch(dim) {
  case 2:
    t = hnk_cartesian_34_create(h,
				kmax);
    break;

  case 3:
    t = hnk_cartesian_78_create(h,
				kmax);
    break;

  default:
    fprintf(stderr, "error: unhandled dimension\n");
    return -1;
  }

  if (t == NULL) {
    fprintf(stderr, "error: failed to create tree\n");
    return -1;
  }

  mpz_init(a);

  for (i = 1; i <= 1000; i ++) {
    if (hnk_get_hnk(t, h, i, a) < 0) {
      fprintf(stderr, "error: failed to get count\n");
      return -1;
    }
    
    gmp_printf("%d: %Zd\n", i, a);
  }
  
  mpz_clear(a);
  hnk_destroy(t);
  
  return 0;
}
    
  
