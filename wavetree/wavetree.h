//
//    Wavetree Library : A library for performed trans-dimensional tree inversion,
//    See
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

#ifndef wavetree_h
#define wavetree_h

typedef enum {

  WT_PERTURB_INVALID = -1,
  WT_PERTURB_NONE = 0,
  WT_PERTURB_BIRTH = 1,
  WT_PERTURB_DEATH = 2,
  WT_PERTURB_VALUE = 3,
  WT_PERTURB_MOVE  = 4,
  WT_PERTURB_HIERARCHICAL = 5,
  WT_PERTURB_PTEXCHANGE = 6,
  WT_PERTURB_PTMODELEXCHANGE = 7

} wavetree_perturb_t;

#endif /* wavetree_h */
