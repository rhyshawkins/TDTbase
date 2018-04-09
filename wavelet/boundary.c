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

#include "boundary.h"

int
wavelet_boundary_periodic(int i, int width)
{
  if (i < 0) {
    return wavelet_boundary_periodic(i + width, width);
  } else {
    return i % width;
  }
}

int
wavelet_boundary_reflect(int i, int width)
{
  if (i < 0) {
    return wavelet_boundary_reflect(-i, width);
  } else if (i >= width) {
    return wavelet_boundary_reflect(2*width - 2 - i, width);
  } else {
    return i;
  }
}



