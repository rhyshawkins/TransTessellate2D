//
//    TransTesselate2D : A general Trans-dimensional Tesselation program
//    for 2D Cartesian problems.
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

#include <stdarg.h>
#include <string.h>
#include <ctype.h>

#include "generalvoronoicartesianutil.hpp"

std::string mkformatstring(const char *fmt, ...)
{
  static char *buffer = nullptr;
  static int buffer_size = -1;

  if (buffer == nullptr) {
    buffer_size = 512;
    buffer = new char[buffer_size];
  }

  va_list ap;
  int size;
  
  va_start(ap, fmt);
  size = vsnprintf(buffer, buffer_size, fmt, ap);
  while (size >= buffer_size) {
    delete [] buffer;
    buffer_size *= 2;
    buffer = new char[buffer_size];
    size = vsnprintf(buffer, buffer_size, fmt, ap);
  }
  va_end(ap);

  return std::string(buffer);
}
