//
//    TransTessellate2D : A general Trans-dimensional Tessellation program
//    for 2D Cartesian problems.
//
//    Copyright (C) 2014 - 2019 Rhys Hawkins
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

#pragma once
#ifndef generalvoronoicartesianutil_hpp
#define generalvoronoicartesianutil_hpp

#include <ostream>
#include <istream>
#include <string>

#include <stdio.h>

//
// Templated IO of values that can be over-ridden
//
template <typename value>
bool value_read(std::istream &s, value &v) {
  s.read((char*)&v, sizeof(v));
    
  return true;
}

template <typename value>
bool value_read(FILE *fp, value &v) {
  if (fread(&v, sizeof(v), 1, fp) != 1) {
    return false;
  }

  return true;
}

template <typename value>
bool value_write(std::ostream &s, const value &v)
{
  s.write((char*)&v, sizeof(v));
  
  return true;
}

template <typename value>
bool value_write(FILE *fp, const value &v) {
  if (fwrite(&v, sizeof(v), 1, fp) != 1) {
    return false;
  }

  return true;
}

template <typename value>
double scalartodouble(const value &v)
{
  return (double)v;
}

template <typename value>
value doubletoscalar(const double &d)
{
  return (value)d;
}

template
<
  typename T
>
bool encode(const T &i, char *buffer, int &offset, int length)
{
  if ((offset + (int)sizeof(T)) > length) {
    return false;
  }

  T *p = (T*)(buffer + offset);
  *p = i;
  offset += sizeof(T);
  return true;
}

template
<
  typename T
>
bool decode(T &i, const char *buffer, int &offset, int length)
{
  if ((offset + (int)sizeof(T)) > length) {
    return false;
  }

  const T *p = (const T*)(buffer + offset);
  i = *p;
  offset += sizeof(T);
  return true;
}

std::string mkformatstring(const char *fmt, ...);

#endif // generalvoronoicartesianutil_hpp
