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
#ifndef generalvoronoicartesianexception_hpp
#define generalvoronoicartesianexception_hpp

#include <exception>

#define GENERALVORONOICARTESIANEXCEPTION(fmt, ...) generalvoronoicartesianexception(__FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)

class generalvoronoicartesianexception : public std::exception {
public:

  
  generalvoronoicartesianexception(const char *srcfile,
			    const char *function,
			    int lineno,
			    const char *fmt, ...);
  ~generalvoronoicartesianexception();
  
};

#endif // generalvoronoicartesianexception_hpp
