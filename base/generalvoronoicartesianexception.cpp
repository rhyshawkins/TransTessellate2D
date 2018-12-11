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

#include <stdio.h>
#include <stdarg.h>

#include "generalvoronoicartesianexception.hpp"

extern "C" {
#include "slog.h"
};

generalvoronoicartesianexception::generalvoronoicartesianexception(const char *srcfile,
						     const char *function,
						     int lineno,
						     const char *fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);

  vslog(SLOG_ERROR, srcfile, function, lineno, fmt, ap);
  
  va_end(ap);

}

generalvoronoicartesianexception::~generalvoronoicartesianexception()
{
}
