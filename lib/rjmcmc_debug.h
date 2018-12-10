//
//    Trans-d Voronoi/Delaunay example code.
//    
//    Copyright (C) 2018 Rhys Hawkins
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

#ifndef rjmcmc_debug_h
#define rjmcmc_debug_h

#include <stdarg.h>

typedef enum {
  RJMCMC_FATAL = 0,
  RJMCMC_ERROR,
  RJMCMC_WARNING,
  RJMCMC_DEBUG
} rjmcmc_debug_level_t;

typedef void (*rjmcmc_debug_function_t)(rjmcmc_debug_level_t level,
					const char *fmt,
					va_list ap);

rjmcmc_debug_function_t 
rjmcmc_debug_set_output(rjmcmc_debug_function_t fn);

void 
rjmcmc_debug_set_default_output(void);

rjmcmc_debug_level_t
rjmcmc_debug_set_level(rjmcmc_debug_level_t new_level);

void rjmcmc_fatal(const char *fmt, ...);
void rjmcmc_error(const char *fmt, ...);
void rjmcmc_warning(const char *fmt, ...);
void rjmcmc_debug(const char *fmt, ...);

#endif /* rjmcmc_debug_h */
