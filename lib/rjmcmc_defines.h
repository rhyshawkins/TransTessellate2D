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

#ifndef RJMCMC_DEFINES_H
#define RJMCMC_DEFINES_H

#include "rjmcmc_debug.h"

/* */
#define RJMCMC_INDEXCHECKVOID(i, len, msg) \
  if (i < 0 || i >= len) { \
  rjmcmc_error(msg); \
  return; \
  }

#define RJMCMC_INDEXCHECKPTR(i, len, msg) \
  if (i < 0 || i >= len) { \
  rjmcmc_error(msg); \
  return NULL; \
  }

#define RJMCMC_NULLCHECKVOID(p, msg) \
  if (p == NULL) { \
  rjmcmc_error(msg); \
  return; \
  }

#define RJMCMC_NULLCHECKINT(p, msg) \
  if (p == NULL) { \
  rjmcmc_error(msg); \
  return -1; \
  }

#define RJMCMC_NULLCHECKPTR(p, msg) \
  if (p == NULL) { \
  rjmcmc_error(msg); \
  return NULL; \
  }

#define RJMCMC_CONDITIONCHECKVOID(cond, msg) \
  if (cond) { \
  rjmcmc_error(msg); \
  return;\
  }

#define RJMCMC_CONDITIONCHECKINT(cond, msg) \
  if (cond) { \
  rjmcmc_error(msg); \
  return -1;\
  }

#define RJMCMC_CONDITIONCHECKPTR(cond, msg) \
  if (cond) { \
  rjmcmc_error(msg); \
  return NULL; \
  }

#define RJMCMC_INTCHECKPTR(i, msg) \
  if (i < 0) { \
  rjmcmc_error(msg); \
  return NULL; \
  }

#define RJMCMC_INTCHECKINT(i, msg) \
  if (i < 0) { \
  rjmcmc_error(msg); \
  return -1; \
  }


#endif /* RJMCMC_DEFINES_H */
