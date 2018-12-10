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

#include "bbox2d.h"

void
bbox2d_set(bbox2d_t *bound,
	   double xmin,
	   double xmax,
	   double ymin,
	   double ymax)
{
  bound->xmin = xmin;
  bound->xmax = xmax;
  bound->ymin = ymin;
  bound->ymax = ymax;
}

void 
bbox2d_initialize(bbox2d_t *bound,
		  double x,
		  double y)
{
  bound->xmin = x;
  bound->xmax = x;

  bound->ymin = y;
  bound->ymax = y;
}

void 
bbox2d_expand(bbox2d_t *bound,
	      double x,
	      double y)
{
  if (x < bound->xmin) {
    bound->xmin = x;
  }

  if (x > bound->xmax) {
    bound->xmax = x;
  }

  if (y < bound->ymin) {
    bound->ymin = y;
  }

  if (y > bound->ymax) {
    bound->ymax = y;
  }
}

void
bbox2d_bound_expand(bbox2d_t *bound,
		    const bbox2d_t *exp_bound)
{
  if (exp_bound->xmin < bound->xmin) {
    bound->xmin = exp_bound->xmin;
  }

  if (exp_bound->xmax > bound->xmax) {
    bound->xmax = exp_bound->xmax;
  }

  if (exp_bound->ymin < bound->ymin) {
    bound->ymin = exp_bound->ymin;
  }

  if (exp_bound->ymax > bound->ymax) {
    bound->ymax = exp_bound->ymax;
  }
}
