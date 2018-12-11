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

#include "global.hpp"

extern "C" {
#include "slog.h"
};

global *
global::current_state = nullptr;

int
global::addobservation(int *npoints,
		       int *modelindices,
		       double *lons,
		       double *lats)
{
  if (current_state != nullptr) {

    for (int i = 0; i < (*npoints); i ++) {

      if (lons[i] < current_state->xmin) {
	ERROR("X out of range: %d %10.6f < %10.6f\n", i, lons[i], current_state->xmin);
	return -1;
      }
      if (lons[i] > current_state->xmax) {
	ERROR("X out of range: %d %10.6f > %10.6f\n", i, lons[i], current_state->xmax);
	return -1;
      }
    }

    for (int i = 0; i < (*npoints); i ++) {
  
      if (lats[i] < current_state->ymin) {
	ERROR("Y out of range: %d %10.6f < %10.6f\n", i, lats[i], current_state->ymin);
	return -1;
      }
      if (lats[i] > current_state->ymax) {
	ERROR("Y out of range: %d %10.6f > %10.6f\n", i, lats[i], current_state->ymax);
	return -1;
      }
    }
    
    current_state->data.add(*npoints, modelindices, lons, lats);
    return 0;
    
  }

  return -1;
}
