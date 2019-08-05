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
#ifndef genericinterface_hpp
#define genericinterface_hpp

#include "synthetic.hpp"

extern "C" {

  //
  // This method must tell the framework how many models and how many
  // hierarchical parameters are needed.
  //
  int gvcart_initialise_(int *nmodels,
			 int *nhierarchical);

  typedef int (gvcart_addobservation_t)(int *npoints,
					int *modelindex,
					double *xs,
					double *ys);
  //
  // Load data and for each observation, call a callback to tell
  // the inversion about the points involved
  //
  int gvcart_loaddata_(int *n,
		       const char *filename,
		       gvcart_addobservation_t addobs);

  //
  // For a single observation, compute the prediction given
  // model values for each point
  //
  int gvcart_compute_prediction_(int *nmodels,
				 int *observation,
				 int *npoints,
				 const double *value,
				 double *unused,
				 double *prediction);
  //
  // For the observations, given predictions, compute residuals, likelihood
  // and norm
  //
  int gvcart_compute_likelihood_(int *nmodels,
				 int *nhierarchical,
				 int *nobservation,
				 double *hierarchical,
				 double *predictions,
				 double *residuals,
				 double *unused,
				 double *like,
				 double *norm);
  
  //
  // Used for making synthetic datasets, save data in correct format
  // using predictions to overwrite observations
  //
  int gvcart_savedata_(int *n,
		       const char *filename,
		       double *hierarchical,
		       int *nobservations,
		       double *predictions);

  //
  // For (optional) post processing, this function can be defined to
  // compute derived mean, stddev etc from one or more models
  //
  double gvcart_compute_derived_(int *nmodels,
				 double *x,
				 double *y,
				 double *values);
  
};

#endif // genericinterface_hpp

		 
