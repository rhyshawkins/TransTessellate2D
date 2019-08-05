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
#ifndef generalvoronoicartesianobservations_hpp
#define generalvoronoicartesianobservations_hpp

#include <vector>

#include "cartesianvoronoimodel.hpp"

//
// A basic container for point based observations
//

class cvcache;

class GeneralVoronoiCartesianObservations {
public:

  struct observation {
    double pred;
    double res;

    std::vector<int> mi;    // Model index
    std::vector<double> xs; // X points
    std::vector<double> ys; // Y points

    std::vector<cvcache*> refs;     // Model references
    
    std::vector<double> values;  // Cached value for forward model
    std::vector<double> weights; // Unused at the moment 
  };

  
  GeneralVoronoiCartesianObservations(int _nmodels) :
    nmodels(_nmodels)
  {
  }

  void add(int npoints, const int *mi, const double *xs, const double *ys)
  {
    int n = obs.size();
    obs.push_back(observation());

    obs[n].pred = 0.0;
    obs[n].res = 0.0;

    obs[n].mi.resize(npoints);
    obs[n].xs.resize(npoints);
    obs[n].ys.resize(npoints);

    obs[n].values.resize(npoints);
    obs[n].weights.resize(npoints);
    
    for (int i = 0; i < npoints; i ++) {

      obs[n].mi[i] = mi[i];
      obs[n].xs[i] = xs[i];
      obs[n].ys[i] = ys[i];
      
    }

  }

  void initialize_model_references(std::vector<cartesianvoronoimodel *> &models)
  {
    for (auto &o: obs) {

      o.refs.clear();
      
      for (int i = 0; i < (int)o.mi.size(); i ++) {

	cvcache *r = models[o.mi[i]]->create_reference(o.xs[i], o.ys[i]);
	if (r == nullptr) {
	  throw GENERALVORONOICARTESIANEXCEPTION("Failed to create point reference");
	}
	
	o.refs.push_back(r);
      }
    }
  }

  int nmodels;
  std::vector<observation> obs;
};

#endif // generalvoronois2observations_hpp
