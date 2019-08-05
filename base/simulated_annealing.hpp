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
#ifndef simulated_annealing_hpp
#define simulated_annealing_hpp

#include "global.hpp"

double
simulated_annealing_optimize(global &global,
			     int iterations,
			     double Tmax,
			     int &fail_count,
			     double &acceptance)
{
  //
  // Auto scaling for efficient jump sizes
  //
  double scale_factor = 1.0e-3;
  
  //
  // Initial likelihood
  //
  double current_N = 0.0;
  double current_P = global.likelihood(current_N, true);

  //
  // Prune unused model parameters
  //
  std::vector<bool> used(global.model->cells.size(), false);
  for (auto &o : global.data.obs) {
    for (auto &mi : o.idx) {
      used[mi] = true;
    }
  }

  int pruned_count = 0;
  for (int j = used.size() - 1; j >= 0; j --) {
    if (!used[j]) {
      global.model->cells.erase(global.model->cells.begin() + j);
      pruned_count ++;
    }
  }

  if (pruned_count > 0) {
    printf("Pre  prune: %12.6f %6d\n", current_P, (int)used.size());
    current_P = global.likelihood(current_N, true);
    printf("Post prune: %12.6f %6d %d\n", current_P, (int)global.model->cells.size(), pruned_count);
  }

  std::vector<double> scale(global.model->cells.size(), 1.0);
  std::vector<int> local_proposed(global.model->cells.size(), 0);
  std::vector<int> local_accepted(global.model->cells.size(), 0);

  fail_count = 0;
  int proposed = 0;
  int accepted = 0;

  for (int i = 0; i < iterations; i ++) {

    double T = temperature_power(i, iterations, Tmax, 2.0);

    //
    // Perturb model parameters
    //
    // If the perturbation violates prior bounds, we just keep going hoping
    // that the next perturbation won't. The number of prior violations is
    // counted which can be used as a diagnostic. (Proposal width is too large)
    //
    int j = 0;
    for (auto &c : global.model->cells) {

      double oldv = c.v;
      double newv;
      double logpriorratio;
      
      if (global.prior->propose(global.random,
				T * scale[j],
				oldv,
				newv,
				logpriorratio)) {
	local_proposed[j] ++;
	proposed ++;
       
	c.v = newv;

	//
	// Compute perturbed likelihood
	//
	double proposed_N = 0.0;
	double proposed_P = global.likelihood(proposed_N, false);
	
	//
	// Accept/reject
	//
	double u = log(global.random.uniform());
	if (u < (current_P - proposed_P)/T) {
	  //
	  // Accept: set likelihood and keep model as is
	  //
	  current_P = proposed_P;
	  accepted ++;
	  local_accepted[j] ++;
	  
	} else {
	  //
	  // Reject: restore old model
	  //
	  c.v = oldv;
	}

	//
	// Update adaptive step size
	//
	if (local_proposed[j] >= 10 && local_proposed[j] % 10 == 0) {
	  
	  double local_acceptance = (double)local_accepted[j]/(double)local_proposed[j];
	  if (local_acceptance < 0.5) {
	    scale[j] *= (1.0 - scale_factor);
	  } else {
	    scale[j] *= (1.0 + scale_factor);
	  }
	}
	
      } else {
	fail_count ++;
      }

    }


    if (proposed > 0) {
      acceptance = (double)accepted/(double)proposed;
    }

    if ((i + 1) % 1000 == 0) {

      INFO("%6d %12.6f %10.6f\n",
	   i + 1,
	   current_P,
	   acceptance);

    }
	
  }

  return current_P;
}
			     
#endif // simulated_annealing_hpp
  
