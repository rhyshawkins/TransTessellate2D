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
#ifndef hierarchical_hpp
#define hierarchical_hpp

#include <string>

#include "global.hpp"
#include "perturbation.hpp"

class Hierarchical : public Perturbation {
public:

  typedef deltaVoronoi delta_t;
  typedef model_deltaVoronoi model_delta_t;

  Hierarchical(int _nhierarchical) :
    undo_index(-1),
    undo_v(0.0),
    nhierarchical(_nhierarchical),
    p(new int[_nhierarchical]),
    a(new int[_nhierarchical])
  {
    for (int i = 0; i < nhierarchical; i ++) {
      p[i] = 0;
      a[i] = 0;
    }
  }
  
  ~Hierarchical()
  {
  }

  virtual bool propose(int nmodels,
		       int nhierarchical,
		       int maxcells,
		       int nobs,
		       Rng &random,
		       std::vector<PriorProposal*> &prior,
		       std::vector<PositionPriorProposal*> &position_prior,
		       std::vector<PriorProposal*> &hierarchical_priors,
		       std::vector<cartesianvoronoimodel*> &models,
		       hierarchical_model &hierarchical,
		       double temperature,
		       double &log_prior_ratio,
		       double &log_proposal_ratio,
		       double &log_extra,
		       delta_t *&perturbation,
		       bool &relocate)
  {
    bool validproposal = false;
    int hindex;
    double newv;

    relocate = false;

    if (this->primary()) {
      
    
      hindex = 0;
      if (hierarchical.get_nhierarchical() > 1) {
	hindex = random.uniform(hierarchical.get_nhierarchical());
      }

      p[hindex] ++;

      double oldv = hierarchical.get(hindex);
      if (hierarchical_priors[hindex]->propose(random, temperature, oldv, newv, log_prior_ratio)) {
	
	perturbation = new hierarchical_deltaVoronoi(1, &hindex, &oldv, &newv);
	validproposal = true;
	
      } else {
	perturbation = new hierarchical_deltaVoronoi(1, &hindex, &oldv, &newv);
      }
    }

    this->communicate(validproposal);

    if (validproposal) {

      this->communicate(hindex);
      this->communicate(newv);
			
      undo_index = hindex;
      undo_v = hierarchical.get(hindex);
	
      log_prior_ratio = 0.0;
      log_proposal_ratio = 0.0;
      log_extra = 0.0;
	
      hierarchical.set(hindex, newv);
    }
    
    return validproposal;
  }
    
  
  virtual void accept()
  {
    if (undo_index == -1) {
      throw GENERALVORONOICARTESIANEXCEPTION("No undo information\n");
    }

    a[undo_index] ++;
    
    undo_index = -1;
    undo_v = 0.0;
  }    

  virtual void reject(std::vector<cartesianvoronoimodel*> &models,
		      hierarchical_model &hierarchical)
  {
    if (undo_index == -1) {
      throw GENERALVORONOICARTESIANEXCEPTION("No undo information\n");
    }
    
    hierarchical.set(undo_index, undo_v);
    
    undo_index = -1;
    undo_v = 0.0;
  }

  virtual int pa_categories() const
  {
    return nhierarchical;
  }
  
  virtual int proposal_count(int hi) const
  {
    return p[hi];
  }
  
  virtual int acceptance_count(int hi) const
  {
    return a[hi];
  }

  virtual const char *displayname() const
  {
    return "Hierarchical";
  }
  
private:

  int undo_index;
  double undo_v;

  int nhierarchical;
  int *p;
  int *a;
  
  
};

#endif // hierarchical_hpp
