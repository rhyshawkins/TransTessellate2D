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
#ifndef value_hpp
#define value_hpp

#include <mpi.h>

#include <string>

#include "global.hpp"
#include "perturbation.hpp"

class Value : public Perturbation {
public:

  typedef deltaVoronoi delta_t;
  typedef model_deltaVoronoi model_delta_t;
  
  Value(int _nmodels) :
    undo_mi(-1),
    undo_cell(-1),
    undo_v(0.0),
    nmodels(_nmodels),
    p(new int[_nmodels]),
    a(new int[_nmodels])
  {
    for (int i = 0; i < nmodels; i ++) {
      p[i] = 0;
      a[i] = 0;
    }
  }
  
  ~Value()
  {
    delete [] p;
    delete [] a;
  }

  virtual bool propose(int nmodels,
		       int nhierarchical,
		       int maxcells,
		       int nobs,
		       Rng &random,
		       std::vector<PriorProposal*> &priors,
		       std::vector<PositionPriorProposal*> &position_priors,
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
    int cell = -1;
    double oldv = 0.0;
    double newv;
    int mi = 0;

    relocate = false;
    
    if (this->primary()) {

      //
      // Primary does generation of new value
      //

      if (nmodels > 1) {
	mi = random.uniform(nmodels);
      }
      
      p[mi] ++;
      
      cell = models[mi]->selectvaluecell(random);
      if (cell < 0) {
	perturbation = model_delta_t::mkvalue(mi, -1, 0.0, 0.0);
	validproposal = false;
      } else {	
      
	//
	// Next get the active node
	//
	oldv = models[mi]->value_at_index(cell);
	
	if (priors[mi]->propose(random, temperature, oldv, newv, log_prior_ratio)) {
	  
	  perturbation = model_delta_t::mkvalue(mi, cell, oldv, newv);
	  
	  validproposal = true;
	  
	} else {
	  
	  perturbation = model_delta_t::mkvalue(mi, cell, oldv, 0.0);
	  
	}
	
      }
    }

    this->communicate(validproposal);

    if (validproposal) {
      this->communicate(mi);
      this->communicate(cell);
      this->communicate(newv);

      undo_mi = mi;
      undo_cell = cell;
      undo_v = models[mi]->value_at_index(cell);

      models[mi]->set_value_at_index(cell, newv);

      log_prior_ratio = 0.0;
      
      log_proposal_ratio = priors[mi]->log_proposal_ratio(random, temperature, undo_v, newv);

      log_extra = 0.0;
    }
    
    return validproposal;
  }

  void accept()
  {
    if (undo_cell == -1) {
      throw GENERALVORONOICARTESIANEXCEPTION("No undo information\n");
    }

    a[undo_mi] ++;
    
    undo_mi = -1;
    undo_cell = -1;
    undo_v = 0.0;
  }
  
  void reject(std::vector<cartesianvoronoimodel*> &models,
	      hierarchical_model &hierarchical)
  {
    if (undo_cell == -1) {
      throw GENERALVORONOICARTESIANEXCEPTION("No undo information\n");
    }

    models[undo_mi]->set_value_at_index(undo_cell, undo_v);

    undo_mi = -1;
    undo_cell = -1;
    undo_v = 0.0;
  }

  virtual int pa_categories() const
  {
    return nmodels;
  }
  
  virtual int proposal_count(int mi) const
  {
    return p[mi];
  }
  
  virtual int acceptance_count(int mi) const
  {
    return a[mi];
  }

  virtual const char *displayname() const
  {
    return "Value";
  }
  
private:

  int undo_mi;
  int undo_cell;
  double undo_v;

  int nmodels;
  int *p;
  int *a;
  
};

#endif // value_hpp
