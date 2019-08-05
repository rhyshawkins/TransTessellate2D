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
#ifndef move_hpp
#define move_hpp

#include "global.hpp"
#include "perturbation.hpp"

class Move : public Perturbation {
public:

  typedef deltaVoronoi delta_t;
  typedef model_deltaVoronoi model_delta_t;

  Move(int _nmodels) :
    undo_mi(-1),
    undo_cell(-1),
    nmodels(_nmodels),
    p(new int[_nmodels]),
    a(new int[_nmodels])
  {
    for (int i = 0; i < nmodels; i ++) {
      p[i] = 0;
      a[i] = 0;
    }
  }
  
  ~Move()
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
    int mi = 0;
    int cell = -1;
    double newx;
    double newy;

    relocate = false;
    
    if (this->primary()) {

      if (nmodels > 1) {
	mi = random.uniform(nmodels);
      }
      
      p[mi] ++;

      cell = models[mi]->selectmobilecell(random);
      if (cell < 0) {
	validproposal = false;
	perturbation = model_deltaVoronoi::mkmove(mi, -1, 0.0, 0.0, 0.0, 0.0);
      } else {
    
	//
	// Next get the active node
      //
	double oldx, oldy;
	
	models[mi]->position_at_index(cell, &oldx, &oldy);
	
	if (position_priors[mi]->propose(random,
					 temperature,
					 oldx,
					 oldy,
					 newx,
					 newy,
					 log_prior_ratio)) {

	  validproposal = true;
	  perturbation = model_deltaVoronoi::mkmove(mi, cell, oldx, oldy, newx, newy);
	  
	} else {
	  
	  perturbation = model_deltaVoronoi::mkmove(mi, cell, oldx, oldy, 0.0, 0.0);

	}
      }
    }

    this->communicate(validproposal);

    if (validproposal) {

      relocate = true;
      
      this->communicate(mi);
      this->communicate(cell);
      this->communicate(newx);
      this->communicate(newy);

      undo_mi = mi;
      undo_cell = cell;
      undo_v = models[mi]->value_at_index(cell);
      models[mi]->position_at_index(cell, &undo_x, &undo_y);
      
      models[mi]->set_position_at_index(cell, newx, newy);

      log_prior_ratio = 0.0;
      
      log_proposal_ratio = position_priors[mi]->log_proposal_ratio(random,
								   temperature,
								   undo_x,
								   undo_y,
								   newx,
								   newy);
      
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
  }
  
  void reject(std::vector<cartesianvoronoimodel *> &models,
	      hierarchical_model &hierarchical)
  {
    if (undo_cell == -1) {
      throw GENERALVORONOICARTESIANEXCEPTION("No undo information\n");
    }

    models[undo_mi]->set_position_at_index(undo_cell, undo_x, undo_y);

    double t = models[undo_mi]->value_at_index(undo_cell);
    if (t != undo_v) {
      throw GENERALVORONOICARTESIANEXCEPTION("Mismatch: %f %f", t, undo_v);
    }

    undo_mi = -1;
    undo_cell = -1;

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
    return "Move";
  }
  
private:
  
  int undo_mi;
  int undo_cell;
  double undo_x;
  double undo_y;
  double undo_v;

  int nmodels;
  int *p;
  int *a;
  
};

#endif // move_hpp
