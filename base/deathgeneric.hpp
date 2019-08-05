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
#ifndef deathgeneric_hpp
#define deathgeneric_hpp

#include "global.hpp"
#include "perturbation.hpp"

class DeathGeneric : public Perturbation {
public:

  typedef deltaVoronoi delta_t;
  typedef model_deltaVoronoi model_delta_t;
  
  DeathGeneric(int _nmodels,
	       std::vector<Proposal *> &_value_proposals,
	       std::vector<PositionProposal *> &_position_proposals) :
    value_proposals(_value_proposals),
    position_proposals(_position_proposals),
    undo_index(-1),
    nmodels(_nmodels),
    p(new int[_nmodels]),
    a(new int[_nmodels])
  {
    for (int i = 0; i < nmodels; i ++) {
      p[i] = 0;
      a[i] = 0;
    }
  }

  ~DeathGeneric()
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
    int cell;
    
    relocate = false;
    
    if (this->primary()) {
      
      if (nmodels > 1) {
	mi = random.uniform(nmodels);
      }
    
      p[mi] ++;

      int k = models[mi]->nmobilecells();
      if (k > 1) {
	
	cell = models[mi]->selectmobilecell(random);
	validproposal = true;
	perturbation = model_deltaVoronoi::mkdeath(mi, cell);

      } else {
	validproposal = false;
	perturbation = model_deltaVoronoi::mkdeath(mi, -1);
      }
    }

    this->communicate(validproposal);

    if (validproposal) {

      relocate = true;
      
      this->communicate(mi);
      this->communicate(cell);

      //
      // Store undo information
      //
      undo_mi = mi;
      undo_index = cell;
      undo_value = models[mi]->value_at_index(cell);
      models[mi]->position_at_index(cell, &undo_x, &undo_y);


      //
      // Compute ratios
      //
      log_prior_ratio = -(priors[mi]->logpdf(undo_value) + position_priors[mi]->logpdf(undo_x,
										       undo_y));
      models[mi]->delete_cell(cell);

      int t0 = 0;
      double new_value = models[mi]->value_at_point(undo_x, undo_y, t0);
      
      log_proposal_ratio =
	position_proposals[mi]->log_proposal(random,
					     temperature,
					     0.0,
					     0.0,
					     undo_x,
					     undo_y) +
	value_proposals[mi]->log_proposal(random,
					  temperature,
					  new_value,
					  undo_value);

      log_extra = 0.0;
    }

    return validproposal;
  }

  void accept()
  {
    if (undo_index < 0) {
      throw GENERALVORONOICARTESIANEXCEPTION("No undo information\n");
    }
    
    a[undo_mi] ++;
    
    undo_index = -1;
  }

  void reject(std::vector<cartesianvoronoimodel*> &models,
	      hierarchical_model &hierarchical)
  {
    if (undo_index < 0) {
      throw GENERALVORONOICARTESIANEXCEPTION("No undo information\n");
    }

    //
    // Re-insert the deleted node
    //
    models[undo_mi]->insert_cell(undo_index, undo_x, undo_y, undo_value);

    undo_index = -1;
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
    return "Death";
  }

private:

  std::vector<Proposal *> &value_proposals;
  std::vector<PositionProposal *> &position_proposals;

  int undo_mi;
  int undo_index;
  double undo_value;
  double undo_x;
  double undo_y;

  int nmodels;
  int *p;
  int *a;

};

#endif // deathgeneric_hpp
