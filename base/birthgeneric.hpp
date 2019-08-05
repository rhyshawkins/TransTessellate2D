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
#ifndef birthgeneric_hpp
#define birthgeneric_hpp

#include "global.hpp"
#include "perturbation.hpp"

class BirthGeneric : public Perturbation {
public:

  typedef deltaVoronoi delta_t;
  typedef model_deltaVoronoi model_delta_t;

  BirthGeneric(int _nmodels,
	       std::vector<Proposal*> &_value_proposals,
	       std::vector<PositionProposal*> &_position_proposals) :
    value_proposals(_value_proposals),
    position_proposals(_position_proposals),
    undo_available(false),
    nmodels(_nmodels),
    p(new int[_nmodels]),
    a(new int[_nmodels])
  {
    for (int i = 0; i < nmodels; i ++) {
      p[i] = 0;
      a[i] = 0;
    }
  }
  
  ~BirthGeneric()
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
    double oldvalue;
    double newvalue;
    double newx;
    double newy;

    int mi = 0;

    relocate = false;
    
    if (undo_available) {
      throw GENERALVORONOICARTESIANEXCEPTION("Proposal in progress\n");
    }

    if (this->primary()) {
      
      if (nmodels > 1) {
	mi = random.uniform(nmodels);
      }
    
      p[mi] ++;

      int k = models[mi]->nvaluecells();
    
      if (k == maxcells) {

	perturbation = model_deltaVoronoi::mkbirth(mi,
						   0.0, 0.0, 
						   0.0);
	validproposal = false;

      } else {

	//
	// Generate new point
	//

	if (!position_proposals[mi]->propose(random,
					     temperature,
					     0.0,
					     0.0,
					     newx,
					     newy,
					     log_prior_ratio)) {
	  perturbation = model_deltaVoronoi::mkbirth(mi,
						     0.0, 0.0,
						     0.0);
	  validproposal = false;
	} else {
	  

	  //
	  // Get existing value at point
	  //
	  int t0 = 0;
	  oldvalue = models[mi]->value_at_point(newx, newy, t0);
	  
	  //
	  // Propose new value
	  //
	  if (!value_proposals[mi]->propose(random,
					    temperature,
					    oldvalue,
					    newvalue,
					    log_prior_ratio)) {
	    perturbation = model_deltaVoronoi::mkbirth(mi,
						       0.0, 0.0,
						       0.0);
	    validproposal = false;
	  } else {
	    
	    perturbation = model_deltaVoronoi::mkbirth(mi, newx, newy,
						       newvalue);
	    
	    validproposal = true;
	    
	  }
	}
      }
    }

    this->communicate(validproposal);

    if (validproposal) {

      relocate = true;
      
      this->communicate(mi);
      this->communicate(newx);
      this->communicate(newy);
      this->communicate(oldvalue);
      this->communicate(newvalue);

      log_prior_ratio =
	position_priors[mi]->logpdf(newx, newy) +
	priors[mi]->logpdf(newvalue);

      log_proposal_ratio =
	-position_proposals[mi]->log_proposal(random,
					      temperature,
					      0.0,
					      0.0,
					      newx,
					      newy)
	-value_proposals[mi]->log_proposal(random,
					   temperature,
					   oldvalue,
					   newvalue);      
      
      //
      // Proposal valid
      //
      models[mi]->add_cell(newx, newy, newvalue);

      undo_available = true;
      undo_mi = mi;

      log_extra = 0.0;
    }

    return validproposal;
  }

  void accept()
  {
    if (!undo_available) {
      throw GENERALVORONOICARTESIANEXCEPTION("No proposal in progress\n");
    }
    
    a[undo_mi] ++;
    undo_available = false;
  }

  void reject(std::vector<cartesianvoronoimodel*> &models,
	      hierarchical_model &hierarchical)
  {
    if (!undo_available) {
      throw GENERALVORONOICARTESIANEXCEPTION("No proposal in progress\n");
    }

    models[undo_mi]->pop();
    
    undo_available = false;
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
    return "Birth";
  }

private:

  std::vector<Proposal *> &value_proposals;
  std::vector<PositionProposal *> &position_proposals;

  bool undo_available;
  int undo_mi;

  int nmodels;
  int* p;
  int* a;

};

#endif // birthgenerics2_hpp
