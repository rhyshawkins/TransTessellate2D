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
#ifndef hmc_hpp
#define hmc_hpp

#include "perturbation.hpp"

class HMC : public Perturbation {
public:

  HMC(global &_state, int _L, double _epsilon) :
    state(_state),
    L(_L),
    epsilon(_epsilon),
    p(0),
    a(0),
    last_adjustment(0)
  {
    cells.resize(state.nmodels);
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
    perturbation = nullptr;
    
    //
    // Initialise (count no. cells and allocate structures)
    //
    int total = 0;
    for (int mi = 0; mi < nmodels; mi ++) {
      int n = models[mi]->nvaluecells();
      cells[mi] = n;
      total += n;
    }

    dU.resize(total);
    momentum.resize(total);
    
    undo_v.resize(total);
    
    //
    // Store current model(s) value(s) for undo
    //
    int offset = 0;
    for (int mi = 0; mi < nmodels; mi ++) {
      int n = cells[mi];
      for (int j = 0; j < n; j ++) {
	undo_v[offset + j] = models[mi]->value_at_index(models[mi]->selectvaluecell(j));
      }
      offset += n;
    }

    //
    // Random momentum
    //
    double currentK = 0.0;
    offset = 0;
    for (int mi = 0; mi < nmodels; mi ++) {
      int n = cells[mi];
      if (this->primary()) {
	for (int j = 0; j < n; j ++) {
	  momentum[offset + j] = random.normal(1.0);
	}
      }

      communicate(n, momentum.data() + offset);

      for (int j = 0; j < n; j ++) {
	currentK += momentum[offset + j]*momentum[offset + j];
      }
      
      offset += n;
    }
    currentK /= 2.0;
    

    //
    // Compute initial likelihood and prior
    //
    double logprior;
    double likelihood;
    double norm;
    double U;
    
    if (!computeU(true, priors, logprior, likelihood, norm, temperature, U, dU)) {
      //
      // Prior violated at initial point?
      //
      if (this->primary()) {
	perturbation = new model_deltaValuesVoronoi(models);
      }
    
      fprintf(stderr, "error: prior violated in initial likelihood\n");
      return false;
    }

    // for (auto &du : dU) {
    //   printf("%6.e3 ", du);
    // }
    // printf("\n");

    p ++;

    //
    // Record the starting prior (Note U only does value priors not position)
    //
    log_prior_ratio = -logprior;
    
    //
    // Initial Half step
    //
    offset = 0;
    for (int mi = 0; mi < nmodels; mi ++) {
      int n = cells[mi];
      for (int j = 0; j < n; j ++) {
	momentum[offset + j] -= epsilon * dU[offset + j]/2.0;
      }
      offset += n;
    }

    //
    // Full Iteration steps
    //
    for (int i = 0; i < L; i ++) {

      //
      // Model updates
      //
      offset = 0;
      for (int mi = 0; mi < nmodels; mi ++) {
	int n = cells[mi];
	double pmin, pmax;
	Prior::bound_t bounded = priors[mi]->get_prior()->bounded(pmin, pmax);
	  
	for (int j = 0; j < n; j ++) {
	  int ci = models[mi]->selectvaluecell(j);
	  
	  double v = models[mi]->value_at_index(ci);
	  double newv = v + epsilon*momentum[offset + j];
	  if (momentum[offset + j] < 0.0 &&
	      (bounded == Prior::LOWBOUND || bounded == Prior::BOUNDED)) {

	    if (newv < pmin) {
	      newv = 2.0*pmin - newv;
	      momentum[offset + j] = -momentum[offset + j];
	    }
	  }

	  if (momentum[offset + j] > 0.0 &&
	      (bounded == Prior::HIGHBOUND || bounded == Prior::BOUNDED)) {

	    if (newv > pmax) {
	      newv = 2.0*pmax - newv;
	      momentum[offset + j] = -momentum[offset + j];
	    }
	  }
	  
	  models[mi]->set_value_at_index(ci, newv);
	}
	offset += n;
      }

      //
      // Likelihood
      //
      if (!computeU(false, priors, logprior, likelihood, norm, temperature, U, dU)) {
	// fprintf(stderr, "error: failed to compute U (prior violation)\n");

	if (this->primary()) {
	  perturbation = new model_deltaValuesVoronoi(models);
	}
	
	//
	// Restore model
	//
	offset = 0;
	for (int mi = 0; mi < nmodels; mi ++) {
	  int n = cells[mi];
	  for (int j = 0; j < n; j ++) {
	    int ci = models[mi]->selectvaluecell(j);
	    
	    models[mi]->set_value_at_index(ci, undo_v[offset + j]);
	  }
	  
	  offset += n;
	}

	// fprintf(stderr, "prior violation\n");
	return false;
      }

      if (i == (L - 1)) {
	//
	// Final half step
	//
	offset = 0;
	for (int mi = 0; mi < nmodels; mi ++) {
	  int n = cells[mi];
	  for (int j = 0; j < n; j ++) {
	    momentum[offset + j] -= epsilon * dU[offset + j]/(2.0);
	    momentum[offset + j] = -momentum[offset + j];
	  }
	  offset += n;
	}

      } else {
	//
	// Full step
	//
	offset = 0;
	for (int mi = 0; mi < nmodels; mi ++) {
	  int n = cells[mi];
	  for (int j = 0; j < n; j ++) {
	    momentum[offset + j] -= epsilon * dU[offset + j];
	  }
	  offset += n;
	}
      }
    }

    //
    // Compute momentum
    //
    double proposedK = 0.0;
    offset = 0;
    for (int mi = 0; mi < nmodels; mi ++) {
      int n = cells[mi];
      for (int j = 0; j < n; j ++) {
	proposedK += momentum[offset + j]*momentum[offset + j];
      }
      offset += n;
    }
    proposedK /= 2.0;

    log_prior_ratio += logprior;
    log_proposal_ratio = 0.0;
    log_extra = currentK - proposedK;
    
    relocate = false;
    if (this->primary()) {
      perturbation = new model_deltaValuesVoronoi(models);
    }
      
    return true;
  }

  bool computeU(bool relocate,
		std::vector<PriorProposal*> &priors,
		double &logprior,
		double &likelihood,
		double &norm,
		double temperature,
		double &U,
		std::vector<double> &dU)
  {
    U = 0.0;
    logprior = 0.0;
    for (auto &d: dU) {
      d = 0.0;
    }

    //
    // Compute prior and gradient w.r.t. prior
    //
    int offset = 0;
    for (int mi = 0; mi < state.nmodels; mi ++) {
      int n = state.models[mi]->nvaluecells();

      for (int j = 0; j < n; j ++) {
	int ci = state.models[mi]->selectvaluecell(j);

	double dlp;
	double v = state.models[mi]->value_at_index(ci);
	double lp = priors[mi]->logpdf(v,
				       dlp);
	if (lp == Prior::INVALID_LOGP) {
	  // Prior violation
	  // printf("Prior violation: %16.9e %16.9e\n", v, undo_v[offset + j]);
	  return false;
	}

	U += lp;
	logprior += lp;
	dU[offset + j] = dlp;
      }

      offset += n;
    }

    //
    // Likelihood (give real likelihood back to likelihood and norm variables, but tempered
    // value for U).
    //
    likelihood = state.likelihood(norm, relocate);

    U += (likelihood + norm)/temperature;

    //
    // Likelihood gradient
    //
    state.likelihood_backproject();

    offset = 0;
    for (int mi = 0; mi < state.nmodels; mi ++) {
      int n = state.models[mi]->nvaluecells();

      for (int j = 0; j < n; j ++) {

	dU[offset + j] += state.dLdm[mi][j]/temperature;

      }
      
      offset += n;
    }

    return true;
  }

  virtual void accept()
  {
    a ++;
  }

  virtual void reject(std::vector<cartesianvoronoimodel*> &models,
                      hierarchical_model &hierarchical)
  {
    int offset = 0;
    int nmodels = models.size();
    for (int mi = 0; mi < nmodels; mi ++) {
      int n = cells[mi];
      for (int j = 0; j < n; j ++) {
	int ci = models[mi]->selectvaluecell(j);
	  
	models[mi]->set_value_at_index(ci, undo_v[offset + j]);
      }

      offset += n;
    }
  }

  virtual int pa_categories() const
  {
    return 1;
  }
  
  virtual int proposal_count(int ci) const
  {
    return p;
  }
  
  virtual int acceptance_count(int ci) const
  {
    return a;
  }
  
  virtual const char *displayname() const
  {
    return "HMC";
  }

  void adjust_epsilon(double target_acceptance, double scale, int steps)
  {
    if (p == 0) {
      return;
    }

    if (p - last_adjustment < steps) {
      return;
    }

    
    double acceptance = (double)a/(double)p;

    if (acceptance < target_acceptance) {

      epsilon *= (1.0 - scale);

    } else if (acceptance > target_acceptance) {

      epsilon *= (1.0 + scale);

    }

    last_adjustment = p;
  }

  double get_epsilon() const
  {
    return epsilon;
  }
  

private:

  global &state;

  int L;
  double epsilon;

  int p;
  int a;
  int last_adjustment;
  
  std::vector<int> cells; // No. value cells in each model
  std::vector<double> momentum; // Hamiltonian momentum
  std::vector<double> dU;       // Gradient of likelihood wrt model parameters
  std::vector<double> undo_v;   // Old values for rejection
};


#endif // hmc_hpp
