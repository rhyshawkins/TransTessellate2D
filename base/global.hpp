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
#ifndef globalspherical_hpp
#define globalspherical_hpp

#include <mpi.h>
#include <math.h>
#include <string.h>

#include "generalvoronoicartesianexception.hpp"

#include "cartesianvoronoimodel.hpp"

#include "generalvoronoicartesianobservations.hpp"
#include "generalvoronoicartesianutil.hpp"
#include "prior.hpp"
#include "positionprior.hpp"

#include "hierarchical_model.hpp"

#include "genericinterface.hpp"

#include "simulated_annealing_scales.hpp"

extern "C" {
  #include "slog.h"
};

class global {
public:

  global(int _nmodels,
	 int _nhierarchical,
	 double _xmin,
	 double _xmax,
	 double _ymin,
	 double _ymax,
	 std::vector<int> &_parameterization_type,
	 const char *input,
	 std::vector<const char *> &initial_models,
	 std::vector<const char *> &prior_files,
	 std::vector<const char *> &position_prior_files,
	 std::vector<const char *> &hierarchical_prior_files,
	 std::vector<const char *> &birthdeath_proposal_files,
	 int _maxcells,
	 std::vector<double> _lambdas,
	 double _temperature,
	 int seed,
	 bool _posterior) :
    nmodels(_nmodels),
    xmin(_xmin),
    xmax(_xmax),
    ymin(_ymin),
    ymax(_ymax),
    communicator(MPI_COMM_NULL),
    mpi_rank(-1),
    mpi_size(-1),
    mpi_counts(nullptr),
    mpi_offsets(nullptr),
    data(nmodels),
    hierarchical(new hierarchical_model(_nhierarchical)),
    temperature(_temperature),
    residual_size(0),
    mean_residual_n(0),
    maxcells(_maxcells),
    random(seed),
    posterior(_posterior)
  {
    //
    // Load value priors (required)
    //
    if ((int)prior_files.size() != nmodels) {
      throw GENERALVORONOICARTESIANEXCEPTION("Incorrect number of value prior files, require %d, have %d",
					     nmodels,
					     (int)prior_files.size());
    } else {
      for (auto &filename : prior_files) {

	PriorProposal *pp = PriorProposal::load(filename);
	if (pp == nullptr) {
	  throw GENERALVORONOICARTESIANEXCEPTION("Failed to load prior from %s", filename);
	}

	priors.push_back(pp);
      }
    }

    //
    // Load position priors (required)
    //
    if ((int)position_prior_files.size() != nmodels) {
      throw GENERALVORONOICARTESIANEXCEPTION("Incorrect number of position prior files");
    } else {
      for (auto &filename : position_prior_files) {

	PositionPriorProposal *pp = PositionPriorProposal::load(filename, xmin, xmax, ymin, ymax);
	if (pp == nullptr) {
	  throw GENERALVORONOICARTESIANEXCEPTION("Failed to load position prior from %s", filename);
	}

	position_priors.push_back(pp);
      }
    }

    //
    // Load hierarchical priors (optional)
    //
    for (int i = 0; i < _nhierarchical; i ++) {

      if (i < (int)hierarchical_prior_files.size()) {
	PriorProposal *pp = PriorProposal::load(hierarchical_prior_files[i]);
	if (pp == nullptr) {
	  throw GENERALVORONOICARTESIANEXCEPTION("Failed to load prior from %s", hierarchical_prior_files[i]);
	}
	hierarchical_priors.push_back(pp);
      } else {
	Prior *p = new UniformPrior(0.1, 5.0);
	Proposal *pp = new GaussianProposal(*p, 0.05);
	hierarchical_priors.push_back(new PriorProposal(p, pp));
      }

      if (i < (int)_lambdas.size()) {
	hierarchical->set(i, _lambdas[i]);
      } else {
	hierarchical->set(i, 1.0);
      }
    }

    //
    // Load birth/death proposal file (optional)
    //
    for (int i = 0; i < nmodels; i ++) {

      if (i < (int)birthdeath_proposal_files.size()) {
	FILE *fp = fopen(birthdeath_proposal_files[i], "r");
	if (fp == NULL) {
	  throw GENERALVORONOICARTESIANEXCEPTION("Failed to open birth/death proposal file %s",
						 birthdeath_proposal_files[i]);
	}

	Proposal *p = Proposal::load(fp, *(priors[i]->get_prior()));
	if (p == nullptr) {
	  throw GENERALVORONOICARTESIANEXCEPTION("Failed to read birth/death proposal file %s",
						 birthdeath_proposal_files[i]);
	}

	PositionProposal *pp = PositionProposal::load(fp, *(position_priors[i]->get_prior()));
	if (pp == nullptr) {
	  throw GENERALVORONOICARTESIANEXCEPTION("Failed to read birth/death proposal file %s",
						 birthdeath_proposal_files[i]);
	}
	fclose(fp);

	birthdeath_proposals.push_back(p);
	birthdeath_position_proposals.push_back(pp);
	
      } else {
	Proposal *p = new PriorSampleProposal(*(priors[i]->get_prior()));
	birthdeath_proposals.push_back(p);

	PositionProposal *pp = new PriorSamplePositionProposal(*(position_priors[i]->get_prior()));
	birthdeath_position_proposals.push_back(pp);
      }
    }
    
    if (!posterior) {

      current_state = this;
      int n = strlen(input);
      if (gvcart_loaddata_(&n, input, addobservation) < 0) {
	throw GENERALVORONOICARTESIANEXCEPTION("Failed to load observations from %s\n", input);
      }
      current_state = nullptr;
      
      residual_size = data.obs.size();
      predictions.resize(residual_size);
      residuals.resize(residual_size);
      mean_residuals.resize(residual_size);
      last_valid_residuals.resize(residual_size);
      mean_residual_n = 0;
      for (int i = 0; i < residual_size; i ++) {
	predictions[i] = 0.0;
	residuals[i] = 0.0;
	mean_residuals[i] = 0.0;
	last_valid_residuals[i] = 0.0;
      }

      modelweights.resize(data.obs.size());

    } else {
      
    }

    for (int i = 0; i < nmodels; i ++) {
      dLdm.push_back(new double[_maxcells + 4]);
    }
      
    for (int i = 0; i < nmodels; i ++) {

      //
      // Default to Voronoi cell but allow user to
      // specify each parameterization independently
      //
      int pt = 0;

      if ((int)_parameterization_type.size() > i) {
	pt = _parameterization_type[i];
      }
	
      models.push_back(new cartesianvoronoimodel(maxcells,
						 xmin, xmax,
						 ymin, ymax,
						 pt));

    }
  }

  void initialize_mpi(MPI_Comm _communicator, double _temperature)
  {
    MPI_Comm_dup(_communicator, &communicator);

    MPI_Comm_rank(communicator, &mpi_rank);
    MPI_Comm_size(communicator, &mpi_size);

    temperature = _temperature;

    int observations = (int)data.obs.size();
    int processes = mpi_size;

    mpi_offsets = new int[mpi_size];
    mpi_counts = new int[mpi_size];
      
    for (int i = 0; i < mpi_size; i ++) {
      mpi_counts[i] = observations/processes;
      observations -= mpi_counts[i];
      processes --;
    }

    mpi_offsets[0] = 0;
    for (int i = 1; i < mpi_size; i ++) {
      mpi_offsets[i] = mpi_offsets[i - 1] + mpi_counts[i - 1];
    }

    for (int i = 0; i < mpi_size; i ++) {
      INFO("MPI %2d %6d %6d", i, mpi_offsets[i], mpi_counts[i]);
    }

    if (mpi_offsets[mpi_size - 1] + mpi_counts[mpi_size - 1] != (int)data.obs.size()) {
      throw GENERALVORONOICARTESIANEXCEPTION("Failed to distribute data points properly");
    }
  }

  void initialize(int initialcells, std::vector<const char *> &initial_models)
  {
    for (int mi = 0; mi < nmodels; mi ++) {

      if (mi < (int)initial_models.size()) {
	throw GENERALVORONOICARTESIANEXCEPTION("Unimplemented");
      } else {

	//
	// Random prior sampled model
	//

	for (int j = 0; j < initialcells; j ++) {

	  double xyz[3];

	  if (communicator == MPI_COMM_NULL || mpi_rank == 0) {
	    xyz[2] = priors[mi]->sample(random);
	    position_priors[mi]->sample(random, xyz[0], xyz[1]);
	  }

	  if (communicator != MPI_COMM_NULL) {
	    MPI_Bcast(xyz, 3, MPI_DOUBLE, 0, communicator);
	  }

	  models[mi]->add_cell(xyz[0],
			       xyz[1],
			       xyz[2]);
	  // INFO("Add %2d : %f %f %f", mi, xyz[0], xyz[1], xyz[2]);
	}
      }
    }

    data.initialize_model_references(models);
  }

  double optimize_sa(int iterations, double Tmax, double &acceptance, int verbose)
  {
    double scale_factor = 1.0e-3;
    
    if (communicator == MPI_COMM_NULL) {

      //
      // Initial likelihood
      //
      double current_N = 0.0;
      double current_P = likelihood(current_N, true);
      
      //
      // Prune unused model parameters
      //
      if (not posterior) {
	for (int mi = 0; mi < nmodels; mi ++) {
	  std::vector<bool> used(models[mi]->ntotalcells(), false);
	  for (auto &o : data.obs) {
	    
	    for (int pi = 0; pi < (int)o.xs.size(); pi ++) {
	    
	      if (o.mi[pi] == mi) {
		int nc;
		const int *indices = o.refs[pi]->cached_indices(nc);
		if (nc > 0) {
		  for (int j = 0; j < nc; j ++) {
		    used[indices[j]] = true;
		  }
		  
		} else {
		  printf("warning: no cached indices\n");
		}
	      }
	    }
	  }
	  
	  int pruned_count = 0;
	  int minpoints = 4;
	  if (models[mi]->type == cartesianvoronoimodel::VORONOI) {
	    minpoints = 5;
	  }
	  for (int j = used.size() - 1; j >= minpoints; j --) {
	    if (!used[j]) {
	      models[mi]->delete_cell(j);
	      pruned_count ++;
	    }
	  }
	  
	  if (pruned_count > 0) {
	    printf("%2d Pre  prune: %12.6f %6d %6d\n", mi, current_P, (int)used.size(), pruned_count);
	    current_P = likelihood(current_N, true);
	    printf("%2d Post prune: %12.6f %6d %d\n", mi, current_P,
		   (int)models[mi]->ntotalcells(), pruned_count);
	  }
	}
      }
      
      std::vector<std::vector<double>> scales(nmodels);
      std::vector<std::vector<int>> local_proposed(nmodels);
      std::vector<std::vector<int>> local_accepted(nmodels);
      
      for (int mi = 0; mi < nmodels; mi ++) {
	
	int nv = models[mi]->nvaluecells();
	scales[mi].resize(nv, 1.0);
	local_proposed[mi].resize(nv, 0);
	local_accepted[mi].resize(nv, 0);
	
      }
      
      int proposed = 0;
      int accepted = 0;
      
      for (int i = 0; i < iterations; i ++) {
        double T = temperature_power(i, iterations, Tmax, 2.0);
	
	for (int mi = 0; mi < nmodels; mi ++) {
	  //
	  // Perturb model parameters
	  //
	  // If the perturbation violates prior bounds, we just keep going hoping
	  // that the next perturbation won't. The number of prior violations is
	  // counted which can be used as a diagnostic. (Proposal width is too large)
	  //
	  int nv = models[mi]->nvaluecells();
	  for (int j = 0; j < nv; j ++)  {
	    
	    int index = models[mi]->selectvaluecell(j);
	    double oldv = models[mi]->value_at_index(index);
	    double newv;
	    double logpriorratio;
	    
	    while (!priors[mi]->propose(random,
					T * scales[mi][j],
					oldv,
					newv,
					logpriorratio)) {
	    }
	    
	    local_proposed[mi][j] ++;
	    proposed ++;
	    
	    models[mi]->set_value_at_index(index, newv);
	    
	    //
	    // Compute perturbed likelihood
	    //
	    double proposed_N = 0.0;
	    double proposed_P = likelihood(proposed_N, false);
	    
	    //
	    // Accept/reject
	    //
	    double u = ::log(random.uniform());
	    if (u < (current_P - proposed_P)/T) {
	      //
	      // Accept: set likelihood and keep model as is
	      //
	      current_P = proposed_P;
	      local_accepted[mi][j] ++;
	      accepted ++;
	      
	    } else {
	      //
	      // Reject: restore old model
	      //
	      models[mi]->set_value_at_index(index, oldv);
	    }
            
	    //
	    // Update adaptive step size
	    //
	    if (local_proposed[mi][j] >= 10 && local_proposed[mi][j] % 10 == 0) {
	      
	      double local_acceptance = (double)local_accepted[mi][j]/(double)local_proposed[mi][j];
	      if (local_acceptance < 0.5) {
		scales[mi][j] *= (1.0 - scale_factor);
	      } else {
		scales[mi][j] *= (1.0 + scale_factor);
	      }
	    }
	    
	  }
        }
        
        if (proposed > 0) {
          acceptance = (double)accepted/(double)proposed;
        }
        
        if (verbose > 0 && (i + 1) % verbose == 0) {
          INFO("%6d %12.6f %10.6f (%8d/%8d)\n",
               i + 1,
               current_P,
               acceptance,
               accepted, proposed);
	  
        }
        
      }

      //
      // Compute mean proposal width for each model
      //
      for (int mi = 0; mi < nmodels; mi ++) {
	if (scales[mi].size() > 0) {
	  double t = 0.0;
	  for (auto &s : scales[mi]) {
	    t += s;
	  }
	  
	  t /= (double)scales[mi].size();
	  INFO("Model %d Mean %10.6f", mi, t);

	}
      }
      
      return current_P;

    } else {
      //
      // Initial likelihood
      //
      double current_N = 0.0;
      double current_P = likelihood(current_N, true);      

      //
      // Prune unused model parameters
      //
      for (int mi = 0; mi < nmodels; mi ++) {
	std::vector<int> used(models[mi]->ntotalcells(), 0);
	for (int oi = 0; oi < mpi_counts[mpi_rank]; oi ++) {

	  auto &o = data.obs[mpi_offsets[mpi_rank] + oi];
	  
	  for (int pi = 0; pi < (int)o.xs.size(); pi ++) {
	    
	    if (o.mi[pi] == mi) {
	      int nc;
	      const int *indices = o.refs[pi]->cached_indices(nc);
	      if (nc > 0) {
		for (int j = 0; j < nc; j ++) {
		  used[indices[j]] = 1;
		}
	      } else {
		printf("warning: no cached indices\n");
	      }
	    }
	  }
	}

	if (mpi_rank == 0) {
	  MPI_Reduce(MPI_IN_PLACE, used.data(), used.size(), MPI_INT, MPI_SUM, 0, communicator);
	} else {
	  MPI_Reduce(used.data(), NULL, used.size(), MPI_INT, MPI_SUM, 0, communicator);
	}
	MPI_Bcast(used.data(), used.size(), MPI_INT, 0, communicator);

	int pruned_count = 0;
	for (int j = used.size() - 1; j >= 4; j --) {
	  if (used[j] == 0) {
	    models[mi]->delete_cell(j);
	    pruned_count ++;
	  }
	}
	
	if (pruned_count > 0) {
	  if (mpi_rank == 0) {
	    INFO("%2d Pre  prune: %12.6f %6d %6d\n", mi, current_P, (int)used.size(), pruned_count);
	  }

	  current_P = likelihood(current_N, true);

	  if (mpi_rank == 0) {
	    INFO("%2d Post prune: %12.6f %6d %6d\n", mi, current_P,
		 (int)models[mi]->ntotalcells(), pruned_count);
	  }
	}
      }

      std::vector<std::vector<double>> scales(nmodels);
      std::vector<std::vector<int>> local_proposed(nmodels);
      std::vector<std::vector<int>> local_accepted(nmodels);
      
      for (int mi = 0; mi < nmodels; mi ++) {
	
	int nv = models[mi]->nvaluecells();
	scales[mi].resize(nv, 1.0);
	local_proposed[mi].resize(nv, 0);
	local_accepted[mi].resize(nv, 0);
	
      }
      
      int proposed = 0;
      int accepted = 0;
      
      for (int i = 0; i < iterations; i ++) {

        double T = temperature_power(i, iterations, Tmax, 2.0);

	for (int mi = 0; mi < nmodels; mi ++) {

	  int nv = models[mi]->nvaluecells();
	  for (int j = 0; j < nv; j ++)  {
	  
	    int index = models[mi]->selectvaluecell(j);
	    double oldv = models[mi]->value_at_index(index);
	    double newv;
	    double logpriorratio;
	    
	    if (mpi_rank == 0) {
	      while (!priors[mi]->propose(random,
					  T * scales[mi][j],
					  oldv,
					  newv,
					  logpriorratio)) {
		// Keep proposing until we get a valid proposal
	      }
                                           
	    }

	    MPI_Bcast(&newv, 1, MPI_DOUBLE, 0, communicator);

	    proposed ++;
	    local_proposed[mi][j] ++;
	    models[mi]->set_value_at_index(index, newv);
	    
	    //
	    // Compute perturbed likelihood
	    //
	    double proposed_N = 0.0;
	    double proposed_P = likelihood(proposed_N, false);
	    
	    int accept;

	    if (mpi_rank == 0) {
	      double u = log(random.uniform());
	      if (u < (current_P - proposed_P)/T) {
		accept = 1;
	      } else {
		accept = 0;
	      }
	    }

	    MPI_Bcast(&accept, 1, MPI_INT, 0, communicator);
	    
	    if (accept == 1) {
	      current_P = proposed_P;
	      accepted ++;
	      local_accepted[mi][j] ++;
	      
	    } else {
	      //
	      // Reject: restore old model
	      //
	      models[mi]->set_value_at_index(index, oldv);
	    }
	    
	    if (mpi_rank == 0) {
	      //
	      // Update adaptive step size
	      //
	      if (local_proposed[mi][j] >= 10 && local_proposed[mi][j] % 10 == 0) {
		
		double local_acceptance = (double)local_accepted[mi][j]/(double)local_proposed[mi][j];
		if (local_acceptance < 0.5) {
		  scales[mi][j] *= (1.0 - scale_factor);
		} else {
		  scales[mi][j] *= (1.0 + scale_factor);
		}
	      }
	      
	    }
	  }	    
        }

        if (proposed > 0) {
          acceptance = (double)accepted/(double)proposed;
        }
        
        if (mpi_rank == 0 && verbose > 0 && (i + 1) % verbose == 0) {
          INFO("%6d %12.6f %10.6f (%8d/%8d)\n",
               i + 1,
               current_P,
               acceptance,
               accepted, proposed);
          
        }
      } 

      return current_P;
    }
  }
  
  double likelihood(double &norm, bool relocate)
  {
    if (data.obs.size() > 0) {
      if (communicator == MPI_COMM_NULL) {

	if (relocate) {
	  for (auto &m : models) {
	    m->recompute();
	  }
	}

	for (int i = 0; i < (int)data.obs.size(); i ++) {
	  int npoints = data.obs[i].xs.size();

	  //
	  // Look up points
	  //
	  for (int j = 0; j < npoints; j ++) {
	    data.obs[i].values[j] = data.obs[i].refs[j]->evaluate(relocate);
	  }

	  //
	  // Custom forwardmodel
	  //
	  if (gvcart_compute_prediction_(&nmodels,
					 &i,
					 &npoints,
					 data.obs[i].values.data(),
					 data.obs[i].weights.data(),
					 &predictions[i]) < 0) {
	    throw GENERALVORONOICARTESIANEXCEPTION("Failed to compute predictions");
	  }
	}

	//
	// Likelihood
	//
	int nobs = data.obs.size();
	double like;
	int nhierarchical = hierarchical->get_nhierarchical();
	
	if (gvcart_compute_likelihood_(&nmodels,
				       &nhierarchical,
				       &nobs,
				       hierarchical->data(),
				       predictions.data(),
				       residuals.data(),
				       modelweights.data(),
				       &like,
				       &norm) < 0) {
	  throw GENERALVORONOICARTESIANEXCEPTION("Failed to compute likelihood");
	}

	return like;
      } else {

	if (relocate) {
	  for (auto &m : models) {
	    m->recompute();
	  }
	}

	for (int i = 0; i < mpi_counts[mpi_rank]; i ++) {
	  int o = mpi_offsets[mpi_rank] + i;
	  int npoints = data.obs[o].xs.size();

	  //
	  // Look up points
	  //
	  for (int j = 0; j < npoints; j ++) {
	    data.obs[o].values[j] = data.obs[o].refs[j]->evaluate(relocate);
	  }

	  //
	  // Custom forwardmodel
	  //
	  
	  if (gvcart_compute_prediction_(&nmodels,
					 &o,
					 &npoints,
					 data.obs[o].values.data(),
					 data.obs[o].weights.data(),
					 &predictions[o]) < 0) {
	    throw GENERALVORONOICARTESIANEXCEPTION("Failed to compute predictions");
	  }
	}
	  
	//
	// Share predictions
	//
	if (mpi_size > 1) {
	  MPI_Gatherv(predictions.data() + mpi_offsets[mpi_rank],
		      mpi_counts[mpi_rank],
		      MPI_DOUBLE,
		      predictions.data(),
		      mpi_counts,
		      mpi_offsets,
		      MPI_DOUBLE,
		      0,
		      communicator);
	}

	
	int nobs = data.obs.size();
	double like;
	int nhierarchical = hierarchical->get_nhierarchical();
	
	if (gvcart_compute_likelihood_(&nmodels,
				       &nhierarchical,
				       &nobs,
				       hierarchical->data(),
				       predictions.data(),
				       residuals.data(),
				       modelweights.data(),
				       &like,
				       &norm) < 0) {
	  throw GENERALVORONOICARTESIANEXCEPTION("Failed to compute likelihood");
	}

	MPI_Bcast(&like, 1, MPI_DOUBLE, 0, communicator);
	MPI_Bcast(&norm, 1, MPI_DOUBLE, 0, communicator);

	return like;
	
      }
    } else {
      norm = 0.0;
      return 1.0;
    }
  }

  void likelihood_backproject()
  {
    if (data.obs.size() > 0) {
      if (communicator == MPI_COMM_NULL) {

	//
	// Zero
	//
	for (int mi = 0; mi < nmodels; mi ++) {
	  int n = models[mi]->nvaluecells();
	  for (int j = 0; j < n; j ++) {
	    dLdm[mi][j] = 0.0;
	  }
	}

	//
	// Loop through by observation
	//
	for (int i = 0; i < (int)data.obs.size(); i ++) {
	  int npoints = data.obs[i].xs.size();

	  double rw = modelweights[i];
	  
	  for (int j = 0; j < npoints; j ++) {

	    int nc;
	    const int *indices = data.obs[i].refs[j]->cached_indices(nc);
	    if (nc < 0) {
	      throw GENERALVORONOICARTESIANEXCEPTION("No cached indices in backproject");
	    }
	    const double *weights = data.obs[i].refs[j]->cached_weights(nc);
	    if (nc < 0) {
	      throw GENERALVORONOICARTESIANEXCEPTION("No cached weights in backproject");
	    }
	      
	    int mi = data.obs[i].mi[j];
	    
	    for (int l = 0; l < nc; l ++) {
	      // printf("%2d %3d %3d %3d %16.9e %16.9e %16.9e\n",
	      // 	     mi,
	      // 	     i,
	      // 	     j,
	      // 	     indices[l],
	      // 	     rw, data.obs[i].weights[j], weights[l]);
	      int ci = models[mi]->ordinal(indices[l]);
	      
	      dLdm[mi][ci] += rw * data.obs[i].weights[j] * weights[l];
	    }
	    
	  }
	}

      } else {

	//
	// Modelweights are only available on root node so Broadcast
	//
	MPI_Bcast(modelweights.data(), data.obs.size(), MPI_DOUBLE, 0, communicator);

	//
	// Zero
	//
	for (int mi = 0; mi < nmodels; mi ++) {
	  int n = models[mi]->nvaluecells();
	  for (int j = 0; j < n; j ++) {
	    dLdm[mi][j] = 0.0;
	  }
	}

	//
	// Loop through by observation
	//
	for (int i = 0; i < mpi_counts[mpi_rank]; i ++) {
	  int oi = mpi_offsets[mpi_rank];
	  int npoints = data.obs[oi].xs.size();

	  double rw = modelweights[oi];
	  
	  for (int j = 0; j < npoints; j ++) {

	    int nc;
	    const int *indices = data.obs[oi].refs[j]->cached_indices(nc);
	    if (nc < 0) {
	      throw GENERALVORONOICARTESIANEXCEPTION("No cached indices in backproject");
	    }
	    const double *weights = data.obs[oi].refs[j]->cached_weights(nc);
	    if (nc < 0) {
	      throw GENERALVORONOICARTESIANEXCEPTION("No cached weights in backproject");
	    }
	      
	    int mi = data.obs[oi].mi[j];
	    
	    for (int l = 0; l < nc; l ++) {
	      // printf("%2d %3d %3d %3d %16.9e %16.9e %16.9e\n",
	      // 	     mi,
	      // 	     i,
	      // 	     j,
	      // 	     indices[l],
	      // 	     rw, data.obs[i].weights[j], weights[l]);
	      int ci = models[mi]->ordinal(indices[l]);
	      
	      dLdm[mi][ci] += rw * data.obs[oi].weights[j] * weights[l];
	    }
	    
	  }
	}

	//
	// Reduce Gradient array(s)
	//
	for (int mi = 0; mi < nmodels; mi ++) {
	  MPI_Allreduce(MPI_IN_PLACE, dLdm[mi], models[mi]->nvaluecells(),
			MPI_DOUBLE, MPI_SUM, communicator);
	}
      }

    } else {
      //
      // Posterior (no data) testing -> Zero gradient
      //
      for (int mi = 0; mi < nmodels; mi ++) {
	int n = models[mi]->nvaluecells();
	for (int j = 0; j < n; j ++) {
	  dLdm[mi][j] = 0.0;
	}
      }
    }
  
  }

  void accept()
  {
    for (int i = 0; i < residual_size; i ++) {
      last_valid_residuals[i] = residuals[i];
    }
    update_mean_residual();
  }

  void reject()
  {
    update_mean_residual();
  }

  void update_mean_residual()
  {
    mean_residual_n ++;

    for (int i = 0; i < residual_size; i ++) {
      double delta = last_valid_residuals[i] - mean_residuals[i];
      mean_residuals[i] += delta/(double)mean_residual_n;
    }
  }

  int encode(char *buffer, int length)
  {
    int o = 0;

    if (!::encode<int>(models.size(), buffer, o, length)) {
      return -1;
    }

    for (auto &m : models) {

      int n = m->ntotalcells();
      if (!::encode<int>(n, buffer, o, length)) {
	return -1;
      }

      for (int i = 0; i < n; i ++) {
	double x, y, z;
	
	z = m->value_at_index(i);
	m->position_at_index(i, &x, &y);

	if (!::encode<double>(x, buffer, o, length) ||
	    !::encode<double>(y, buffer, o, length) ||
	    !::encode<double>(z, buffer, o, length)) {
	  return -1;
	}
      }
    }

    int nh = hierarchical->get_nhierarchical();
    if (!::encode<int>(nh, buffer, o, length)) {
      return -1;
    }

    for (int i = 0; i < nh; i ++) {
      double h = hierarchical->get(i);
      if (!::encode<double>(h, buffer, o, length)) {
	return -1;
      }
    }

    return o;
  }

  int decode(const char *buffer, int length)
  {
    int o = 0;

    int nm;
    if (!::decode<int>(nm, buffer, o, length)) {
      ERROR("Failed to decode nmodels %d", length);
      return -1;
    }

    if (nm != (int)models.size()) {
      ERROR("nmodels mismatch %d %d", nm, (int)models.size());
      return -1;
    }

    for (auto &m : models) {

      m->reset();
      
      int n;
      if (!::decode<int>(n, buffer, o, length)) {
	ERROR("Failed to decode ncells");
	return -1;
      }

      for (int j = 0; j < n; j ++) {

	double x, y, z;
	if (!::decode<double>(x, buffer, o, length) ||
	    !::decode<double>(y, buffer, o, length) ||
	    !::decode<double>(z, buffer, o, length)) {
	  ERROR("Failed to decode point");
	  return -1;
	}
	  
	if (j < 4) {
	  m->set_value_at_index(j, z);
	} else {
	  m->add_cell(x, y, z);
	}
      }
    }
    
    int nh;
    if (!::decode<int>(nh, buffer, o, length)) {
      ERROR("Failed to decode nhierarchical");
      return -1;
    }

    if (nh != hierarchical->get_nhierarchical()) {
      ERROR("nhierarchical mismatch");
      return -1;
    }

    for (int i = 0; i < nh; i ++) {
      double h;
      if (!::decode<double>(h, buffer, o, length)) {
	ERROR("Failed to decode hierarchical value");
	return -1;
      }

      hierarchical->set(i, h);   
    }
    
    return o;
  }

  static global *current_state;
  static int addobservation(int *npoints, int *modelindices, double *lons, double *lats);

  int nmodels;
  int nhierarchical;
  
  double xmin;
  double xmax;
  double ymin;
  double ymax;

  MPI_Comm communicator;
  int mpi_rank;
  int mpi_size;
  int *mpi_counts;
  int *mpi_offsets;

  GeneralVoronoiCartesianObservations data;
  std::vector<cartesianvoronoimodel *> models;

  std::vector<PriorProposal *> priors;
  std::vector<PositionPriorProposal *> position_priors;
  std::vector<PriorProposal *> hierarchical_priors;

  std::vector<Proposal *> birthdeath_proposals;
  std::vector<PositionProposal *> birthdeath_position_proposals;

  hierarchical_model *hierarchical;
  double temperature;

  int residual_size;
  int mean_residual_n;

  std::vector<double> mean_residuals;
  std::vector<double> predictions;
  std::vector<double> residuals;
  std::vector<double> last_valid_residuals;

  std::vector<double> modelweights;
  std::vector<double*> dLdm;
  
  int maxcells;
  
  Rng random;

  bool posterior;

};

#endif // globalspherical_hpp
