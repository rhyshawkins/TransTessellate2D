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
#ifndef perturbationcollection_hpp
#define perturbationcollection_hpp

#include "global.hpp"
#include "chainhistoryVoronoi.hpp"
#include "perturbation.hpp"

class PerturbationCollection {
public:

  typedef deltaVoronoi delta_t;
  
  
  PerturbationCollection() :
    communicator(MPI_COMM_NULL),
    rank(-1),
    size(-1),
    weight_sum(0.0),
    active(-1)
  {
  }
  
  ~PerturbationCollection()
  {
    for (auto &i : perturbations) {
      delete i.p;
    }
  }

  void initialize_mpi(MPI_Comm _communicator)
  {
    MPI_Comm_dup(_communicator, &communicator);
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);

    for (auto &i : perturbations) {
      i.p->initialize_mpi(_communicator);
    }
  }
  
  void add(Perturbation *p, double weight)
  {
    if (weight <= 0.0) {
      throw GENERALVORONOICARTESIANEXCEPTION("Invalid weight: %f\n", weight);
    }
    
    weight_sum += weight;
    perturbations.push_back(WeightedPerturbation(p, weight));
  }

  bool propose(global &g,
	       double &log_prior_ratio,
	       double &log_proposal_ratio,
	       double &log_extra,
	       delta_t *&perturbation,
	       bool &relocate)
  {
    if (active >= 0) {
      throw GENERALVORONOICARTESIANEXCEPTION("Already have active perturbation\n");
    }
    
    if (primary()) {
      double u = g.random.uniform();
      active = 0;
      for (auto &wp: perturbations) {
	double w = wp.w/weight_sum;
	
	if (u < w) {
	  break;
	} else {
	  active ++;
	  u -= w;
	}
      }
    
      if (active >= (int)perturbations.size()) {
	throw GENERALVORONOICARTESIANEXCEPTION("Failed to determine active perturbation\n");
      }
    }

    communicate(active);
    
    WeightedPerturbation &wp = perturbations[active];
    
    int nobs = g.data.obs.size();

    bool r = wp.p->propose(g.nmodels,
			   g.nhierarchical,
			   g.maxcells,
			   nobs,
			   g.random,
			   g.priors,
			   g.position_priors,
			   g.hierarchical_priors,
			   g.models,
			   *g.hierarchical,
			   g.temperature,
			   log_prior_ratio,
			   log_proposal_ratio,
			   log_extra,
			   perturbation,
			   relocate);
    if (!r) {
      active = -1;
    }
    
    return r;
  }

  void accept(global &g)
  {
    if (active < 0 || active >= (int)perturbations.size()) {
      throw GENERALVORONOICARTESIANEXCEPTION("Invalid active perturbation %d (%d)\n", active, (int)perturbations.size());
    }
    
    WeightedPerturbation &wp = perturbations[active];
    wp.p->accept();
    
    active = -1;
  }
  

  void reject(global &g)
  {
    if (active < 0 || active >= (int)perturbations.size()) {
      throw GENERALVORONOICARTESIANEXCEPTION("Invalid active perturbation %d\n", active);
    }
    
    WeightedPerturbation &wp = perturbations[active];
    wp.p->reject(g.models, *(g.hierarchical));
    
    active = -1;
  }
  

  void writeacceptancereport(FILE *fp)
  {
    fprintf(fp, "%s", generateacceptancereport().c_str());
  }

  std::string generateacceptancereport()
  {
    std::string s;
    char linebuffer[1024];
    
    for (auto &wp: perturbations) {

      sprintf(linebuffer, "  %12s:", wp.p->displayname());
      s += linebuffer;

      for (int i = 0; i < wp.p->pa_categories(); i ++) {
	
	int p = wp.p->proposal_count(i);
	int a = wp.p->acceptance_count(i);

	double f = 0.0;
	if (p > 0) {
	  f = (double)a/(double)p * 100.0;
	}
	
	sprintf(linebuffer, "|%8d/%8d : %6.2f |", a, p, f);
	s += linebuffer;
      }
      
      s += "\n";
    }

    return s;
  }

  int get_active() const
  {
    return active;
  }
  
private:

  bool primary()
  {
    return (communicator == MPI_COMM_NULL) || (rank == 0);
  }

  void communicate(int &i)
  {
    if (communicator != MPI_COMM_NULL) {
      MPI_Bcast(&i, 1, MPI_INT, 0, communicator);
    }
  }

  void communicate(double &d)
  {
    if (communicator != MPI_COMM_NULL) {
      MPI_Bcast(&d, 1, MPI_DOUBLE, 0, communicator);
    }
  }
  

  struct WeightedPerturbation {
    WeightedPerturbation(Perturbation *_p, double _w) :
      p(_p),
      w(_w)
    {
    }
    
    Perturbation *p;
    double w;
  };

  MPI_Comm communicator;
  int rank;
  int size;
  
  double weight_sum;
  std::vector<WeightedPerturbation> perturbations;
  int active;

};

#endif // perturbationcollection_hpp
