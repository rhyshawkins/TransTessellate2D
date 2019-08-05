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
#ifndef perturbation_hpp
#define perturbation_hpp

#include <mpi.h>

#include "rng.hpp"
#include "prior.hpp"
#include "positionprior.hpp"
#include "hierarchical_model.hpp"

#include "cartesianvoronoimodel.hpp"

#include "chainhistoryVoronoi.hpp"

class deltaVoronoi;

class Perturbation {
public:

  typedef deltaVoronoi delta_t;
  
  Perturbation() :
    communicator(MPI_COMM_NULL),
    rank(-1),
    size(-1)
  {
  }
  
  virtual ~Perturbation()
  {
  }

  virtual void initialize_mpi(MPI_Comm _communicator)
  {
    MPI_Comm_dup(_communicator, &communicator);
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);
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
		       bool &relocate) = 0;

  virtual void accept() = 0;

  virtual void reject(std::vector<cartesianvoronoimodel*> &models,
		      hierarchical_model &hierarchical) = 0;

  virtual int pa_categories() const = 0;
  
  virtual int proposal_count(int ci) const = 0;
  
  virtual int acceptance_count(int ci) const = 0;

  virtual const char *displayname() const = 0;

protected:

  bool primary()
  {
    return (communicator == MPI_COMM_NULL) || (rank == 0);
  }

  bool secondary()
  {
    return (communicator != MPI_COMM_NULL) && (rank != 0);
  }

  void communicate(bool &b)
  {
    if (communicator != MPI_COMM_NULL) {
      int i = (int)b;
      MPI_Bcast(&i, 1, MPI_INT, 0, communicator);
      b = (bool)i;
    }
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

  void communicate(int n, double *d)
  {
    if (communicator != MPI_COMM_NULL) {
      MPI_Bcast(d, n, MPI_DOUBLE, 0, communicator);
    }
  }

private:

  MPI_Comm communicator;
  int rank;
  int size;
  

};

#endif // perturbation_hpp
