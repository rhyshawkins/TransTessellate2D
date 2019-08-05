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
#ifndef positionprior_hpp
#define positionprior_hpp

#include <map>
#include <string>

#include "rng.hpp"

class PositionPrior {
public:

  PositionPrior();
  virtual ~PositionPrior();

  virtual void sample(Rng &rng, double &x, double &y) = 0;

  virtual bool valid(double x, double y) = 0;

  virtual double pdf(double x, double y) = 0;

  virtual double logpdf(double x, double y) = 0;

  static PositionPrior *load(FILE *fp,
			     double xmin, double xmax,
			     double ymin, double ymax);

  typedef PositionPrior* (*position_prior_reader_f)(FILE *fp,
						    double xmin, double xmax,
						    double ymin, double ymax);
  
  static bool register_prior(const char *name, position_prior_reader_f reader);
  
private:

  static std::map<std::string, position_prior_reader_f> readers;
};
  
//
// Abstract proposal
//
class PositionProposal {
public:

  PositionProposal(PositionPrior &prior);
  virtual ~PositionProposal();

  virtual PositionPrior *get_prior();
  
  virtual bool propose(Rng &rng,
		       double temperature,
		       double oldx,
		       double oldy,
		       double &newx,
		       double &newy,
		       double &logpriorratio) = 0;

  virtual double log_proposal(Rng &rng, double temperature,
			      double oldx,
			      double oldy,
			      double newx,
			      double newy) = 0;
  
  virtual double log_proposal_ratio(Rng &rng, double temperature,
				    double oldx,
				    double oldy,
				    double newx,
				    double newy) = 0;

  static PositionProposal* load(FILE *fp, PositionPrior &prior);

  typedef PositionProposal* (*position_proposal_reader_f)(FILE *fp, PositionPrior &prior);
  
  static bool register_position_proposal(const char *name, position_proposal_reader_f reader);

protected:

  PositionPrior &prior;

private:

  static std::map<std::string, position_proposal_reader_f> readers;

};

class PositionPriorProposal {
public:

  PositionPriorProposal(PositionPrior *prior, PositionProposal *proposal);
  ~PositionPriorProposal();
  
  void sample(Rng &rng, double &x, double &y);

  double pdf(double x, double y);

  double logpdf(double x, double y);

  bool propose(Rng &rng,
	       double temperature,
	       double oldx, double oldy,
	       double &newx, double &newy,
	       double &logpriorratio);
  
  double log_proposal_ratio(Rng &rng, double temperature,
			    double oldx, double oldy,
			    double newx, double newy);

  double log_proposal(Rng &rng, double temperature,
		      double oldx, double oldy,
		      double newx, double newy);

  PositionPrior *get_prior()
  {
    return prior;
  }
  
  PositionProposal *get_proposal()
  {
    return proposal;
  }

  static PositionPriorProposal* load(const char *filename,
				     double xmin, double xmax,
				     double ymin, double ymax);

private:

  PositionPrior *prior;
  PositionProposal *proposal;

};

//
// Standard Priors
//

class UniformPositionPrior : public PositionPrior { 
public:

  UniformPositionPrior(double xmin, double xmax,
		       double ymin, double ymax);
  ~UniformPositionPrior();

  virtual bool valid(double x, double y);
  
  virtual void sample(Rng &rng, double &x, double &y);

  virtual double pdf(double x, double y);

  virtual double logpdf(double x, double y);

  static PositionPrior* reader(FILE *fp,
			     double xmin, double xmax,
			     double ymin, double ymax);
  
private:

  static const bool REGISTRATION;
  
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  
};

//
// Standard Proposals
//

class PriorSamplePositionProposal : public PositionProposal {
public:
  
  PriorSamplePositionProposal(PositionPrior &prior);
  ~PriorSamplePositionProposal();

  virtual bool propose(Rng &rng,
		       double temperature,
		       double oldx,
		       double oldy,
		       double &newx,
		       double &newy,
		       double &logpriorratio);

  virtual double log_proposal(Rng &rng, double temperature,
			      double oldx,
			      double oldy,
			      double newx,
			      double newy);
  
  virtual double log_proposal_ratio(Rng &rng,
				    double temperature,
				    double oldx,
				    double oldy,
				    double newx,
				    double newy);

  static PositionProposal *read(FILE *fp, PositionPrior &prior);

private:

  static const bool REGISTRATION;

};

class GaussianPositionProposal : public PositionProposal {
public:

  GaussianPositionProposal(PositionPrior &prior, double sigma);
  ~GaussianPositionProposal();

  virtual bool propose(Rng &rng,
		       double temperature,
		       double oldx,
		       double oldy,
		       double &newx,
		       double &newy,
		       double &logpriorratio);

  virtual double log_proposal(Rng &rng, double temperature,
			      double oldx,
			      double oldy,
			      double newx,
			      double newy);

  virtual double log_proposal_ratio(Rng &rng,
				    double temperature,
				    double oldx,
				    double oldy,
				    double newx,
				    double newy);

  static PositionProposal *read(FILE *fp, PositionPrior &prior);

private:

  static const bool REGISTRATION;
  
  double sigma;

};
  

#endif // positionprior_hpp
