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
#ifndef prior_hpp
#define prior_hpp

#include <map>
#include <string>

#include "rng.hpp"

class Prior {
public:

  static const double INVALID_LOGP;

  typedef enum {
    NOBOUNDS = 0,
    LOWBOUND = 1,
    HIGHBOUND = 2,
    BOUNDED = 3
  } bound_t;
  
  Prior();
  virtual ~Prior();

  virtual double sample(Rng &rng) = 0;

  virtual bool valid(double v) = 0;

  virtual double pdf(double v) = 0;

  virtual bound_t bounded(double &_vmin, double &_vmax) const
  {
    return NOBOUNDS;
  }

  virtual double logpdf(double v) = 0;

  virtual double logpdf(double v, double &ddv) = 0;
  
  static Prior *load(FILE *fp);

  typedef Prior* (*prior_reader_f)(FILE *fp);

  static bool register_prior(const char *name, prior_reader_f reader);
  
private:

  static std::map<std::string, prior_reader_f> readers;
};

//
// Abstract proposal
//
class Proposal {
public:

  Proposal(Prior &prior);
  virtual ~Proposal();

  virtual Prior *get_prior();
  
  virtual bool propose(Rng &rng,
		       double temperature,
		       double oldv,
		       double &newv,
		       double &logpriorratio) = 0;

  virtual double log_proposal(Rng &rng, double temperature,
			      double oldv,
			      double newv) = 0;
  
  virtual double log_proposal_ratio(Rng &rng, double temperature,
				    double oldv,
				    double newv) = 0;

  static Proposal* load(FILE *fp, Prior &prior);

  typedef Proposal* (*proposal_reader_f)(FILE *fp, Prior &prior);
  
  static bool register_proposal(const char *name, proposal_reader_f reader);

protected:

  Prior &prior;

private:

  static std::map<std::string, proposal_reader_f> readers;

};

class PriorProposal {
public:

  PriorProposal(Prior *prior, Proposal *proposal);
  ~PriorProposal();
  
  double sample(Rng &rng);

  double pdf(double v);

  double logpdf(double v);

  double logpdf(double v, double &ddv);
  
  bool propose(Rng &rng, double temperature, double oldv, double &newv, double &logpriorratio);
  
  double log_proposal(Rng &rng, double temperature,
		      double oldv,
		      double newv);
  
  double log_proposal_ratio(Rng &rng, double temperature,
			    double oldv,
			    double newv);

  Prior *get_prior()
  {
    return prior;
  }
  
  Proposal *get_proposal()
  {
    return proposal;
  }

  static PriorProposal* load(const char *filename);

private:

  Prior *prior;
  Proposal *proposal;

};

//
// Standard Priors
//

class UniformPrior : public Prior { 
public:

  UniformPrior(double vmin, double vmax);
  ~UniformPrior();

  virtual bool valid(double v);
  
  virtual double sample(Rng &rng);

  virtual double pdf(double v);

  virtual bound_t bounded(double &_vmin, double &_vmax) const
  {
    _vmin = vmin;
    _vmax = vmax;
    
    return Prior::BOUNDED;
  }
  
  virtual double logpdf(double v);

  virtual double logpdf(double v, double &ddv);
  
  static Prior* reader(FILE *fp);
  
private:

  static const bool REGISTRATION;
  
  double vmin;
  double vmax;
  
};

class CosinePrior : public Prior {
public:
  CosinePrior();
  ~CosinePrior();

  virtual bool valid(double v);
  
  virtual double sample(Rng &rng);

  virtual double pdf(double v);

  virtual double logpdf(double v);

  virtual double logpdf(double v, double &ddv);

  static Prior* reader(FILE *fp);
  
private:

  static const bool REGISTRATION;
  
};
  
class GaussianPrior : public Prior {
public:

  GaussianPrior(double mu, double sigma);
  ~GaussianPrior();

  virtual bool valid(double v);
  
  virtual double sample(Rng &rng);

  virtual double pdf(double v);

  virtual double logpdf(double v);

  virtual double logpdf(double v, double &ddv);

  static Prior* reader(FILE *fp);
  
private:

  static const bool REGISTRATION;

  double mu;
  double sigma;
  
};

class LogNormalPrior : public Prior {
public:
  
  LogNormalPrior(double muhat, double sigmahat);
  ~LogNormalPrior();

  virtual bool valid(double v);
  
  virtual double sample(Rng &rng);

  virtual double pdf(double v);

  virtual double logpdf(double v);

  virtual double logpdf(double v, double &ddv);

  static Prior* reader(FILE *fp);
  
private:

  static const bool REGISTRATION;

  double muhat;
  double sigmahat;
  
};

//
// Truncated Jeffreys prior. Basically a Jeffreys prior with limits so that
// it can be sampled from and normalized (i.e. proper).
//
class JeffreysPrior : public Prior {
public:

  JeffreysPrior(double vmin, double vmax);
  ~JeffreysPrior();

  virtual bool valid(double v);
  
  virtual double sample(Rng &rng);

  virtual double pdf(double v);

  virtual double logpdf(double v);

  virtual double logpdf(double v, double &ddv);

  static Prior* reader(FILE *fp);
  
private:

  static const bool REGISTRATION;

  double vmin;
  double vmax;
  double C;
  

};

//
// Standard Proposals
//

class PriorSampleProposal : public Proposal {
public:
  
  PriorSampleProposal(Prior &prior);
  ~PriorSampleProposal();

  virtual bool propose(Rng &rng,
		       double temperature,
		       double oldv,
		       double &newv,
		       double &logpriorratio);

  virtual double log_proposal(Rng &rng, double temperature,
			      double oldv,
			      double newv);
  
  virtual double log_proposal_ratio(Rng &rng, double temperature,
				    double oldv,
				    double newv);

  static Proposal *read(FILE *fp, Prior &prior);

private:

  static const bool REGISTRATION;
  
};

class GaussianProposal : public Proposal {
public:

  GaussianProposal(Prior &prior, double std);
  ~GaussianProposal();

  virtual bool propose(Rng &rng,
		       double temperature,
		       double oldv,
		       double &newv,
		       double &logpriorratio);

  virtual double log_proposal(Rng &rng, double temperature,
			      double oldv,
			      double newv);
  
  virtual double log_proposal_ratio(Rng &rng, double temperature,
				    double oldv,
				    double newv);

  static Proposal *read(FILE *fp, Prior &prior);

private:

  static const bool REGISTRATION;
  
  double std;

};

class LogGaussianProposal : public Proposal {
public:

  LogGaussianProposal(Prior &prior, double std);
  ~LogGaussianProposal();

  virtual bool propose(Rng &rng,
                       double temperature,
                       double oldv,
                       double &newv,
                       double &logpriorratio);

  virtual double log_proposal(Rng &rng, double temperature,
                              double oldv,
                              double newv);
  
  virtual double log_proposal_ratio(Rng &rng, double temperature,
                                    double oldv,
                                    double newv);

  static Proposal *read(FILE *fp, Prior &prior);

private:

  static const bool REGISTRATION;
  
  double std;

};
  

#endif // prior_hpp
