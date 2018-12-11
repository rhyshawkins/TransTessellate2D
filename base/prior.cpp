//
//    TransTesselate2D : A general Trans-dimensional Tesselation program
//    for 2D Cartesian problems.
//
//    Copyright (C) 2014 - 2018 Rhys Hawkins
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

#include <string.h>
#include <math.h>

#include "prior.hpp"

#include "generalvoronoicartesianexception.hpp"

const double Prior::INVALID_LOGP = 1.0e99;

std::map<std::string, Prior::prior_reader_f> Prior::readers;

static double
newtonsolve(double (*f)(double x),
	    double (*fprime)(double x),
	    double xmin,
	    double xmax,
	    double y,
	    double threshold,
	    int maxsteps);

Prior::Prior()
{
}

Prior::~Prior()
{
}

Prior *
Prior::load(FILE *fp)
{
  char type[256];
  char name[256];

  if (fscanf(fp, "%s %s\n", type, name) != 2) {
    return nullptr;
  }

  if (strcmp(type, "prior") != 0) {
    fprintf(stderr, "Prior::load: exptected 'prior', got '%s'\n", type);
    return nullptr;
  }

  std::map<std::string, Prior::prior_reader_f>::iterator i = readers.find(name);
  if (i == readers.end()) {
    fprintf(stderr, "Prior::load: no prior name '%s'\n", name);
    return nullptr;
  } else {

    return (i->second)(fp);

  }
}

bool
Prior::register_prior(const char *name, prior_reader_f reader)
{
  std::map<std::string, prior_reader_f>::iterator i = readers.find(name);
  if (i == readers.end()) {
    readers[name] = reader;
    return true;
  }

  throw GENERALVORONOICARTESIANEXCEPTION("Duplicate prior registered: %s\n", name);
}


std::map<std::string, Proposal::proposal_reader_f> Proposal::readers;

Proposal::Proposal(Prior &_prior) :
  prior(_prior)
{
}

Proposal::~Proposal()
{
}

Prior *
Proposal::get_prior()
{
  return &prior;
}

Proposal*
Proposal::load(FILE *fp, Prior &prior)
{
  char type[256];
  char name[256];

  if (fscanf(fp, "%s %s\n", type, name) != 2) {
    return nullptr;
  }

  if (strcmp(type, "proposal") != 0) {
    fprintf(stderr, "Proposal::load: exptected 'proposal', got '%s'\n", type);
    return nullptr;
  }

  std::map<std::string, proposal_reader_f>::iterator i = readers.find(name);
  if (i == readers.end()) {
    fprintf(stderr, "Proposal::load: no proposal name '%s'\n", name);
    return nullptr;
  } else {

    return (i->second)(fp, prior);

  }
}

typedef Proposal* (*proposal_reader_f)(FILE *fp, Prior &prior);
  
bool
Proposal::register_proposal(const char *name, proposal_reader_f reader)
{
  std::map<std::string, proposal_reader_f>::iterator i = readers.find(name);
  if (i == readers.end()) {
    readers[name] = reader;
    return true;
  }

  throw GENERALVORONOICARTESIANEXCEPTION("Duplicate proposal registered: %s\n", name);
}


PriorProposal::PriorProposal(Prior *_prior, Proposal *_proposal) :
  prior(_prior),
  proposal(_proposal)
{
}

PriorProposal::~PriorProposal()
{
  delete prior;
  delete proposal;
}
  
double
PriorProposal::sample(Rng &rng)
{
  return prior->sample(rng);
}

double
PriorProposal::pdf(double v)
{
  return prior->pdf(v);
}

double
PriorProposal::logpdf(double v)
{
  return prior->logpdf(v);
}

double
PriorProposal::logpdf(double v, double &ddv)
{
  return prior->logpdf(v, ddv);
}

bool
PriorProposal::propose(Rng &rng, double temperature, double oldv, double &newv, double &logpriorratio)
{
  return proposal->propose(rng, temperature, oldv, newv, logpriorratio);
}
  
double
PriorProposal::log_proposal_ratio(Rng &rng, double temperature,
				  double oldv,
				  double newv)
{
  return proposal->log_proposal_ratio(rng, temperature, oldv, newv);
}

PriorProposal*
PriorProposal::load(const char *filename)
{
  FILE *fp;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stderr, "PriorProposal::load: failed to open file for reading\n");
    return nullptr;
  }

  Prior *prior = Prior::load(fp);
  if (prior == nullptr) {
    return nullptr;
  }

  Proposal *proposal = Proposal::load(fp, *prior);
  if (proposal == nullptr) {
    return nullptr;
  }

  fclose(fp);

  return new PriorProposal(prior, proposal);
}

//
// Uniform prior
//

const bool UniformPrior::REGISTRATION = Prior::register_prior("Uniform", UniformPrior::reader);

Prior*
UniformPrior::reader(FILE *fp)
{
  double vmin, vmax;

  if (fscanf(fp, "%lf %lf\n", &vmin, &vmax) != 2) {
    fprintf(stderr, "UniformPrior::reader: failed to read min/max for prior\n");
    return nullptr;
  }

  return new UniformPrior(vmin, vmax);
}
  
  

UniformPrior::UniformPrior(double _vmin, double _vmax) :
  vmin(_vmin),
  vmax(_vmax)
{
}

UniformPrior::~UniformPrior()
{
}

bool
UniformPrior::valid(double v)
{
  return (v >= vmin && v <= vmax);
}

double
UniformPrior::sample(Rng &rng)
{
  return rng.uniform() * (vmax - vmin) + vmin;
}

double
UniformPrior::pdf(double v)
{
  if (v < vmin || v > vmax) {
    return 0.0;
  }

  return 1.0/(vmax - vmin);
}

double
UniformPrior::logpdf(double v)
{
  if (v < vmin || v > vmax) {
    return INVALID_LOGP;
  }
  
  return log(pdf(v));
}

double
UniformPrior::logpdf(double v, double &ddv)
{
  if (v < vmin || v > vmax) {
    return INVALID_LOGP;
  }

  ddv = 0.0;
  return log(pdf(v));
}

//
// Cosine Prior
//

const bool CosinePrior::REGISTRATION = Prior::register_prior("Cosine", CosinePrior::reader);

CosinePrior::CosinePrior()
{
}

CosinePrior::~CosinePrior()
{
}

bool
CosinePrior::valid(double v)
{
  return (v > -1.0) && (v < 1.0);
}
  
double
CosinePrior::sample(Rng &rng)
{
  //
  // CDF is 0.5*(1/M_PI * sin(M_PI * t) + t + 1.0)
  // Gradient is 0.5*(cos(M_PI * t) + 1.0)
  //
  // Can generate uniform t and then use newton method to find root of CDF - t
  //
  auto cdf = [](double x) -> double { return 0.5 * (sin(M_PI * x)/M_PI + x + 1.0); };
  auto pdf = [](double x) -> double { return 0.5 * (cos(M_PI * x) + 1.0); };

  double u = rng.uniform();

  return newtonsolve(cdf,
		     pdf,
		     -1.0,
		     1.0,
		     u,
		     1.0e-9,
		     1000);
}

double
CosinePrior::pdf(double v)
{
  if ((v > -1.0) && (v < 1.0)) {
    return 0.5 * (cos(M_PI * v) + 1.0);
  }
  
  return 0.0;
}

double
CosinePrior::logpdf(double v)
{
  double p = pdf(v);
  if (p > 0.0) {
    return log(p);
  }

  return INVALID_LOGP;
}

double
CosinePrior::logpdf(double v, double &ddv)
{
  if ((v > -1.0) && (v < 1.0)) {

    double theta = M_PI*v;
    double ctheta = cos(theta);
    double stheta = sin(theta);
    
      
    ddv = (-M_PI * stheta)/(ctheta + 1.0);
    return log(0.5 * (ctheta + 1.0));
  }

  return INVALID_LOGP;
}

Prior*
CosinePrior::reader(FILE *fp)
{
  return new CosinePrior();
}
  
//
// Jeffreys Prior
//

const bool JeffreysPrior::REGISTRATION = Prior::register_prior("Jeffreys", JeffreysPrior::reader);

Prior*
JeffreysPrior::reader(FILE *fp)
{
  double vmin, vmax;

  if (fscanf(fp, "%lf %lf\n", &vmin, &vmax) != 2) {
    fprintf(stderr, "JeffreysPrior::reader: failed to read min/max for prior\n");
    return nullptr;
  }

  return new JeffreysPrior(vmin, vmax);
}


JeffreysPrior::JeffreysPrior(double _vmin, double _vmax) :
  vmin(_vmin),
  vmax(_vmax)
{
  C = 1.0/(log(vmax) - log(vmin));
}

JeffreysPrior::~JeffreysPrior()
{
}

bool
JeffreysPrior::valid(double v)
{
  return (v >= vmin && v <= vmax);
}
  
double
JeffreysPrior::sample(Rng &rng)
{
  return exp(rng.uniform()/C + log(vmin));
}

double
JeffreysPrior::pdf(double v)
{
  if (v >= vmin && v <= vmax) {
    return C/v;
  }

  return 0.0;
}

double
JeffreysPrior::logpdf(double v)
{
  if (v >= vmin && v <= vmax) {
    return log(C) - log(v);
  }

  return INVALID_LOGP;
}

double
JeffreysPrior::logpdf(double v, double &ddv)
{
  if (v >= vmin && v <= vmax) {
    ddv = 1.0/v;
    return log(C) - log(v);
  }

  return INVALID_LOGP;
}
//
// Gaussian Prior
//

const bool GaussianPrior::REGISTRATION = Prior::register_prior("Gaussian", GaussianPrior::reader);

GaussianPrior::GaussianPrior(double _mu, double _sigma) :
  mu(_mu),
  sigma(_sigma)
{
}

GaussianPrior::~GaussianPrior()
{
}

bool
GaussianPrior::valid(double v)
{
  return true;
}
  
double
GaussianPrior::sample(Rng &rng)
{
  return mu + rng.normal(sigma);
}

double
GaussianPrior::pdf(double v)
{
  return Rng::pdf_normal(v, mu, sigma);
}

double
GaussianPrior::logpdf(double v)
{
  return log(pdf(v));
}

double
GaussianPrior::logpdf(double v, double &ddv)
{
  ddv = (v - mu)/(sigma*sigma);
  return log(pdf(v));
}

Prior*
GaussianPrior::reader(FILE *fp)
{
  double mu;
  double sigma;

  if (fscanf(fp, "%lf %lf\n", &mu, &sigma) != 2) {
    fprintf(stderr, "GaussianPrior::reader: failed to read mu, sigma for prior\n");
    return nullptr;
  }

  return new GaussianPrior(mu, sigma);
}

//
// Log Normal Prior
//
const bool LogNormalPrior::REGISTRATION = Prior::register_prior("LogNormal", LogNormalPrior::reader);

LogNormalPrior::LogNormalPrior(double _muhat, double _sigmahat) :
  muhat(_muhat),
  sigmahat(_sigmahat)
{
}

LogNormalPrior::~LogNormalPrior()
{
}

bool
LogNormalPrior::valid(double v)
{
  return v > 0.0;
}
  
double
LogNormalPrior::sample(Rng &rng)
{
  return exp(muhat + rng.normal(sigmahat));
}

double
LogNormalPrior::pdf(double v)
{
  if (v > 0.0) {

    double dx = (log(v) - muhat);
    return exp(-dx*dx/(2.0 * sigmahat * sigmahat))/(v * sigmahat * sqrt(2.0 * M_PI));

  } else {

    return 0.0;

  }
}

double
LogNormalPrior::logpdf(double v)
{
  if (v > 0.0) {
    return log(pdf(v));
  } else {
    return INVALID_LOGP;
  }
}
  
double
LogNormalPrior::logpdf(double v, double &ddv)
{
  if (v > 0.0) {

    ddv = (log(v) - muhat)/(sigmahat * sigmahat * v);
    return log(pdf(v));
    
  } else {
    return INVALID_LOGP;
  }
}

Prior*
LogNormalPrior::reader(FILE *fp)
{
  double muhat, sigmahat;

  if (fscanf(fp, "%lf %lf\n", &muhat, &sigmahat) != 2) {
    fprintf(stderr, "LogNormalPrior::reader: failed to parameters for prior\n");
    return nullptr;
  }

  return new LogNormalPrior(muhat, sigmahat);
}
		       
//
// Prior proposal
//
Proposal *
PriorSampleProposal::read(FILE *fp, Prior &prior)
{
  return new PriorSampleProposal(prior);
}

const bool PriorSampleProposal::REGISTRATION = Proposal::register_proposal("PriorSample", PriorSampleProposal::read);

PriorSampleProposal::PriorSampleProposal(Prior &prior) :
  Proposal(prior)
{
}

PriorSampleProposal::~PriorSampleProposal()
{
}

bool
PriorSampleProposal::propose(Rng &rng,
			     double temperate,
			     double oldv,
			     double &newv,
			     double &logpriorratio)
{
  newv = prior.sample(rng);
  logpriorratio = prior.logpdf(newv) - prior.logpdf(oldv);

  return true;
}

double
PriorSampleProposal::log_proposal(Rng &rng,
					double temperature,
					double oldv,
					double newv)
{
  return prior.logpdf(newv);
}

double
PriorSampleProposal::log_proposal_ratio(Rng &rng,
					double temperature,
					double oldv,
					double newv)
{
  return 0.0;
}

//
// Gaussian proposal
//
Proposal *
GaussianProposal::read(FILE *fp, Prior &prior)
{
  double std;

  if (fscanf(fp, "%lf\n", &std) != 1) {
    fprintf(stderr, "GaussianProposal::read: failed to read std dev\n");
    return nullptr;
  }

  return new GaussianProposal(prior, std);
}

const bool GaussianProposal::REGISTRATION = Proposal::register_proposal("Gaussian", GaussianProposal::read);

GaussianProposal::GaussianProposal(Prior &prior, double _std) :
  Proposal(prior),
  std(_std)
{
}
  
GaussianProposal::~GaussianProposal()
{
}

bool
GaussianProposal::propose(Rng &rng,
			  double temperature,
			  double oldv,
			  double &newv,
			  double &logpriorratio)
{
  newv = oldv + rng.normal(std * sqrt(temperature));

  if (prior.valid(newv)) {
    logpriorratio = prior.logpdf(newv) - prior.logpdf(oldv);
    return true;
  }

  return false;
}

double
GaussianProposal::log_proposal(Rng &rng,
			       double temperature,
			       double oldv,
			       double newv)
{
  return log(rng.pdf_normal(newv, oldv, std * sqrt(temperature)));
}

double
GaussianProposal::log_proposal_ratio(Rng &rng,
				     double temperature,
				     double oldv,
				     double newv)
{
  return 0.0;
}

//
// Log Gaussian proposal
//
Proposal *
LogGaussianProposal::read(FILE *fp, Prior &prior)
{
  double std;

  if (fscanf(fp, "%lf\n", &std) != 1) {
    fprintf(stderr, "LogGaussianProposal::read: failed to read std dev\n");
    return nullptr;
  }

  return new LogGaussianProposal(prior, std);
}

const bool LogGaussianProposal::REGISTRATION = Proposal::register_proposal("LogGaussian", LogGaussianProposal::read);

LogGaussianProposal::LogGaussianProposal(Prior &prior, double _std) :
  Proposal(prior),
  std(_std)
{
}
  
LogGaussianProposal::~LogGaussianProposal()
{
}

bool
LogGaussianProposal::propose(Rng &rng,
                             double temperature,
                             double oldv,
                             double &newv,
                             double &logpriorratio)
{
  newv = exp(log(oldv) + rng.normal(std * sqrt(temperature)));

  if (prior.valid(newv)) {
    logpriorratio = prior.logpdf(newv) - prior.logpdf(oldv);
    return true;
  }

  return false;
}

double
LogGaussianProposal::log_proposal(Rng &rng,
                                  double temperature,
                                  double oldv,
                                  double newv)
{
  return log(rng.pdf_normal(log(newv), log(oldv), std * sqrt(temperature)));
}

double
LogGaussianProposal::log_proposal_ratio(Rng &rng,
                                     double temperature,
                                     double oldv,
                                     double newv)
{
  return 0.0;
}




static double
newtonsolve(double (*f)(double x),
	    double (*fprime)(double x),
	    double xmin,
	    double xmax,
	    double y,
	    double threshold,
	    int maxsteps)
{
  double x = (xmin + xmax)/2.0;
  double yt = f(x);
  int i;

  i = 0;
  while ((i < maxsteps) && (fabs(yt - y) > threshold)) {

    double g = fprime(x);
    double dx;
    
    if (g == 0.0) {
      //
      // Since the cdf's should be monotonicly increasing, this should only happen at the boundaries
      //
      if (y < yt) {
	dx = 1.0e-6;
      } else {
	dx = -1.0e-6;
      }
    } else {
      dx = (y - yt)/g;
    }

    x += dx;
    if (x < xmin) {
      x = xmin + threshold;
    }
    if (x > xmax) {
      x = xmax - threshold;
    }

    yt = f(x);

    i ++;
  }

  return x;
}


