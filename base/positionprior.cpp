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

#include "positionprior.hpp"

#include "generalvoronoicartesianexception.hpp"

extern "C" {
  #include "slog.h"
};

std::map<std::string, PositionPrior::position_prior_reader_f> PositionPrior::readers;

static double distance(double x0, double y0,
		       double x1, double y1);

PositionPrior::PositionPrior()
{
}

PositionPrior::~PositionPrior()
{
}

PositionPrior *
PositionPrior::load(FILE *fp,
		    double xmin, double xmax,
		    double ymin, double ymax)
{
  char type[256];
  char name[256];

  if (fscanf(fp, "%s %s\n", type, name) != 2) {
    return nullptr;
  }

  if (strcmp(type, "positionprior") != 0) {
    fprintf(stderr, "PositionPrior::load: exptected 'positionprior', got '%s'\n", type);
    return nullptr;
  }

  std::map<std::string, PositionPrior::position_prior_reader_f>::iterator i = readers.find(name);
  if (i == readers.end()) {
    fprintf(stderr, "PositionPrior::load: no prior name '%s'\n", name);
    return nullptr;
  } else {

    return (i->second)(fp, xmin, xmax, ymin, ymax);

  }
}

bool
PositionPrior::register_prior(const char *name, position_prior_reader_f reader)
{
  std::map<std::string, position_prior_reader_f>::iterator i = readers.find(name);
  if (i == readers.end()) {
    readers[name] = reader;
    return true;
  }

  throw GENERALVORONOICARTESIANEXCEPTION("Duplicate prior registered: %s\n", name);
}


std::map<std::string, PositionProposal::position_proposal_reader_f> PositionProposal::readers;

PositionProposal::PositionProposal(PositionPrior &_prior) :
  prior(_prior)
{
}

PositionProposal::~PositionProposal()
{
}

PositionPrior *
PositionProposal::get_prior()
{
  return &prior;
}

PositionProposal*
PositionProposal::load(FILE *fp, PositionPrior &prior)
{
  char type[256];
  char name[256];

  if (fscanf(fp, "%s %s\n", type, name) != 2) {
    return nullptr;
  }

  if (strcmp(type, "positionproposal") != 0) {
    fprintf(stderr, "PositionProposal::load: exptected 'proposal', got '%s'\n", type);
    return nullptr;
  }

  std::map<std::string, position_proposal_reader_f>::iterator i = readers.find(name);
  if (i == readers.end()) {
    fprintf(stderr, "PositionProposal::load: no proposal name '%s'\n", name);
    return nullptr;
  } else {

    return (i->second)(fp, prior);

  }
}

bool
PositionProposal::register_position_proposal(const char *name, position_proposal_reader_f reader)
{
  std::map<std::string, position_proposal_reader_f>::iterator i = readers.find(name);
  if (i == readers.end()) {
    readers[name] = reader;
    return true;
  }

  throw GENERALVORONOICARTESIANEXCEPTION("Duplicate proposal registered: %s\n", name);
}


PositionPriorProposal::PositionPriorProposal(PositionPrior *_prior, PositionProposal *_proposal) :
  prior(_prior),
  proposal(_proposal)
{
}

PositionPriorProposal::~PositionPriorProposal()
{
  delete prior;
  delete proposal;
}
  
void
PositionPriorProposal::sample(Rng &rng, double &x, double &y)
{
  
  prior->sample(rng, x, y);
}

double
PositionPriorProposal::pdf(double x, double y)
{
  return prior->pdf(x, y);
}

double
PositionPriorProposal::logpdf(double x, double y)
{
  return prior->logpdf(x, y);
}

bool
PositionPriorProposal::propose(Rng &rng,
				double temperature,
				double oldx,
				double oldy,
				double &newx,
				double &newy,
				double &logpriorratio)
{
  return proposal->propose(rng, temperature, oldx, oldy, newx, newy, logpriorratio);
}
  
double
PositionPriorProposal::log_proposal(Rng &rng,
					   double temperature,
					   double oldx,
					   double oldy,
					   double newx,
					   double newy)
{
  return proposal->log_proposal(rng, temperature, oldx, oldy, newx, newy);
}

double
PositionPriorProposal::log_proposal_ratio(Rng &rng,
					   double temperature,
					   double oldx,
					   double oldy,
					   double newx,
					   double newy)
{
  return proposal->log_proposal_ratio(rng, temperature, oldx, oldy, newx, newy);
}

PositionPriorProposal*
PositionPriorProposal::load(const char *filename,
			    double xmin, double xmax,
			    double ymin, double ymax)
{
  FILE *fp;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stderr, "PositionPriorProposal::load: failed to open file for reading\n");
    return nullptr;
  }

  PositionPrior *prior = PositionPrior::load(fp,
					     xmin, xmax,
					     ymin, ymax);
  if (prior == nullptr) {
    return nullptr;
  }

  PositionProposal *proposal = PositionProposal::load(fp, *prior);
  if (proposal == nullptr) {
    return nullptr;
  }

  fclose(fp);

  return new PositionPriorProposal(prior, proposal);
}

//
// Uniform Prior
//

const bool UniformPositionPrior::REGISTRATION =
  PositionPrior::register_prior("UniformPosition", UniformPositionPrior::reader);

UniformPositionPrior::UniformPositionPrior(double _xmin, double _xmax,
					   double _ymin, double _ymax) :
  xmin(_xmin),
  xmax(_xmax),
  ymin(_ymin),
  ymax(_ymax)
{
}

UniformPositionPrior::~UniformPositionPrior()
{
}

bool
UniformPositionPrior::valid(double x, double y)
{
  if (x < xmin || x > xmax) {
    return false;
  }

  if (y < ymin || y > ymax) {
    return false;
  }
  
  return true;
}
  
void
UniformPositionPrior::sample(Rng &rng, double &x, double &y)
{
  y = rng.uniform() * (ymax - ymin) + ymin;
  x = rng.uniform() * (xmax - xmin) + xmin;
}

double
UniformPositionPrior::pdf(double x, double y)
{
  if (valid(x, y)) {
    return 1.0/((xmax - xmin) * (ymax - ymin));
  } else {
    return 0.0;
  }
}
  

double
UniformPositionPrior::logpdf(double x, double y)
{
  if (valid(x, y)) {
    return -log((xmax - xmin) * (ymax - ymin));
  } else {
    return 0.0;
  }
}

PositionPrior*
UniformPositionPrior::reader(FILE *fp,
			     double xmin, double xmax,
			     double ymin, double ymax)
{
  return new UniformPositionPrior(xmin, xmax,
				  ymin, ymax);
}

//
//
//
const bool PriorSamplePositionProposal::REGISTRATION = PositionProposal::register_position_proposal("PriorSamplePosition", PriorSamplePositionProposal::read);

PriorSamplePositionProposal::PriorSamplePositionProposal(PositionPrior &_prior) :
  PositionProposal(_prior)
{
}

PriorSamplePositionProposal::~PriorSamplePositionProposal()
{
}

bool
PriorSamplePositionProposal::propose(Rng &rng,
				      double temperature,
				      double oldx,
				      double oldy,
				      double &newx,
				      double &newy,
				      double &logpriorratio)
{
  prior.sample(rng, newx, newy);
  logpriorratio = prior.logpdf(newx, newy) - prior.logpdf(oldx, oldy);
  
  return true;
}

double
PriorSamplePositionProposal::log_proposal(Rng &rng,
					   double temperature,
					   double oldx,
					   double oldy,
					   double newx,
					   double newy)
{
  return prior.logpdf(newx, newy);
}

double
PriorSamplePositionProposal::log_proposal_ratio(Rng &rng,
						 double temperature,
						 double oldx,
						 double oldy,
						 double newx,
						 double newy)
{
  return prior.logpdf(newx, newy) - prior.logpdf(oldx, oldy);
}

PositionProposal *
PriorSamplePositionProposal::read(FILE *fp, PositionPrior &prior)
{
  return new PriorSamplePositionProposal(prior);
}

//
// von Mises Position Proposal
//

const bool GaussianPositionProposal::REGISTRATION = PositionProposal::register_position_proposal("GaussianPosition", GaussianPositionProposal::read);

GaussianPositionProposal::GaussianPositionProposal(PositionPrior &prior, double _sigma) :
  PositionProposal(prior),
  sigma(_sigma)
{
}

GaussianPositionProposal::~GaussianPositionProposal()
{
}

bool
GaussianPositionProposal::propose(Rng &rng,
				  double temperature,
				  double oldx,
				  double oldy,
				  double &newx,
				  double &newy,
				  double &logpriorratio)
{
  double theta = 2.0 * M_PI * rng.uniform();
  double r = rng.normal(sigma);

  newx = oldx + r*cos(theta);
  newy = oldy - r*sin(theta);

  if (prior.valid(newx, newy)) {
  
    logpriorratio = prior.logpdf(newx, newy) - prior.logpdf(oldx, oldy);
    return true;

  } else {

    return false;
  }
}

double
GaussianPositionProposal::log_proposal(Rng &rng, double temperature,
				       double oldx,
				       double oldy,
				       double newx,
				       double newy)
{
  double r = distance(oldx, oldy, newx, newy);

  return rng.pdf_normal(r, 0.0, sigma)/(2.0 * M_PI);
}

double GaussianPositionProposal::log_proposal_ratio(Rng &rng,
						    double temperature,
						    double oldx,
						    double oldy,
						    double newx,
						    double newy)
{
  return 0.0;
}

PositionProposal *
GaussianPositionProposal::read(FILE *fp, PositionPrior &prior)
{
  double sigma;

  if (fscanf(fp, "%lf", &sigma) != 1) {
    ERROR("Failed to read parameters\n");
    return nullptr;
  }

  return new GaussianPositionProposal(prior, sigma);
}

//
// Helper functions
//

static double distance(double x0, double y0,
		       double x1, double y1)
{
  double dx = x1 - x0;
  double dy = y1 - y0;

  return sqrt(dx*dx + dy*dy);
}




