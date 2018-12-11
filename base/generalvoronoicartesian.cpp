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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <getopt.h>

#include "global.hpp"

#include "perturbationcollection.hpp"
#include "value.hpp"
#include "birthgeneric.hpp"
#include "deathgeneric.hpp"
#include "move.hpp"
#include "hierarchical.hpp"

#include "pathutil.hpp"

static char short_options[] = "i:I:o:x:X:y:Y:A:P:H:M:B:T:S:t:l:v:b:pC:z:Z:h";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"initial", required_argument, 0, 'I'},
  {"output", required_argument, 0, 'o'},

  {"xmin", required_argument, 0, 'x'},
  {"xmax", required_argument, 0, 'X'},
  {"ymin", required_argument, 0, 'y'},
  {"ymax", required_argument, 0, 'Y'},
  {"parameterization", required_argument, 0, 'A'},
  

  {"prior", required_argument, 0, 'P'},
  {"hierarchical-prior", required_argument, 0, 'H'},
  {"move-prior", required_argument, 0, 'M'},
  {"birth-death-prior", required_argument, 0, 'B'},

  {"max-cells", required_argument, 0, 'T'},
  {"seed", required_argument, 0, 'S'},
  
  {"total", required_argument, 0, 't'},
  {"lambda", required_argument, 0, 'l'},
  
  {"verbosity", required_argument, 0, 'v'},

  {"birth-probability", required_argument, 0, 'b'},
  {"posterior", no_argument, 0, 'p'},

  {"initial-cells", required_argument, 0, 'C'},
  {"optimize-iterations", required_argument, 0, 'z'},
  {"optimize-temperature", required_argument, 0, 'Z'},
  
  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};

static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;

  //
  // Configuration
  //
  int nmodels;
  int nhierarchical;
    
  char *input;
  std::vector<const char *> initial_models;
  char *output;

  double xmin;
  double xmax;
  double ymin;
  double ymax;

  std::vector<int> parameterization_type;
  
  std::vector<const char *> priors;
  std::vector<const char *> hierarchicalpriors;
  std::vector<const char *> positionpriors;
  std::vector<const char *> birthdeathpriors;

  int maxcells;
  int initialcells;

  int optimizeiterations;
  double optimizetemperature;
  
  std::vector<double> lambdas;

  int seed;

  bool posterior;

  int total;
  int verbosity;

  double Pb;
  double Pm;

  //
  // State
  //
  global *global_state;
  
  double current_likelihood;
  double current_norm;
  
  //
  // Misc
  //
  char filename[1024];

  //
  // Initialize defaults
  //
  input = nullptr;
  output = nullptr;

  xmin = -1.0;
  xmax = 1.0;
  ymin = -1.0;
  ymax = 1.0;

  maxcells = 1000;
  initialcells = 1;

  optimizeiterations = 10000;
  optimizetemperature = 100.0;
  
  seed = 983;

  posterior = false;

  total = 1000;
  verbosity = 1000;

  Pb = 0.05;
  Pm = 0.25;

  if (gvcart_initialise_(&nmodels, &nhierarchical) < 0) {
    fprintf(stderr, "error: failed to initialise\n");
    return -1;
  }
  if (nmodels < 1 || nhierarchical < 1) {
    fprintf(stderr, "error: invalid initialisation parameters %d %d\n",
	    nmodels, nhierarchical);
    return -1;
  }
  
  option_index = 0;
  while (1) {

    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {

    case 'i':
      input = optarg;
      break;

    case 'I':
      initial_models.push_back(optarg);
      break;

    case 'x':
      xmin = atof(optarg);
      break;

    case 'X':
      xmax = atof(optarg);
      break;

    case 'y':
      ymin = atof(optarg);
      break;

    case 'Y':
      ymax = atof(optarg);
      break;

    case 'A':
      {
	int pt = atoi(optarg);
	if (pt < 0 || pt > 2) {
	  fprintf(stderr, "error: parameterization type must be between 0 .. 2\n");
	  return -1;
	}

	parameterization_type.push_back(pt);
      }
      break;
	
    case 'o':
      output = optarg;
      break;

    case 'P':
      priors.push_back(optarg);
      break;

    case 'H':
      hierarchicalpriors.push_back(optarg);
      break;

    case 'M':
      positionpriors.push_back(optarg);
      break;

    case 'B':
      birthdeathpriors.push_back(optarg);
      break;

    case 'T':
      maxcells = atoi(optarg);
      if (maxcells < 1) {
	fprintf(stderr, "error: max cells must be 1 or greater\n");
	return -1;
      }
      break;

    case 'S':
      seed = atoi(optarg);
      break;
	

    case 't':
      total = atoi(optarg);
      if (total < 1) {
	fprintf(stderr, "error: total must be greater than 0\n");
	return -1;
      }
      break;

    case 'l':
      {
	double l = atof(optarg);
	if (l <= 0.0) {
	  fprintf(stderr, "error: lambda must be greater than 0\n");
	  return -1;
	}

	lambdas.push_back(l);
      }
      break;

    case 'v':
      verbosity = atoi(optarg);
      if (verbosity < 0) {
	fprintf(stderr, "error: verbosity must be greater than 0\n");
	return -1;
      }
      break;

    case 'b':
      Pb = atof(optarg);
      if (Pb < 0.0 || Pb >= 0.5) {
	fprintf(stderr, "error: Pb must be between 0 and 0.5\n");
	return -1;
      }
      break;
      
    case 'p':
      posterior = true;
      break;

    case 'C':
      initialcells = atoi(optarg);
      if (initialcells < 1) {
	fprintf(stderr, "error: need at least one initial cell\n");
	return -1;
      }
      break;

    case 'z':
      optimizeiterations = atoi(optarg);
      if (optimizeiterations < 0) {
	fprintf(stderr, "error: optimization iterations must be positive\n");
	return -1;
      }
      break;

    case 'Z':
      optimizetemperature = atof(optarg);
      if (optimizetemperature < 1.0) {
	fprintf(stderr, "error: optimization temperature must be 1 or greater\n");
	return -1;
      }
      break;
      
    default:
      fprintf(stderr, "invalid option %c\n", c);
    case 'h':
      usage(argv[0]);
      return -1;

    }
  }
  
  if (input == nullptr) {
    fprintf(stderr, "error: required parameter input not set\n");
    return -1;
  }

  global_state = new global(nmodels,
			    nhierarchical,
			    xmin,
			    xmax,
			    ymin,
			    ymax,
			    parameterization_type,
			    input,
			    initial_models,
			    priors,
			    positionpriors,
			    hierarchicalpriors,
			    birthdeathpriors,
			    maxcells,
			    lambdas,
			    1.0, // temperature
			    seed,
			    posterior);

  global_state->initialize(initialcells, initial_models);

  current_norm = 0.0;
  current_likelihood = global_state->likelihood(current_norm, true);

  printf(" Initial likelihood: %10.6f %10.6f\n", current_likelihood, current_norm);

  if (optimizeiterations > 0) {

    double acceptance = 0.0;
    current_likelihood = global_state->optimize_sa(optimizeiterations,
						   optimizetemperature,
						   acceptance,
						   verbosity);
    
    printf("Optimized likelihood: %10.6f %10.6f\n", current_likelihood, current_norm);
    printf("Optimization AR     : %10.6f\n", acceptance);
  }
  

  int *khistogram = new int[nmodels * (global_state->maxcells + 1)];
  for (int i = 0; i < (nmodels * (global_state->maxcells + 1)); i ++) {
    khistogram[i] = 0;
  }

  PerturbationCollection pc;
  mkpath(output, "ch.dat", filename);
  chainhistorywriterVoronoi *history = new chainhistorywriterVoronoi(filename,
								     global_state->models,
								     *(global_state->hierarchical),
								     current_likelihood,
								     current_norm);

  if (posterior) {
    Pb = 0.45;
    
    pc.add(new Value(nmodels), (1.0 - (2.0*Pb)) * (1.0 - Pm));
    pc.add(new Move(nmodels), (1.0 - (2.0*Pb)) * Pm);

    pc.add(new BirthGeneric(nmodels,
			    global_state->birthdeath_proposals,
			    global_state->birthdeath_position_proposals), Pb);
    pc.add(new DeathGeneric(nmodels,
			    global_state->birthdeath_proposals,
			    global_state->birthdeath_position_proposals), Pb);
  } else {
    pc.add(new Value(nmodels), (1.0 - (2.0*Pb)) * (1.0 - Pm));
    pc.add(new Move(nmodels), (1.0 - (2.0*Pb)) * Pm);
    
    if (Pb > 0.0) {
      pc.add(new BirthGeneric(nmodels,
			      global_state->birthdeath_proposals,
      			      global_state->birthdeath_position_proposals), Pb);
      pc.add(new DeathGeneric(nmodels,
			      global_state->birthdeath_proposals,
      			      global_state->birthdeath_position_proposals), Pb);
    }

    
    if (hierarchicalpriors.size() > 0) {
      pc.add(new Hierarchical(nhierarchical), 0.5);
    }
  }

  bool last_relocate = false;
  
  for (int i = 0; i < total; i ++) {
    
    double log_prior_ratio;
    double log_proposal_ratio;
    double log_extra;
    Perturbation::delta_t *perturbation = nullptr;

    bool propose_relocate;
    if (pc.propose(*global_state, log_prior_ratio, log_proposal_ratio, log_extra, perturbation, propose_relocate)) {

      if (perturbation == nullptr) {
	throw GENERALVORONOICARTESIANEXCEPTION("Valid proposal has null perturbation\n");
      }

      double u = log(global_state->random.uniform());

      double proposed_norm = 0.0;
      double proposed_likelihood = global_state->likelihood(proposed_norm, last_relocate || propose_relocate);
      perturbation->set_proposed_likelihood(proposed_likelihood, proposed_norm);

      // printf("%3d %2d Proposed: %10.6f %10.6f (%10.6f %10.6f %10.6f) %10.6f %10.6f\n",
      // 	     i, pc.get_active(),
      // 	     proposed_likelihood, proposed_norm,
      // 	     log_prior_ratio, log_proposal_ratio, log_extra,
      // 	     current_likelihood, current_norm);
      
      // if (!relocate) {
      // 	printf("Proposed: %10.6f %10.6f (%10.6f %10.6f)\n",
      // 	       proposed_likelihood, proposed_norm,
      // 	       current_likelihood, current_norm);
      // }
	
      if (u < ((current_likelihood + current_norm) -
	       (proposed_likelihood + proposed_norm) +
	       log_prior_ratio + log_proposal_ratio)) {

	pc.accept(*global_state);
	perturbation->accept();
	global_state->accept();
	
	current_likelihood = proposed_likelihood;
	current_norm = proposed_norm;

	last_relocate = false;
      } else {
	pc.reject(*global_state);
	perturbation->reject(); 
	global_state->reject();

	last_relocate = propose_relocate;
      }
    }

    if (verbosity > 0 && (i + 1) % verbosity == 0) {

      printf("%5d: Likelihood %10.6f (%10.6f) Lambda(s) ",
             i + 1,
             current_likelihood,
	     current_norm);

      for (int hi = 0; hi < global_state->hierarchical->get_nhierarchical(); hi ++) {
	printf("%10.6f ", global_state->hierarchical->get(hi));
      }

      printf("Cell ");
      for (auto &m : global_state->models) {
	printf("%3d ", m->nvaluecells());
      }
      printf("\n");

      pc.writeacceptancereport(stdout);
    }

    for (int mi = 0; mi < nmodels; mi ++) {
      int k = global_state->models[mi]->nvaluecells();
      if (k < 1 || k > maxcells) {
	throw GENERALVORONOICARTESIANEXCEPTION("k out of range: %d (%d)\n", k, maxcells);
      }
      khistogram[mi * (maxcells + 1) + k] ++;
    }
    
    history->add(perturbation);
  }

  //
  // Save khistogram
  //
  mkpath(output, "khistogram.txt", filename);
  FILE *fp = fopen(filename, "w");
  for (int i = 0; i <= global_state->maxcells; i ++) {
    fprintf(fp, "%d ", i);
    for (int mi = 0; mi < nmodels; mi ++) {
      fprintf(fp, "%d ", khistogram[mi * (maxcells + 1) + i]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  delete [] khistogram;

  //
  // Save residuals
  //
  mkpath(output, "residuals.txt", filename);
  fp = fopen(filename, "w");
  for (int i = 0; i < global_state->residual_size; i ++) {
    fprintf(fp, "%.9g\n", global_state->mean_residuals[i]);
  }
  fclose(fp);

  //
  // Save final model
  //
  if (nmodels == 1) {
    mkpath(output, "finalmodel.txt", filename);
    
    if (!global_state->models[0]->save(filename)) {
      throw GENERALVORONOICARTESIANEXCEPTION("Failed to save final model\n");
    }
  } else {
    for (int mi = 0; mi < nmodels; mi ++) {
      mkmodelpath(mi, output, "finalmodel.txt", filename);
      
      if (!global_state->models[mi]->save(filename)) {
	throw GENERALVORONOICARTESIANEXCEPTION("Failed to save final model\n");
      }
    }
  }
    
  history->flush();

  delete history;

  return 0;
}

static void usage(const char *pname)
{
}
