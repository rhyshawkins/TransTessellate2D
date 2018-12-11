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

#include <mpi.h>

#include "global.hpp"

#include "perturbationcollection.hpp"
#include "value.hpp"
#include "birthgeneric.hpp"
#include "deathgeneric.hpp"
#include "move.hpp"
#include "hierarchical.hpp"
#include "ptexchange.hpp"

#include "pathutil.hpp"

static char short_options[] = "i:I:o:x:X:y:Y:A:P:H:M:B:T:S:t:l:v:b:pC:z:Z:c:K:m:e:h";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"initial", required_argument, 0, 'I'},
  {"output", required_argument, 0, 'o'},

  {"xmin", required_argument, 0, 'x'},
  {"xmax", required_argument, 0, 'X'},
  {"ymin", required_argument, 0, 'y'},
  {"ymax", required_argument, 0, 'Y'},

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
  
  {"chains", required_argument, 0, 'c'},
  {"temperatures", required_argument, 0, 'K'},
  {"max-temperature", required_argument, 0, 'm'},
  {"exchange-rate", required_argument, 0, 'e'},
  
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

  int seed_base;
  int seed_mult;

  bool posterior;

  int total;
  int verbosity;

  double Pb;
  double Pm;

  int chains;
  int temperatures;
  double max_temperature;
  int exchange_rate;
  
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
  int flush_rate;

  int mpi_size;
  int mpi_rank;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

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

  seed_base = 983;
  seed_mult = 101;

  posterior = false;

  total = 1000;
  verbosity = 1000;

  Pb = 0.05;
  Pm = 0.25;

  flush_rate = 100000;

  chains = 1;
  temperatures = 1;
  max_temperature = 10.0;
  exchange_rate = 10;

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

    case 'c':
      chains = atoi(optarg);
      if (chains <= 0) {
	fprintf(stderr, "error: no. chains must be greater than 0\n");
	return -1;
      }
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

    case 'K':
      temperatures = atoi(optarg);
      if (temperatures < 1) {
	fprintf(stderr, "error: there must be at least one temperature\n");
	return -1;
      }
      break;

    case 'm':
      max_temperature = atof(optarg);
      if (max_temperature < 1.0) {
	fprintf(stderr, "error: max temperature must be greater than 1\n");
	return -1;
      }
      break;

    case 'e':
      exchange_rate = atoi(optarg);
      if (exchange_rate < 1) {
	fprintf(stderr, "error: exchange rate must be 1 or more\n");
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

  int ntotalchains = temperatures * chains;
  if (ntotalchains == 0 || mpi_size % ntotalchains != 0) {
    fprintf(stderr, "error: no. chains (%d) must be a divisor of MPI processes (%d)\n",
	    ntotalchains,
	    mpi_size);
    return -1;
  }

  if (temperatures > 1 && ntotalchains % 2 != 0) {
    fprintf(stderr, "error: no. chains * no. temperatures must be even\n");
    return -1;
  }
  
  mkrankpath(mpi_rank, output, "log.txt", filename);
  if (slog_set_output_file(filename,
                           SLOG_FLAGS_CLEAR) < 0) {
    fprintf(stderr, "error: failed to redirect log file\n");
    return -1;
  }

  MPI_Comm chain_communicator;
  MPI_Comm temperature_communicator;

  double temperature;
  int chain_rank;
  int chain_id;
  
  if (chains == 1 && temperatures == 1) {
    chain_id = 0;
    temperature = 1.0;
    chain_communicator = MPI_COMM_WORLD;
  } else {
    int processesperchain = mpi_size/ntotalchains;
    chain_id = mpi_rank/processesperchain;

    if (temperatures > 1) {
      int temperature_id = chain_id / chains;
      temperature = pow(10.0, log10(max_temperature) * (double)temperature_id/(double)(temperatures - 1));
    } else {
      temperature = 1.0;
    }
    MPI_Comm_split(MPI_COMM_WORLD, chain_id, mpi_rank, &chain_communicator);
  }

  MPI_Comm_rank(chain_communicator, &chain_rank);

  // const char *initial_model_ptr = nullptr;
  // char initial_model_filename[1024];
  // if (initial != nullptr) {
  //   mkrankpath(chain_id, initial, "finalmodel.txt", initial_model_filename);
  //   initial_model_ptr = initial_model_filename;
  // }

  try {
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
			      temperature,
			      seed_base + seed_mult * mpi_rank,
			      posterior);
  } catch (...) {
    fprintf(stderr, "Exception in process %d, check log files\n", mpi_rank);
    return -1;
  }
  

  Value *value = new Value(nmodels);
  Move *move = new Move(nmodels);
  BirthGeneric *birth = new BirthGeneric(nmodels,
					 global_state->birthdeath_proposals,
					 global_state->birthdeath_position_proposals);
  DeathGeneric *death = new DeathGeneric(nmodels,
					 global_state->birthdeath_proposals,
					 global_state->birthdeath_position_proposals);
  
  global_state->initialize_mpi(chain_communicator, temperature);

  
  global_state->initialize(initialcells, initial_models);

  current_norm = 0.0;
  current_likelihood = global_state->likelihood(current_norm, true);
  if (chain_rank == 0) {
    INFO("Chain %03d: Initial likelihood: %10.6f (%10.6f) T = %6.3f\n", chain_id, current_likelihood, current_norm, temperature);
  }
  global_state->accept();

  if (optimizeiterations > 0) {

    double acceptance = 0.0;
    current_likelihood = global_state->optimize_sa(optimizeiterations,
						   optimizetemperature,
						   acceptance,
						   verbosity);

    INFO("Optimized likelihood: %10.6f %10.6f", current_likelihood, current_norm);
    INFO("Optimization AR     : %10.6f", acceptance);
  }
  
  int *khistogram = nullptr;
  chainhistorywriterVoronoi *history = nullptr;
  
  if (chain_rank == 0) {
    khistogram = new int[nmodels * (global_state->maxcells + 1)];
    for (int i = 0; i < (nmodels * (global_state->maxcells + 1)); i ++) {
      khistogram[i] = 0;
    }

    mkrankpath(chain_id, output, "ch.dat", filename);

    history = new chainhistorywriterVoronoi(filename,
					    global_state->models,
					    *(global_state->hierarchical),
					    current_likelihood,
					    current_norm);
  }

  //
  //
  //
  PerturbationCollection pc;
  
  if (posterior) {
    Pb = 0.45;
  }
  
  pc.add(value, (1.0 - (2.0 * Pb)) * (1.0 - Pm));
  pc.add(move, (1.0 - (2.0 * Pb)) * Pm);
  
  pc.add(birth, Pb);
  pc.add(death, Pb);


  if (hierarchicalpriors.size() > 0) {
    //
    // Add Hierarchical Sampling if required
    //
    Hierarchical *hierarchical = new Hierarchical(nhierarchical);

    pc.add(hierarchical, 0.5);
  }

  pc.initialize_mpi(chain_communicator);

  PTExchange *ptexchange = nullptr;
  if (temperatures > 1) {
    //
    // Add PT Exchanges if required
    //
    ptexchange = new PTExchange(*global_state);

    MPI_Comm_split(MPI_COMM_WORLD, chain_rank == 0, mpi_rank, &temperature_communicator);

    ptexchange->initialize_mpi(MPI_COMM_WORLD,
			       temperature_communicator,
			       chain_communicator,
			       temperatures);
  }
  
  bool last_relocate = false;
  
  for (int i = 0; i < total; i ++) {
    
    double log_prior_ratio;
    double log_proposal_ratio;
    double log_extra;
    Perturbation::delta_t *perturbation = nullptr;
    bool accepted;

    bool propose_relocate;
    if (pc.propose(*global_state, log_prior_ratio, log_proposal_ratio, log_extra, perturbation, propose_relocate)) {

      double proposed_norm = 0.0;
      double proposed_likelihood = global_state->likelihood(proposed_norm,
							    last_relocate || propose_relocate);
      accepted = false;
      
      // printf("%3d %2d Proposed: %10.6f %10.6f (%10.6f %10.6f %10.6f) %10.6f %10.6f\n",
      // 	     i, pc.get_active(),
      // 	     proposed_likelihood, proposed_norm,
      // 	     log_prior_ratio, log_proposal_ratio, log_extra,
      // 	     current_likelihood, current_norm);

      if (chain_rank == 0) {
	
	if (perturbation == nullptr) {
	  throw GENERALVORONOICARTESIANEXCEPTION("Valid proposal has null perturbation\n");
	}
	
	double u = log(global_state->random.uniform());
	
	perturbation->set_proposed_likelihood(proposed_likelihood, proposed_norm);
	
	accepted = u < (((current_likelihood + current_norm) -
			 (proposed_likelihood + proposed_norm))/global_state->temperature +
			log_prior_ratio + log_proposal_ratio);

	if (accepted) {
	  perturbation->accept();
	} else {
	  perturbation->reject();
	}
      }

      int t = accepted;
      MPI_Bcast(&t, 1, MPI_INT, 0, chain_communicator);
      accepted = t;

      if (accepted) {
	pc.accept(*global_state);
	current_likelihood = proposed_likelihood;
	current_norm = proposed_norm;
	global_state->accept();

	last_relocate = false;
      } else {
	pc.reject(*global_state);
	global_state->reject();

	last_relocate = propose_relocate;
      }
    }

   if (ptexchange != nullptr && exchange_rate > 0 && ((i + 1) % exchange_rate == 0)) {

      int exchanged = ptexchange->step(current_likelihood, current_norm);
      if (exchanged < 0) {
        ERROR("Failed to do PT exchange\n");
        return -1;
      }

      if (exchanged) {
	// Completely new model
	last_relocate = true;

	if (chain_rank == 0) {
	  history->add(perturbation);

	  perturbation = new model_initializationVoronoi(global_state->models,
							 *global_state->hierarchical,
							 current_likelihood,
							 current_norm);
	}
      }

    }

    if (chain_rank == 0) {
      if (verbosity > 0 && (i + 1) % verbosity == 0) {

	INFO_LARGE_START("\n%5d: Likelihood %10.6f (%10.6f)\nLambda(s) ",
	     i + 1,
	     current_likelihood,
	     current_norm);

	for (int hi = 0; hi < nhierarchical; hi ++) {
	  INFO_LARGE_WRITE("%7.5f ", global_state->hierarchical->get(hi));
	}
	INFO_LARGE_WRITE("Cell(s) ");
	for (int mi = 0; mi < nmodels; mi ++) {
	  INFO_LARGE_WRITE("%3d ", global_state->models[mi]->nvaluecells());
	}
	std::string report = pc.generateacceptancereport();
	INFO_LARGE_WRITE("\n%s", report.c_str());

	if (ptexchange != nullptr) {
	  
	  std::string pt = ptexchange->write_short_stats();
	  INFO_LARGE_WRITE("%s", pt.c_str());
	  
	}
	
	INFO_LARGE_END();
      }

      for (int mi = 0; mi < nmodels; mi ++) {
	int k = global_state->models[mi]->nvaluecells();
	if (k < 1 || k > maxcells) {
	  throw GENERALVORONOICARTESIANEXCEPTION("k out of range: %d (%d)\n", k, maxcells);
	}
	khistogram[mi * (maxcells + 1) + k] ++;
      }
    
      history->add(perturbation);

      if (i > 0 && i % flush_rate == 0) {

	history->flush();
	
      }
    }
  }

  if (chain_rank == 0) {
    //
    // Save khistogram
    //
    mkrankpath(chain_id, output, "khistogram.txt", filename);
    FILE *fp = fopen(filename, "w");
    for (int i = 0; i <= global_state->maxcells; i ++) {
      fprintf(fp, "%d ", i);
      for (int mi = 0; mi < nmodels; mi ++) {
	fprintf(fp, "%d ", khistogram[mi * (global_state->maxcells + 1) + i]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
    delete [] khistogram;
    
    history->flush();
    
    delete history;

    if (nmodels == 1) {
      mkrankpath(chain_id, output, "finalmodel.txt", filename);
      if (!global_state->models[0]->save(filename)) {
	throw GENERALVORONOICARTESIANEXCEPTION("Failed to save final model\n");
      }
    } else {
      for (int mi = 0; mi < nmodels; mi ++) {
	mkmodelrankpath(mi, chain_id, output, "finalmodel.txt", filename);
	if (!global_state->models[mi]->save(filename)) {
	  throw GENERALVORONOICARTESIANEXCEPTION("Failed to save final model\n");
	}
      }
    }
	

    //
    // Save residuals
    //
    mkrankpath(chain_id, output, "residuals.txt", filename);
    fp = fopen(filename, "w");
    for (int i = 0; i < global_state->residual_size; i ++) {
      fprintf(fp, "%.9g\n", global_state->mean_residuals[i]);
    }
    fclose(fp);
      
  }

  MPI_Finalize();

  return 0;
}

static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  " -i|--input <filename>                   Observations input file\n"
	  " -I|--initial <filename>                 Initial model input file\n"
	  " -o|--output <path>                      Path/prefix for output files\n"
	  "\n"
	  " -P|--prior <filename>                   Prior input file\n"
	  " -H|--hierarchical-prior <filename>      Hierarchical prior input file\n"
	  " -M|--move-prior <filename>              Move prior input file\n"
	  " -B|--birth-death-prior <filename>       Birth/Death proposal file\n"
	  "\n"
	  " -t|--total <int>                        Total number of iterations\n"
	  " -v|--verbosity <int>                    Number of iterations between status updates (0 = none)\n"
	  "\n"
	  " -l|--lambda <float>                     Initial/fixed lambda parameter\n"
	  "\n"
	  " -T|--max-cells <int>                    Max no. Voronoi cells\n"
	  " -b|--birth-probability <float>          Relative probability of birth\n"
	  " -p|--posterior                          Posterior test\n"
	  "\n"
	  " -c|--chains <int>                       No. of chains to run\n"
	  " -K|--temperatures <int>                 No. of temperatures to run\n"
	  "\n"
	  " -h|--help                               Usage information\n"
	  "\n",
	  pname);

}
