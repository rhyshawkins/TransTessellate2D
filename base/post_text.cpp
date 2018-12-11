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

#include <getopt.h>

#include "cartesianvoronoimodel.hpp"
#include "chainhistoryVoronoi.hpp"
#include "genericinterface.hpp"

static char short_options[] = "i:I:x:X:y:Y:A:T:o:t:s:h";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"modelindex", required_argument, 0, 'I'},
  
  {"xmin", required_argument, 0, 'x'},
  {"xmax", required_argument, 0, 'X'},
  {"ymin", required_argument, 0, 'y'},
  {"ymax", required_argument, 0, 'Y'},
  {"parameterization", required_argument, 0, 'A'},
  {"maxcells", required_argument, 0, 'T'},
  
  {"output", required_argument, 0, 'o'},
  
  {"thin", required_argument, 0, 't'},
  {"skip", required_argument, 0, 's'},

  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
 
};

static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;
  
  int maxcells;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  std::vector<int> parameterization_type;

  //
  // Chain processing
  //
  int skip;
  int thin;

  //
  // Input files
  //
  char *input;
  int modelindex;

  //
  // Output Files
  //
  char *output;

  int nmodels;
  int nhierarchical;
  
  //
  // Defaults
  //

  input = nullptr;
  output = nullptr;
  
  maxcells = 1000;
  xmin = -1.0;
  xmax = 1.0;
  ymin = -1.0;
  ymax = 1.0;

  skip = 0;
  thin = 0;
  
  modelindex = 0;
  
  if (gvcart_initialise_(&nmodels, &nhierarchical) < 0) {
    fprintf(stderr, "error: failed to initialise\n");
    return -1;
  }

  if (nmodels < 1 || nhierarchical < 1) {
    fprintf(stderr, "error: invalid no. models/hierarchical\n");
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
      modelindex = atoi(optarg);
      if (modelindex < 0 || modelindex >= nmodels) {
	fprintf(stderr, "error: model index out of range\n");
	return -1;
      }
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

    case 'T':
      maxcells = atoi(optarg);
      if (maxcells < 1) {
	fprintf(stderr, "error: max cells must be 1 or greater\n");
	return -1;
      }
      break;

    case 'o':
      output = optarg;
      break;

    case 't':
      thin = atoi(optarg);
      if (thin < 0) {
	fprintf(stderr, "error: thin must be 0 or positive\n");
	return -1;
      }
      break;

    case 's':
      skip = atoi(optarg);
      if (skip < 0) {
	fprintf(stderr, "error: skip must be 0 or greater\n");
	return -1;
      }
      break;

    default:
      fprintf(stderr, "error: invalid option '%c'\n", c);
      
    case 'h':
      usage(argv[0]);
      return -1;
    }
  }

  if (input == nullptr) {
    fprintf(stderr, "error: required input file parameter missing\n");
    return -1;
  }

  if (output == nullptr) {
    fprintf(stderr, "error: required output file parameter missing\n");
    return -1;
  }

  chainhistoryreaderVoronoi reader(input);

  std::vector<cartesianvoronoimodel *> models;

  for (int i = 0; i < nmodels; i ++) {
    int pt = 0;
    if (i < (int)parameterization_type.size()) {
      pt = parameterization_type[i];
    }
    models.push_back(new cartesianvoronoimodel(maxcells,
					       xmin, xmax,
					       ymin, ymax,
					       pt));
  }
  
  hierarchical_model hierarchical(nhierarchical);

  double likelihood;
  double norm;
    
  int status = reader.step(models, hierarchical, likelihood, norm);
  int step = 0;

  FILE *fp = fopen(output, "w");
  if (fp == NULL) {
    fprintf(stderr, "error: failed to create text file\n");
    return -1;
  }

  while (status > 0) {
      
    if (step >= skip) {
	
      if (thin <= 1 || (step - skip) % thin == 0) {
	
	fprintf(fp, "%6d %15.9f %2d",
		step, likelihood, hierarchical.get_nhierarchical());

	for (int i = 0; i < hierarchical.get_nhierarchical(); i ++) {
	  fprintf(fp, "%15.9f ", hierarchical.get(i));
	}

	fprintf(fp, "%2d ", models[modelindex]->ntotalcells());

	for (int i = 0; i < models[modelindex]->ntotalcells(); i ++) {

	  double x, y, z;

	  models[modelindex]->position_at_index(i, &x, &y);
	  z = models[modelindex]->value_at_index(i);
	    
	  fprintf(fp, "%15.9f %15.9f %15.9f ",
		  x, y, z);
	}

	fprintf(fp, "\n");
      }
    }

    status = reader.step(models, hierarchical, likelihood, norm);
    step ++;
    
    if ((step % 100000) == 0) {
      printf("%d\n", step);
    }
  }
  
  if (status < 0) {
    fprintf(stderr, "error: failed to step through chain history\n");
    return -1;
  }

  fclose(fp);

  return 0;
}

static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  " -i|--input <filename>        Input observations file (required)\n"
	  " -o|--output <filename>       Output text file (required)\n"
	  "\n"
	  " -t|--thin <int>             Only use every nth model\n"
	  " -s|--skip <int>             Skip first n models\n"
	  "\n"
	  " -h|--help                   Usage\n"
	  "\n",
	  pname);
}
