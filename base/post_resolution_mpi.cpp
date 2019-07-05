//
//    GeneralVoronoiCartesian : A general Trans-dimensional Voronoic Cell program
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

#include <mpi.h>

#include <getopt.h>

#include "cartesianvoronoimodel.hpp"
#include "chainhistoryVoronoi.hpp"
#include "genericinterface.hpp"

#include "pathutil.hpp"

static char short_options[] = "i:I:x:X:y:Y:A:T:o:m:M:B:D:V:e:E:g:b:z:Z:W:H:t:s:h";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"modelindex", required_argument, 0, 'I'},
  
  {"output", required_argument, 0, 'o'},
  {"stddev", required_argument, 0, 'D'},
  
  {"credinterval", required_argument, 0, 'I'},

  {"lonsamples", required_argument, 0, 'W'},
  {"latsamples", required_argument, 0, 'H'},
  
  {"thin", required_argument, 0, 't'},
  {"skip", required_argument, 0, 's'},

  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
 
};

static void usage(const char *pname);

static int saveimage(const char *filename, double *image, int width, int height);

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

  int lonsamples;
  int latsamples;
  
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
  char *output; // mean
  char *stddev_file;

  //
  // State
  //
  int imagesize;
  double *image;

  int meann;
  double delta;
  double *mean;
  double *variance;
  double *workspace;

  int mpi_size;
  int mpi_rank;

  int nmodels;
  int nhierarchical;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  //
  // Defaults
  //

  lonsamples = 16;
  latsamples = 16;
  
  maxcells = 1000;
  xmin = -1.0;
  xmax = 1.0;
  ymin = -1.0;
  ymax = 1.0;

  input = nullptr;
  output = nullptr;
  stddev_file = nullptr;

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

    case 'D':
      stddev_file = optarg;
      break;

    case 'W':
      lonsamples = atoi(optarg);
      if (lonsamples < 1) {
	fprintf(stderr, "error: xsamples must be 1 or greater\n");
	return -1;
      }
      break;

    case 'H':
      latsamples = atoi(optarg);
      if (latsamples < 1) {
	fprintf(stderr, "error: latsamples must be 1 or greater\n");
	return -1;
      }
      break;

    case 's':
      skip = atoi(optarg);
      break;

    case 't':
      thin = atoi(optarg);
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

  //
  // Initialize state
  //

  meann = 0;
  imagesize = lonsamples * latsamples;
  mean = new double[imagesize];
  image = new double[imagesize];
  variance = new double[imagesize];
  workspace = new double[imagesize];

  for (int i = 0; i < imagesize; i ++) {
    image[i] = 0.0;
    mean[i] = 0.0;
    variance[i] = 0.0;
  }

  char rankfile[1024];
  mkrankpath(mpi_rank, nullptr, input, rankfile);
  
  chainhistoryreaderVoronoi reader(rankfile);
  
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
  while (status > 0) {
    
    if (step >= skip) {
      
      if (thin <= 1 || (step - skip) % thin == 0) {

	for (auto &m : models) {
	  m->recompute();
	}
	
	//
	// Compute the sub-sampled image
	//
	for (int j = 0; j < latsamples; j ++) {
	  
	  double y = ((double)j + 0.5)/(double)latsamples * (ymax - ymin) + ymin;
	  
	  for (int i = 0; i < lonsamples; i ++) {
	    
	    double x = ((double)i + 0.5)/(double)lonsamples * (xmax - xmin) + xmin;
	    
	    image[j * lonsamples + i] = models[modelindex]->incircle_radius_at_point(x, y);
	    
	  }
	}
	  
	//
	// Update mean/variance and hist counts
	//
	meann ++;
	for (int i = 0; i < imagesize; i ++) {
	  
	  delta = image[i] - mean[i];
	  mean[i] += delta/(double)meann;
	  variance[i] += delta * (image[i] - mean[i]);
	  
	}

      }
    }
    
    status = reader.step(models, hierarchical, likelihood, norm);
    step ++;
    
    // if (mpi_rank == 0 && (step % 100000) == 0) {
    //   printf("%d\n", step);
    // }
  }
  
  if (status < 0) {
    fprintf(stderr, "error: failed to step through chain history\n");
    return -1;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  //
  // Update local variances first
  //
  for (int i = 0; i < imagesize; i ++) {
    variance[i] /= (double)(meann - 1);
    variance[i] += (mean[i] * mean[i]);
  }

  //
  // Aggregate means
  //
  MPI_Reduce(mean, workspace, imagesize, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (mpi_rank == 0) {

    for (int i = 0; i < imagesize; i ++) {
      mean[i] = workspace[i]/(double)mpi_size;
    }
  
    //
    // Always save the mean
    //
    if (saveimage(output, mean, lonsamples, latsamples) < 0) {
      fprintf(stderr, "error: failed to save mean\n");
      return -1;
    }
  }

  //
  // Aggregate variance
  //
  MPI_Reduce(variance, workspace, imagesize, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (mpi_rank == 0) {

    for (int i = 0; i < imagesize; i ++) {
      variance[i] = workspace[i]/(double)mpi_size - mean[i]*mean[i];
    }

    if (stddev_file != nullptr) {
      for (int i = 0; i < imagesize; i ++) {
	variance[i] = sqrt(variance[i]);
      }
      
      if (saveimage(stddev_file, variance, lonsamples, latsamples) < 0) {
	fprintf(stderr, "error: failed to save std dev\n");
	return -1;
      }
    }
  }

  
  delete [] image;
  delete [] mean;
  delete [] variance;
  delete [] workspace;

  MPI_Finalize();
  
  return 0;
}

static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  " -i|--input <filename>       Input observations file (required)\n"
	  " -o|--output <filename>      Output likelihood file (required)\n"
	  "\n"
	  " -D|--stddev <filename>      Std dev. output\n"
	  "\n"
	  " -W|--lonsamples <int>       No. samples in longitude direction\n"
	  " -H|--latsamples <int>       No. samples in latitude direction\n"
	  "\n"
	  " -t|--thin <int>             Only use every nth model\n"
	  " -s|--skip <int>             Skip first n models\n"
	  "\n"
	  " -h|--help                   Usage\n"
	  "\n",
	  pname);
}


static int saveimage(const char *filename, double *image, int width, int height)
{

  FILE *fp;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    return -1;
  }

  for (int j = 0; j < height; j ++) {
    for (int i = 0; i < width; i ++) {

      fprintf(fp, "%10.6f ", image[width * j + i]);

    }
    fprintf(fp, "\n");
  }

  fclose(fp);
  return 0;
}
