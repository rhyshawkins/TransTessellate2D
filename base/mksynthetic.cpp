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

#include <math.h>

#include <map>
#include <string>

#include <getopt.h>
#include <string.h>

#include "generalvoronoicartesianobservations.hpp"
#include "rng.hpp"

#include "genericinterface.hpp"

#include "synthetic.hpp"
#include "pathutil.hpp"

static char short_options[] = "i:x:X:y:Y:o:O:I:m:ln:S:W:H:s:f:h";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"xmin", required_argument, 0, 'x'},
  {"xmax", required_argument, 0, 'X'},
  {"ymin", required_argument, 0, 'y'},
  {"ymax", required_argument, 0, 'Y'},
  
  {"output", required_argument, 0, 'o'},
  {"output-true", required_argument, 0, 'O'},

  {"model", required_argument, 0, 'm'},
  {"list-models", no_argument, 0, 'l'},
  
  {"noise", required_argument, 0, 'n'},
  {"seed", required_argument, 0, 'S'},

  {"image-output", required_argument, 0, 'I'},
  {"image-width", required_argument, 0, 'W'},
  {"image-height", required_argument, 0, 'H'},

  {"scale", required_argument, 0, 's'},
  {"offset", required_argument, 0, 'f'},
  
  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};

static GeneralVoronoiCartesianObservations *observations = nullptr;

static int addobservation(int *npoints,
			  int *modelindices,
			  double *lons,
			  double *lats)
{
  observations->add(*npoints, modelindices, lons, lats);

  return 0;
}


static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;

  int nmodels;
  int nhierarchical;
  
  char *source_data;

  double xmin;
  double xmax;
  double ymin;
  double ymax;
  
  char *output_file;
  char *output_true;
  int seed;
  double noise_sigma;
  std::vector<char*> model_names;

  char *image_output;
  int image_width;
  int image_height;

  double scale;
  double offset;
    
  //
  // Defaults
  //
  source_data = nullptr;

  xmin = -1.0;
  xmax = 1.0;
  ymin = -1.0;
  ymax = 1.0;
    
  output_file = nullptr;
  output_true = nullptr;
  noise_sigma = 0.1;
  seed = 983;

  image_output = nullptr;
  image_width = 128;
  image_height = 64;

  scale = 1.0;
  offset = 0.0;

  if (gvcart_initialise_(&nmodels, &nhierarchical) < 0) {
    fprintf(stderr, "error: gvcart initialization failure\n");
    return -1;
  }

  if (nmodels <= 0) {
    fprintf(stderr, "error: invalid number of models: %d\n", nmodels);
    return -1;
  }

  if (nhierarchical <= 0) {
    fprintf(stderr, "error: invalid number of hierarchical parameters: %d\n", nhierarchical);
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
      source_data = optarg;
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

    case 'o':
      output_file = optarg;
      break;

    case 'O':
      output_true = optarg;
      break;
      
    case 'I':
      image_output = optarg;
      break;

    case 'm':
      model_names.push_back(optarg);
      break;

    case 'l':
      {
	fprintf(stderr, "Models:\n");
	for (auto &mn : synthetic_models) {
	  fprintf(stderr, "  `%s'\n", mn.first.c_str());
	}
	return -1;
      }
      break;

    case 'n':
      noise_sigma = atof(optarg);
      if (noise_sigma <= 0.0) {
	fprintf(stderr, "error: noise must be greater than 0\n");
	return -1;
      }
      break;

    case 'S':
      seed = atoi(optarg);
      break;

    case 'W':
      image_width = atoi(optarg);
      if (image_width < 16) {
	fprintf(stderr, "error: image width must be 16 or greater\n");
	return -1;
      }
      break;

    case 'H':
      image_height = atoi(optarg);
      if (image_height < 16) {
	fprintf(stderr, "error: image_height must be 16 or greater\n");
	return -1;
      }
      break;

    case 's':
      scale = atof(optarg);
      if (scale <= 0.0) {
	fprintf(stderr, "error: scale must be greater than 0\n");
	return -1;
      }
      break;

    case 'f':
      offset = atof(optarg);
      break;

    case 'h':
    default:
      usage(argv[0]);
      return -1;

    }
  }
  
  if (source_data == nullptr) {
    fprintf(stderr, "error: required parameter input not set\n");
    return -1;
  }

  if (output_file == nullptr) {
    fprintf(stderr, "error: required parameter output not set\n");
    return -1;
  }

  std::vector<synthetic_model_f> models;
  
  for (int mi = 0; mi < nmodels; mi ++) {

    if (mi < (int)model_names.size()) {
      std::map<std::string, synthetic_model_f>::const_iterator i = synthetic_models.find(model_names[mi]);
      
      if (i == synthetic_models.end()) {
	fprintf(stderr, "error: invalid model name: %s\n", model_names[mi]);
	return -1;
      }

      printf("Model %2d : %s\n", mi, model_names[mi]);
      models.push_back(i->second);
    } else {
      std::map<std::string, synthetic_model_f>::const_iterator i = synthetic_models.find("Constant");
      
      if (i == synthetic_models.end()) {
	fprintf(stderr, "error: invalid model name Constant");
	return -1;
      }
      
      printf("Model %2d : Constant\n", mi);
      models.push_back(i->second);
    }
  }

  //
  // Load data
  //
  int n = strlen(source_data);

  observations = new GeneralVoronoiCartesianObservations(nmodels);
  
  if (gvcart_loaddata_(&n, source_data, addobservation) < 0) {
    fprintf(stderr, "error: failed to load data from %s\n", source_data);
    return -1;
  }

  printf("%d observations\n", (int)observations->obs.size());

  //
  // Compute predictions
  //
  int oi = 0;
  std::vector<double> predictions;
  for (auto &o : observations->obs) {

    //
    // Look up synthetic model values
    //
    int npoints = o.xs.size();


    for (int pi = 0; pi < npoints; pi ++) {
      int mi = o.mi[pi];

      if (o.xs[pi] < xmin || o.xs[pi] > xmax) {
	fprintf(stderr, "error: x coord out of range: %10.6f (%10.6f .. %10.6f)\n",
		o.xs[pi], xmin, xmax);
	return -1;
      }

      if (o.ys[pi] < ymin || o.ys[pi] > ymax) {
	fprintf(stderr, "error: y coord out of range: %10.6f (%10.6f .. %10.6f)\n",
		o.ys[pi], ymin, ymax);
	return -1;
      }
	
      double nx = (o.xs[pi] - xmin)/(xmax - xmin)*2.0 - 1.0;
      double ny = (o.ys[pi] - ymin)/(ymax - ymin)*2.0 - 1.0;
	
      o.values[pi] = models[mi](nx, ny) * scale + offset;
    }
    
    //
    // Compute prediction
    //
    if (gvcart_compute_prediction_(&nmodels,
				   &oi,
				   &npoints,
				   o.values.data(),
				   o.weights.data(),
				   &o.pred)) {
      fprintf(stderr, "error: received error from compute prediction\n");
      return -1;
    }

    predictions.push_back(o.pred);
    oi ++;
  }				

  //
  // Save true observations
  //
  if (output_true) {

    double noise = 0.0;
    int nobs = observations->obs.size();
    int n = strlen(output_true);
    
    if (gvcart_savedata_(&n,
			 output_true,
			 &noise,
			 &nobs,
			 predictions.data()) < 0) {
      fprintf(stderr, "error: failed to save synthetic true observations\n");
      return -1;
    }
  }

  //
  // Add noise
  //
  Rng random(seed);
  
  for (auto &p : predictions) {
    p += random.normal(noise_sigma);
  }

  //
  // Save noisy observations
  //
  int nobs = observations->obs.size();

  n = strlen(output_file);
  if (gvcart_savedata_(&n,
		       output_file,
		       &noise_sigma,
		       &nobs,
		       predictions.data()) < 0) {
    fprintf(stderr, "error: failed to save synthetic observations\n");
    return -1;
  }
  
  if (image_output != nullptr) {

    if (nmodels == 1) {

      FILE *fp_image = fopen(image_output, "w");
      if (fp_image == NULL) {
	fprintf(stderr, "error: failed to create image output file\n");
	return -1;
      }
      
      for (int j = 0; j < image_height; j ++) {
	
	double ny = 2.0*((double)j + 0.5)/(double)(image_height) - 1.0;
	
	for (int i = 0; i < image_width; i ++) {
	  
	  double nx = 2.0*((double)i + 0.5)/(double)(image_width) - 1.0;
	  
	  fprintf(fp_image, "%15.9f ", models[0](nx, ny) * scale + offset);
	  
	}
	
	fprintf(fp_image, "\n");
      }
    } else {
      
      for (int mi = 0; mi < nmodels; mi ++) {
	
	char filename[1024];
	mkmodelpath(mi, NULL, image_output, filename);
	
	FILE *fp_image = fopen(filename, "w");
	if (fp_image == NULL) {
	  fprintf(stderr, "error: failed to create image output file\n");
	  return -1;
	}
	
	for (int j = 0; j < image_height; j ++) {
	  
	  double ny = 2.0*((double)j + 0.5)/(double)(image_height) - 1.0;
	  
	  for (int i = 0; i < image_width; i ++) {
	    
	    double nx = 2.0*((double)i + 0.5)/(double)(image_width) - 1.0;
	    
	    fprintf(fp_image, "%15.9f ", models[mi](nx, ny) * scale + offset);
	    
	  }
	  
	  fprintf(fp_image, "\n");
	}

	fclose(fp_image);
      }
    }
  }

  delete observations;
  
  return 0;
}

static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of\n"
	  "\n"
	  " -i | --input <filename>           Input data to base recompute with synthetic model\n"
	  " -o | --output <filename>          Output file to write synthetic noisy observations to\n"
	  " -O | --output-true <filename>     Output file to write synthetic true observations to\n"
	  "\n"
	  " -m | --model <name>               Synthetic model to use\n"
	  " -l | --list-models                List available synthetic models and exit\n"
	  "\n"
	  " -n | --noise <float>              Std dev of gaussian noise to add to observations\n"
	  " -S | --seed <int>                 Random seed\n"
	  "\n"
	  " -I | --image-output <filename>    Write synthetic model image\n"
	  " -W | --image-width <int>          Image width\n"
	  " -H | --image-height <int>         Image height\n"
	  "\n"
	  " -h | --help                       Usage\n"
	  "\n",
	  pname);
}
