#include <vector>
#include <stdio.h>
#include <math.h>

#include "genericinterface.hpp"

#include "tidesynthetic.hpp"

typedef enum {
  SEA_MODEL = 0,
  LAND_MODEL = 1
} sealand_model_t;

typedef enum {
  TIDE_OBS = 1,
  LAND_OBS = 2,
  SEA_OBS = 3
} sealand_obs_t;

struct observation {
  double lon;
  double lat;
  int type;
  double value;
  double sigma;
  double logsigma;
};

static std::vector<struct observation> obs;

extern "C" {

  int gvcart_initialise_(int *nmodels,
			 int *nhierarchical)
  {
    *nmodels = 2;
    *nhierarchical = 3;

    synthetic_add("HorizontalCosine1", tidesynthetic_horizontal_cosine1);
    synthetic_add("HorizontalCosine2", tidesynthetic_horizontal_cosine2);
    synthetic_add("HorizontalCosine3", tidesynthetic_horizontal_cosine3);

    synthetic_add("VerticalCosine1", tidesynthetic_vertical_cosine1);
    synthetic_add("VerticalCosine2", tidesynthetic_vertical_cosine2);
    synthetic_add("VerticalCosine3", tidesynthetic_vertical_cosine3);

    synthetic_add("TasSea", tidesynthetic_tas_sea);
    synthetic_add("TasLand", tidesynthetic_tas_land);
      

    return 0;
  }

  int gvcart_loaddata_(int *filename_len,
		       const char *filename,
		       gvcart_addobservation_t addobs)
  {
    int n;
    int mi[2];
    double xs[2];
    double ys[2];
    
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
      fprintf(stderr, "error: failed to open file %s for reading\n", filename);
      return -1;
    }

    int nobs;
    if (fscanf(fp, "%d\n", &nobs) != 1) {
      fprintf(stderr, "error: failed to read no. entries\n");
      return -1;
    }

    obs.resize(nobs);
    for (int i = 0; i < nobs; i ++) {
      if (fscanf(fp, "%lf %lf %d %lf %lf\n",
		 &obs[i].lon,
		 &obs[i].lat,
		 &obs[i].type,
		 &obs[i].value,
		 &obs[i].sigma) != 5) {
	fprintf(stderr, "error: failed to parse observation %d/%d\n", i, nobs);
	return -1;
      }

      switch(obs[i].type) {
      case TIDE_OBS:
	//
	// For a prediction we need the value of the sea and land model
	// at a point.
	//
	n = 2;
	mi[0] = SEA_MODEL;
	mi[1] = LAND_MODEL;

	xs[0] = obs[i].lon;
	xs[1] = obs[i].lon;

	ys[0] = obs[i].lat;
	ys[1] = obs[i].lat;
	break;

      case LAND_OBS:
	n = 1;
	mi[0] = LAND_MODEL;

	xs[0] = obs[i].lon;
	ys[0] = obs[i].lat;
	break;

      case SEA_OBS:
	n = 1;
	
	mi[0] = SEA_MODEL;

	xs[0] = obs[i].lon;
	ys[0] = obs[i].lat;
	break;
	
      default:
	fprintf(stderr, "error: invalid type %d\n", obs[i].type);
	return -1;
      }

      if (addobs(&n, mi, xs, ys) < 0) {
	fprintf(stderr, "error: failed to add observation\n");
	return -1;
      }

      if (obs[i].sigma > 0.0) {
	obs[i].logsigma = log(obs[i].sigma);
      } else {
	fprintf(stderr, "error: invalid sigma %f\n", obs[i].sigma);
	return -1;
      }
	  
     }

    fclose(fp);
    return 0;
  }

  int gvcart_compute_prediction_(int *nmodels,
				 int *observation,
				 int *npoints,
				 const double *value,
				 double *weight,
				 double *prediction)
  {
    if ((*observation) < 0 || (*observation) >= (int)obs.size()) {
      return -1;
    }

    switch (obs[*observation].type) {
    case TIDE_OBS: // Tide Gauge
      if ((*npoints) != 2) {
	fprintf(stderr, "error: invalid number of points for tide observation %d\n",
		(*npoints));
	return -1;
      }

      // Sea model is value[0] and land model is value[1] (see load data)
      prediction[0] = value[0] - value[1];
      weight[0] = 1.0;
      weight[1] = -1.0;
      break;
      
    case LAND_OBS: // Land GPS
      prediction[0] = value[0];
      weight[0] = 1.0;
      break;

    case SEA_OBS: // Sea level (Satellite)
      prediction[0] = value[0];
      weight[0] = 1.0;
      break;

    default:
      return -1;
    }
    return 0;
  }

  int gvcart_compute_likelihood_(int *nmodel,
				 int *nhierarchical,
				 int *nobservation,
				 double *hierarchical,
				 double *predictions,
				 double *residuals,
				 double *weight,
				 double *_like,
				 double *_norm)

  {
    constexpr double NORMSCALE = 0.9189385332046727; //0.5*log(2.0*M_PI);
    double sum = 0.0;
    double norm = 0.0;
    double loghierarchical[3];
    for (int i = 0; i < *nhierarchical; i ++) {
      loghierarchical[i] = log(hierarchical[i]);
    }

    for (int i = 0; i < (*nobservation); i ++) {
      double res = predictions[i] - obs[i].value;
      residuals[i] = res;

      double n = obs[i].sigma;
      double ln = NORMSCALE + obs[i].logsigma;
      
      switch (obs[i].type) {
      case TIDE_OBS:
	n *= hierarchical[0];
	ln += loghierarchical[0];
	break;
	
      case LAND_OBS:
	n *= hierarchical[1];
	ln += loghierarchical[1];
	break;

      case SEA_OBS:
	n *= hierarchical[2];
	ln += loghierarchical[2];
	break;

      default:
	fprintf(stderr, "gvcart_compute_likelihood: unknown observation type\n");
	return -1;
      }	
	
      weight[i] = res/(n*n);
      
      sum += (res*res)/(2.0*n*n);
      norm += ln;
    }

    *_like = sum;
    *_norm = norm;

    return 0;
  }

  int gvcart_savedata_(int *n,
		       const char *filename,
		       double *noiselevel,
		       int *nobservations,
		       double *predictions)
  {
    FILE *fp;

    fp = fopen(filename, "w");
    if (fp == NULL) {
      return -1;
    }

    fprintf(fp, "%d\n", (int)obs.size());
    int i = 0;
    for (auto &o : obs) {
      fprintf(fp, "%16.9f %16.9f %d  %16.9f %16.9f\n",
	      o.lon, o.lat, o.type, predictions[i], noiselevel[0]);
      i ++;
    }

    fclose(fp);
    return 0;
  }

  //
  // Derived variable is predicted tide gauge, i.e. Sea rate model - Land rate model
  //
  double gvcart_compute_derived_(int *nmodels,
				 double *x, double *y,
				 double *values)
  {
    return values[SEA_MODEL] - values[LAND_MODEL];
  }
  
}

