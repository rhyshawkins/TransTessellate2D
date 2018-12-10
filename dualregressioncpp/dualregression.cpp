#include <vector>
#include <stdio.h>
#include <math.h>

#include "genericinterface.hpp"

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
    *nhierarchical = 2;

    return 0;
  }

  int gvcart_loaddata_(int *n,
		       const char *filename,
		       gvcart_addobservation_t addobs)
  {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
      return -1;
    }

    int nobs;
    if (fscanf(fp, "%d\n", &nobs) != 1) {
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
	return -1;
      }

      if (obs[i].type < 0 || obs[i].type > 1) {
	fprintf(stderr, "error: invalid data type %d (must be 0 or 1)\n", obs[i].type);
	return -1;
      }

      int n;
      n = 1;
      
      if (addobs(&n, &obs[i].type, &obs[i].lon, &obs[i].lat) < 0) {
	return -1;
      }

      if (obs[i].sigma > 0.0) {
	obs[i].logsigma = log(obs[i].sigma);
      } else {
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
				 double *unused,
				 double *prediction)
  {
    if ((*observation) < 0 || (*observation) >= (int)obs.size()) {
      return -1;
    }

    prediction[0] = value[0];

    return 0;
  }

  int gvcart_compute_likelihood_(int *nmodel,
				 int *nhierarchical,
				 int *nobservation,
				 double *hierarchical,
				 double *predictions,
				 double *residuals,
				 double *unused,
				 double *_like,
				 double *_norm)

  {
    constexpr double NORMSCALE = 0.9189385332046727; //0.5*log(2.0*M_PI);
    double sum = 0.0;
    double norm = 0.0;
    double loghierarchical[2];
    for (int i = 0; i < *nhierarchical; i ++) {
      loghierarchical[i] = log(hierarchical[i]);
    }

    for (int i = 0; i < (*nobservation); i ++) {
      double res = predictions[i] - obs[i].value;
      residuals[i] = res;

      double n = obs[i].sigma;
      double ln = NORMSCALE + obs[i].logsigma;

      if (obs[i].type == 0) {
	n *= hierarchical[0];
	ln += loghierarchical[0];
      } else {
	n *= hierarchical[0];
	ln += loghierarchical[0];
      }	
	
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

  
}

