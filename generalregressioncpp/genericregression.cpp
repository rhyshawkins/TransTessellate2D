#include <vector>
#include <stdio.h>
#include <math.h>

#include "genericinterface.hpp"

struct observation {
  double x;
  double y;
  double value;
  double sigma;
};

static std::vector<struct observation> obs;

extern "C" {

  int gvcart_initialise_(int *nmodels,
			 int *nhierarchical)
  {
    *nmodels = 1;
    *nhierarchical = 1;

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
      if (fscanf(fp, "%lf %lf %lf %lf\n",
		 &obs[i].x,
		 &obs[i].y,
		 &obs[i].value,
		 &obs[i].sigma) != 4) {
	return -1;
      }

      int mi;
      int n;
      
      n = 1;
      mi = 0;
      if (addobs(&n, &mi, &obs[i].x, &obs[i].y) < 0) {
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
    prediction[0] = value[0];
    weight[0] = 1.0;
    return 0;
  }

  int gvcart_compute_likelihood_(int *nmodels,
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

    for (int i = 0; i < (*nobservation); i ++) {
      double res = predictions[i] - obs[i].value;
      residuals[i] = res;
      double n = hierarchical[0] * obs[i].sigma;

      weight[i] = res/(n*n);
      
      sum += (res*res)/(2.0*n*n);
      norm += log(NORMSCALE*n);
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
      fprintf(fp, "%16.9f %16.9f %16.9f %16.9f\n",
	      o.x, o.y, predictions[i], noiselevel[0]);
      i ++;
    }

    fclose(fp);
    return 0;
  }

  
}

