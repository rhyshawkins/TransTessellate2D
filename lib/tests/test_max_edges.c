#include <stdio.h>
#include <stdlib.h>

#include "delaunay2d.h"

static delaunay2d_t *mkrandom(int npoints);

int main(int argc, char *argv[])
{
  delaunay2d_t *d;
  int i;
  int e;
  int te;
  
  int maxe = 0;
  int maxte = 0;
  
  for (i = 0; i < 10000; i ++) {
    d = mkrandom(10000);

    e = delaunay2d_max_edges(d);
    te = delaunay2d_max_triangle_edges(d);
    
    printf("%6d Edges %6d Triangle Edges %6d\n",
	   i, 
	   e, te);

    if (e > maxe) {
      maxe = e;
    }
    if (te > maxte) {
      maxte = te;
    }
    
    delaunay2d_destroy(d);
  }

  printf("---\n");
  printf("Total Edges %6d Triangle Edges %6d\n",
	 maxe, maxte);
  
  return 0;
}


static delaunay2d_t *mkrandom(int npoints)
{
  int i;
  delaunay2d_t *d;
  double x;
  double y;
  bbox2d_t bound;

  d = delaunay2d_create(npoints + 4,
			-1.0, 1.0,
			-1.0, 1.0);
  
  if (d == NULL) {
    return NULL;
  }

  for (i = 0; i < npoints; i ++) {

    x = (double)random()/(double)RAND_MAX * 2.0 - 1.0;
    y = (double)random()/(double)RAND_MAX * 2.0 - 1.0;
    
    if (delaunay2d_add(d, x, y, 0.0, &bound) < 0) {
      fprintf(stderr, "error: failed to add point\n");
      return NULL;
    }

  }

  return d;
}
