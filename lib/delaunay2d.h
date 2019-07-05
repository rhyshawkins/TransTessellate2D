//
//    Trans-d Voronoi/Delaunay example code.
//    
//    Copyright (C) 2018 Rhys Hawkins
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

#ifndef delaunay_h
#define delaunay_h

#include "bbox2d.h"

#define DELAUNAY2D_MAXTRIANGLEEDGES 64

typedef struct _delaunay2d delaunay2d_t;

delaunay2d_t *
delaunay2d_create(int maxpoints,
		  double xmin,
		  double xmax,
		  double ymin,
		  double ymax);

void
delaunay2d_reset(delaunay2d_t *d);

void 
delaunay2d_destroy(delaunay2d_t *d);

delaunay2d_t *
delaunay2d_load(const char *filename);

int
delaunay2d_save(const delaunay2d_t *d,
		const char *filename);

int
delaunay2d_clone(const delaunay2d_t *src,
		 delaunay2d_t *dst);

int 
delaunay2d_add(delaunay2d_t *d, 
	       double x,
	       double y,
	       double z,
	       bbox2d_t *bound);

int 
delaunay2d_delete(delaunay2d_t *d,
		  int pi,
		  bbox2d_t *bound);

int 
delaunay2d_shift_replace(delaunay2d_t *d,
			 int pi);

int 
delaunay2d_find_enclosing_triangle(const delaunay2d_t *d,
				   int t0,
				   double px,
				   double py,
				   int *pa,
				   int *pb, 
				   int *pc,
				   double *ba,
				   double *bb,
				   double *bc);

int
delaunay2d_nearest_update(delaunay2d_t *d);

int
delaunay2d_linear_update(delaunay2d_t *d);

int
delaunay2d_nearest_value_at(const delaunay2d_t *d,
			    int t0,
			    double px,
			    double py,
			    double *z);

int
delaunay2d_linear_value_at(const delaunay2d_t *d,
			   int t0,
			   double px,
			   double py,
			   double *z);

int
delaunay2d_ct_update(delaunay2d_t *d);

int 
delaunay2d_ct_value_at(const delaunay2d_t *d,
		       int t0,
		       double px,
		       double py,
		       double *z);

int 
delaunay2d_ct_value_at_gradient(const delaunay2d_t *d,
				int t0,
				double px,
				double py,
				double *z,
				int *npoints,
				int *indices,
				double *weights);

int 
delaunay2d_nearest(const delaunay2d_t *d,
		   int include_corners,
		   double px,
		   double py);

int 
delaunay2d_nearest_from(const delaunay2d_t *d,
			int i,
			int include_corners,
			double px,
			double py);

int 
delaunay2d_point_of_index(const delaunay2d_t *d,
			  int i,
			  double *px,
			  double *py);

int
delaunay2d_value_of_index(const delaunay2d_t *d,
			  int i,
			  double *z);

int
delaunay2d_set_value_of_index(const delaunay2d_t *d,
			      int i,
			      double z);


int 
delaunay2d_index_of_point(const delaunay2d_t *d,
			  double x,
			  double y);

int
delaunay2d_npoints(const delaunay2d_t *d);

int
delaunay2d_ntriangles(const delaunay2d_t *d);

int
delaunay2d_polygon_bound(const delaunay2d_t *d,
			 int i,
			 bbox2d_t *bound);

int 
delaunay2d_validate_circumcircles(const delaunay2d_t *d);

int
delaunay2d_validate_delaunay(const delaunay2d_t *d);

int
delaunay2d_validate_neighbours(const delaunay2d_t *d);

int
delaunay2d_validate_nonintersecting(const delaunay2d_t *d);

int
delaunay2d_validate_edges(const delaunay2d_t *d);

int 
delaunay2d_validate_delaunay(const delaunay2d_t *d);

void
delaunay2d_print_triangles(const delaunay2d_t *d);

void
delaunay2d_print_points(const delaunay2d_t *d);

void
delaunay2d_print_edges(const delaunay2d_t *d);

int
delaunay2d_save_geo(const delaunay2d_t *d, const char *filename);

int
delaunay2d_save_cc_geo(const delaunay2d_t *d, const char *filename);

int
triangle_circumcircle(double x1, double y1,
		      double x2, double y2,
		      double x3, double y3,
		      double *cx, double *cy,
		      double *cr2);

int
triangle_incircle_radius(double x1, double y1,
			 double x2, double y2,
			 double x3, double y3,
			 double *ir);

int 
point_in_triangle(double px,
		  double py,
		  double x1, double y1,
		  double x2, double y2,
		  double x3, double y3);

int delaunay2d_max_edges(const delaunay2d_t *d);

int delaunay2d_max_triangle_edges(const delaunay2d_t *d);

double delaunay2d_circumcircle_radius_at_point(const delaunay2d_t *d, double x, double y);

double delaunay2d_incircle_radius_at_point(const delaunay2d_t *d, double x, double y);

#endif /* delaunay2d_h */
