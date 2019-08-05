//
//    TransTessellate2D : A general Trans-dimensional Tessellation program
//    for 2D Cartesian problems.
//
//    Copyright (C) 2014 - 2019 Rhys Hawkins
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

#pragma once
#ifndef cartesianvoronoimodel_hpp
#define cartesianvoronoimodel_hpp

#include <vector>

extern "C" {
  #include "delaunay2d.h"
};

#include "rng.hpp"
#include "generalvoronoicartesianexception.hpp"

class cvcache
{
public:
  cvcache(delaunay2d_t *_model, double _x, double _y) :
    model(_model),
    x(_x),
    y(_y)
  {
  }

  virtual ~cvcache()
  {
  }

  virtual double evaluate(bool force) = 0;

  virtual const int *cached_indices(int &nindices) = 0;
  virtual const double *cached_weights(int &nweights) = 0;

protected:

  delaunay2d_t *model;
  double x, y;

};

class cvcachevoronoi : public cvcache
{
public:
  cvcachevoronoi(delaunay2d_t *model, double x, double y) :
    cvcache(model, x, y),
    cached(false),
    index(-1)
  {
  }

  virtual ~cvcachevoronoi()
  {
  }

  virtual double evaluate(bool force)
  {
    if (force || !cached) {
      index = delaunay2d_nearest(model,
				 0, // include corners (no for voronoi)
				 x, y);
      if (index < 0) {
	delaunay2d_print_points(model);
	throw GENERALVORONOICARTESIANEXCEPTION("Failed to get nearest cell");
      }

      cached = true;
    }

    double v;
    if (delaunay2d_value_of_index(model, index, &v) < 0) {
      int tindex = delaunay2d_nearest(model,
				      0, // include corners (no for voronoi)
				      x, y);

      throw GENERALVORONOICARTESIANEXCEPTION("Failed to get value of index %d (%d)", index, tindex);
    }
    
    return v;
  }

  virtual const int* cached_indices(int &nindices)
  {
    if (cached) {
      nindices = 1;
      return &index;
    }

    nindices = -1;
    return nullptr;
  }

  virtual const double *cached_weights(int &nweights)
  {
    if (cached) {
      nweights = 1;
      weight = 1.0;
      return &weight;
    }

    nweights = -1;
    return nullptr;
  }

private:

  bool cached;
  int index;
  double weight;
    
};

class cvcachedelaunay : public cvcache
{
public:
  cvcachedelaunay(delaunay2d_t *model, double x, double y) :
    cvcache(model, x, y),
    cached(false)
  {
  }

  virtual ~cvcachedelaunay()
  {
  }

  virtual double evaluate(bool force)
  {
    if (force || !cached) {

      if (delaunay2d_find_enclosing_triangle(model,
					     0, // Starting triangle
					     x,
					     y,
					     &indices[0],
					     &indices[1], 
					     &indices[2],
					     &weights[0],
					     &weights[1],
					     &weights[2]) < 0) {
	throw GENERALVORONOICARTESIANEXCEPTION("Failed to get enclosing triangle");
      }

      cached = true;
    }

    double v = 0.0;
    
    for (int i = 0; i < 3; i ++) {
      double t;
      if (delaunay2d_value_of_index(model, indices[i], &t) < 0) {
	throw GENERALVORONOICARTESIANEXCEPTION("Failed to get value of index");
      }
      v += weights[i] * t;
    }

    return v;
  }

  virtual const int* cached_indices(int &nindices)
  {
    if (cached) {
      nindices = 3;
      return indices;
    }

    nindices = -1;
    return nullptr;
  }
  
  virtual const double* cached_weights(int &nweights)
  {
    if (cached) {
      nweights = 3;
      return weights;
    }

    nweights = -1;
    return nullptr;
  }

private:

  bool cached;
  int indices[3];
  double weights[3];
  
};

class cvcachecloughtocher : public cvcache
{
public:
  cvcachecloughtocher(delaunay2d_t *model, double x, double y) :
    cvcache(model, x, y),
    cached(false)
  {
  }

  virtual ~cvcachecloughtocher()
  {
  }

  virtual double evaluate(bool force)
  {
    if (true || force || !cached) {

      double z;
      
      if (delaunay2d_ct_value_at_gradient(model,
					  0, // Starting triangle
					  x,
					  y,
					  &z,
					  &npoints,
					  indices,
					  weights) < 0) {
	throw GENERALVORONOICARTESIANEXCEPTION("Failed to get enclosing triangle");
      }

      cached = true;
    }

    double v = 0.0;
    
    for (int i = 0; i < npoints; i ++) {
      double t;
      if (delaunay2d_value_of_index(model, indices[i], &t) < 0) {
	throw GENERALVORONOICARTESIANEXCEPTION("Failed to get value of index");
      }
      v += weights[i] * t;
    }

    return v;
  }

  virtual const int* cached_indices(int &nindices)
  {
    if (cached) {
      nindices = npoints;
      return indices;
    }

    nindices = -1;
    return nullptr;
  }
  
  virtual const double* cached_weights(int &nweights)
  {
    if (cached) {
      nweights = npoints;
      return weights;
    }

    nweights = -1;
    return nullptr;
  }

private:

  bool cached;

  int npoints;
  int indices[DELAUNAY2D_MAXTRIANGLEEDGES];
  double weights[DELAUNAY2D_MAXTRIANGLEEDGES];
  
};


class cartesianvoronoimodel {
public:

  

  typedef enum {
    VORONOI = 0,
    DELAUNAY = 1,
    CLOUGHTOCHER = 2,
  } parameterization_t;
  
  cartesianvoronoimodel(int _maxpoints,
			double _xmin, double _xmax,
			double _ymin, double _ymax,
			int _type) :
    maxpoints(_maxpoints),
    xmin(_xmin),
    xmax(_xmax),
    ymin(_ymin),
    ymax(_ymax),
    type(VORONOI),
    model(delaunay2d_create(_maxpoints + 4,
			    _xmin,
			    _xmax,
			    _ymin,
			    _ymax))
  {
    switch(_type) {
    case 0:
      type = VORONOI;
      break;

    case 1:
      type = DELAUNAY;
      break;

    case 2:
      type = CLOUGHTOCHER;
      break;

    default:
      throw GENERALVORONOICARTESIANEXCEPTION("Unknown paramerization type");
    }
  }
  
  ~cartesianvoronoimodel()
  {
    delaunay2d_destroy(model);
  }

  int parameterization() const
  {
    return (int)type;
  }
  
  cvcache *create_reference(double x, double y)
  {
    switch (type) {
    case VORONOI:
      return new cvcachevoronoi(model, x, y);
      
    case DELAUNAY:
      return new cvcachedelaunay(model, x, y);
      
    case CLOUGHTOCHER:
      return new cvcachecloughtocher(model, x, y);

    default:
      throw GENERALVORONOICARTESIANEXCEPTION("Invalid parameterization");
    }
  }

  void recompute()
  {
    switch (type) {
    case VORONOI:
      if (delaunay2d_nearest_update(model) < 0) {
	throw GENERALVORONOICARTESIANEXCEPTION("Failed to update nearest\n");
      }
      break;

    case DELAUNAY:
      if (delaunay2d_linear_update(model) < 0) {
	throw GENERALVORONOICARTESIANEXCEPTION("Failed to update linear\n");
      }
      break;

    case CLOUGHTOCHER:
      if (delaunay2d_ct_update(model) < 0) {
	throw GENERALVORONOICARTESIANEXCEPTION("Failed to update C-T\n");
      }
      break;

    default:
      throw GENERALVORONOICARTESIANEXCEPTION("Invalid parameterization");
    }
  }

  void reset()
  {
    delaunay2d_reset(model);
  }

  int ntotalcells() const
  {
    return delaunay2d_npoints(model);
  }

  int nvaluecells() const
  {
    int n = delaunay2d_npoints(model);
    if (type == VORONOI) {
      return n - 4;
    } else {
      return n;
    }
  }

  int nmobilecells() const
  {
    int n = delaunay2d_npoints(model);
    return n - 4;
  }

  //
  // Due to supporting both Delaunay and Voronoi, when model indices
  // are used within the code we need to use the following two
  // functions. The first maps 0 .. n to 4 .. n + 4 for Voronoi
  // parameterizations and 0 .. n to 0 .. n for Delaunay.
  // The difference is caused by the Voronoi parmeterization not
  // using the first four corner points
  //
  int selectvaluecell(int index)
  {
    if (type == VORONOI) {
      return index + 4;
    } else {
      return index;
    }
  }

  //
  // See description of selectvaluecell above. This functions is
  // its inverse.
  //
  int ordinal(int index)
  {
    if (type == VORONOI) {
      if (index < 4) {
	throw GENERALVORONOICARTESIANEXCEPTION("Invalid index");
      }
      return index - 4;
    } else {
      return index;
    }
  }
  
  int selectvaluecell(Rng &random)
  {
    int n = delaunay2d_npoints(model);
    if (type == VORONOI) {
      if (n == 4) {
	return -1;
      }
      return 4 + random.uniform(n - 4);
    } else {
      return random.uniform(n);
    }
  }

  
  
  int selectmobilecell(Rng &random)
  {
    int n = delaunay2d_npoints(model);
    if (n == 4) {
      return -1;
    }
    
    return 4 + random.uniform(n - 4);
  }

  void add_cell(double x, double y, double z)
  {
    bbox2d_t bound;
    if (delaunay2d_add(model,
		       x, y, z, &bound) < 0) {
      throw GENERALVORONOICARTESIANEXCEPTION("Failed to add cell");
    }
  }

  void pop()
  {
    int n = delaunay2d_npoints(model);
    delete_cell(n - 1);
  }

  void delete_cell(int index)
  {
    bbox2d_t bound;
    if (delaunay2d_delete(model, index, &bound) < 0) {
      throw GENERALVORONOICARTESIANEXCEPTION("Failed to remove cell");
    }
  }

  void insert_cell(int index, double x, double y, double z)
  {
    bbox2d_t bound;
    
    if (delaunay2d_add(model, x, y, z, &bound) < 0) {
      throw GENERALVORONOICARTESIANEXCEPTION("Failed to add point");
    }

    if (delaunay2d_shift_replace(model, index) < 0) {
      throw GENERALVORONOICARTESIANEXCEPTION("Failed to reorder points");
    }
  }

  void nearest(double x, double y,
	       double *cx, double *cy,
	       int *cell_index, double *cell_value) const
  {
    int include_corners = !(type == VORONOI);
    
    int i = delaunay2d_nearest(model,
			       include_corners,
			       x, y);
    if (i < 0) {
      throw GENERALVORONOICARTESIANEXCEPTION("Failed to find nearest");
    }

    delaunay2d_point_of_index(model, i, cx, cy);
    delaunay2d_value_of_index(model, i, cell_value);
    *cell_value = i;
  }

  double value_at_point(double x, double y, int &t0) const
  {
    double z;
    
    switch (type) {
    case VORONOI:
      if (t0 < 0 || t0 >= delaunay2d_npoints(model)) {
	// Reset silently.
	t0 = 0;
      }
      t0 = delaunay2d_nearest_value_at(model, t0,
				       x, y,
				       &z);
      if (t0 < 0) {
	throw GENERALVORONOICARTESIANEXCEPTION("Failed to find nearest point %10.6f %10.6f %d",
					       x, y, 
					       delaunay2d_npoints(model));
      }
      break;
      
    case DELAUNAY:
      if (t0 < 0 || t0 >= delaunay2d_ntriangles(model)) {
	t0 = 0;
      }
      t0 = delaunay2d_linear_value_at(model, t0,
				      x, y,
				      &z);
      if (t0 < 0) {
	throw GENERALVORONOICARTESIANEXCEPTION("Failed to find linear value point %10.6f %10.6f %d",
					       x, y,
					       delaunay2d_npoints(model));
      }
      
      break;
      
    case CLOUGHTOCHER:
      if (t0 < 0 || t0 >= delaunay2d_ntriangles(model)) {
	t0 = 0;
      }
      t0 = delaunay2d_ct_value_at(model, t0,
				  x, y,
				  &z);
      if (t0 < 0) {
	throw GENERALVORONOICARTESIANEXCEPTION("Failed to find ct value point %10.6f %10.6f %d",
					       x, y,
					       delaunay2d_npoints(model));
      }
      break;

    default:
      throw GENERALVORONOICARTESIANEXCEPTION("Invalid parameterization type");

    }

    return z;
  }

  double circumcircle_radius_at_point(double x, double y) const
  {
    return delaunay2d_circumcircle_radius_at_point(model, x, y);
  }

  double incircle_radius_at_point(double x, double y) const
  {
    return delaunay2d_incircle_radius_at_point(model, x, y);
  }

  double value_at_index(int index) const
  {
    double z;

    if (delaunay2d_value_of_index(model, index, &z) < 0) {
      throw GENERALVORONOICARTESIANEXCEPTION("Failed to get value");
    }

    return z;
  }

  void range(double &vmin, double &vmax)
  {
    int start = 0;
    int n = delaunay2d_npoints(model);
    if (type == VORONOI) {
      start = 4;
    }
    
    for (int i = start; i < n; i ++) {

      double z;
      if (delaunay2d_value_of_index(model, i, &z) < 0) {
	throw GENERALVORONOICARTESIANEXCEPTION("Failed to get value");
      }

      if (z < vmin) {
	vmin = z;
      }

      if (z > vmax) {
	vmax = z;
      }
    }
  }

  void set_value_at_index(int index, double z)
  {
    if (delaunay2d_set_value_of_index(model, index, z) < 0) {
      throw GENERALVORONOICARTESIANEXCEPTION("Failed to set value %d", index);
    }
  }

  void position_at_index(int index, double *x, double *y)
  {
    if (delaunay2d_point_of_index(model, index, x, y) < 0) {
      throw GENERALVORONOICARTESIANEXCEPTION("Failed to get position");
    }
  }

  void set_position_at_index(int index, double x, double y)
  {
    double oldv;
    bbox2d_t bound;
    
    if (delaunay2d_value_of_index(model, index, &oldv) < 0) {
      throw GENERALVORONOICARTESIANEXCEPTION("Failed to get value");
    }

    if (delaunay2d_delete(model, index, &bound) < 0) {
      throw GENERALVORONOICARTESIANEXCEPTION("Failed to remove point");
    }

    if (delaunay2d_add(model, x, y, oldv, &bound) < 0) {
      throw GENERALVORONOICARTESIANEXCEPTION("Failed to add point");
    }

    if (delaunay2d_shift_replace(model, index) < 0) {
      throw GENERALVORONOICARTESIANEXCEPTION("Failed to reorder points");
    }
  }

  bool save(const char *filename)
  {
    return delaunay2d_save(model, filename) == 0;
  }

  bool load(const char *filename)
  {
    delaunay2d_t *new_model = delaunay2d_load(filename);
    if (new_model == NULL) {
      return false;
    }

    if (model != NULL) {
      delaunay2d_destroy(model);
    }
    model = new_model;
    return true;
  }

  int maxpoints;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  parameterization_t type;
  
  delaunay2d_t *model;

};

#endif // cartesianvoronoimodel_hpp
  
  
