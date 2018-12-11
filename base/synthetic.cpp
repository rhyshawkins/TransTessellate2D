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

#include "synthetic.hpp"

double synthetic_constant(double nx, double ny)
{
  return 0.5;
}

double synthetic_eastwest(double nx, double ny)
{
  if (nx < 0.0) {
    return 0.25;
  } else {
    return 0.75;
  }
}

double synthetic_northsouth(double nx, double ny)
{
  if (ny < 0.0) {
    return 0.33;
  } else {
    return 0.66;
  }
}

double synthetic_2x2(double nx, double ny)
{
  if (nx < 0.0) {
    if (ny < 0.0) {
      return 0.25;
    } else {
      return 0.625;
    }
  } else {
    if (ny < 0.0) {
      return 0.5;
    } else {
      return 0.375;
    }
  }
}

double synthetic_gaussian(double nx, double ny)
{
  double r2 = nx*nx + ny*ny;
  return exp(-r2/0.5);
}

double synthetic_franke(double nx, double ny)
{
  //
  // Franke function from
  //
  //@TechReport{Mann:1998:A,
  //  author = 	 {Mann S.},
  //  title = 	 {Cubic precision {C}lough-{T}ocher interpolation},
  //  institution =  {Computer Science Department, University of Waterloo},
  //  year = 	 {1998},
  //  key = 	 {CS-98-15},
  //}
  //

  double xi = (nx + 1.0)/2.0;
  double eta = (ny + 1.0)/2.0;
    
  double dx;
  double dy;

  double v = 0.0;

  dx = 9.0*xi - 2.0;
  dy = 9.0*eta - 2.0;

  v = 0.75 * exp(-(dx*dx + dy*dy)/4.0);

  dx = 9.0*xi + 1.0;
  dy = 9.0*eta + 1.0;
  
  v += 0.75 * exp(-(dx*dx)/49.0 - (dy*dy)/10.0);

  dx = 9.0*xi - 7.0;
  dy = 9.0*eta - 3.0;
  
  v += 0.5 * exp(-(dx*dx + dy*dy)/4.0);

  dx = 9.0*xi - 4.0;
  dy = 9.0*eta - 7.0;
  
  v -= 0.2 * exp(-(dx*dx) - (dy*dy));
  
  return v;
}

double synthetic_pseudofranke(double nx, double ny)
{
  //
  // Franke function from 
  //
  //@TechReport{Mann:1998:A,
  //  author = 	 {Mann S.},
  //  title = 	 {Cubic precision {C}lough-{T}ocher interpolation},
  //  institution =  {Computer Science Department, University of Waterloo},
  //  year = 	 {1998},
  //  key = 	 {CS-98-15},
  //}
  //

  double v = synthetic_franke(nx, ny);

  if (v < 0.0) {
    return 0.0;
  } else if (v < 0.1) {
    return 0.1;
  } else if (v < 0.3) {
    return 0.3;
  } else if (v < 0.6) {
    return 0.6;
  } else {
    return 0.9;
  }
}

double synthetic_tesselation(double nx, double ny)
{
  double horizontal_threshold_c = -0.25*nx;
  double horizontal_threshold_l = -0.25*nx - 0.25;

  double horizontal_threshold_wu = -0.25*nx + 0.5;

  double vertical_threshold_l = 0.25*ny - 0.25;
  double vertical_threshold_u = 0.25*ny + 0.25;

  double vertical_threshold_wl = 0.25*ny - 0.5;
  double vertical_threshold_wu = 0.25*ny + 0.5;

  if (ny > horizontal_threshold_wu) {
    if (nx > vertical_threshold_wl && nx < vertical_threshold_wu) {
      return 0.0;
    }
  } else if (ny < horizontal_threshold_c) {

    if (nx > vertical_threshold_u) {
      return 0.7;
    } else {

      if (ny < horizontal_threshold_l &&
	  nx < vertical_threshold_l) {

	return 1.0;
      }
    }
  }

  return 0.3;
}



void synthetic_add(const char *name,
		   synthetic_model_f model_f)
{
  synthetic_models.insert(std::pair<std::string, synthetic_model_f>(name, model_f));
}

std::map<std::string, synthetic_model_f> synthetic_models = { {"Constant", synthetic_constant},
							      {"EastWest", synthetic_eastwest},
							      {"NorthSouth", synthetic_northsouth},
							      {"2x2", synthetic_2x2},
							      {"Gaussian", synthetic_gaussian},
							      {"Franke", synthetic_franke},
							      {"PseudoFranke", synthetic_pseudofranke},
							      {"Tesselation", synthetic_tesselation}};

