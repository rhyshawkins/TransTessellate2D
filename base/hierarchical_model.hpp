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
#ifndef hierarchical_model_hpp
#define hierarchical_model_hpp

#include <ostream>
#include <istream>

class hierarchical_model {
public:

  hierarchical_model(int nhierarchical);
  virtual ~hierarchical_model();

  virtual int get_nhierarchical() const;
  virtual void set(int i, double v);
  virtual double get(int i) const;

  virtual bool write(std::ostream &s) const;
  virtual bool read(std::istream &s);

  virtual double *data();

protected:

  int nhierarchical;
  double *hierarchical;
  
};

class singlescaling_hierarchical_model : public hierarchical_model {
public:

  singlescaling_hierarchical_model(double lambda = 1.0);
  ~singlescaling_hierarchical_model();

};

#endif // hierarchical_model_hpp
  
