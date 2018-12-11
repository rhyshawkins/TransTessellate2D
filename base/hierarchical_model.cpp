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

#include "hierarchical_model.hpp"

#include "generalvoronoicartesianexception.hpp"

hierarchical_model::hierarchical_model(int _nhierarchical) :
  nhierarchical(_nhierarchical),
  hierarchical(new double[_nhierarchical])
{
}

hierarchical_model::~hierarchical_model()
{
  delete [] hierarchical;
}

int
hierarchical_model::get_nhierarchical() const
{
  return nhierarchical;
}

void
hierarchical_model::set(int i, double v)
{
  if (i < 0 || i >= nhierarchical) {
    throw GENERALVORONOICARTESIANEXCEPTION("Index out of range %d (%d)\n", i, nhierarchical);
  }

  hierarchical[i] = v;
}

double
hierarchical_model::get(int i) const
{
  if (i < 0 || i >= nhierarchical) {
    throw GENERALVORONOICARTESIANEXCEPTION("Index out of range %d (%d)\n", i, nhierarchical);
  }

  return hierarchical[i];
}

bool
hierarchical_model::write(std::ostream &s) const
{
  s.write((char*)&nhierarchical, sizeof(int));
  
  for (int i = 0; i < nhierarchical; i ++) {
    s.write((char*)(&hierarchical[i]), sizeof(double));
  }

  return true;
}

bool
hierarchical_model::read(std::istream &s)
{
  int n;

  s.read((char*)&n, sizeof(int));
  if (n != nhierarchical) {
    return false;
  }

  for (int i = 0; i < nhierarchical; i ++) {
    s.read((char*)(&hierarchical[i]), sizeof(double));
  }

  return true;
}

double *
hierarchical_model::data()
{
  return hierarchical;
}

singlescaling_hierarchical_model::singlescaling_hierarchical_model(double lambda) :
  hierarchical_model(1)
{
  set(0, lambda);
}

singlescaling_hierarchical_model::~singlescaling_hierarchical_model()
{
}
