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
#ifndef chainhistoryVoronoi_hpp
#define chainhistoryVoronoi_hpp

#include <vector>

#include "hierarchical_model.hpp"
#include "generalvoronoicartesianexception.hpp"

#include "cartesianvoronoimodel.hpp"

extern "C" {
  #include "slog.h"
};

class deltaVoronoi;

class model_initializationVoronoi;
  
class model_deltaVoronoi;

class model_deltaValuesVoronoi;

class hierarchical_deltaVoronoi;

class deltaVoronoi {
public:

  enum {
    DELTA_INITIALIZATION = 0,
    DELTA_DELTA = 1,
    DELTA_HIERARCHICAL = 2,
    DELTA_VALUES = 3
  };

  deltaVoronoi(int _id) :
    id(_id),
    proposed_like(0.0),
    proposed_norm(0.0),
    accepted(false)
  {
  }
  
  virtual ~deltaVoronoi()
  {
  }

  virtual int write_header(FILE *fp)
  {
    if (fwrite(&id, sizeof(int), 1, fp) != 1) {
      return -1;
    }
    
    if (fwrite(&proposed_like, sizeof(double), 1, fp) != 1) {
      return -1;
    }
    
    if (fwrite(&proposed_norm, sizeof(double), 1, fp) != 1) {
      return -1;
    }

    int a = (int)accepted;
    if (fwrite(&a, sizeof(int), 1, fp) != 1) {
      return -1;
    }
    
    return 0;
  }

  static int read_header(FILE *fp, double &like, double &norm, bool &accepted)
  {
    if (fread(&like, sizeof(double), 1, fp) != 1) {
      return -1;
    }
    
    if (fread(&norm, sizeof(double), 1, fp) != 1) {
      return -1;
    }

    int a;
    if (fread(&a, sizeof(int), 1, fp) != 1) {
      return -1;
    }
    
    accepted = (bool)a;
    
    return 0;
  }

  virtual int write(FILE *fp) = 0;

  virtual int apply(std::vector<cartesianvoronoimodel*> &model, hierarchical_model &hierarchical) = 0;

  virtual void accept()
  {
    accepted = true;
  }
    

  virtual void reject()
  {
    // Nothing to do: assumed rejected on creation
  }
    

  virtual void set_proposed_likelihood(double like, double norm)
  {
    proposed_like = like;
    proposed_norm = norm;
  }

  virtual bool isaccepted() const
  {
    return accepted;
  }
  
  virtual double get_proposed_likelihood(double &norm) const
  {
    norm = proposed_norm;
    return proposed_like;
  }

  static deltaVoronoi *read(FILE *fp);
    
private:

  int id;
  double proposed_like;
  double proposed_norm;
  bool accepted;
  
  static std::vector<deltaVoronoi* (*)(FILE *fp)> readers;
  
};

std::vector<deltaVoronoi* (*)(FILE *fp)> deltaVoronoi::readers;


class model_initializationVoronoi : public deltaVoronoi {
public:
  
  model_initializationVoronoi(std::vector<cartesianvoronoimodel*> &_models,
			      hierarchical_model &_hierarchical,
			      double _likelihood,
			      double _norm) :
    deltaVoronoi(deltaVoronoi::DELTA_INITIALIZATION)
  {
    // Initialization always accepted
    deltaVoronoi::accept();
    
    deltaVoronoi::set_proposed_likelihood(_likelihood, _norm);

    initial_models.resize(_models.size());
    int mi = 0;
    for (auto &m : _models) {

      int ncells = m->ntotalcells();
      
      initial_models[mi].parameterization = m->parameterization();
      initial_models[mi].xs.resize(ncells);
      initial_models[mi].ys.resize(ncells);
      initial_models[mi].zs.resize(ncells);
      
      for (int i = 0; i < ncells; i ++) {

	double x, y;
	m->position_at_index(i, &x, &y);

	double z = m->value_at_index(i);

	initial_models[mi].xs[i] = x;
	initial_models[mi].ys[i] = y;
	initial_models[mi].zs[i] = z;
      }

      mi ++;
    }
    
    int nh = _hierarchical.get_nhierarchical();
    for (int i = 0; i < nh; i ++) {
      hierarchical.push_back(_hierarchical.get(i));
    }
  }
    
  virtual ~model_initializationVoronoi()
  {
  }

  virtual int write(FILE *fp)
  {
    if (deltaVoronoi::write_header(fp) < 0) {
      return -1;
    }

    int nmodels = initial_models.size();

    if (fwrite(&nmodels, sizeof(int), 1, fp) != 1) {
      return -1;
    }

    for (auto &m : initial_models) {

      int ncells = m.xs.size();

      if (fwrite(&m.parameterization, sizeof(int), 1, fp) != 1) {
	return -1;
      }
      
      if (fwrite(&ncells, sizeof(int), 1, fp) != 1) {
	return -1;
      }
      
      for (int j = 0; j < ncells; j ++) {
	double xyz[3];

	xyz[0] = m.xs[j];
	xyz[1] = m.ys[j];
	xyz[2] = m.zs[j];

	if (fwrite(xyz, 3*sizeof(double), 1, fp) != 1) {
	  return -1;
	}
      }
    }
    
    int nh = (int)hierarchical.size();
    
    if (fwrite(&nh, sizeof(int), 1, fp) != 1) {
      return -1;
    }
    
    for (auto &h : hierarchical) {
      if (fwrite(&h, sizeof(double), 1, fp) != 1) {
	return -1;
      }
    }
    
    return 0;
  }

  virtual int apply(std::vector<cartesianvoronoimodel*> &_models, hierarchical_model &_hierarchical)
  {
    if (_models.size() != initial_models.size()) {
      fprintf(stderr, "model_initialization: no. models mismatch\n");
      return -1;
    }

    for (int i = 0; i < (int)_models.size(); i ++) {

      if (_models[i]->parameterization() != initial_models[i].parameterization) {
	fprintf(stderr, "model_initialization: parameterization mismatch Configured %d File %d\n",
		_models[i]->parameterization(), initial_models[i].parameterization);
	return -1;
      }
      
      _models[i]->reset();

      for (int j = 0; j < (int)initial_models[i].xs.size(); j ++) {

	if (j < 4) {
	  _models[i]->set_value_at_index(j, initial_models[i].zs[j]);
	} else {
	  _models[i]->add_cell(initial_models[i].xs[j],
			       initial_models[i].ys[j],
			       initial_models[i].zs[j]);
	}
      }

    }


    int i = 0;
    for (auto &h : hierarchical) {
      _hierarchical.set(i, h);
      i ++;
    }
    
    return 0;
  }

  static deltaVoronoi *read(FILE *fp)
  {
    double like;
    double norm;
    bool accepted;
    
    if (deltaVoronoi::read_header(fp, like, norm, accepted) < 0) {
      fprintf(stderr, "model_initialization::read: failed to read header\n");
      return nullptr;
    }
    
    model_initializationVoronoi *r = new model_initializationVoronoi();
    
    r->accept();
    r->set_proposed_likelihood(like, norm);

    int nmodels;

    if (fread(&nmodels, sizeof(int), 1, fp) != 1) {
      fprintf(stderr, "model_initialization::read: failed to read no. models\n");
      return nullptr;
    }

    r->initial_models.resize(nmodels);
    
    for (int i = 0; i < nmodels; i ++) {

      if (fread(&r->initial_models[i].parameterization, sizeof(int), 1, fp) != 1) {
	fprintf(stderr, "model_initialization::read: failed to read parameterization\n");
	return nullptr;
      }

      int ncells;

      if (fread(&ncells, sizeof(int), 1, fp) != 1) {
	fprintf(stderr, "model_initialization::read: failed to read no. cells\n");
	return nullptr;
      }

      r->initial_models[i].xs.resize(ncells);
      r->initial_models[i].ys.resize(ncells);
      r->initial_models[i].zs.resize(ncells);
    
      for (int j = 0; j < ncells; j ++) {

	double xyz[3];
	if (fread(xyz, 3*sizeof(double), 1, fp) != 1) {
	  fprintf(stderr, "model_initialization::read: failed to read xyz\n");
	  return nullptr;
	}

	r->initial_models[i].xs[j] = xyz[0];
	r->initial_models[i].ys[j] = xyz[1];
	r->initial_models[i].zs[j] = xyz[2];
      }
    }

    int nh;

    if (fread(&nh, sizeof(int), 1, fp) != 1) {
      fprintf(stderr, "model_initialization::read: failed to no. hierarchical\n");
      return nullptr;
    }
    
    for (int i = 0; i < nh; i ++) {
      double h;
      if (fread(&h, sizeof(double), 1, fp) != 1) {
	fprintf(stderr, "model_initialization::read: failed to hierarchical %d\n", i);
	return nullptr;
      }
      
      r->hierarchical.push_back(h);
    }
    
    return r;
  }
  
private:

  model_initializationVoronoi() :
    deltaVoronoi(deltaVoronoi::DELTA_INITIALIZATION)
  {
  }

  struct model_init {
    int parameterization;
    std::vector<double> xs;
    std::vector<double> ys;
    std::vector<double> zs;
  };

  std::vector<model_init> initial_models;
  std::vector<double> hierarchical;
};

//
// For global change to model, eg HMC proposals
//
class model_deltaValuesVoronoi : public deltaVoronoi {
public:

  model_deltaValuesVoronoi(std::vector<cartesianvoronoimodel *> &models) :
    deltaVoronoi(deltaVoronoi::DELTA_VALUES)
  {
    int nmodels = models.size();
    values.resize(nmodels);
    for (int mi = 0; mi < nmodels; mi ++) {

      int n = models[mi]->nvaluecells();
      values[mi].resize(n);
      for (int i = 0; i < n; i ++) {
	int ci = models[mi]->selectvaluecell(i);
	values[mi][i] = models[mi]->value_at_index(ci);
      }
    }
  }

  virtual int write(FILE *fp)
  {
    if (deltaVoronoi::write_header(fp) < 0) {
      return -1;
    }

    int nmodels = values.size();
    if (fwrite(&nmodels, sizeof(int), 1, fp) != 1) {
      return -1;
    }

    for (int mi = 0; mi < nmodels; mi ++) {

      int n = values[mi].size();
      if (fwrite(&n, sizeof(int), 1, fp) != 1) {
	return -1;
      }

      for (int i = 0; i < n; i ++) {
	if (fwrite(&(values[mi][i]), sizeof(double), 1, fp) != 1) {
	  return -1;
	}
      }
    }    

    return 0;
  }

  virtual int apply(std::vector<cartesianvoronoimodel*> &models, hierarchical_model &hierarchical)
  {
    if (models.size() != values.size()) {
      fprintf(stderr, "model_deltaValues: no. models mismatch %d %d\n",
	      (int)models.size(),
	      (int)values.size());
      return -1;
    }

    if (isaccepted()) {
      for (int mi = 0; mi < (int)models.size(); mi ++) {
	
	int n = models[mi]->nvaluecells();
	if (n != (int)values[mi].size()) {
	  fprintf(stderr, "model_deltaValues: no. cells mismatch %d %d\n",
		  n, (int)values[mi].size());
	  return -1;
	}
	
	for (int i = 0; i < n; i ++) {
	  int ci = models[mi]->selectvaluecell(i);
	  
	  models[mi]->set_value_at_index(ci, values[mi][i]);
	}
      }
    }
    
    return 0;
  }

  static deltaVoronoi *read(FILE *fp)
  {
    double like;
    double norm;
    bool accepted;
    
    if (deltaVoronoi::read_header(fp, like, norm, accepted) < 0) {
      fprintf(stderr, "model_initialization::read: failed to read header\n");
      return nullptr;
    }

    model_deltaValuesVoronoi *r = new model_deltaValuesVoronoi();

    if (accepted) {
      r->accept();
    }
    r->set_proposed_likelihood(like, norm);

    int nmodels;
    if (fread(&nmodels, sizeof(int), 1, fp) != 1) {
      fprintf(stderr, "model_deltaValuesVoronoi::read: failed to read no. models\n");
      return nullptr;
    }

    r->values.resize(nmodels);

    for (int mi = 0; mi < nmodels; mi ++) {

      int n;
      if (fread(&n, sizeof(int), 1, fp) != 1) {
	fprintf(stderr, "model_deltaValuesVoronoi::read: failed to read no. values\n");
	return nullptr;
      }

      r->values[mi].resize(n);

      for (int i = 0; i < n; i ++) {

	if (fread(&(r->values[mi][i]), sizeof(double), 1, fp) != 1) {
	  fprintf(stderr, "model_deltaValuesVoronoi::read: failed to read values\n");
	  return nullptr;
	}
      }
    }

    return r;
  }

private:

  model_deltaValuesVoronoi() :
    deltaVoronoi(deltaVoronoi::DELTA_VALUES)
  {
  }
  
  std::vector<std::vector<double>> values;
  
};

class model_deltaVoronoi : public deltaVoronoi {
public:

  //
  // Change value
  //
  static model_deltaVoronoi *mkvalue(int modelindex,
				     int cellindex,
				     double oldvalue,
				     double newvalue)
  {
    return new model_deltaVoronoi(VALUE,
				  modelindex,
				  cellindex,
				  oldvalue,
				  newvalue);
  }

  //
  // Move
  //
  static model_deltaVoronoi *mkmove(int modelindex,
				    int cellindex,
				    double oldx,
				    double oldy,
				    double newx,
				    double newy)
  {
    return new model_deltaVoronoi(MOVE,
				  modelindex,
				  cellindex,
				  oldx, oldy,
				  newx, newy);
  }

  //
  // Birth
  //
  static model_deltaVoronoi *mkbirth(int modelindex,
				     double newx,
				     double newy,
				     double newvalue)
  {
    return new model_deltaVoronoi(BIRTH,
				  modelindex,
				  newx,
				  newy,
				  newvalue);
  }

  static model_deltaVoronoi *mkdeath(int modelindex,
				     int cellindex)
  {
    return new model_deltaVoronoi(DEATH,
				  modelindex,
				  cellindex);
  }
			      
  ~model_deltaVoronoi()
  {
  }
  
  virtual int write(FILE *fp)
  {
    if (deltaVoronoi::write_header(fp) < 0) {
      return -1;
    }

    int itype = (int)type;
    if (fwrite(&itype, sizeof(int), 1, fp) != 1) {
      ERROR("Failed to write type");
      return -1;
    }

    if (fwrite(&modelindex, sizeof(int), 1, fp) != 1) {
      ERROR("Failed to write model index\n");
      return -1;
    }

    if (fwrite(&cellindex, sizeof(int), 1, fp) != 1) {
      ERROR("Failed to write cell index\n");
      return -1;
    }

    if (fwrite(&oldvalue, sizeof(double), 1, fp) != 1) {
      ERROR("Failed to write old value\n");
      return -1;
    }

    if (fwrite(&newvalue, sizeof(double), 1, fp) != 1) {
      ERROR("Failed to write old value\n");
      return -1;
    }
    
    if (fwrite(&oldx, sizeof(double), 1, fp) != 1) {
      ERROR("Failed to write old x\n");
      return -1;
    }

    if (fwrite(&oldy, sizeof(double), 1, fp) != 1) {
      ERROR("Failed to write old y\n");
      return -1;
    }

    if (fwrite(&newx, sizeof(double), 1, fp) != 1) {
      ERROR("Failed to write new x\n");
      return -1;
    }

    if (fwrite(&newy, sizeof(double), 1, fp) != 1) {
      ERROR("Failed to write new y\n");
      return -1;
    }

    return 0;
  }
  
    
  virtual int apply(std::vector<cartesianvoronoimodel*> &models, hierarchical_model &hierarchical)
  {
    if (deltaVoronoi::isaccepted()) {

      switch (type) {
      case BIRTH:
	{
	  models[modelindex]->add_cell(newx, newy, newvalue);
	}
	break;
	
      case DEATH:
	{
	  models[modelindex]->delete_cell(cellindex);
	}
	break;

      case MOVE:
	{
	  models[modelindex]->set_position_at_index(cellindex, newx, newy);
	}
	break;

      case VALUE:
	{
	  models[modelindex]->set_value_at_index(cellindex, newvalue);
	}
	break;
	  
      default:
	throw GENERALVORONOICARTESIANEXCEPTION("Unhandled delta type\n");
      }
    }
    return 0;
  }
  
  static deltaVoronoi *read(FILE *fp)
  {
    double like;
    double norm;
    bool accepted;
    if (deltaVoronoi::read_header(fp, like, norm, accepted) < 0) {
      fprintf(stderr, "model_deltaVoronoi::read: failed to read header\n");
      return nullptr;
    }
    
    model_deltaVoronoi *r = new model_deltaVoronoi();
    
    if (accepted) {
      r->accept();
    }
    r->set_proposed_likelihood(like, norm);

    int itype;
    if (fread(&itype, sizeof(int), 1, fp) != 1) {
      ERROR("Failed to read type");
      return nullptr;
    }
    if (itype < BIRTH || itype > MOVE) {
      ERROR("Invalid type: %d\n", itype);
    }
    r->type = (delta_t)itype;

    if (fread(&r->modelindex, sizeof(int), 1, fp) != 1) {
      ERROR("failed to read model index\n");
      return nullptr;
    }

    if (fread(&r->cellindex, sizeof(int), 1, fp) != 1) {
      ERROR("Failed to read cell index\n");
      return nullptr;
    }

    if (fread(&r->oldvalue, sizeof(double), 1, fp) != 1) {
      ERROR("Failed to read old value\n");
      return nullptr;
    }

    if (fread(&r->newvalue, sizeof(double), 1, fp) != 1) {
      ERROR("Failed to read old value\n");
      return nullptr;
    }

    
    if (fread(&r->oldx, sizeof(double), 1, fp) != 1) {
      ERROR("Failed to read old x\n");
      return nullptr;
    }

    if (fread(&r->oldy, sizeof(double), 1, fp) != 1) {
      ERROR("Failed to read old y\n");
      return nullptr;
    }

    if (fread(&r->newx, sizeof(double), 1, fp) != 1) {
      ERROR("Failed to read new x\n");
      return nullptr;
    }

    if (fread(&r->newy, sizeof(double), 1, fp) != 1) {
      ERROR("Failed to read new y\n");
      return nullptr;
    }
    
    return r;
  }

private:

  typedef enum {
    UNKNOWN = -1,
    BIRTH = 0,
    DEATH = 1,
    VALUE = 2,
    MOVE = 3
  } delta_t;

  model_deltaVoronoi() :
    deltaVoronoi(deltaVoronoi::DELTA_DELTA)
  {
  }
  
  //
  // Birth Cell
  //
  model_deltaVoronoi(delta_t _type,
		     int _modelindex,
		     double _newx,
		     double _newy,
		     double _newvalue) :
    deltaVoronoi(deltaVoronoi::DELTA_DELTA),
    type(_type),
    modelindex(_modelindex),
    cellindex(-1),
    oldvalue(0.0),
    newvalue(_newvalue),
    newx(_newx),
    newy(_newy)
  {
  }

  //
  // Death Cell
  model_deltaVoronoi(delta_t _type,
		     int _modelindex,
		     int _cellindex) :
    deltaVoronoi(deltaVoronoi::DELTA_DELTA),
    type(_type),
    modelindex(_modelindex),
    cellindex(_cellindex),
    oldvalue(0.0),
    newvalue(0.0)
  {
  }
    

  //
  // Value
  //
  model_deltaVoronoi(delta_t _type,
		     int _modelindex,
		     int _cellindex,
		     double _oldvalue,
		     double _newvalue) :
    deltaVoronoi(deltaVoronoi::DELTA_DELTA),
    type(_type),
    modelindex(_modelindex),
    cellindex(_cellindex),
    oldvalue(_oldvalue),
    newvalue(_newvalue)
  {
  }

  //
  // Move
  //
  model_deltaVoronoi(delta_t _type,
		     int _modelindex,
		     int _cellindex,
		     double _oldx,
		     double _oldy,
		     double _newx,
		     double _newy) :
    deltaVoronoi(deltaVoronoi::DELTA_DELTA),
    type(_type),
    modelindex(_modelindex),
    cellindex(_cellindex),
    oldvalue(0.0),
    newvalue(0.0),
    oldx(_oldx),
    oldy(_oldy),
    newx(_newx),
    newy(_newy)
  {
  }

  delta_t type;

  int modelindex;
  int cellindex;
  double oldvalue;
  double newvalue;
  double oldx;
  double oldy;
  double newx;
  double newy;
  
};

class hierarchical_deltaVoronoi : public deltaVoronoi {
public:

  hierarchical_deltaVoronoi(int nhierarchical, int *indices, double *old_value, double *new_value) :
    deltaVoronoi(deltaVoronoi::DELTA_HIERARCHICAL)
  {
    for (int i = 0; i < nhierarchical; i ++) {
      
      hierarchical_term_delta d;
      
      d.index = indices[i];
      d.old_value = old_value[i];
      d.new_value = new_value[i];
      
      hierarchical.push_back(d);
    }
  }
  ~hierarchical_deltaVoronoi()
  {
  }
  
  virtual int write(FILE *fp)
  {
    if (deltaVoronoi::write_header(fp) < 0) {
      return -1;
    }
    
    int nh = (int)hierarchical.size();
    if (fwrite(&nh, sizeof(int), 1, fp) != 1) {
      return -1;
    }
    
    for (auto &dh : hierarchical) {
      if (fwrite(&(dh.index), sizeof(int), 1, fp) != 1) {
	return -1;
      }
      if (fwrite(&(dh.old_value), sizeof(double), 1, fp) != 1) {
	return -1;
      }
      if (fwrite(&(dh.new_value), sizeof(double), 1, fp) != 1) {
	return -1;
      }
    }
    
    return 0;
  }

  virtual int apply(std::vector<cartesianvoronoimodel*> &models, hierarchical_model &_hierarchical)
  {
    if (deltaVoronoi::isaccepted()) {
      for (auto &dh : hierarchical) {
	
	if (dh.index < 0 || dh.index >= _hierarchical.get_nhierarchical()) {
	  throw GENERALVORONOICARTESIANEXCEPTION("Hierarchical index out of range: %d (%d)\n",
						 dh.index, _hierarchical.get_nhierarchical());
	}
	
	_hierarchical.set(dh.index, dh.new_value);
      }
    }
  
    return 0;
  }
  
  static deltaVoronoi *read(FILE *fp)
  {
    double like;
    double norm;
    bool accepted;
    if (deltaVoronoi::read_header(fp, like, norm, accepted) < 0) {
      return nullptr;
    }
    
    hierarchical_deltaVoronoi *r = new hierarchical_deltaVoronoi(0, nullptr, nullptr, nullptr);
    
    if (accepted) {
      r->accept();
    }
    r->set_proposed_likelihood(like, norm);
    
    int nh;
    if (fread(&nh, sizeof(int), 1, fp) != 1) {
      return nullptr;
    }
    
    for (int i = 0; i < nh; i ++) {
      
      hierarchical_term_delta dh;
      
      if (fread(&(dh.index), sizeof(int), 1, fp) != 1) {
	return nullptr;
      }

      if (fread(&(dh.old_value), sizeof(double), 1, fp) != 1) {
	return nullptr;
      }
      if (fread(&(dh.new_value), sizeof(double), 1, fp) != 1) {
	return nullptr;
      }
      
      r->hierarchical.push_back(dh);
    }
    
    return r;
  }
  
private:
  
  struct hierarchical_term_delta {
    int index;
    double old_value;
    double new_value;
  };
  std::vector<hierarchical_term_delta> hierarchical;
};

class chainhistorywriterVoronoi {
public:

  chainhistorywriterVoronoi(const char *_filename,
			    std::vector<cartesianvoronoimodel*> &_initial_models,
			    hierarchical_model &_hierarchical,
			    double _likelihood,
			    double _norm) :
    fp(fopen(_filename, "w"))
  {
    if (fp == NULL) {
      throw GENERALVORONOICARTESIANEXCEPTION("Failed to create chain history file: %s\n", _filename);
    }
    
    steps.push_back(new model_initializationVoronoi(_initial_models, _hierarchical, _likelihood, _norm));
    
    flush();
  }

  ~chainhistorywriterVoronoi()
  {
    flush();
    fclose(fp);
  }


  void add(deltaVoronoi *d)
  {
    steps.push_back(d);
  }

  void flush()
  {
    for (auto &s : steps) {
      s->write(fp);
      delete s;
    }
    
    steps.clear();
  }

private:

  FILE *fp;
  std::vector<deltaVoronoi*> steps;
};

class chainhistoryreaderVoronoi {
public:
  
  chainhistoryreaderVoronoi(const char *filename) :
    fp(fopen(filename, "r")),
    current_likelihood(-1.0),
    current_norm(-1.0)
  {
    if (fp == NULL) {
      throw GENERALVORONOICARTESIANEXCEPTION("Failed to open file for reading: %s\n", filename);
    }
  }
  ~chainhistoryreaderVoronoi()
  {
    fclose(fp);
  }

  int step(std::vector<cartesianvoronoimodel*> &models, hierarchical_model &hierarchical, double &likelihood, double &norm)
  {
    deltaVoronoi *d = deltaVoronoi::read(fp);

    if (d == nullptr) {
      if (feof(fp)) {
	return 0;
      } else {
	fprintf(stderr, "chainhistoryreaderVoronoi::step: failed to read next step\n");
	return -1;
      }
    }
    
    if (d->apply(models, hierarchical) < 0) {
      fprintf(stderr, "chainhistoryreaderVoronoi::step: failed to apply step to model/hierarchical/likelihood\n");
      return -1;
    }
    
    if (d->isaccepted()) {
      current_likelihood = d->get_proposed_likelihood(current_norm);
    }
    
    likelihood = current_likelihood;
    norm = current_norm;
    
    delete d;
    
    return 1;
  }
  
private:

  FILE *fp;
  double current_likelihood;
  double current_norm;
  
};

deltaVoronoi *
deltaVoronoi::read(FILE *fp)
{
  int id;
  
  if (fread(&id, sizeof(int), 1, fp) != 1) {
    return nullptr;
  }
  
  switch (id) {
    case DELTA_INITIALIZATION:
      return model_initializationVoronoi::read(fp);
      
  case DELTA_DELTA:
    return model_deltaVoronoi::read(fp);
    
  case DELTA_HIERARCHICAL:
    return hierarchical_deltaVoronoi::read(fp);

  case DELTA_VALUES:
    return model_deltaValuesVoronoi::read(fp);
    
  default:
    fprintf(stderr, "deltaVoronoi::read: invalid unique index: %d (%d)\n", id, (int)readers.size());
    return nullptr;
  }
}


#endif // chainhistoryVoronoi_hpp
  
