
INCLUDES = \
	-I../lib \
	$(shell mpicc -showme:compile)

EXTRA_LIBS =


CXX ?= g++
CXXFLAGS = -c -g -Wall --std=c++11 $(INCLUDES)

CXXFLAGS += -O3

LIBS = $(EXTRA_LIBS) \
	-lm \
	$(shell gsl-config --libs) \
	$(shell mpicxx -showme:link)

OBJS = rng.o \
	prior.o \
	positionprior.o \
	hierarchical_model.o \
	pathutil.o \
	generalvoronoicartesianexception.o \
	global.o \
	simulated_annealing_scales.o \
	post_likelihood.o \
	post_mean.o \
	post_mean_derived.o \
	post_mean_mpi.o \
	post_mean_derived_mpi.o \
	post_text.o \
	generalvoronoicartesian.o \
	generalvoronoicartesianpt.o \
	generalvoronoicartesianhmc.o \
	generalvoronoicartesianhmcpt.o \
	generalvoronoicartesianutil.o \
	synthetic.o \
	mksynthetic.o \
	ptexchange.o

all : $(OBJS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -o $*.o $*.cpp

generalvoronoicartesianexception.o : generalvoronoicartesianexception.hpp

generalvoronoicartesianhmc.o : global.hpp \
	generalvoronoicartesianexception.hpp \
	cartesianvoronoimodel.hpp \
	generalvoronoicartesianobservations.hpp \
	prior.hpp \
	positionprior.hpp \
	hierarchical_model.hpp \
	genericinterface.hpp \
	simulated_annealing_scales.hpp \
	perturbationcollection.hpp \
	birthgeneric.hpp \
	deathgeneric.hpp \
	move.hpp \
	hierarchical.hpp \
	pathutil.hpp

generalvoronoicartesian.o : global.hpp \
	generalvoronoicartesianexception.hpp \
	cartesianvoronoimodel.hpp \
	generalvoronoicartesianobservations.hpp \
	prior.hpp \
	positionprior.hpp \
	hierarchical_model.hpp \
	genericinterface.hpp \
	simulated_annealing_scales.hpp \
	perturbationcollection.hpp \
	birthgeneric.hpp \
	deathgeneric.hpp \
	move.hpp \
	hierarchical.hpp \
	pathutil.hpp

generalvoronoicartesianpt.o : global.hpp \
	generalvoronoicartesianexception.hpp \
	cartesianvoronoimodel.hpp \
	generalvoronoicartesianobservations.hpp \
	prior.hpp \
	positionprior.hpp \
	hierarchical_model.hpp \
	genericinterface.hpp \
	simulated_annealing_scales.hpp \
	perturbationcollection.hpp \
	birthgeneric.hpp \
	deathgeneric.hpp \
	move.hpp \
	hierarchical.hpp \
	pathutil.hpp

global.hpp : generalvoronoicartesianexception.hpp \
	cartesianvoronoimodel.hpp \
	generalvoronoicartesianobservations.hpp \
	prior.hpp \
	positionprior.hpp \
	hierarchical_model.hpp \
	genericinterface.hpp \
	simulated_annealing_scales.hpp

rng.o : rng.hpp



clean :
	rm -f *.o




