
CXX = mpicxx

#
# For intel compilers
#
#CXX =mpiicpc

CXXFLAGS=-c -g -Wall --std=c++11 -I../base

CXXFLAGS+=-O3

BASEOBJS = \
	../base/generalvoronoicartesianexception.o \
	../base/hierarchical_model.o \
	../base/pathutil.o \
	../base/prior.o \
	../base/rng.o \
	../base/positionprior.o \
	../base/global.o \
	../base/simulated_annealing_scales.o \
	../base/synthetic.o \
	../base/ptexchange.o \
	../base/generalvoronoicartesianutil.o

LIBS = -L../lib -ldelaunay \
	$(shell gsl-config --libs)

TARGETOBJS = genericregression.o

TARGETS = mksynthetic \
	regression \
	regressionpt \
	hmc \
	hmcpt \
	post_mean \
	post_likelihood \
	post_likelihood_mpi \
	post_mean_mpi \
	post_resolution \
	post_resolution_mpi \
	post_text


all : $(TARGETS)

genericregression.o : genericregression.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

../base/%.o : ../base/%.cpp
	make -C ../base  $*.o

mksynthetic : $(TARGETOBJS) $(BASEOBJS) ../base/mksynthetic.o
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/mksynthetic.o $(LIBS)

regression : $(TARGETOBJS) $(BASEOBJS) ../base/generalvoronoicartesian.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/generalvoronoicartesian.o $(LIBS)

hmc : $(TARGETOBJS) $(BASEOBJS) ../base/generalvoronoicartesianhmc.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/generalvoronoicartesianhmc.o $(LIBS)

hmcpt : $(TARGETOBJS) $(BASEOBJS) ../base/generalvoronoicartesianhmcpt.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/generalvoronoicartesianhmcpt.o $(LIBS)

regressionpt : $(TARGETOBJS) $(BASEOBJS) ../base/generalvoronoicartesianpt.o 
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/generalvoronoicartesianpt.o $(LIBS)

post_mean : $(TARGETOBJS) $(BASEOBJS) ../base/post_mean.o
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/post_mean.o $(LIBS)

post_mean_mpi : $(TARGETOBJS) $(BASEOBJS) ../base/post_mean_mpi.o
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/post_mean_mpi.o $(LIBS)

post_likelihood : $(TARGETOBJS) $(BASEOBJS) ../base/post_likelihood.o
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/post_likelihood.o $(LIBS)

post_likelihood_mpi : $(TARGETOBJS) $(BASEOBJS) ../base/post_likelihood_mpi.o
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/post_likelihood_mpi.o $(LIBS)

post_resolution : $(TARGETOBJS) $(BASEOBJS) ../base/post_resolution.o
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/post_resolution.o $(LIBS)

post_resolution_mpi : $(TARGETOBJS) $(BASEOBJS) ../base/post_resolution_mpi.o
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/post_resolution_mpi.o $(LIBS)

post_text : $(TARGETOBJS) $(BASEOBJS) ../base/post_text.o
	$(CXX) -o $@ $(TARGETOBJS) $(BASEOBJS) ../base/post_text.o $(LIBS)

clean :
	rm -f $(TARGETS) *.o
