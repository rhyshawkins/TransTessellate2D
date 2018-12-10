
TARGETS = \
	pointtemplate_50_50_50.txt \
	syntheticobs_50_50_50_tas.txt \
	results_50_50_50_vv/mean_sea.txt \
	results_50_50_50_vv/mean_land.txt \
	results_50_50_50_vv/mean_tide.txt \
	results_50_50_50_dv/mean_sea.txt \
	results_50_50_50_dv/mean_land.txt \
	results_50_50_50_dv/mean_tide.txt

PROCESSORS = 4
CHAINS = 1

all : $(TARGETS)

MINLON=144.0
MAXLON=149.0
MINLAT=-44.0
MAXLAT=-40.0

#
# Point templates
#

pointtemplate_50_50_50.txt : sdf_tas_fine.txt
	python2 ../scripts/generatesdftemplatepoints.py -T 50 -L 50 -S 50 \
	--sdf sdf_tas_fine.txt \
	--xmin $(MINLON) --xmax $(MAXLON) \
	--ymin $(MINLAT) --ymax $(MAXLAT) \
	--land-boundary 0.5 \
	--tide-boundary 0.5 \
	-o $@

#
# Synthetic observations Sea model is first, the Land model
#

syntheticobs_50_50_50_tas.txt : pointtemplate_50_50_50.txt
	../mksynthetic -i $< -o $@ \
	-m TasSea \
	-m TasLand \
	--xmin $(MINLON) --xmax $(MAXLON) \
	--ymin $(MINLAT) --ymax $(MAXLAT) \
	-n 1.0 \
	-I syntheticobs_50_50_50_tas

#
# Priors
#
prior_land.txt :
	echo "prior Uniform" > $@
	echo "-10.0 10.0" >> $@
	echo "proposal Gaussian" >> $@
	echo "0.2" >> $@

prior_sea.txt :
	echo "prior Uniform" > $@
	echo "-5.0 5.0" >> $@
	echo "proposal Gaussian" >> $@
	echo "0.2" >> $@

positionprior_land.txt :
	echo "positionprior UniformPosition" > $@
	echo "positionproposal GaussianPosition" >> $@
	echo "0.2">> $@

positionprior_sea.txt :
	echo "positionprior UniformPosition" > $@
	echo "positionproposal GaussianPosition" >> $@
	echo "0.2">> $@

hierarchical_land.txt :
	echo "prior Uniform" > $@
	echo "0.001 2.0" >> $@
	echo "proposal Gaussian" >> $@
	echo "0.05">> $@

hierarchical_sea.txt :
	echo "prior Uniform" > $@
	echo "0.001 2.0" >> $@
	echo "proposal Gaussian" >> $@
	echo "0.05">> $@

#
# Synthetic inversion using Voronoi/Voronoi parameterization
#
results_50_50_50_vv/ch.dat : syntheticobs_50_50_50_tas.txt \
	prior_land.txt \
	prior_sea.txt \
	positionprior_land.txt \
	positionprior_sea.txt \
	hierarchical_land.txt \
	hierarchical_sea.txt
	mkdir -p results_50_50_50_vv
	mpirun -np $(PROCESSORS) ../tideshmcpt -i syntheticobs_50_50_50_tas.txt \
	-o results_50_50_50_vv/ \
	-P prior_sea.txt -P prior_land.txt \
	-M positionprior_sea.txt -M positionprior_land.txt \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-A 0 -A 0 \
	-C 10 \
	-T 100 \
	-t 1000000 \
	-v 10000 \
	-c $(CHAINS) \
	-K 4 -m 5.0 -e 10

results_50_50_50_vv/mean_sea.txt : results_50_50_50_vv/ch.dat
	mpirun -np $(CHAINS) ../post_mean_mpi -i $< -o $@ \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-T 100 \
	-A 0 -A 0 \
	-W 50 -H 50 \
	-I 0 \
	-s 500000 -t 10

results_50_50_50_vv/mean_land.txt : results_50_50_50_vv/ch.dat
	mpirun -np $(CHAINS) ../post_mean_mpi -i $< -o $@ \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-T 100 \
	-A 0 -A 0 \
	-W 50 -H 50 \
	-I 1 \
	-s 500000 -t 10

results_50_50_50_vv/mean_tide.txt : results_50_50_50_vv/ch.dat
	mpirun -np $(CHAINS) ../post_mean_tide_mpi -i $< -o $@ \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-T 100 \
	-A 0 -A 0 \
	-W 50 -H 50 \
	-s 500000 -t 10

#
# Synthetic inversion using Delaunay/Voronoi parameterization
#
results_50_50_50_dv/ch.dat : syntheticobs_50_50_50_tas.txt \
	prior_land.txt \
	prior_sea.txt \
	positionprior_land.txt \
	positionprior_sea.txt \
	hierarchical_land.txt \
	hierarchical_sea.txt
	mkdir -p results_50_50_50_dv
	mpirun -np $(PROCESSORS) ../tideshmcpt -i syntheticobs_50_50_50_tas.txt \
	-o results_50_50_50_dv/ \
	-P prior_sea.txt -P prior_land.txt \
	-M positionprior_sea.txt -M positionprior_land.txt \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-A 1 -A 0 \
	-C 10 \
	-T 100 \
	-t 1000000 \
	-v 10000 \
	-c $(CHAINS) \
	-K 4 -m 5.0 -e 10

results_50_50_50_dv/mean_sea.txt : results_50_50_50_dv/ch.dat
	mpirun -np $(CHAINS) ../post_mean_mpi -i $< -o $@ \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-T 100 \
	-A 1 -A 0 \
	-W 50 -H 50 \
	-I 0 \
	-s 500000 -t 10

results_50_50_50_dv/mean_land.txt : results_50_50_50_dv/ch.dat
	mpirun -np $(CHAINS) ../post_mean_mpi -i $< -o $@ \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-T 100 \
	-A 1 -A 0 \
	-W 50 -H 50 \
	-I 1 \
	-s 500000 -t 10

results_50_50_50_dv/mean_tide.txt : results_50_50_50_dv/ch.dat
	mpirun -np $(CHAINS) ../post_mean_tide_mpi -i $< -o $@ \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-T 100 \
	-A 1 -A 0 \
	-W 50 -H 50 \
	-s 500000 -t 10

