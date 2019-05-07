
TARGETS = datatemplate.txt \
	syntheticobs_franke.txt \
	results_voronoi/mean.txt \
	results_delaunay/mean.txt \
	results_ct/mean.txt

NOISE = 0.05

PROCESSORS = 4
CHAINS = 1

all : $(TARGETS)

datatemplate.txt :
	python2 ../../scripts/generatetemplatepoints.py -N 100 -o $@ \
	--seed 983 \
	--xmin 0.0 --xmax 1.0 \
	--ymin 0.0 --ymax 1.0 

syntheticobs_franke.txt : datatemplate.txt 
	../mksynthetic -m Franke \
	-i datatemplate.txt \
	-x 0.0 -X 1.0 \
	-y 0.0 -Y 1.0 \
	-o $@ -n $(NOISE) \
	-O $@.true \
	-I syntheticobs_franke.img \
	-W 1024 -H 1024

priorproposal.txt :
	echo "prior Uniform" > $@
	echo "-0.5 2.0" >> $@
	echo "proposal Gaussian" >> $@
	echo "0.05" >> $@

positionpriorproposal.txt :
	echo "positionprior UniformPosition" > $@
	echo "positionproposal GaussianPosition" >> $@
	echo "0.1" >> $@

hierarchicalpriorproposal.txt :
	echo "prior Uniform" > $@
	echo "0.5 2.0" >> $@
	echo "proposal Gaussian" >> $@
	echo "0.07" >> $@

results_voronoi/ch.dat-000 : syntheticobs_franke.txt \
	priorproposal.txt \
	positionpriorproposal.txt \
	hierarchicalpriorproposal.txt 
	mkdir -p results_voronoi
	mpirun --oversubscribe -np $(PROCESSORS) ../hmcpt -i syntheticobs_franke.txt \
	-P priorproposal.txt \
	-M positionpriorproposal.txt \
	-H hierarchicalpriorproposal.txt \
	-o results_voronoi/ \
	-v 500 \
	-x 0.0 -X 1.0 \
	-y 0.0 -Y 1.0 \
	-t 1000000 \
	-v 10000 \
	-C 50 \
	-c $(CHAINS) \
	-K 4 -m 5.0 -e 10

results_voronoi/mean.txt : results_voronoi/ch.dat-000
	mpirun -np $(CHAINS) ../post_mean_mpi \
	-x 0.0 -X 1.0 -y 0.0 -Y 1.0 \
	-i results_voronoi/ch.dat -o $@ \
	-t 50 \
	-s 500000 -W 50 -H 50 \
	-D results_voronoi/stddev.txt

results_voronoi/like.txt : results_voronoi/ch.dat-000
	../post_likelihood \
	-x 0.0 -X 1.0 -y 0.0 -Y 1.0 \
	-i results_voronoi/ch.dat-000 -o $@ \
	-H results_voronoi/hierarchical.txt \
	-K results_voronoi/khistory.txt \
	-t 50 \
	-s 500000

results_delaunay/ch.dat-000 : syntheticobs_franke.txt \
	priorproposal.txt \
	positionpriorproposal.txt \
	hierarchicalpriorproposal.txt 
	mkdir -p results_delaunay
	mpirun -np $(PROCESSORS) ../hmcpt -i syntheticobs_franke.txt \
	-P priorproposal.txt \
	-M positionpriorproposal.txt \
	-H hierarchicalpriorproposal.txt \
	-o results_delaunay/ \
	-v 500 \
	-x 0.0 -X 1.0 \
	-y 0.0 -Y 1.0 \
	-t 1000000 \
	-v 10000 \
	-C 50 \
	-A 1 \
	-c $(CHAINS) \
	-K 4 -m 5.0 -e 10

results_delaunay/mean.txt : results_delaunay/ch.dat-000
	mpirun -np $(CHAINS) ../post_mean_mpi \
	-x 0.0 -X 1.0 -y 0.0 -Y 1.0 \
	-i results_delaunay/ch.dat -o $@ \
	-t 50 -s 500000 -W 50 -H 50 \
	-D results_delaunay/stddev.txt \
	-A 1

results_ct/ch.dat-000 : syntheticobs_franke.txt \
	priorproposal.txt \
	positionpriorproposal.txt \
	hierarchicalpriorproposal.txt 
	mkdir -p results_ct
	mpirun -np $(PROCESSORS) ../hmcpt -i syntheticobs_franke.txt \
	-P priorproposal.txt \
	-M positionpriorproposal.txt \
	-H hierarchicalpriorproposal.txt \
	-o results_ct/ \
	-v 500 \
	-x 0.0 -X 1.0 \
	-y 0.0 -Y 1.0 \
	-t 1000000 \
	-v 10000 \
	-C 50 \
	-A 2 \
	-c $(CHAINS) \
	-K 4 -m 5.0 -e 10

results_ct/mean.txt : results_ct/ch.dat-000
	mpirun -np $(CHAINS) ../post_mean_mpi \
	-x 0.0 -X 1.0 -y 0.0 -Y 1.0 \
	-i results_ct/ch.dat -o $@ \
	-t 50 -s 500000 -W 50 -H 50 \
	-D results_ct/stddev.txt \
	-A 2


