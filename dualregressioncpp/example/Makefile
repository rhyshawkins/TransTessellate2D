
TARGETS = datatemplate.txt \
	syntheticobs_constant.txt \
	syntheticobs_gaussian.txt \
	results_constant/ch.dat 

all : $(TARGETS)

datatemplate.txt :
	python2 ../../scripts/generatedualtemplatepoints.py -A 50 -B 50 -o $@ \
	--xmin -10.0 --xmax 10.0 \
	--ymin -10.0 --ymax 10.0 

syntheticobs_constant.txt : datatemplate.txt 
	../mksynthetic -m Constant -m Constant \
	-i datatemplate.txt \
	-x -10.0 -X 10.0 \
	-y -10.0 -Y 10.0 \
	-o $@ -n 0.05

syntheticobs_gaussian.txt : datatemplate.txt 
	../mksynthetic -m Gaussian -m Gaussian \
	-i datatemplate.txt \
	-x -10.0 -X 10.0 \
	-y -10.0 -Y 10.0 \
	-o $@ \
	-O $@.true \
	-n 0.05

priorproposal.txt :
	echo "prior Uniform" > $@
	echo "0.0 1.0" >> $@
	echo "proposal Gaussian" >> $@
	echo "0.05" >> $@

positionpriorproposal.txt :
	echo "positionprior UniformPosition" > $@
	echo "positionproposal GaussianPosition" >> $@
	echo "0.2" >> $@

hierarchicalpriorproposal.txt :
	echo "prior Uniform" > $@
	echo "0.5 1.5" >> $@
	echo "proposal Gaussian" >> $@
	echo "0.05" >> $@


#
# Voronoi Constant True model
#
results_constant/ch.dat : syntheticobs_constant.txt \
	priorproposal.txt \
	positionpriorproposal.txt
	mkdir -p results_constant
	../regressionvoronoicartesian -i syntheticobs_constant.txt \
	-P priorproposal.txt \
	-M positionpriorproposal.txt \
	-o results_constant/ \
	-x -10.0 -X 10.0 \
	-y -10.0 -Y 10.0 \
	-t 50000 \
	-v 500

results_constant/mean.txt : results_constant/ch.dat
	../post_mean -x -10.0 -X 10.0 -y -10.0 -Y 10.0 -i $< -o $@ -t 5 -s 0

#
# Delaunay Constant True model
#
results_constant_delaunay/ch.dat : syntheticobs_constant.txt \
	priorproposal.txt \
	positionpriorproposal.txt
	mkdir -p results_constant_delaunay
	../regressionvoronoicartesian -i syntheticobs_constant.txt \
	-A 1 \
	-P priorproposal.txt -P priorproposal.txt \
	-M positionpriorproposal.txt -M positionpriorprosal.txt \
	-o results_constant_delaunay/ \
	-x -10.0 -X 10.0 \
	-y -10.0 -Y 10.0 \
	-t 50000 \
	-v 500

results_constant_delaunay/mean.txt : results_constant_delaunay/ch.dat
	../post_mean -A 1 -x -10.0 -X 10.0 -y -10.0 -Y 10.0 -i $< -o $@ -t 5 -s 0

#
# Serial Gaussian Hill True model
#
results_gaussian/ch.dat : syntheticobs_gaussian.txt \
	priorproposal.txt \
	positionpriorproposal.txt \
	hierarchicalpriorproposal.txt
	mkdir -p results_gaussian
	../regressionvoronoicartesian -i syntheticobs_gaussian.txt \
	-P priorproposal.txt -P priorproposal.txt \
	-M positionpriorproposal.txt -M positionpriorproposal.txt \
	-H hierarchicalpriorproposal.txt -H hierarchicalpriorproposal.txt \
	-o results_gaussian/ \
	-v 500 \
	-x -10.0 -X 10.0 \
	-y -10.0 -Y 10.0 \
	-t 1000000 \
	-v 10000 \
	-C 50

results_gaussian/meanA.txt : results_gaussian/ch.dat
	../post_mean -I 0 -x -10.0 -X 10.0 -y -10.0 -Y 10.0 -i $< -o $@ -t 50 -s 500000 -W 50 -H 50

results_gaussian/meanB.txt : results_gaussian/ch.dat
	../post_mean -I 1 -x -10.0 -X 10.0 -y -10.0 -Y 10.0 -i $< -o $@ -t 50 -s 500000 -W 50 -H 50


results_gaussian_opt/ch.dat : syntheticobs_gaussian.txt \
	priorproposal.txt \
	positionpriorproposal.txt \
	hierarchicalpriorproposal.txt
	mkdir -p results_gaussian_opt
	../regressionvoronoicartesian -i syntheticobs_gaussian.txt \
	-P priorproposal.txt -P priorproposal.txt \
	-M positionpriorproposal.txt -M positionpriorproposal.txt \
	-H hierarchicalpriorproposal.txt -H hierarchicalpriorproposal.txt \
	-o results_gaussian_opt/ \
	-v 500 \
	-x -10.0 -X 10.0 \
	-y -10.0 -Y 10.0 \
	-t 1 \
	-v 1000 \
	-C 50 \
	-z 10000 \
	-Z 5.0

results_gaussian_opt/meanA.txt : results_gaussian_opt/ch.dat
	../post_mean -I 0 -x -10.0 -X 10.0 -y -10.0 -Y 10.0 -i $< -o $@ -t 0 -s 0 -W 50 -H 50

results_gaussian_opt/meanB.txt : results_gaussian_opt/ch.dat
	../post_mean -I 1 -x -10.0 -X 10.0 -y -10.0 -Y 10.0 -i $< -o $@ -t 0 -s 0 -W 50 -H 50

results_gaussian_opt_pt/ch.dat-000 : syntheticobs_gaussian.txt \
	priorproposal.txt \
	positionpriorproposal.txt \
	hierarchicalpriorproposal.txt
	mkdir -p results_gaussian_opt_pt
	mpirun -np 2 ../regressionvoronoicartesianpt -i syntheticobs_gaussian.txt \
	-P priorproposal.txt -P priorproposal.txt \
	-M positionpriorproposal.txt -M positionpriorproposal.txt \
	-H hierarchicalpriorproposal.txt -H hierarchicalpriorproposal.txt \
	-o results_gaussian_opt_pt/ \
	-v 500 \
	-x -10.0 -X 10.0 \
	-y -10.0 -Y 10.0 \
	-t 1 \
	-v 1000 \
	-C 50 \
	-z 10000 \
	-Z 5.0

results_gaussian_opt_pt/meanA.txt : results_gaussian_opt_pt/ch.dat-000
	../post_mean -I 0 -x -10.0 -X 10.0 -y -10.0 -Y 10.0 -i $< -o $@ -t 0 -s 0 -W 50 -H 50

results_gaussian_opt_pt/meanB.txt : results_gaussian_opt_pt/ch.dat-000
	../post_mean -I 1 -x -10.0 -X 10.0 -y -10.0 -Y 10.0 -i $< -o $@ -t 0 -s 0 -W 50 -H 50

#
# Serial Hierarchical Prior Constant
#
results_constant_hp/ch.dat : syntheticobs_constant.txt
	mkdir -p results_constant_hp
	../regressiontdtwave2d -i syntheticobs_constant.txt \
	-M priorproposal.txt \
	-o results_constant_hp/ \
	-x 5 -y 5 -t 100000 -w 4 \
	-p 0.1

results_constant_hp/mean.txt : results_constant_hp/ch.dat
	../postprocess_mean -w 4 -x 5 -y 5 -i $< -o $@ -t 5 -s 50000

#
# Serial Hierarchical Prior Gaussian
#
results_gaussian_hp/ch.dat : syntheticobs_gaussian.txt
	mkdir -p results_gaussian_hp
	../regressiontdtwave2d -i syntheticobs_gaussian.txt \
	-M priorproposal.txt \
	-o results_gaussian_hp/ \
	-x 5 -y 5 -t 1000000 -w 4 \
	-p 0.1

results_gaussian_hp/mean.txt : results_gaussian_hp/ch.dat
	../postprocess_mean -w 4 -x 5 -y 5 -i $< -o $@ -t 5 -s 500000

#
# Serial Hierarchical Noise Constant
#
results_constant_hn/ch.dat : syntheticobs_constant.txt
	mkdir -p results_constant_hn
	../regressiontdtwave2d -i syntheticobs_constant.txt \
	-M priorproposal.txt \
	-o results_constant_hn/ \
	-x 5 -y 5 -t 100000 -w 4 \
	-L 0.01

results_constant_hn/mean.txt : results_constant_hn/ch.dat
	../postprocess_mean -w 4 -x 5 -y 5 -i $< -o $@ -t 5 -s 50000

#
# Serial Hierarchical Noise Gaussian
#
results_gaussian_hn/ch.dat : syntheticobs_gaussian.txt
	mkdir -p results_gaussian_hn
	../regressiontdtwave2d -i syntheticobs_gaussian.txt \
	-M priorproposal.txt \
	-o results_gaussian_hn/ \
	-x 5 -y 5 -t 1000000 -w 4 \
	-v 10000 \
	-L 0.01

results_gaussian_hn/mean.txt : results_gaussian_hn/ch.dat
	../postprocess_mean -w 4 -x 5 -y 5 -i $< -o $@ -t 5 -s 500000


#
# Parallel Constant
#
results_constant_mpi/ch.dat-000 : syntheticobs_constant.txt
	mkdir -p results_constant_mpi
	mpirun -np 4 ../regressiontdtwave2dpt -c 4 -T 1 \
	-i syntheticobs_constant.txt \
	-M priorproposal.txt \
	-o results_constant_mpi/ \
	-x 5 -y 5 -t 100000 -w 4 

results_constant_mpi/mean.txt : results_constant_mpi/ch.dat-000
	mpirun -np 4 ../postprocess_mean_mpi -i results_constant_mpi/ch.dat -o $@ \
	-w 4 -x 5 -y 5 -t 5 -s 50000

#
# Parallel Gaussian
#
results_gaussian_mpi/ch.dat-000 : syntheticobs_gaussian.txt
	mkdir -p results_gaussian_mpi
	mpirun --oversubscribe -np 4 ../regressiontdtwave2dpt -c 4 -T 1 \
	-i syntheticobs_gaussian.txt \
	-M priorproposal.txt \
	-o results_gaussian_mpi/ \
	-x 5 -y 5 -t 1000000 -v 10000 -w 4 

results_gaussian_mpi/mean.txt : results_gaussian_mpi/ch.dat-000
	mpirun --oversubscribe -np 4 ../postprocess_mean_mpi -i results_gaussian_mpi/ch.dat -o $@ \
	-w 4 -x 5 -y 5 -t 20 -s 500000

results_gaussian_mpi_cont/ch.dat-000 : results_gaussian_mpi/ch.dat-000
	mkdir -p results_gaussian_mpi_cont
	mpirun --oversubscribe -np 4 ../regressiontdtwave2dpt -c 4 -T 1 \
	-i syntheticobs_gaussian.txt \
	-I results_gaussian_mpi/ \
	-M priorproposal.txt \
	-o results_gaussian_mpi_cont/ \
	-x 5 -y 5 -t 1000000 -v 10000 -w 4 

results_gaussian_mpi_cont/mean.txt : results_gaussian_mpi_cont/ch.dat-000
	mpirun --oversubscribe -np 4 ../postprocess_mean_mpi -i results_gaussian_mpi_cont/ch.dat -o $@ \
	-w 4 -x 5 -y 5 -t 20 -s 500000
