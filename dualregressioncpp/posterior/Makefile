
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

priorproposal.txt :
	echo "prior Uniform" > $@
	echo "0.0 1.0" >> $@
	echo "proposal Gaussian" >> $@
	echo "0.05" >> $@

positionpriorproposal.txt :
	echo "positionprior UniformPosition" > $@
	echo "positionproposal GaussianPosition" >> $@
	echo "0.2" >> $@

#
# Voronoi Constant True model
#
results_constant/ch.dat : syntheticobs_constant.txt \
	priorproposal.txt \
	positionpriorproposal.txt
	mkdir -p results_constant
	../regressionvoronoicartesian -i syntheticobs_constant.txt \
	-P priorproposal.txt -P priorproposal.txt \
	-M positionpriorproposal.txt -M positionpriorproposal.txt \
	-o results_constant/ \
	-x -10.0 -X 10.0 \
	-y -10.0 -Y 10.0 \
	-t 5000000 \
	-T 30 \
	-v 500 \
	-p

