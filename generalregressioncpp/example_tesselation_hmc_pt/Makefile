
TARGETS = datatemplate.txt \
	syntheticobs_tesselation.txt  

NOISE = 0.05

all : $(TARGETS)

datatemplate.txt :
	python2 ../../scripts/generatetemplatepoints.py -N 100 -o $@ \
	--seed 983 \
	--xmin 0.0 --xmax 1.0 \
	--ymin 0.0 --ymax 1.0 

syntheticobs_tesselation.txt : datatemplate.txt 
	../mksynthetic -m Tesselation \
	-i datatemplate.txt \
	-x 0.0 -X 1.0 \
	-y 0.0 -Y 1.0 \
	-o $@ -n $(NOISE) \
	-O $@.true \
	-I syntheticobs_tesselation.img \
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

