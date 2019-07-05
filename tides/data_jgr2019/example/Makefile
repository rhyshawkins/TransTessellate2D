
TARGETS = results_australia/ch.dat

all : $(TARGETS)

prior_sea.txt :
	echo "prior Uniform" > $@
	echo "-10.0 10.0" >> $@
	echo "proposal Gaussian" >> $@
	echo "0.2" >> $@

prior_land.txt :
	echo "prior Uniform" > $@
	echo "-10.0 10.0" >> $@
	echo "proposal Gaussian" >> $@
	echo "0.2" >> $@

positionprior_sea.txt :
	echo "positionprior UniformPosition" > $@
	echo "positionproposal GaussianPosition" >> $@
	echo "0.5" >> $@

positionprior_land.txt :
	echo "positionprior UniformPosition" > $@
	echo "positionproposal GaussianPosition" >> $@
	echo "0.5" >> $@

hierarchical_prior_tide.txt :
	echo "prior Jeffreys" > $@
	echo "0.01 10.0" >> $@
	echo "proposal LogGaussian" >> $@
	echo "0.1" >> $@

hierarchical_prior_land.txt :
	echo "prior Jeffreys" > $@
	echo "0.01 10.0" >> $@
	echo "proposal LogGaussian" >> $@
	echo "0.025" >> $@

hierarchical_prior_sea.txt :
	echo "prior Jeffreys" > $@
	echo "0.01 10.0" >> $@
	echo "proposal LogGaussian" >> $@
	echo "0.1" >> $@

results_australia/ch.dat : prior_sea.txt \
	prior_land.txt \
	positionprior_sea.txt \
	positionprior_land.txt \
	hierarchical_prior_tide.txt \
	hierarchical_prior_land.txt \
	hierarchical_prior_sea.txt \
	../data/australia_combined.txt
	mkdir -p results_australia
	../../tideshmc \
	-o results_australia/ \
	-P prior_sea.txt -P prior_land.txt \
	-M positionprior_sea.txt -M positionprior_land.txt \
	-H hierarchical_prior_tide.txt -H hierarchical_prior_land.txt -H hierarchical_prior_sea.txt \
	-i ../data/australia_combined.txt \
	-x 110.0 -X 160.0 \
        -y -50.0 -Y -5.0 \
        -C 100 \
        -T 500 \
        -t 100000 \
        -v 1000 \
        -z 0 \
        -Z 5.0 \
        -A 1 -A 1 \
        -S 1.0e-2

results_australia/mean_sea.txt : results_australia/ch.dat
	../../post_mean -i $< -o $@ \
	-x 110.0 -X 160.0 \
        -y -50.0 -Y -5.0 \
	-s 50000 -t 10 \
	-A 1 -A 1

