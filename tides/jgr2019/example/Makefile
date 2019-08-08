
#
# The important output files to create in this demonstration inversion
#

TARGETS = results_australia/ch.dat \
	results_australia/mean_sea.txt \
	results_australia/mean_land.txt \
	results_australia/mean_tide.txt \
	mean_sea.pdf \
	mean_land.pdf \
	mean_tide.pdf

#
# The longitude/latitude extents for the Australian region
#
EXTENTS = -x 110.0 -X 160.0 \
        -y -50.0 -Y -5.0


all : $(TARGETS)

#
# Priors for the sea (GSL) and land (VLM) models. The important part is the
# second line in each. Here we use a +- 10 mm/year uniform prior
#
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

#
# Positional priors. Most important is the last line (0.5) which specifies
# the size (std. dev.) of the random proposal.
#
positionprior_sea.txt :
	echo "positionprior UniformPosition" > $@
	echo "positionproposal GaussianPosition" >> $@
	echo "0.5" >> $@

positionprior_land.txt :
	echo "positionprior UniformPosition" > $@
	echo "positionproposal GaussianPosition" >> $@
	echo "0.5" >> $@

#
# Hierarchical priors
#
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

#
# The actual inversion command that does the main inversion, important parameters are
#  -C initial number of nodes
#  -T maximum number of nodes
#  -t number of iterations
#  -v number of iterations between statistics printout
#  -z 0 Disable initial optimization
#  -A 1 Set parameterization to Delaunay (0 = Voronoi, 2 = Clough-Tocher).
#
# This inversion takes approximation 1 - 2 minutes on a modern desktop
#
# The ch.dat file created is a binary chain history of the ensemble from which
# mean/stddev etc files can be extracted. See the 3 commands following below
# for how to do this
#

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
	$(EXTENTS) \
        -C 100 \
        -T 500 \
        -t 100000 \
        -v 1000 \
        -z 0 \
        -A 1 -A 1 

#
# Extract the sea model mean of the ensemble. Important parameters are:
# -s skip this many iterations (burnin)
# -t use every ith model (thinning)
# -W/-H width and height in pixels of mean etc
# -z/Z range of histogram file
# -b number of bins in histogram
# -A parameterization (must match inversion)
# -I model to use (0 = sea, 1 = land)
#
# The next two commands do the same for land, and tide or RSL
#
results_australia/mean_sea.txt : results_australia/ch.dat
	../../post_mean -i $< -o $@ \
	$(EXTENTS) \
	-s 50000 -t 10 \
	-W 200 -H 200 \
        -D results_australia/stddev_sea.txt \
        -m results_australia/median_sea.txt \
        -g results_australia/histogram_sea.txt \
        -z -10.0 -Z 10.0 -b 400 \
	-A 1 -A 1 \
	-I 0

results_australia/mean_land.txt : results_australia/ch.dat
	../../post_mean -i $< -o $@ \
	$(EXTENTS) \
	-s 50000 -t 10 \
	-W 200 -H 200 \
        -D results_australia/stddev_land.txt \
        -m results_australia/median_land.txt \
        -g results_australia/histogram_land.txt \
        -z -10.0 -Z 10.0 -b 400 \
	-A 1 -A 1 \
	-I 1

results_australia/mean_tide.txt : results_australia/ch.dat
	../../post_mean_tide -i $< -o $@ \
	$(EXTENTS) \
	-s 50000 -t 10 \
	-W 200 -H 200 \
        -D results_australia/stddev_tide.txt \
        -m results_australia/median_tide.txt \
        -g results_australia/histogram_tide.txt \
        -z -10.0 -Z 10.0 -b 400 \
	-A 1 -A 1

#
# Plot images of the mean models using a python script
#
mean_land.pdf : results_australia/mean_land.txt
	python3 ../scripts/plot_image_map.py \
	$(EXTENTS) \
	-i $< --pdf $@ --noshow \
	--intermediate \
	--vmin -5.0 --vmax 5.0 \
	--cm Spectral 

mean_sea.pdf : results_australia/mean_sea.txt
	python3 ../scripts/plot_image_map.py \
	$(EXTENTS) \
	-i $< --pdf $@ --noshow \
	--intermediate \
	--vmin -5.0 --vmax 5.0 \
	--cm Spectral 

mean_tide.pdf : results_australia/mean_tide.txt
	python3 ../scripts/plot_image_map.py \
	$(EXTENTS) \
	-i $< --pdf $@ --noshow \
	--intermediate \
	--vmin -5.0 --vmax 5.0 \
	--cm Spectral 
