
TARGETS = \
	sdf_tas_fine.txt \
	pointtemplate_50_50_50.txt \
	syntheticobs_50_50_50_tas.txt 


all : $(TARGETS)

MINLON=144.0
MAXLON=149.0
MINLAT=-44.0
MAXLAT=-40.0

#
# SDF for region
#
sdf_tas_coarse.txt : /data/rhys/coastlines/GSHHS_shp/f/GSHHS_f_L1.shp
	$(HOME)/PostDocLyon/sealevel/sdfshp/mksdfregion \
	-n $(MINLON) -N $(MAXLON) \
	-a $(MINLAT) -A $(MAXLAT) \
	-i $< -o $@ \
	-W 50 -H 40

sdf_tas_intermediate.txt : /data/rhys/coastlines/GSHHS_shp/f/GSHHS_f_L1.shp
	$(HOME)/PostDocLyon/sealevel/sdfshp/mksdfregion \
	-n $(MINLON) -N $(MAXLON) \
	-a $(MINLAT) -A $(MAXLAT) \
	-i $< -o $@ \
	-W 250 -H 200

sdf_tas_fine.txt : /data/rhys/coastlines/GSHHS_shp/f/GSHHS_f_L1.shp
	$(HOME)/PostDocLyon/sealevel/sdfshp/mksdfregion \
	-n $(MINLON) -N $(MAXLON) \
	-a $(MINLAT) -A $(MAXLAT) \
	-i $< -o $@ \
	-W 500 -H 400

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
# Synthetic optimisation test
#
results_50_50_50_opt/ch.dat : syntheticobs_50_50_50_tas.txt \
	prior_land.txt \
	prior_sea.txt \
	positionprior_land.txt \
	positionprior_sea.txt \
	hierarchical_land.txt \
	hierarchical_sea.txt
	mkdir -p results_50_50_50_opt
	../tides -i syntheticobs_50_50_50_tas.txt -o results_50_50_50_opt/ \
	-P prior_sea.txt -P prior_land.txt \
	-M positionprior_sea.txt -M positionprior_land.txt \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-A 1 -A 0 \
	-C 10 \
	-T 100 \
	-t 1 \
	-v 1000 \
	-z 10000 \
	-Z 5.0

results_50_50_50_opt/mean_sea.txt : results_50_50_50_opt/ch.dat
	../post_mean -i $< -o $@ \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-T 100 \
	-A 1 -A 0 \
	-W 50 -H 50 \
	-I 0 \
	-s 0 -t 0

results_50_50_50_opt/mean_land.txt : results_50_50_50_opt/ch.dat
	../post_mean -i $< -o $@ \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-T 100 \
	-A 1 -A 0 \
	-W 50 -H 50 \
	-I 1 \
	-s 0 -t 0

#
# Synthetic mcmc test Delaunay/Voronoi
#
results_50_50_50_dv_mcmc/ch.dat : syntheticobs_50_50_50_tas.txt \
	prior_land.txt \
	prior_sea.txt \
	positionprior_land.txt \
	positionprior_sea.txt \
	hierarchical_land.txt \
	hierarchical_sea.txt
	mkdir -p results_50_50_50_dv_mcmc
	../tides -i syntheticobs_50_50_50_tas.txt -o results_50_50_50_dv_mcmc/ \
	-P prior_sea.txt -P prior_land.txt \
	-M positionprior_sea.txt -M positionprior_land.txt \
	-H hierarchical_sea.txt -H hierarchical_land.txt \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-A 1 -A 0 \
	-C 10 \
	-T 100 \
	-t 1000000 \
	-v 1000 \
	-z 10000 \
	-Z 5.0

results_50_50_50_dv_mcmc/mean_sea.txt : results_50_50_50_dv_mcmc/ch.dat
	../post_mean -i $< -o $@ \
	-D results_50_50_50_dv_mcmc/stddev_sea.txt \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-T 100 \
	-A 1 -A 0 \
	-W 50 -H 50 \
	-I 0 \
	-s 500000 -t 10

results_50_50_50_dv_mcmc/mean_land.txt : results_50_50_50_dv_mcmc/ch.dat
	../post_mean -i $< -o $@ \
	-D results_50_50_50_dv_mcmc/stddev_sea.txt \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-T 100 \
	-A 1 -A 0 \
	-W 50 -H 50 \
	-I 1 \
	-s 500000 -t 10

results_50_50_50_dv_mcmc/mean_tide.txt : results_50_50_50_dv_mcmc/ch.dat
	../post_mean_tide -i $< -o $@ \
	-D results_50_50_50_dv_mcmc/stddev_sea.txt \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-T 100 \
	-A 1 -A 0 \
	-W 50 -H 50 \
	-s 500000 -t 10

#
# Synthetic mcmc test Clough-Tocher/Voronoi
#
results_50_50_50_cv_mcmc/ch.dat : syntheticobs_50_50_50_tas.txt \
	prior_land.txt \
	prior_sea.txt \
	positionprior_land.txt \
	positionprior_sea.txt \
	hierarchical_land.txt \
	hierarchical_sea.txt
	mkdir -p results_50_50_50_cv_mcmc
	../tides -i syntheticobs_50_50_50_tas.txt -o results_50_50_50_cv_mcmc/ \
	-P prior_sea.txt -P prior_land.txt \
	-M positionprior_sea.txt -M positionprior_land.txt \
	-H hierarchical_sea.txt -H hierarchical_land.txt \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-A 2 -A 0 \
	-C 10 \
	-T 100 \
	-t 1000000 \
	-v 1000 \
	-z 10000 \
	-Z 5.0

results_50_50_50_cv_mcmc/mean_sea.txt : results_50_50_50_cv_mcmc/ch.dat
	../post_mean -i $< -o $@ \
	-D results_50_50_50_cv_mcmc/stddev_sea.txt \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-T 100 \
	-A 2 -A 0 \
	-W 50 -H 50 \
	-I 0 \
	-s 500000 -t 10

results_50_50_50_cv_mcmc/mean_land.txt : results_50_50_50_cv_mcmc/ch.dat
	../post_mean -i $< -o $@ \
	-D results_50_50_50_cv_mcmc/stddev_land.txt \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-T 100 \
	-A 2 -A 0 \
	-W 50 -H 50 \
	-I 1 \
	-s 500000 -t 10

results_50_50_50_cv_mcmc/mean_tide.txt : results_50_50_50_cv_mcmc/ch.dat
	../post_mean_tide -i $< -o $@ \
	-D results_50_50_50_cv_mcmc/stddev_tide.txt \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-T 100 \
	-A 2 -A 0 \
	-W 50 -H 50 \
	-s 500000 -t 10

#
# Synthetic mcmc test Voronoi/Voronoi
#
results_50_50_50_vv_mcmc/ch.dat : syntheticobs_50_50_50_tas.txt \
	prior_land.txt \
	prior_sea.txt \
	positionprior_land.txt \
	positionprior_sea.txt \
	hierarchical_land.txt \
	hierarchical_sea.txt
	mkdir -p results_50_50_50_vv_mcmc
	../tides -i syntheticobs_50_50_50_tas.txt -o results_50_50_50_vv_mcmc/ \
	-P prior_sea.txt -P prior_land.txt \
	-M positionprior_sea.txt -M positionprior_land.txt \
	-H hierarchical_sea.txt -H hierarchical_land.txt \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-A 0 -A 0 \
	-C 10 \
	-T 100 \
	-t 1000000 \
	-v 1000 \
	-z 10000 \
	-Z 5.0

results_50_50_50_vv_mcmc/mean_sea.txt : results_50_50_50_vv_mcmc/ch.dat
	../post_mean -i $< -o $@ \
	-D results_50_50_50_vv_mcmc/stddev_sea.txt \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-T 100 \
	-A 0 -A 0 \
	-W 50 -H 50 \
	-I 0 \
	-s 500000 -t 10

results_50_50_50_vv_mcmc/mean_land.txt : results_50_50_50_vv_mcmc/ch.dat
	../post_mean -i $< -o $@ \
	-D results_50_50_50_vv_mcmc/stddev_land.txt \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-T 100 \
	-A 0 -A 0 \
	-W 50 -H 50 \
	-I 1 \
	-s 500000 -t 10

results_50_50_50_vv_mcmc/mean_tide.txt : results_50_50_50_vv_mcmc/ch.dat
	../post_mean_tide -i $< -o $@ \
	-D results_50_50_50_vv_mcmc/stddev_tide.txt \
	-x $(MINLON) -X $(MAXLON) \
	-y $(MINLAT) -Y $(MAXLAT) \
	-T 100 \
	-A 0 -A 0 \
	-W 50 -H 50 \
	-s 500000 -t 10
