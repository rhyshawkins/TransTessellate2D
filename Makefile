

LIBSRC= $(wildcard lib/*.c) \
	$(wildcard lib/*.h) \
	lib/Makefile 

BASESRC = base/Makefile \
	$(wildcard base/*.?pp)

GENERALREGRESSIONSRC=generalregressioncpp/genericregression.cpp \
	generalregressioncpp/Makefile \
	generalregressioncpp/example/Makefile \
	generalregressioncpp/posterior/Makefile \
	generalregressioncpp/example_franke_hmc_pt/syntheticobs_franke.txt \
	generalregressioncpp/example_pseudo_franke_hmc_pt/syntheticobs_pseudofranke.txt \
	generalregressioncpp/example_tesselation_hmc_pt/syntheticobs_tesselation.txt \
	$(wildcard generalregressioncpp/transcale/*.slurm) \
	$(wildcard generalregressioncpp/transcale/tidegauge/*.slurm) 

SCRIPTS=scripts/generatedualtemplatepoints.py \
	scripts/generatetemplatepoints.py

TIDESRC=tides/Makefile \
	tides/tides.cpp \
	tides/tidesynthetic.cpp \
	tides/tidesynthetic.hpp \
	tides/scripts/generateislandtemplatepoints.py \
	tides/synthetic_constant/Makefile \
	tides/synthetic_cosine1/Makefile \
	tides/synthetic_cosine2/Makefile \
	tides/synthetic_eastwest/Makefile \
	$(wildcard tides/transcale/asia/*.slurm) \
	$(wildcard tides/transcale/europe/*.slurm) \
	$(wildcard tides/transcale/northamerica/*.slurm) \
	$(wildcard tides/transcale/africa/*.slurm) \
	$(wildcard tides/transcale/australia/*.slurm) \
	$(wildcard tides/transcale/southamerica/*.slurm) \
	tides/tas_synthetic/Makefile \
	tides/tas_synthetic/syntheticobs_50_50_50_tas.txt \
	$(wildcard tides/tas_synthetic/transcale/*.slurm)

SRCS=Makefile \
	LICENSE \
	README.md \
	$(LIBSRC) \
	$(BASESRC) \
	$(GENERALREGRESSIONSRC) \
	$(SCRIPTS) \
	$(TIDESRC)

INSTALL = install
INSTALLFLAGS = -D
DATE = $(shell date +"%Y%m%d%H%M")
DIR = TransTessellate2D
TGZ = $(DIR).tar.gz

all :
	make -C lib
	make -C base
	make -C generalregressioncpp
	make -C tides

clean :
	make -C lib clean
	make -C base clean
	make -C generalregressioncpp clean
	make -C tides clean

dist : 
	mkdir -p $(DIR)
	echo $(DATE) > $(DIR)/Version
	for f in Makefile $(SRCS); do \
	    $(INSTALL) $(INSTALLFLAGS) $$f $(DIR)/$$f ; \
	done
	tar -czf $(TGZ) $(DIR)/*
	rm -rf $(DIR)

