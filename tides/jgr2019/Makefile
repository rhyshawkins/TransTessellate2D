

DATA = $(wildcard data/*.txt) \
	$(wildcard points/*.txt) \
	$(wildcard results/*/*.txt) \
	$(wildcard results/*.sh) \
	$(wildcard scripts/*.py) \
	$(wildcard slurm/*/*.slurm)

INSTALL = install
INSTALLFLAGS = -D
DATE = $(shell date +"%Y%m%d%H%M")
DIR = VirtualTideGaugeJGR2019
TGZ = $(DIR).tar.gz

default :
	@echo "make dist to create the data bundle"

dist : 
	mkdir -p $(DIR)
	echo $(DATE) > $(DIR)/Version
	for f in $(DATA); do \
	    $(INSTALL) $(INSTALLFLAGS) $$f $(DIR)/$$f ; \
	done
	$(INSTALL) $(INSTALLFLAGS) ../README.md $(DIR)/README.md
	tar -czf $(TGZ) $(DIR)/*
	rm -rf $(DIR)

