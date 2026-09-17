all: hmmTMB.so 

# Every header is a prerequisite: the template is one translation unit, so a
# header-only edit would otherwise be silently ignored by make.
HEADERS = hmmTMB.hpp dist.hpp dist_def.hpp added_dists.hpp likelihood.hpp

hmmTMB.so: hmmTMB.cpp $(HEADERS)
	$(R_HOME)/bin$(R_ARCH_BIN)/Rscript --vanilla -e "invisible(TMB::compile('hmmTMB.cpp', safebounds = FALSE, safeunload = FALSE))"
	
clean:
	rm -rf *o
