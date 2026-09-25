all: hmmTMB.so 

# TMBad rather than TMB's own CppAD default. A Gaussian field has to
# factorise a parameter-dependent sparse precision inside the likelihood,
# which CppAD records operation by operation onto the tape. Set
# HMMTMB_AD_FRAMEWORK=CppAD before installing to fall back.

# Every header is a prerequisite: the template is one translation unit, so a
# header-only edit would otherwise be silently ignored by make.
HEADERS = hmmTMB.hpp dist.hpp dist_def.hpp added_dists.hpp likelihood.hpp

hmmTMB.so: hmmTMB.cpp $(HEADERS)
	$(R_HOME)/bin$(R_ARCH_BIN)/Rscript --vanilla -e "invisible(TMB::compile('hmmTMB.cpp', safebounds = FALSE, safeunload = FALSE, framework = Sys.getenv('HMMTMB_AD_FRAMEWORK', 'TMBad')))"
	
clean:
	rm -rf *o
