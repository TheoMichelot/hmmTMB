all: hmmTMB.so 

hmmTMB.so: hmmTMB.cpp
	$(R_HOME)/bin$(R_ARCH_BIN)/Rscript --vanilla -e "invisible(TMB::compile('hmmTMB.cpp', safebounds = FALSE, safeunload = FALSE))"
	
clean:
	rm -rf *o