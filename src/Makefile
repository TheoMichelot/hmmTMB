all: hmmTMB.so 

hmmTMB.so: hmmTMB.cpp
	$(R_HOME)/bin/Rscript --vanilla -e "invisible(TMB::compile('hmmTMB.cpp', safebounds = FALSE, safeunload = FALSE))"
	
clean:
	rm -rf *o