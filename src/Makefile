all: hmm_tmb.so 

hmm_tmb.so: hmm_tmb.cpp
	Rscript --vanilla -e "TMB::compile('hmm_tmb.cpp','-O0 -g -Wno-ignored-attributes -I ../inst/include/', safebounds=FALSE)"
	
clean:
	rm -rf *o