EDTDIR := $(shell pwd)/../../data/bca/
TABLES3DDIR :=  $(shell pwd)/../../data/bca/
IIS := $(shell pwd)/../../Obj_bca_g++/iis
export EDTDIR
export IIS
export TABLES3DDIR

all: si_ sic_ gaas_

si_ :
	cd Si;\
	make 

sic_ :
	cd SiC;\
	make

gaas_ :
	cd GaAs;\
	make

plot:
	cd Si;\
	make plot
	cd SiC;\
	make plot
	cd GaAs;\
	make plot

clean:
	cd Si;\
	make clean
	cd SiC;\
	make clean
	cd GaAs;\
	make clean
