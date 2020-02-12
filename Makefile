# ---------------------------------------------------------------------
# Platform independent Makefile
#
# Author: Benoit Sklenard <benoit.sklenard@cea.fr>
#
# October 15, 2014 (v1)
# December 20, 2014 (v2)
#

SHELL     := /bin/bash
MAKE      := make
ROOT      := mmonca

MC_BIN    := $(ROOT)

MACHINES  := $(shell find ./sysmakes/ -name "Makefile.*" | cut -d '.' -f3-)
VERSION   := $(shell cat src/version.h | cut -d '"' -f 2)
DNAME     := mmonca-$(VERSION)

SILENT    ?= yes

ifeq ($(SILENT),yes)
	MAKESILENT = -s
endif

help:
	@echo ''
	@echo 'Usage: '
	@echo '  make dist                  create $(DNAME).tar.gz'
	@echo '  make clean-all             delete all object files'
	@echo '  make clean-<machine>       delete object files from MMonCa compilation for <machine>'
	@echo '  make <machine>             builds MMonCa for <machine>'
	@echo ''
	@echo '  <machine> can be one of these in src/sysmakes:'
	@$(foreach m,$(MACHINES),echo '  - $(m)';)

dist: 
	@mkdir $(DNAME)	
	@cp -r src/ test/ config/ examples/ sysmakes/ doc/ Makefile* $(DNAME)
	@cp AUTHORS README LICENSE NEWS NOTICE $(DNAME)
	@rm -f `find $(DNAME)/test -name test.mc.log`
	@rm -f `find $(DNAME)/test -name time`
	@rm -f `find $(DNAME)/test -name errors`
	@rm -f `find $(DNAME)/test -name *nodist* `
	@rm -f `find $(DNAME)/test -name test.result `
	@rm -f `find $(DNAME)/test -name Inel*`
	@rm -f `find $(DNAME)/test -name EDT_*`
	@rm -f `find $(DNAME)/test -name Zel*`	
	@rm -f  $(DNAME)/examples/*/*.dump
	@rm -f $(DNAME)/doc/images/triptico*
	@tar cvzf $(DNAME).tar.gz \
	$(DNAME)/AUTHORS $(DNAME)/README $(DNAME)/LICENSE $(DNAME)/NEWS $(DNAME)/NOTICE \
	$(DNAME)/src/*.cpp $(DNAME)/src/*.h $(DNAME)/src/*/*.cpp $(DNAME)/src/*/*.h \
	$(DNAME)/test/standard/ $(DNAME)/test/*.sh \
	$(DNAME)/config/*AmorphousSilicon*/ $(DNAME)/config/Silicon/ $(DNAME)/config/Gas/ $(DNAME)/config/Nitride/ \
	$(DNAME)/config/Nitride_SiO2/ $(DNAME)/config/Gas_Silicon/ $(DNAME)/config/S_Iron/ $(DNAME)/config/Gas_S_Iron/ \
	$(DNAME)/config/MC/ $(DNAME)/config/SiO2/ $(DNAME)/config/SiO2_Silicon/ $(DNAME)/config/Mechanics/ \
	$(DNAME)/sysmakes/ $(DNAME)/Makefile $(DNAME)/Makefile.mmonca \
	$(DNAME)/examples/NatureMaterials2004/ \
	$(DNAME)/doc/document.sh $(DNAME)/doc/*.tex $(DNAME)/doc/*/*.tex $(DNAME)/doc/images/ $(DNAME)/doc/MMonCa.pdf

dist-static: 
	@mkdir $(DNAME)	
	@cp -r test/ config/ examples/ doc/ data/ $(DNAME)
	@cp AUTHORS README LICENSE NEWS NOTICE $(DNAME)
	@cp Obj_static/mmonca  $(DNAME)/mmonca
	@rm -f `find $(DNAME)/test -name test.mc.log`
	@rm -f `find $(DNAME)/test -name time`
	@rm -f `find $(DNAME)/test -name errors`
	@rm -f `find $(DNAME)/test -name *nodist* `
	@rm -f `find $(DNAME)/test -name test.result `
	@rm -f `find $(DNAME)/test -name Inel*`
	@rm -f `find $(DNAME)/test -name EDT_*`
	@rm -f `find $(DNAME)/test -name Zel*`	
	@rm -f $(DNAME)/examples/*/*.dump
	@rm -f $(DNAME)/doc/images/triptico*
	@tar cvzf $(DNAME)-static.tar.gz \
	$(DNAME)/AUTHORS $(DNAME)/README $(DNAME)/LICENSE $(DNAME)/NEWS $(DNAME)/NOTICE \
	$(DNAME)/test/standard/ $(DNAME)/test/*.sh \
	$(DNAME)/config/*AmorphousSilicon*/ $(DNAME)/config/Silicon/ $(DNAME)/config/Gas/ $(DNAME)/config/Nitride/ \
	$(DNAME)/config/Nitride_SiO2/ $(DNAME)/config/Gas_Silicon/ $(DNAME)/config/S_Iron/ $(DNAME)/config/Gas_S_Iron/ \
	$(DNAME)/config/MC/ $(DNAME)/config/SiO2/ $(DNAME)/config/SiO2_Silicon/ $(DNAME)/config/Mechanics/ \
	$(DNAME)/examples/NatureMaterials2004/ \
	$(DNAME)/doc/MMonCa.pdf $(DNAME)/mmonca

.DEFAULT:
	@if [ -f ./sysmakes/Makefile.$@ ]; \
	then \
		mkdir -p Obj_$@ ; \
		cp ./sysmakes/Makefile.$@ Obj_$@/Makefile ; \
	else \
		echo 'Machine $@ not found in ./sysmakes/Makefile/*' ; \
		false ; \
	fi ;
	@echo "Compiling MMonCa for machine $@..."
	@$(MAKE) $(MAKESILENT) -C Obj_$@ "EXE = $(MC_BIN)" "TARGET = mmonca"
	@echo "Built target $(shell pwd)/Obj_$@/$(MC_BIN)"

clean-all:
	rm -rf Obj_*

clean-%:
	rm -rf Obj_$(@:clean-%=%)

