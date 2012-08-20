BIN ?= $(CURDIR)
DIR := $(notdir $(CURDIR))

CONFIGFILE := ./mexopts.sh

ifeq ($(64BITS),true)
	MATLABWRAPPERFILES := matlabwrapper64.f checkderm64.f fparamm64.f
	FLAGS              := -largeArrayDims
else
	MATLABWRAPPERFILES := matlabwrapper32.f checkderm32.f fparamm32.f
endif


all: algencan

algencan: $(MATLABWRAPPERFILES) \
          $(ALGSRC)/algencan.o  $(ALGSRC)/auglag.o    $(ALGSRC)/backtrack.o \
          $(ALGSRC)/betra.o     $(ALGSRC)/calc.o      $(ALGSRC)/cgm.o \
          $(ALGSRC)/dogleg.o    $(ALGSRC)/drand.o     $(ALGSRC)/extrap.o \
          $(ALGSRC)/evalal.o    $(ALGSRC)/gencan.o    $(ALGSRC)/happrox.o \
          $(ALGSRC)/lss.o       $(ALGSRC)/moresor.o   $(ALGSRC)/newtd.o \
          $(ALGSRC)/newtonkkt.o $(ALGSRC)/scl.o       $(ALGSRC)/setpoint.o \
          $(ALGSRC)/sevalus.o   $(ALGSRC)/spgls.o     $(ALGSRC)/tevalus.o \
          $(ALGSRC)/tnls.o      $(ALGSRC)/uevalus.o   $(ALGSRC)/vevalus.o
	mex -f $(CONFIGFILE) -I$(ALGSRC) $^ $(FLAGS) -output $@

install:
	mkdir -p               $(BIN)/$(DIR)
	cp    -f algencan.mex* $(BIN)/$(DIR)

clean:
	rm -f algencan.mex*

distclean: clean
	rm -f $(BIN)/$(DIR)/algencan.mex*

.PHONY: all install clean distclean
