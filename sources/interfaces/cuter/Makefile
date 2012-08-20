BIN ?= $(CURDIR)
DIR := $(notdir $(CURDIR))
SIF := ELFUN EXTER GROUP RANGE

all: sifdec algencan

sifdec:
	./call_sifdec $(PROBNAME) $(MASTSIF) $(SIFDEC)
	$(FC) -c $(FFLAGS) $(addsuffix .f,$(SIF))

algencan: cuterwrapper.o $(addsuffix .o,$(SIF)) \
          $(ALGSRC)/algencanma.o $(ALGSRC)/algencan.o $(ALGSRC)/auglag.o \
          $(ALGSRC)/backtrack.o  $(ALGSRC)/betra.o    $(ALGSRC)/calc.o \
          $(ALGSRC)/cgm.o        $(ALGSRC)/checkder.o $(ALGSRC)/dogleg.o \
          $(ALGSRC)/drand.o      $(ALGSRC)/evalal.o   $(ALGSRC)/extrap.o \
          $(ALGSRC)/fparam.o     $(ALGSRC)/gencan.o   $(ALGSRC)/happrox.o \
          $(ALGSRC)/lss.o        $(ALGSRC)/moresor.o  $(ALGSRC)/newtd.o \
          $(ALGSRC)/newtonkkt.o  $(ALGSRC)/scl.o      $(ALGSRC)/setpoint.o \
          $(ALGSRC)/sevalus.o    $(ALGSRC)/spgls.o    $(ALGSRC)/tevalus.o \
          $(ALGSRC)/tnls.o       $(ALGSRC)/uevalus.o  $(ALGSRC)/vevalus.o
	$(FC) -i_dynamic $^ $(CUTER)/double/bin/linpac.o \
	      -L$(CUTER)/double/lib -lcuter -o $@

cuterwrapper.o: cuterwrapper.f cuter.com ../../algencan/dim.par

install:
	mkdir -p        $(BIN)/$(DIR)
	mv -f OUTSDIF.d $(BIN)/$(DIR)
	mv -f algencan  $(BIN)/$(DIR)

clean:
	rm -f *.o
	rm -f $(addsuffix .f,$(SIF))

distclean: clean
	rm -f $(BIN)/$(DIR)/OUTSDIF.d
	rm -f $(BIN)/$(DIR)/algencan

.PHONY: all sifdec install clean distclean
