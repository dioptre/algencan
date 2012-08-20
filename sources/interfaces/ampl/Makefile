BIN ?= $(CURDIR)
DIR := $(notdir $(CURDIR))

all: algencan

algencan: amplwrapper.o \
          $(AMPL)/solvers/amplsolver.a \
          $(ALGSRC)/algencan.o $(ALGSRC)/auglag.o    $(ALGSRC)/backtrack.o \
          $(ALGSRC)/betra.o    $(ALGSRC)/calc.o      $(ALGSRC)/cgm.o \
          $(ALGSRC)/checkder.o $(ALGSRC)/dogleg.o    $(ALGSRC)/drand.o \
          $(ALGSRC)/extrap.o   $(ALGSRC)/evalal.o    $(ALGSRC)/fparam.o \
          $(ALGSRC)/gencan.o   $(ALGSRC)/happrox.o   $(ALGSRC)/lss.o \
          $(ALGSRC)/moresor.o  $(ALGSRC)/newtd.o     $(ALGSRC)/newtonkkt.o \
          $(ALGSRC)/scl.o      $(ALGSRC)/sevalus.o   $(ALGSRC)/spgls.o \
          $(ALGSRC)/tevalus.o  $(ALGSRC)/tnls.o      $(ALGSRC)/uevalus.o \
          $(ALGSRC)/vevalus.o
	$(FC) -ldl $^ -o $@

amplwrapper.o: amplwrapper.c amplwrapper.h ../c/cfortran.h
	$(CC) $(CFLAGS) -I$(AMPL)/solvers -c amplwrapper.c

install:
	mkdir -p       $(BIN)/$(DIR)
	mv -f algencan $(BIN)/$(DIR)

clean:
	rm -f amplwrapper.o

distclean: clean
	rm -f $(BIN)/$(DIR)/algencan

.PHONY: all install clean distclean
