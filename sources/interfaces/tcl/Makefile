BIN ?= $(CURDIR)
DIR := $(notdir $(CURDIR))

LIBNAME := tclwrapper.so

all: $(LIBNAME)

$(LIBNAME): tclwrapper.o \
          $(ALGSRC)/algencan.o $(ALGSRC)/auglag.o    $(ALGSRC)/backtrack.o \
          $(ALGSRC)/betra.o    $(ALGSRC)/calc.o      $(ALGSRC)/cgm.o \
          $(ALGSRC)/checkder.o $(ALGSRC)/dogleg.o    $(ALGSRC)/drand.o \
          $(ALGSRC)/extrap.o   $(ALGSRC)/evalal.o    $(ALGSRC)/fparam.o \
          $(ALGSRC)/gencan.o   $(ALGSRC)/happrox.o   $(ALGSRC)/lss.o \
          $(ALGSRC)/moresor.o  $(ALGSRC)/newtd.o     $(ALGSRC)/newtonkkt.o \
          $(ALGSRC)/scl.o      $(ALGSRC)/setpoint.o  $(ALGSRC)/sevalus.o \
          $(ALGSRC)/spgls.o    $(ALGSRC)/tevalus.o   $(ALGSRC)/tnls.o \
          $(ALGSRC)/uevalus.o  $(ALGSRC)/vevalus.o
	$(CC) -shared $^ -o $@ -lgfortran -DUSE_TCL_STUBS -I$(TCLINC) \
              -L$(TCLLIB) -ltclstub8.5 -lm -pedantic

tclwrapper.o: tclwrapper.c ../c/cfortran.h
	$(CC) $(CFLAGS) -c tclwrapper.c ../c/cfortran.h -DUSE_TCL_STUBS \
              -I$(TCLINC) -L$(TCLLIB) 

install:
	mkdir -p         $(BIN)/$(DIR)
	cp -f $(LIBNAME) $(BIN)/$(DIR)

clean:
	rm -f tclwrapper.o $(LIBNAME)

distclean: clean
	rm -f $(BIN)/$(DIR)/$(LIBNAME)

.PHONY: all sifdec install clean distclean
