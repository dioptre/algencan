BIN ?= $(CURDIR)
DIR := $(notdir $(CURDIR))

all: pywrapper.so

pywrapper.so: pywrapper.o \
          $(ALGSRC)/algencan.o $(ALGSRC)/auglag.o   $(ALGSRC)/backtrack.o \
          $(ALGSRC)/betra.o    $(ALGSRC)/calc.o     $(ALGSRC)/cgm.o \
          $(ALGSRC)/checkder.o $(ALGSRC)/dogleg.o   $(ALGSRC)/drand.o \
          $(ALGSRC)/extrap.o   $(ALGSRC)/evalal.o   $(ALGSRC)/fparam.o \
          $(ALGSRC)/gencan.o   $(ALGSRC)/happrox.o  $(ALGSRC)/lss.o \
          $(ALGSRC)/moresor.o  $(ALGSRC)/newtd.o    $(ALGSRC)/newtonkkt.o \
          $(ALGSRC)/scl.o      $(ALGSRC)/setpoint.o $(ALGSRC)/sevalus.o \
          $(ALGSRC)/spgls.o    $(ALGSRC)/tevalus.o  $(ALGSRC)/tnls.o \
          $(ALGSRC)/uevalus.o  $(ALGSRC)/vevalus.o
	$(FC) $(LDFLAGS) $^ -o $@

pywrapper.o: pywrapper.c pywrapper.h ../c/cfortran.h
	$(CC) $(CFLAGS) -I$(PYTHONINC) \
	                -I$(PYTHONLIB)/site-packages/numpy/core/include \
	                -c pywrapper.c

install:
	mkdir -p           $(BIN)/$(DIR)
	cp -f algencan.py  $(BIN)/$(DIR)
	mv -f pywrapper.so $(BIN)/$(DIR)

clean:
	rm -f pywrapper.o

distclean: clean
	rm -f $(BIN)/$(DIR)/algencan.py*
	rm -f $(BIN)/$(DIR)/pywrapper.so

.PHONY: all install clean distclean
