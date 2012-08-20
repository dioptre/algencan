BIN ?= $(CURDIR)
DIR := $(notdir $(CURDIR))

all: algencan.so

algencan.so: rwrapper.o \
          $(ALGSRC)/algencan.o $(ALGSRC)/auglag.o    $(ALGSRC)/backtrack.o \
          $(ALGSRC)/betra.o    $(ALGSRC)/calc.o      $(ALGSRC)/cgm.o \
          $(ALGSRC)/checkder.o $(ALGSRC)/dogleg.o    $(ALGSRC)/drand.o \
          $(ALGSRC)/extrap.o   $(ALGSRC)/evalal.o    $(ALGSRC)/fparam.o \
          $(ALGSRC)/gencan.o   $(ALGSRC)/happrox.o   $(ALGSRC)/lss.o \
          $(ALGSRC)/moresor.o  $(ALGSRC)/newtd.o     $(ALGSRC)/newtonkkt.o \
          $(ALGSRC)/scl.o      $(ALGSRC)/setpoint.o  $(ALGSRC)/sevalus.o \
          $(ALGSRC)/spgls.o    $(ALGSRC)/tevalus.o   $(ALGSRC)/tnls.o \
          $(ALGSRC)/uevalus.o  $(ALGSRC)/vevalus.o
	R CMD SHLIB -shared -lgfortran $^ -o $@

rwrapper.o: rwrapper.c rwrapper.h cfortran.h
	$(CC) $(CFLAGS) -I$(RINC) -c rwrapper.c

install:
	mkdir -p           $(BIN)/$(DIR)
	cp -f algencanma.r $(BIN)/$(DIR)
	cp -f toyprob.r    $(BIN)/$(DIR)
	mv -f algencan.so  $(BIN)/$(DIR)

clean:
	rm -f rwrapper.o

distclean: clean
	rm -f $(BIN)/$(DIR)/algencanma.r
	rm -f $(BIN)/$(DIR)/toyprob.r
	rm -f $(BIN)/$(DIR)/algencan.so

.PHONY: all sifdec install clean distclean
