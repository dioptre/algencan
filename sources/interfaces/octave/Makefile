BIN ?= $(CURDIR)
DIR := $(notdir $(CURDIR))

all: algencan.oct

algencan.oct: octavewrapper.o \
              $(ALGSRC)/algencan.o $(ALGSRC)/auglag.o   $(ALGSRC)/backtrack.o \
              $(ALGSRC)/betra.o    $(ALGSRC)/calc.o     $(ALGSRC)/cgm.o \
              $(ALGSRC)/checkder.o $(ALGSRC)/dogleg.o   $(ALGSRC)/drand.o \
              $(ALGSRC)/extrap.o   $(ALGSRC)/evalal.o   $(ALGSRC)/fparam.o \
              $(ALGSRC)/gencan.o   $(ALGSRC)/happrox.o  $(ALGSRC)/lss.o \
              $(ALGSRC)/moresor.o  $(ALGSRC)/newtd.o    $(ALGSRC)/newtonkkt.o \
              $(ALGSRC)/scl.o      $(ALGSRC)/setpoint.o $(ALGSRC)/sevalus.o \
              $(ALGSRC)/spgls.o    $(ALGSRC)/tevalus.o  $(ALGSRC)/tnls.o \
              $(ALGSRC)/uevalus.o  $(ALGSRC)/vevalus.o
	$(CC) -shared -loctinterp -loctave -lcruft -L$(OCTLIB) $^ -o algencan.oct

octavewrapper.o: octavewrapper.c
	g++ $(CFLAGS) -I$(OCTINC) -I$(OCTINC)/octave -c $? -o $@

install:
	mkdir -p           $(BIN)/$(DIR)
	cp -f algencanma.m $(BIN)/$(DIR)
	cp -f toyprob.m    $(BIN)/$(DIR)
	mv -f algencan.oct $(BIN)/$(DIR)

clean:
	rm -f octavewrapper.o

distclean: clean
	rm -f $(BIN)/$(DIR)/algencanma.m
	rm -f $(BIN)/$(DIR)/toyprob.m
	rm -f $(BIN)/$(DIR)/algencan.oct

.PHONY: all install clean distclean
