# Set variables FC and CC with the Fortran and C/C++ compilers of your choice.
# g77 and gfortran are already tested valid options for FC. gcc and g++ are
# also already tested valid options for CC. Setting CC is not necessary to run
# the stand-alone Fortran version of ALGENCAN. Leave it blank if this is your
# case.
FC := gfortran
CC := gcc

FFLAGS  := -O4 -xf77-cpp-input -fPIC
CFLAGS  := -O4 -Df2cFortran -fPIC
LDFLAGS := -O4 -shared

# Set variable ALGENCAN with the absolute path of your ALGENCAN installation
# directory. The value shown expands to the path of current working directory.
ALGENCAN := $(CURDIR)

BIN      := $(ALGENCAN)/bin
SRC      := $(ALGENCAN)/sources
ALGSRC   := $(SRC)/algencan
HSLSRC   := $(SRC)/hsl
PROBSRC  := $(SRC)/problems
INTERSRC := $(SRC)/interfaces
INTERFCS := $(notdir $(wildcard $(INTERSRC)/*))

# Set the variables below with the absolute paths of your tools. The values
# shown are mere examples. Leave them blank if you feel that you will not use
# them. None of them is mandatory to use the stand-alone Fortran version of
# ALGENCAN.
AMPL     := $(HOME)/misc/ampl

MASTSIF  := $(HOME)/CUTEr/MastSIF/mastsif
SIFDEC   := $(HOME)/CUTEr/SifDec/SifDec.custom.pc.lnx.gfo
CUTER    := $(HOME)/CUTEr/CUTEr/cuter-export/CUTEr.custom.pc.lnx.gfo

PYTHONINC := /usr/include/python2.5
PYTHONLIB := /usr/lib/python2.5

RINC := /usr/share/R/include

OCTINC := /usr/include/octave-3.0.0
OCTLIB := /usr/lib/octave-3.0.0

JAVAINC := /usr/lib/jvm/java-6-sun-1.6.0.16/include

TCLINC := /usr/include/tcl8.5
TCLLIB := /usr/lib

# Stop your modifications here.

export

all: algencan

algencan: algencan-objects
	$(MAKE) -C $(PROBSRC) all install

algencan-%: algencan-objects
	$(MAKE) -C $(INTERSRC)/$* all install

algencan-objects:
	$(MAKE) -C $(ALGSRC)

clean:
	$(MAKE) -C $(ALGSRC)  clean
	$(MAKE) -C $(PROBSRC) clean
	$(foreach i,$(INTERFCS),$(MAKE) -C $(INTERSRC)/$(i) clean;)

distclean: clean
	$(MAKE) -C $(PROBSRC) distclean
	$(foreach i,$(INTERFCS),$(MAKE) -C $(INTERSRC)/$(i) distclean;)

.PHONY: all clean distclean algencan*
