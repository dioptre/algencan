BIN ?= $(CURDIR)
DIR := $(notdir $(CURDIR))

LIBNAME := libAlgencan.so

BASE := tango/algencan/base
UTIL := tango/algencan/util

all: $(LIBNAME) algencanma.jar algencan.jar

$(LIBNAME): javawrapper.o \
          $(ALGSRC)/algencan.o $(ALGSRC)/auglag.o   $(ALGSRC)/backtrack.o \
          $(ALGSRC)/betra.o    $(ALGSRC)/calc.o     $(ALGSRC)/cgm.o \
          $(ALGSRC)/checkder.o $(ALGSRC)/dogleg.o   $(ALGSRC)/drand.o \
          $(ALGSRC)/extrap.o   $(ALGSRC)/evalal.o   $(ALGSRC)/fparam.o \
          $(ALGSRC)/gencan.o   $(ALGSRC)/happrox.o  $(ALGSRC)/lss.o \
          $(ALGSRC)/moresor.o  $(ALGSRC)/newtd.o    $(ALGSRC)/newtonkkt.o \
          $(ALGSRC)/scl.o      $(ALGSRC)/setpoint.o $(ALGSRC)/sevalus.o \
          $(ALGSRC)/spgls.o    $(ALGSRC)/tevalus.o  $(ALGSRC)/tnls.o \
          $(ALGSRC)/uevalus.o  $(ALGSRC)/vevalus.o
	$(FC) $(LDFLAGS) -rdynamic -Wl,-soname,$(LIBNAME) $^ -o $(LIBNAME)

algencanma.jar: Manifest.txt Algencanma.class \
              $(PROBNAME)Characteristics.class \
              $(PROBNAME)Combined.class \
              $(PROBNAME)Constr.class \
              $(PROBNAME)ObjFunc.class
	jar cfm algencanma.jar $^

algencan.jar: $(BASE)/ALGENCANException.class \
              $(BASE)/ALGENCANProblem.class \
              $(BASE)/Characteristics.class \
              $(BASE)/Combined.class \
              $(BASE)/Constraints.class \
              $(BASE)/GradientAndJacobian.class \
              $(BASE)/HLTimesVector.class \
              $(BASE)/Hessian.class \
              $(BASE)/Jacobian.class \
              $(BASE)/NotImplementedException.class \
              $(BASE)/ObjectiveAndConstraints.class \
              $(BASE)/ObjectiveFunction.class \
              $(BASE)/Status.class \
              $(UTIL)/ALGENCANHessian.class \
              $(UTIL)/ALGENCANJacobian.class \
              $(UTIL)/FC.class \
              $(UTIL)/GJAC.class \
              $(UTIL)/HP.class \
              $(UTIL)/SparseStructure.class
	jar cf algencan.jar $^

%.class: %.java
	javac $<

javawrapper.o: javawrapper.c javawrapper.h
	$(CC) $(CFLAGS) -fPIC -c -I$(JAVAINC) -I$(JAVAINC)/linux javawrapper.c

install:
	mkdir -p             $(BIN)/$(DIR)
	cp -f algencanma.jar $(BIN)/$(DIR)
	cp -f algencan.jar   $(BIN)/$(DIR)
	cp -f $(LIBNAME)     $(BIN)/$(DIR)

clean:
	rm -f algencanma.jar
	rm -f algencan.jar
	rm -f $(LIBNAME)
	rm -f javawrapper.o
	rm -f         *.class 
	rm -f $(BASE)/*.class
	rm -f $(UTIL)/*.class

distclean: clean
	rm -f $(BIN)/$(DIR)/algencanma.jar
	rm -f $(BIN)/$(DIR)/algencan.jar
	rm -f $(BIN)/$(DIR)/$(LIBNAME)

.PHONY: all install clean distclean
