MA27  := $(shell if [ -e $(HSLSRC)/ma27ad.f ]; then echo true; fi)
MC30  := $(shell if [ -e $(HSLSRC)/mc30ad.f ]; then echo true; fi)
MA57  := $(shell if [ -e $(HSLSRC)/ma57ad.f ]; then echo true; fi)
MC77  := $(shell if [ -e $(HSLSRC)/mc77ad.f ]; then echo true; fi)

all: algencanma.o algencan.o auglag.o backtrack.o betra.o calc.o cgm.o \
     checkder.o dogleg.o drand.o evalal.o extrap.o fparam.o gencan.o \
     happrox.o lss.o moresor.o newtd.o newtonkkt.o scl.o setpoint.o \
     sevalus.o spgls.o tevalus.o tnls.o uevalus.o vevalus.o

algencanma.o: algencanma.f dim.par

algencan.o: algencan.f algconst.par algparam.com counters.com dim.par \
            fixvar.com machconst.com outtyp.com scaling.com slacks.com

auglag.o: auglag.f algconst.par algparam.com dim.par graddat.com machconst.com \
          outtyp.com scaling.com

backtrack.o: backtrack.f algconst.par counters.com dim.par machconst.com \
             outtyp.com

betra.o: betra.f algconst.par counters.com dim.par hessdat.com itetyp.com \
         machconst.com outtyp.com sydat.com

calc.o: calc.f dim.par rspace.com

cgm.o: cgm.f algconst.par algparam.com dim.par machconst.com outtyp.com

checkder.o: checkder.f algparam.com dim.par machconst.com

dogleg.o: dogleg.f dim.par machconst.com outtyp.com

drand.o: drand.f

evalal.o: evalal.f algparam.com dim.par graddat.com machconst.com

extrap.o: extrap.f counters.com dim.par machconst.com outtyp.com

fparam.o: fparam.f algparam.com dim.par fixvar.com outtyp.com scaling.com \
          slacks.com

gencan.o: gencan.f algconst.par algparam.com probdata.com counters.com \
          dim.par itetyp.com machconst.com outtyp.com rspace.com sydat.com

happrox.o: happrox.f algconst.par dim.par graddat.com happdat.com hpredat.com \
           itetyp.com machconst.com sydat.com

moresor.o: moresor.f dim.par machconst.com outtyp.com

newtd.o: newtd.f algparam.com dim.par graddat.com itetyp.com machconst.com \
         outtyp.com rspace.com

newtonkkt.o: newtonkkt.f algparam.com dim.par graddat.com hessdat.com \
             machconst.com outtyp.com

setpoint.o: setpoint.f

sevalus.o: sevalus.f algparam.com probdata.com dim.par outtyp.com scaling.com

spgls.o: spgls.f counters.com dim.par outtyp.com

tevalus.o: tevalus.f algparam.com dim.par outtyp.com slacks.com

tnls.o: tnls.f algconst.par counters.com dim.par machconst.com outtyp.com

uevalus.o: uevalus.f dim.par probdata.com fixvar.com outtyp.com

vevalus.o: vevalus.f algparam.com counters.com dim.par graddat.com hessdat.com \
           machconst.com outtyp.com

ifeq ($(MA57),true)

lss.o: lssma57.f $(HSLSRC)/lma57ad.f dim.par ma57dat.com machconst.com \
       outtyp.com
	$(FC) $(FFLAGS)  -c -o $@ $< -I$(HSLSRC)
scl.o: scl.f

$(HSLSRC)/lma57ad.f: $(HSLSRC)/ma57ad.f
	$(HSLSRC)/fix-ma57.sh $<

else ifeq ($(MA27),true)

lss.o: lssma27.f $(HSLSRC)/lma27ad.f dim.par ma27dat.com machconst.com \
       outtyp.com
	$(FC) $(FFLAGS)  -c -o $@ $< -I$(HSLSRC)

$(HSLSRC)/lma27ad.f: $(HSLSRC)/ma27ad.f
	$(HSLSRC)/fix-ma27.sh $<

  ifeq ($(MC77),true)

scl.o: sclmc77.f dim.par mc77dat.com outtyp.com
	$(FC) $(FFLAGS)  -c -o $@ $< -I$(HSLSRC)

  else ifeq ($(MC30),true)

scl.o: sclmc30.f dim.par mc30dat.com outtyp.com
	$(FC) $(FFLAGS)  -c -o $@ $< -I$(HSLSRC)

  else

scl.o: scl.f

  endif

else

lss.o: lss.f
scl.o: scl.f

endif

clean:
	rm -f *.o
	rm -f $(HSLSRC)/lma27ad.f $(HSLSRC)/lma57ad.f

.PHONY: all clean
