# depends on FORTRAN compilers, command to invloke the fortran compiler
FC = ifort

#please check your compiler for these options
FFLAGS = -O3 
       


2d  =   calangle.o calconstrvector.o calconserforce.o calenergy.o \
         calfps.o calvelocity.o calwallgnsub.o choldc.o choldc_origin.o\
         deltat.o erfc.o farfldres.o forces_random.o netforces.o\
	     gasdev.o grndrm.o invert.o lubksb.o ludcmp.o\
	     main.o ppexct.o ppinf.o pwff.o pwexct.o ran2.o random.o readcn.o\
	     rescale.o writcn_es.o wallgreenfun.o calrodop.o permutator.o getbodyframevector.o \
	   
           $(libs)

.f.o:;  $(FC) $(FFLAGS) -c $<

2dmc: $(2d)
	  $(FC) $(FFLAGS) -o $@ $(2d)


clean:   
	rm -f core* 2dmc *.o *.il

