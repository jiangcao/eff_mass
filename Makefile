FC90 = gfortran  
CFLAG = -fbacktrace -c -g -ffree-line-length-500 -fcheck='all'
CCFLAG= -g -fbacktrace -fcheck='all'
CLIB = -llapack -lblas 

# Modules directory
MODDIR = compiled


# Source directory
SRCDIR = src

# Search directories
vpath %.f90 $(SRCDIR)
vpath %.o $(MODDIR)


all: main

main: main.f90 types.o constants.o input.o eff_mass.o 
	    $(FC90) -o $@ $< $(MODDIR)/*.o -I$(MODDIR) $(CCFLAG) $(CLIB)


.PHONY : clean;

clean :
	rm $(MODDIR)/*.o $(MODDIR)/*.mod

# implicit rules

%.o : %.f90
	$(FC90) -o $(MODDIR)/$@ $< -c $(CFLAG) -I$(MODDIR) -J$(MODDIR)
