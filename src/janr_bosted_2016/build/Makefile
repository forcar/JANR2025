.PHONY: all clean

CERN_ROOT = /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el7/cernlib/2005

FC = gfortran
#FCFLAGS = -std=legacy -ffixed-line-length-none -fbounds-check  -fno-automatic
#FCFLAGS = -g -ffixed-line-length-none -fbounds-check  -fno-automatic

FCFLAGS =  -O -fno-automatic -fallow-argument-mismatch  -ffixed-line-length-none -fno-second-underscore -funroll-loops -fomit-frame-pointer -fPIC 
#FCFLAGS =  -g -fbounds-check  -fno-automatic
#FCFLAGS =  -g -fbounds-check  -fno-automatic -ffixed-line-length-none
#FCFLAGS = -ffixed-line-length-none -fbounds-check  -fno-automatic

SRCS = $(wildcard *.F)
OBJS_TEMP = $(patsubst %.F,%.o,$(SRCS))
OBJS = $(filter-out janr_main.o,$(OBJS_TEMP))
OBJS += abort_rut.o

all: janr 

lib: $(OBJS)
		ar rv libjanr.a $(OBJS) 

clean:
		rm -fv $(OBJS) janr

%.o: %.F
		$(FC) $(FCFLAGS) -c $^ 

%.o: %.c
		gcc -c $^

janr:
		$(FC) janr_bosted_2016.F -L./ -ljanr  -L$(CERN_ROOT)/lib -lpawlib -lgraflib -lgrafX11 -lpacklib -lphtools -lmathlib -lkernlib -ljetset74 -L/opt/X11/lib/ -lX11 -o $@

display:
	  @echo "$(SRCS)"	
		@echo "$(OBJS)"
