SHELL=/bin/sh

# To modify this makefile for another code, you should only need to
#  change the entries in the next section and possibly the compiler
#  and options.
# Comments lines begin with "#".

# The name of this makefile
MAKEFILE= Makefile

# The command you type to run the program (executable name)
#COMMAND=  MatrixElement1.x
COMMAND_MAIN = main.x

# Put here all of the C++ or C or Fortran or ... files to be compiled
#  If you get a "missing separator" error pointed to a line here, 
#  make sure that each \ has NO spaces following it.


PATHMAIN=./
_OBJSMAIN= \
   ./mpi_modules.o  ./metropolis.o ./pre_deut_wave.o ./deut_wave.o  ./operations.o  ./main.o
OBJSMAIN = $(patsubst %,$(PATHMAIN)/%,$(_OBJSMAIN))

# Every entry here matches one in SRCS, but with a ".o" ending 
#  Uncomment one of the OBJS= lines according to your preference
#  Put the .o files in the same directory as the corresponding source file
#OBJS= $(addsuffix .o, $(basename $(SRCS)))
#OBJS2= $(addsuffix .o, $(basename $(SRCS2)))
#  Put the .o files in the makefile directory
#OBJS= $(addsuffix .o, $(notdir $(SRCS)))

# Put here any local header files used
HDRS= \


###########################################################################
# Commads and options for different compilers
COMPILER=GNU

DEBUG=FALSE
#
# Compiler parameters
#
# CXX           Name of the C++ compiler to use
# CFLAGS	Flags to the C++ compiler
# F90		Name of the fortran compiler to use 
# FFLAGS	Flags to the fortran compiler 
# LDFLAGS	Flags to the loader
# LIBS		A list of libraries 
#

# Gnu compiler options
ifeq ($(COMPILER),GNU)
  CXX= g++
  F90= mpif90 -mcmodel=medium
  OMP_FLAG=-fopenmp
  LIBS=     -llapack -lblas -lgsl

endif

# Intel compiler options
ifeq ($(COMPILER),INTEL)
  CXX= icpc
  F90= ifort -mcmodel=medium -shared-intel
  OMP_FLAG=-qopenmp
endif

#cray compiler options
ifeq ($(COMPILER),CRAY)
  CXX= cc
  F90= ftn -mcmodel=medium
  OMP_FLAG=-fopenmp
endif	

ifeq ($(DEBUG),TRUE)
  FFLAGS= -O0 -traceback -warn all -g
endif

ifeq ($(DEBUG),FALSE)
  FFLAGS= -O3 -fallow-argument-mismatch
endif

LIBANG=-I/home/agnech/HH/lib -L/home/agnech/HH/lib  -lang -lint -lfun

GSL= -lgsl
###########################################################################
# Instructions to compile and link
########################################################################### 
all:	$(COMMAND_MAIN)

.SUFFIXES:
.SUFFIXES: .o .mod .f90 .f .c

%.o:	%.mod 

# This is the command to link all of the object files together 
$(COMMAND_MAIN): $(OBJSMAIN) $(MAKEFILE) 
	$(F90) $(FFLAGS) $(OMP_FLAG) -o $(COMMAND_MAIN) $(OBJSMAIN) $(LIBANG) $(GSL)

.f90.mod:
	$(F90) -c $(FFLAGS) $(OMP_FLAG) -o $@ $< $(LIBANG)
 
.f90.o:	
	$(F90) -c $(FFLAGS) $(OMP_FLAG) -o $@ $< $(LIBANG)
 
.f.o:	
	$(F90) -c $(FFLAGS) $(OMP_FLAG) -o $@ $< $(LIBANG)

.f.mod:
	$(F90) -c $(FFLAGS) $(OMP_FLAG) -o $@ $< $(LIBANG) 

.c.o:
	$(CXX) -c $(CFLAGS) -o $@ $<

##########################################################################
# Additional tasks      
##########################################################################
      
# Delete the program and the object files (and any module files)
clean:
	/bin/rm -f $(COMMAND_MAIN) $(OBJSMAIN)
	/bin/rm -f $(OBJS1) $(OBJS2) 
	/bin/rm -f *~ *.mod
 
# Pack up the code in a compressed gnu tar file 
#tarz:
#	tar cfvz $(TARFILE) $(README) $(MAKEFILE) $(SRCS)   
 
