#
# makefile for dynCore 
#

F90 = gfortran
#F90 = ifort
EXE = main

SOURCE  = ../source
INCLUDE = ../include

# compiler flags
#CFLAGS  = -check all
CFLAGS = -fbounds-check -fcheck-array-temporaries -ffpe-trap=invalid,zero,overflow #-Wall -Wextra

# linker flags
LFLAGS  =

vpath %.f90 $(SOURCE) 
vpath %.h   $(INCLUDE)

# Link
$(EXE): $(EXE).o $(OBJECTS) 
	$(F90) $(EXE).o $(OBJECTS) $(LFLAGS) -o $(EXE)

include $(SOURCE)/Makefile.Database
include $(SOURCE)/Makefile.Eos
include $(SOURCE)/Makefile.Driver
include $(SOURCE)/Makefile.Grid
include $(SOURCE)/Makefile.Hydro
include $(SOURCE)/Makefile.RuntimeParameters
include $(SOURCE)/Makefile.Simulation
include $(SOURCE)/Makefile.Io
include $(SOURCE)/Makefile.PhysicalConstants

OBJECTS = $(Database) $(Driver) $(Grid) $(Hydro) \
          $(RuntimeParameters) $(Simulation) $(Io) \
          $(Eos) $(PhysicalConstants) 

# Compile
$(EXE).o : $(OBJECTS)

%.o: %.f90 
	$(F90) -c $(CFLAGS) $<

clean:
	rm $(OBJECTS) $(EXE) $(EXE).o *.mod 
