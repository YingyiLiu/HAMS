# Sources in specific dependency order
SRCS = WavDynMods.f90 PatclVelct.f90 BodyIntgr.f90 BodyIntgr_irr.f90 AssbMatx.f90 AssbMatx_irr.f90 \
SingularIntgr.f90 InfGreen_Appr.f90 FinGrnExtSubs.f90 FinGreen3D.f90 CalGreenFunc.f90 NormalProcess.f90 ReadPanelMesh.f90 \
HydroStatic.f90 ImplementSubs.f90 InputFiles.f90 PotentWavForce.f90 \
PressureElevation.f90 PrintOutput.f90 SolveMotion.f90 WavDynSubs.f90 HAMS_Prog.f90 

# Find all source files, create a list of corresponding object files
#SRCS=$(wildcard *.f90)

OBJS=$(patsubst %.f90,%.o,$(SRCS))

# Ditto for mods (They will be in both lists)
MODS=$(wildcard mod*.f90)
MOD_OBJS=$(patsubst %.f90,%.o,$(MODS))

# Compiler/Linker settings
FC ?= gfortran
ifeq ($(FC), f77)
    FC = gfortran 
endif

ifeq ($(FC), ifort)
	FLFLAGS = -g -xHost -fopenmp -qmkl
	FCFLAGS = -g -xHost -c -O3 -fPIC -check bounds,pointers -warn all -fopenmp -qmkl
	LDCOMP  = -lmkl_rt
else ifeq ($(FC), ifx)
	FLFLAGS = -g -xHost -fopenmp -qmkl
	FCFLAGS = -g -xHost -c -O3 -fPIC -check bounds,pointers -warn all -fopenmp -qmkl
	LDCOMP  = -lmkl_rt
else
	FLFLAGS = -g -march=native -mtune=native -fopenmp
	FCFLAGS = -g -march=native -mtune=native -c -O3 -fPIC -fopenmp -fdec-math -fcheck=bounds,mem,pointer -Wall
	LDCOMP  = -llapack
endif
PROGRAM = hams
PRG_OBJ = $(PROGRAM).o

ifeq ($(OS),Windows_NT)
    LIB := libhams.dll
    LDFLAGS=-shared
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
	LIB := libhams.so
	LDFLAGS=-shared -Wl,-soname,$(LIB)
    endif
    ifeq ($(UNAME_S),Darwin)
	LIB := libhams.so
	LDFLAGS=-dynamiclib
    endif
endif

# make without parameters will make first target found.
default : $(PROGRAM) shared

# Compiler steps for all objects
$(OBJS) : %.o : %.f90
	$(FC) $(FCFLAGS) -o $@ $<

# Linker
$(PROGRAM) : $(OBJS)
	$(FC) $(FLFLAGS) -o $@ $^ $(LDCOMP)

shared : $(OBJS)
	$(FC) $(FLFLAGS) $(LDFLAGS) -o $(LIB) $^ $(LDCOMP)

# If something doesn't work right, have a 'make debug' to 
# show what each variable contains.
debug:
	@echo "FC = $(FC)"
	@echo "SRCS = $(SRCS)"
	@echo "OBJS = $(OBJS)"
	@echo "MODS = $(MODS)"
	@echo "MOD_OBJS = $(MOD_OBJS)"
	@echo "PROGRAM = $(PROGRAM)"
	@echo "SHARED = $(SHARED)"
	@echo "PRG_OBJ = $(PRG_OBJ)"

clean:
	rm -rf $(OBJS) $(PROGRAM) $(LIB) $(patsubst %.o,%.mod,$(MOD_OBJS)) *.mod

.PHONY: debug default clean

# Dependencies

# Main program depends on all modules
$(PRG_OBJ) : $(MOD_OBJS)

# Blocks and allocations depends on shared
mod_blocks.o mod_allocations.o : mod_shared.o
