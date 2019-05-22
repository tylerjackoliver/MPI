# Makefile for compiling on ARCHER only.
# Changed from the original makefile as
# a result of Cray peculiarities

MF = Makefile

FC = ftn

# Select runtime here. Options: PARALLEL2d, serial.

RUNTIME = PARALLEL2d

FFLAGS = -eZ -D PARALLEL2d -O3 -m 2

EXE = edge_to_image

SRC = \
    pgmio.f90 \
    problem_constants.f90 \
    neighbour_indexes.f90 \
	utility_functions.f90 \
    wrappers_2d.f90 \
    main.f90

.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ= $(SRC:.f90=.o)
MOD= $(SRC:.f90=.mod)

.f90.o: 
		$(FC) $(FFLAGS) -c $<

all:	$(EXE)

$(EXE): $(OBJ)
		$(FC) $(FFLAGS) -o ../$@ $(OBJ)

$(OBJ): $(MF)

clean:
		rm -f $(OBJ) $(EXE) $(MOD) core

