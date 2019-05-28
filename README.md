# 2-D Image Re-Composition

The program included here performs a two-dimensional image recomposition from a target `.pgm` file, and continues with the re-construction until convergence to within a user-defined tolerance is achieved.

This program has two main modes of operation: a serial run-time, and a parallel run-time. The compilation of the correct mode of operation is achieved using a compiler flag, set in the respective Makefile.

Two make-files are provided with the source code: to compile using mpif90 on a generic machine, use the Makefile in the case\_study/ directory. To compile using ftn on ARCHER, use the Makefile in the case\_study/src/ directory. Different Makefiles were generated to navigate some of the intricacies behind the Cray compiler.

### Serial Version

To compile in serial, edit the Makefile in either directory, such that the RUNTIME flag is set to SERIAL. By doing so, the serial wrapping library (src/serial.f90) will be loaded, and much of the MPI functionality removed. However, the MPI library will still initialise, so as to be able to use the `MPI_WTIME()` function for performance monitoring purposes.

To run the serial version of the file, type the following at your shell of choice:

`./edge\_to\_image`

### Parallel Version

To compile in parallel using a two-dimensional Cartesian topology, edit the Makefile in either directory such that the RUNTIME flag is set to PARALLEL2d. By doing so, MPI functionality will be activated and the MPI wrapping library (`src/wrappers_2d.f90`) will be loaded.

To run the parallel version of the file, it is necessary to specify the number of processors the executable is being run on, as thus:

`mpirun -np ${PROCNUM} ./edge_to_image ${PROCNUM}`

Should the two values of PROCNUM not match, an `ERROR STOP` will be triggered, and the program will promptly exit violently.

## Program Peculiarities

Some issues may be encountered with the Makefile: since the Makefile works by compiling all files found in a target directory to objects and linking only afterward, it is possible (but rare) that the correct module file will not be available at time of compilation. Should that occur, please manually compile the corresponding source code into an object file, before re-issuing the `make` command.

## Contact

If any features are unclear, or you cannot get the code to run or compile, please email me at <jack.tyler@soton.ac.uk>.
