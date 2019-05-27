#!/bin/bash

# This script runs all of the tests available for the edge reconstruction routines.
# In particular, it checks for the following:
#
#       - Checks that the results from the 2-D cartesian topology agree with the sample
#         solution from the MPI course for all files
#
#       - Checks that the program falls over somewhat gracefully when improper problem
#         constants are used (wrong M, N; P; .)
#
# In the ref_soln folder is the four outputs for the sample solution code; we can run our
# implementation on the same images and check the image is correct to within some tolerance.
#

# First copy the current version of problem_constants.f90 into the working folder

cp ../src/problem_constants.f90 .

## Test 1 - check sample solution for 192 x 168

# Move the correct problem_constants.f90 file (config'd for 192 x 168) into the src/ directory and make

cp prob_files/problem_constants_192x128.f90 ../src/problem_constants.f90

# Move into the case_study folder and remake the solutions

cd ../ && make clean
mpif90 src/utility_functions.f90 -c
make

# Run the file on edge192x168.pgm - writes into output.pgm

mpirun -np 2 ./edge_to_image 2

# Now move this file back into our ref_solns folder

mv output.pgm tests/ref_soln/output.pgm

# Move into the tests folder

cd tests

# Compute the differences between the average values of the two files

python check_average.py ref_soln/ref_192x128.pgm ref_soln/output.pgm

# Wait to read if there has been an assertion error...

sleep 3

# Delete temporary files

rm ref_soln/output.pgm ../src/problem_constants.f90

# ==========================================================================================================

## Test 2 - check sample solution for 256 x 192

# Move the correct problem_constants.f90 file (config'd for 256 x 192) into the src/ directory and make

cp prob_files/problem_constants_256x192.f90 ../src/problem_constants.f90

# Move into the case_study folder and remake the solutions

cd ../ && make clean
mpif90 src/utility_functions.f90 -c
make

# Run the file on edge256x192.pgm - writes into output.pgm

mpirun -np 2 ./edge_to_image 2

# Now move this file back into our ref_solns folder

mv output.pgm tests/ref_soln/output.pgm

# Move into the tests folder

cd tests

# Compute the differences between the average values of the two files

python check_average.py ref_soln/ref_256x192.pgm ref_soln/output.pgm

# Wait to read if there has been an assertion error...

sleep 3

# Delete temporary files

rm ref_soln/output.pgm ../src/problem_constants.f90

# ============================================================================================================

## Test 3 - check sample solution for 512 x 384

cp prob_files/problem_constants_512x384.f90 ../src/problem_constants.f90

# Move into the case_study folder and remake the solutions

cd ../ && make clean
mpif90 src/utility_functions.f90 -c
make

# Run the file on edge512x384.pgm - writes into output.pgm

mpirun -np 2 ./edge_to_image 2

# Now move this file back into our ref_solns folder

mv output.pgm tests/ref_soln/output.pgm

# Move into the tests folder

cd tests

# Compute the differences between the average values of the two files

python check_average.py ref_soln/ref_512x384.pgm ref_soln/output.pgm

# Wait to read if there has been an assertion error...

sleep 3

# Delete temporary files

rm ref_soln/output.pgm ../src/problem_constants.f90

# ============================================================================================================

## Test 4 - check sample solution for 768 x 768

cp prob_files/problem_constants_768x768.f90 ../src/problem_constants.f90

# Move into the case_study folder and remake the solutions

cd ../ && make clean
mpif90 src/utility_functions.f90 -c
make

# Run the file on edge768x768.pgm - writes into output.pgm

mpirun -np 2 ./edge_to_image 2

# Now move this file back into our ref_solns folder

mv output.pgm tests/ref_soln/output.pgm

# Move into the tests folder

cd tests

# Compute the differences between the average values of the two files

python check_average.py ref_soln/ref_768x768.pgm ref_soln/output.pgm

# Wait to read if there has been an assertion error...

sleep 3

# Delete temporary files

rm ref_soln/output.pgm ../src/problem_constants.f90

# =============================================================================================================

echo "=========================================================================================="

## Test the program fails gracefully for an incorrect number of cores

cd .. && mpirun -np 2 ./edge_to_image 4

sleep 3

cd tests/

## Test the program fails gracefully for an incorrect number of cores

cp prob_files/wrong_size.f90 ../src/problem_constants.f90

cd .. && make clean && mpif90 src/utility_functions.f90 -c

make

mpirun -np 2 ./edge_to_image 2

sleep 3

cd tests/

# Copy the old file back

mv problem_constants.f90 ../src/.

