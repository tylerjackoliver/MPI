import sys # Provides argv
import cv2 # Provides imread
import numpy # Provides arrays

print("Checking average of " + sys.argv[1] + ' and ' + sys.argv[2] + '\n')

# Open the first file; get mean value

first_file = cv2.imread(sys.argv[1], -1)

# Second file

second_file = cv2.imread(sys.argv[2], -1)

# Compute the differences in the means and print

diff = abs(numpy.mean(first_file) - numpy.mean(second_file))

# Compute target amount: 5% of average value

tol = numpy.mean(first_file) * 0.05

try:

    assert(numpy.allclose(diff, 0, rtol=tol, atol=tol))
    print ("Test successful!")
  

except AssertionError:

    print("Test failed!")

