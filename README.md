# inviscidFluxSchemes
This code tests various inviscid flux schemes used to solve the compressible Euler equations
The test runs the 1D Sod shock tube

# to compile
go to "bld" folder and type make
you will need a working gfortran compiler

# to run the shocktube test
copy the executable from bld (stube.exe)
do ./stube.exe
The stube.inp file can be used to change the various schemes,order of accuracy and flux limiters
the exact solution to Sod shock tube is stored in the "exact_soln" folder in the root directory
There is also an stube.html file which is javascript that you can run on a browser
