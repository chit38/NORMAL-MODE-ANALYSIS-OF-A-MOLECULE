#1. You need to install blas and lapack in your computer
#https://coral.ise.lehigh.edu/jild13/2016/07/27/install-lapack-and-blas-on-linux-based-systems/

#Compile the BFGS Code
gfortran -o bfgs.x kinds.f90 system.f90 distance.f90 energy.f90 bfgs.f90 backtracking.f90 coordinates.f90 main.f90 -lblas

