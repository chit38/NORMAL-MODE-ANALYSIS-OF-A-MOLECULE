PROGRAM main
 USE kinds
 USE system
 USE function_derivative 
 USE bfgs_mod
 USE backtracking
 USE coordinates
!
 IMPLICIT NONE

 REAL(KIND=dp), POINTER :: x(:), g(:), b(:,:), xold(:), gold(:)

 REAL(KIND=dp) :: energy, gnorm, dum, energy_old
 INTEGER :: i, j, iopt
 LOGICAL :: converged

 REAL(KIND=dp), EXTERNAL :: DDOT

!Paramters
 REAL(KIND=dp), PARAMETER :: gnorm_cut=5.E-5, step_max=0.5_dp
 INTEGER, PARAMETER :: maxiopt=100


  OPEN(1,FILE="input.xyz")
  READ(1,*)n  !Reading number of atoms from the input; saved in system module
  REWIND(1) !rewinding as the full XYZ file is read later by read_xyz routine

  PRINT "(A,I5)", "  Number of Atoms =", n

  ALLOCATE(al(n))   !system module


  ALLOCATE(x(3*n),xold(3*n))
  ALLOCATE(g(3*n),b(3*n,3*n),gold(3*n))

  !read input coordinates in XYZ format
  CALL read_xyz(1,x)
  CLOSE(1)
  PRINT "(A,I5)", "  Input coordinates read from input.xyz"


  !Initialize Hessian Inverse
  b(:,:)=0._dp
  DO i=1,3*n
    b(i,i)=1._dp
  END DO


  !open log files
  OPEN(1,FILE="optimize.xyz")
  OPEN(2,FILE="opt.log")

  !starting the counter
  iopt=0

!Quasi-NR loop starts here
  opt_loop: DO 
    IF(iopt>maxiopt)EXIT opt_loop

    energy=func(x)
    CALL grad(x,g)

!   BFGS update of inverse Hessian
    IF(iopt>0) CALL bfgs(x,xold,g,gold,b)

!   Find the magnitude of gradient vector
    gnorm=DSQRT(DDOT(3*n,g,1,g,1))/DFLOAT(n)   !normalized with number of atoms

    PRINT *, "OPT| ", iopt, energy,gnorm       !print some info on screen
    WRITE(2,"(i8,2f16.6)") iopt, energy, gnorm !pring some info in opt.log file

    IF(gnorm<gnorm_cut)THEN
      PRINT *, " "
      PRINT *, "*** Optimization Completed ***" 
      EXIT opt_loop
    END IF

    !setting the search direction for linesearch
    xold(:)=x(:)
    gold(:)=g(:)

    energy_old=energy
    !Linesearch
    CALL line_search_backtrack(xold,energy_old,g,b,step_max,x,converged)

    !Write coordinates of atoms in XYZ format (append)
    CALL write_xyz(1,iopt,x)

    iopt=iopt+1 !Increment optimization number
  END DO opt_loop
!Quasi-NR loop ends here
  CLOSE(1)
  CLOSE(2)

END PROGRAM main
