MODULE coordinates

  CONTAINS

  SUBROUTINE read_xyz(iunit,x)
  USE kinds
  USE system
  IMPLICIT NONE
  REAL(KIND=dp), POINTER :: x(:)
  INTEGER :: iunit
  INTEGER :: i

  READ(1,*) !skipping number of atoms  entry
  READ(1,*) !skipping blank line 
  DO i=1,n
     READ(1,*)al(i), x((i-1)*3+1:3*i)
  END DO
  END SUBROUTINE read_xyz

  SUBROUTINE write_xyz(iunit,iopt,x)
  USE kinds
  USE system
  IMPLICIT NONE
  REAL(KIND=dp), POINTER :: x(:)
  INTEGER :: iunit, iopt
  
  INTEGER :: i

  WRITE(iunit,"(I10)")n
  WRITE(iunit,"(I10)")iopt
  DO i=1,n
    WRITE(iunit,"(A2,3F16.6)")al(i), x((i-1)*3+1:3*i)
  END DO
  END SUBROUTINE write_xyz

END MODULE coordinates
