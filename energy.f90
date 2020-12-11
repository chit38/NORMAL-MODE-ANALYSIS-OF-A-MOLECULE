MODULE function_derivative

CONTAINS
!This is a function which calculates total energy
  FUNCTION func(x) RESULT(f)
    USE kinds
    USE system
    USE distance

    IMPLICIT NONE
    REAL(KIND=dp) :: f
    REAL(KIND=dp), DIMENSION(:), POINTER :: x
    real(KIND = dp) :: d, d2, d6, d12
    integer :: i, j
    f = 0.d0
    do i = 1, n-1
      do j = i+1, n
        d = x((i-1)*3+1:3*i) .distance. x((j-1)*3+1:3*j)
        d2=1.d0/(d**2)
        d6=d2*d2*d2
        d12=d6*d6
        f=f+d12-d6
      end do
    end do      
    f = f*4.d0
  END FUNCTION func

!This is a subroutine which calculates gradient vector "g"
!x: Input Coordinates. It should be an array of rank 1 of size 3N
!g: Output Gradients.  It is an array of rank 1 of size 3N
  SUBROUTINE grad(x,g) 
    USE kinds
    USE system
    USE distance
    IMPLICIT NONE
    REAL(KIND=dp), POINTER, DIMENSION(:) :: g, x
    real(kind=dp) :: d, d2, d6, d12, dum(3)
    integer :: i, j

    g(1:3*n) = 0.d0
!WRITE IT YOURSELF
    do i = 1, n
      do j = 1, n
        if(i .ne. j) then
        d = x((i-1)*3+1:3*i) .distance. x((j-1)*3+1:3*j)
        dum(1:3) = x((i-1)*3+1:3*i) - x((j-1)*3+1:3*j)
        d2=1.d0/(d**2)
        d6=d2*d2*d2
        d12=d6*d6
        g((i-1)*3+1:3*i) = g((i-1)*3+1:3*i) + (-2.D0*d12+d6)*dum(:)*d2
        dum(:) = 0
        end if
      end do
    end do
  END SUBROUTINE grad

END MODULE