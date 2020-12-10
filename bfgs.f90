MODULE bfgs_mod

  CONTAINS

  SUBROUTINE bfgs(x,xold,g,gold,b)
    USE kinds
    USE system
    IMPLICIT NONE
    
    REAL(KIND=dp), POINTER :: x(:),xold(:), g(:), gold(:), b(:,:)

    REAL(KIND=dp), ALLOCATABLE :: s(:), y(:), scratch(:,:), bb(:,:), by(:)
    REAL(KIND=dp) :: sy,fac
    REAL(KIND=dp), EXTERNAL :: DDOT

!    return
    ALLOCATE(s(3*n),y(3*n),scratch(3*n,3*n),bb(3*n,3*n),by(3*n))

    s(:)=x(:)-xold(:)
    y(:)=g(:)-gold(:)

!copy of B matrix
    bb(:,:)=b(:,:)
!s^Ty
    sy=DDOT(3*n,s,1,y,1)
    sy=1._dp/sy
!sy^T/s^Ty
    CALL DGEMM('N','T',3*n,3*n,1  ,sy     ,s      ,3*n,y,3*n,0._dp,scratch,3*n)
!B_k = B_k- sy^T B_k/s^Ty
    CALL DGEMM('N','N',3*n,3*n,3*n,-1._dp  ,scratch,3*n,bb,3*n,1._dp,b     ,3*n)
!ys^T/s^Ty
    CALL DGEMM('N','T',3*n,3*n,1  ,sy     ,y      ,3*n,s,3*n,0._dp,scratch,3*n)
!B_k = B_k - B_k ys^T/s^Ty
    CALL DGEMM('N','N',3*n,3*n,3*n,-1._dp  ,bb,3*n,scratch,3*n,1._dp,b     ,3*n)
!B_k y/s^Ty
    CALL DGEMV('N',3*n,3*n,sy,bb,3*n,y,1,0._dp,by,1)
!1+y^T B_k y/s^Ty
    fac=DDOT(3*n,y,1,by,1)
    fac=1._dp+fac
!ss^T/s^Ty
    CALL DGEMM('N','T',3*n,3*n,1  ,sy     ,s      ,3*n,s,3*n,0._dp,scratch,3*n)
!Bk+(1+y^T B_k y/S^Ty)ss^T/s^Ty
    CALL DAXPY(3*n*3*n,fac,scratch,1,b,1)
    DEALLOCATE(s,y,bb,scratch,by)
  END SUBROUTINE bfgs
  

END MODULE bfgs_mod
