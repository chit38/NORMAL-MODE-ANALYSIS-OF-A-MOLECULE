MODULE qr

CONTAINS 
  SUBROUTINE qrd(a,q,r)
    IMPLICIT NONE
    REAL(KIND=8), POINTER, DIMENSION(:,:) :: a,q,r
 
    INTEGER :: n, k, i
    REAL(KIND=8), ALLOCATABLE :: v(:),d(:), pr(:,:), tmp(:,:)
    REAL(KIND=8):: dnorm, p 


    n=SIZE(a,1)

    ALLOCATE(v(n))
    ALLOCATE(d(n))
    ALLOCATE(pr(n,n))
    ALLOCATE(tmp(n,n))

! Initialize q
    q(:,:)=0.d0
    DO k=1,n
     q(k,k)=1.d0
    END DO
! Initialize r
    r(:,:)=a(:,:)

    DO k=1,n-1
!Copying k th column of the A matrix
      d(1:n)=r(1:n,k)
!Compute norm of column k of A matrix
      dnorm=DSQRT(DOT_PRODUCT(d(1:n),d(1:n)))
!Compute "d" matrix
      d(:)=d(:)/dnorm

      dnorm=DSQRT(DOT_PRODUCT(d(k:n),d(k:n)) )
      IF(d(k) >0.d0) dnorm=-dnorm

!Compute v_j matrix elements for j=1,..k-1 (setting them to zero)
      v(1:k-1)=0.d0
!Compute v_k element of the v matrix 
      v(k)=DSQRT(0.5d0*(1.d0-d(k)/dnorm))

!Compute p
      p=-dnorm*v(k)

!Compute v_j matrix elements for j=k+1,...,n
      v(k+1:n)=d(k+1:n)/2.d0/p

      pr(:,:)=0.d0
      DO i=1,n
       pr(i,i)=1.d0
      END DO

!     Compute Pk
      CALL compute_pk(n,v,pr)
!
      DO i=1,n
       !PRINT *, "P=", pr(i,1:n)
      END DO

!    update R
      tmp(:,:)=r(:,:)
      r=MATMUL(pr,tmp)
!
      DO i=1,n
       !PRINT *, "R=", r(i,1:n)
      END DO
     
!     update Q
      tmp(:,:)=q(:,:)
      q=MATMUL(pr,q)
!
      DO i=1,n
       !PRINT *, "Q^T=", q(i,1:n)
      END DO
    END DO

   !Q^T is now copied to Q
    tmp(:,:)=q(:,:)
    DO i=1,n
      DO k=1,n
        q(i,k)=tmp(k,i)
      END DO
    END DO

    DEALLOCATE(v)
    DEALLOCATE(d)
    DEALLOCATE(pr)


  END SUBROUTINE qrd

  SUBROUTINE compute_pk(n,v,pk)
!
! Compute:
! P_k = I - 2 VV^T
!
    IMPLICIT NONE
    INTEGER :: n
    REAL*8 :: v(n), pk(n,n)
    REAL*8, ALLOCATABLE :: tmp(:,:), tmpt(:,:)
    INTEGER :: i

    ALLOCATE(tmp(n,1),tmpt(1,n))
    tmp(:,1)=v(:)
    tmpt=TRANSPOSE(tmp) 
    pk=-2.d0*MATMUL(tmp,tmpt)
    DO i=1,n
     pk(i,i)=1.d0+pk(i,i)
    END DO 
    DEALLOCATE(tmp,tmpt)
  END SUBROUTINE compute_pk

END MODULE qr
