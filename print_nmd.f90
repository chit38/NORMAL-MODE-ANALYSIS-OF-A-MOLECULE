MODULE nmd
CONTAINS
SUBROUTINE print_nmd(al,x,eigen_values,eigen_vectors)
! al => atom labels (size equal to the number of atoms)
! x  => atom coordinates (size equal to 3 x number of atoms)

 IMPLICIT NONE
 CHARACTER(LEN=2), POINTER :: al(:)
 REAL*8, POINTER :: x(:), eigen_values(:), eigen_vectors(:,:)

 INTEGER :: nn, n, i,k
 CHARACTER(LEN=25) :: myfmt, myfmt2

 nn=SIZE(eigen_values)
 n=SIZE(al)

! Get printing format
  IF(n<10)THEN
    WRITE(myfmt,"(A4,i1,A3)") "(A5,",n,"A4)"
    WRITE(myfmt2,"(A5,i2,A6)") "(A11,",3*n,"F12.6)"
  ELSE IF(n<100)THEN
    WRITE(myfmt,"(A4,i2,A3)") "(A5,",n,"A4)"  
    WRITE(myfmt2,"(A5,i3,A6)") "(A11,",3*n,"F12.6)" 
  ELSE
!   This may not work
    myfmt="*"
    myfmt2="*"
  END IF
  print*,myfmt2

  open(32, file = "vib.nmd")
  WRITE(32,myfmt) "names", (al(i),i=1,n)

  WRITE(32,myfmt2)"coordinates", (x(i),i=1,3*n)


! Get printing format
  IF(nn+1<10)THEN
    WRITE(myfmt,"(A10,i1,A6)") "(a4,2x,i2,",nn+1,"F10.3)"
  ELSE IF(nn+1<100)THEN
    WRITE(myfmt,"(A10,i2,A6)") "(a4,2x,i2,",nn+1,"F15.3)"
  ELSE
!   This may not work
    myfmt="*"
  END IF

  DO i = 1, nn
     WRITE(32,myfmt) "mode", i, eigen_values(i), &
                               (eigen_vectors(k,i),k=1,nn)
  ENDDO 
  CLOSE(32)

END SUBROUTINE print_nmd
END MODULE nmd
