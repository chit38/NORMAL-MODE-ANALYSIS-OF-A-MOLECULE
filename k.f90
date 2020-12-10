program project3
use kinds
use system
use function_derivative
use qr
implicit none

integer :: i, j, k
real(KIND = dp), pointer :: x(:), H(:,:), gplus(:), gminus(:), xdum(:), A(:,:), q(:,:), r(:,:), u(:,:)
character(len=2) :: at

open(3, file = "optimize.xyz")
read(3,*)n
rewind(3)
allocate(x(3*n), H(3*n,3*n), gplus(3*n), gminus(3*n), xdum(3*n), A(3*n,3*n), q(3*n,3*n), r(3*n,3*n), u(3*n,3*n))

do i = 1, 36
    read(3,*)
    read(3,*)
    do j = 1, n
        read(3,*)at, x(3*(j-1)+1:3*j)
    end do
end do

H(:,:) = 0.0d0
xdum(:) = x(:)
do i = 1, 3*n
    do j = 1, 3*n
        xdum(i) = x(i) + 0.01d0
        call grad(xdum, gplus)
        xdum(i) = x(i) - 0.01d0
        call grad(xdum, gminus)
        H(i,j) = (gplus(j)-gminus(j))/(2*0.01d0)
        xdum(:) = x(:)
    end do
end do

print *,"Hessian Matrix = "
do i = 1, 3*n
    print "(18F10.6)", H(i,1:3*n) 
end do

A(:,:) = H(:,:)
DO k=1,20
    CALL qrd(a,q,r)
    a=MATMUL(r,q)
    IF(k>1)THEN
       r(:,:)=u(:,:)
       u=MATMUL(r,q) 
    ELSE
      u(:,:)=q(:,:)
    END IF
END DO

PRINT *, "********************eigenvalues*****************"
PRINT "(18F10.6)", (A(k,k), k=1,3*n)
PRINT *, "********************eigenvectors*****************"
do i=1,3*n
    PRINT "(18F10.6)", u(i,1:3*n)
end do


end program project3