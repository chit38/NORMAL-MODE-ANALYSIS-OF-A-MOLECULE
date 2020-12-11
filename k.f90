program project3
use kinds
use system
use function_derivative
use qr
use molden_output
use nmd
implicit none

integer :: i, j, k
real(KIND = dp), pointer :: x(:), H(:,:), gplus(:), gminus(:), xdum(:), A(:,:), q(:,:), r(:,:), u(:,:), frq(:)
real(KIND = dp), pointer :: xc(:), yc(:), zc(:), mass(:)
character(len=2) :: at

open(3, file = "optimize.xyz")
read(3,*)n
rewind(3)
print *, "Number of atoms: ", n
allocate(x(3*n), H(3*n,3*n), gplus(3*n), gminus(3*n), xdum(3*n), A(3*n,3*n), q(3*n,3*n), r(3*n,3*n), u(3*n,3*n), frq(3*n))

do i = 1, 36
    read(3,*)
    read(3,*)
    do j = 1, n
        read(3,*)at, x(3*(j-1)+1:3*j)
    end do
end do
close(3)

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

print *, "--------------------------------------------------------------------------------------------------------"

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

print *, "--------------------------------------------------------------------------------------------"


do i = 1, 3*n
    if(A(i,i) .lt. 0.d0)then
        A(i,i) = -1*A(i,i)
        frq(i) = -1*dsqrt(A(i,i))/(2*pi)
    else
        frq(i) = dsqrt(A(i,i))/(2*pi)
    end if
end do

frq(:) = frq(:)/(3.4E-8)

print *, "Frequencies: "
do i = 1, 3*n
print "(2ES16.6)", frq(i)
end do
!---------------------------------------------------------------------------------------------------

allocate(al(n), xc(n), yc(n), zc(n), mass(n))

al(1:n) = "Ar"
do i = 1, n
    xc(i) = x(3*i - 2)
    yc(i) = x(3*i -1)
    zc(i) = x(3*i)
end do
mass(1:n) = 39.948d0

call create_molden_vib_output(frq,u,al,xc,yc,zc,mass)

call print_nmd(al, x, frq, u)

end program project3