MODULE molden_output
CONTAINS
      SUBROUTINE create_molden_vib_output(vib_freq,normal_mode,al,xc,yc,zc,mass)
!
!     vib_freq: Eigen Values vibrational frequencies (Array of 3N size)
!     normal_mode: Eigen vectors for the frequency  (Array of 3Nx3N size)
!     al : atom labels  (Array of N size)
!     xc : array of x coordiantes (Array of size N)
!     yc : array of y coordiantes (Array of size N)
!     zc : array of z coordiantes (Array of size N)
!     mass : array of masses of atoms in AMU (array of size N)
!     ==--------------------------------------------------------------==

!
      IMPLICIT NONE
!     Arguments
      REAL*8, POINTER :: vib_freq(:),normal_mode(:,:),xc(:),yc(:),zc(:)
      REAL*8, POINTER :: mass(:),an(:)
      CHARACTER(LEN=2), POINTER :: al(:)

!     ==--------------------------------------------------------------==
      CHARACTER :: typ*2,a*8,ck(3)*7,fr*17,rm(5)*17
      INTEGER ::   i,j,k,l,n1,is,ia,at,he,m
      REAL*8  ::  gz,xcm,ycm,zcm
      INTEGER  :: natoms,nc,nh
!     ==--------------------------------------------------------------==
      
      natoms=SIZE(al)

      ALLOCATE(an(natoms))

      DO i = 1, natoms
         nc = index(al(i),"C")
         nh = index(al(i),"H")

         IF(nc /= 0) THEN
            an(i) = 6.d0
         ELSEIF(nh /= 0) THEN
             an(i) = 1.d0
         ELSE
             an(i)   = 18.d0
         ENDIF

      ENDDO
           
!     write gaussianformated output into vib.log
!     which are readable in molden/molekel to visualise the vibrations.
      fr=' Frequencies --  '
      rm(1)=' Red. masses --  ' 
      rm(2)=' Frc consts  --  ' 
      rm(3)=' IR Inten    --  ' 
      rm(4)=' Raman Activ --  ' 
      rm(5)=' Depolar     --  '
      typ='?A'
      a=' Atom AN'
      ck(1)='      X'
      ck(2)='      Y'
      ck(3)='      Z'
      at=0
      gz=0.d0
! ---- read coordinates and atomtyps ------------
      OPEN(UNIT=15,FILE="vib.log")

! ---- write some lines needed for molden

         WRITE (15,*)'Entering Gaussian System'
         WRITE (15,*)'this file is generated from the MLDGAU',&
                       ' subroutine in the file secder.F'
         WRITE (15,*)'Please note, that this is a "faked" output;'
         WRITE (15,*)'there are no intensities computed in CPMD.'
         WRITE (15,*)'Standard orientation:'
         WRITE (15,*)'---------------------------------------',&
                       '------------------------------'
         WRITE (15,'(A,2(5X,A),14X,A)') &
                   'Center','Atomic','Atomic','Coordinates (Angstroms)'
         WRITE (15,'(2(A,5X),1X,A,3X,3(11X,A))') &
                       'Number','Number','Type','X','Y','Z'
         WRITE (15,*)'---------------------------------------',&
                       '------------------------------'


! ---- make sure that the center of mass is the origin or move the atoms
! ---- from the center of the box to the origin and printout
             
           xcm = SUM(mass*xc)/SUM(mass)
           ycm = SUM(mass*yc)/SUM(mass)
           zcm = SUM(mass*zc)/SUM(mass)
             
         DO i=1,natoms
            xc(i)=xcm - xc(i)
            yc(i)=ycm - yc(i)
            zc(i)=zcm - zc(i)
            WRITE (15,22)  i, INT(an(i)), at, xc(i), yc(i), zc(i)
         ENDDO

! ---- write some lines for molden -----------------------------
         WRITE(15,*)'--------------------------------------------',&
                      '-------------------------'
         WRITE(15,*)'      basis functions          primitive ', &
                      'gaussians'
         WRITE(15,*)'      alpha electrons          beta electrons'
         WRITE(15,*)'********************************************',&
                      '**************************'
         WRITE(15,*)
         WRITE(15,*)'Harmonic frequencies (cm**-1), IR intensities ',&
                      '(KM/Mole),'
         WRITE(15,*)'Raman scattering activities (A**4/AMU), Raman ',&
                      'depolarization ratios,'
         WRITE(15,*)'reduced masses (AMU), force constants ',&
                      '(mDyne/A) and normal coordinates:'
     

! ---- write eigenvalues and eigenvectors in both files
      DO 15 i=1,(3*natoms),3
        WRITE(15,23) i, i+1, i+2
        WRITE(15,24) typ, typ, typ
        WRITE(15,25) fr, (vib_freq(l),l=i,i+2)
        DO n1=1,5
           WRITE(15,25) rm(n1), gz, gz, gz
        ENDDO
        WRITE(15,26) a,(ck(n1),n1=1,3),(ck(n1),n1=1,3),(ck(n1),n1=1,3)
        DO 16 j=1,3*natoms,3
          he=(j-1)/3+1
          WRITE(15,27) he,INT(an(he)),&
                      (normal_mode(j,m),normal_mode(j+1,m),normal_mode(j+2,m),m=i,i+2) 
 16     CONTINUE
 15   CONTINUE
      WRITE(15,*) 'Normal termination of Gaussian 98.'

      CLOSE(15)

 22   FORMAT(i5,i11,i14,4x,3(3x,f9.6))
 23   FORMAT(i22,2i23)
 24   FORMAT(20x,a2,2(21x,a2))
 25   FORMAT(a17,f15.4,2f23.4)
 26   FORMAT(a8,3a7,2(2x,3a7))
 27   FORMAT(2i4,3(f9.2,2f7.2))
      DEALLOCATE(an)
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE create_molden_vib_output  
END MODULE
