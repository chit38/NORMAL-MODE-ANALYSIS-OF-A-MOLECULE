MODULE system
  INTEGER :: n                         !Number of atoms
  SAVE :: n                            !

  CHARACTER(LEN=2), POINTER :: al(:)   !Array pointer for storing Atom Labels 
  SAVE :: al                           !Store the values of "al" arrays in the memory

  REAL*8, PARAMETER :: pi=4.d0*ATAN(1.d0)
END MODULE system
