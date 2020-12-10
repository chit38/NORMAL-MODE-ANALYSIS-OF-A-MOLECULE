MODULE distance
USE kinds 

  INTERFACE OPERATOR(.distance.)
    MODULE PROCEDURE calculate_distance
  END INTERFACE

CONTAINS
   FUNCTION calculate_distance(v1,v2) RESULT(d)
     IMPLICIT NONE
     REAL (KIND=dp) :: d
     REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: v1, v2
  
     REAL(KIND=dp) :: buffer(3)

     buffer(1:3) = v1(1:3)-v2(1:3)
     d=SQRT(DOT_PRODUCT(buffer,buffer))     
   END FUNCTION calculate_distance 

END MODULE distance
