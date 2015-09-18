      SUBROUTINE CLOSEZR
      !
      IMPLICIT NONE
      !
      INCLUDE "org.h"
      INCLUDE "corg.h"
      INCLUDE "unitzr.h"
      !
      ! Local Variables
      !
      INTEGER :: I
      !
      DO I=1,INUMZRM
         CLOSE(NUEDAT(I))
      ENDDO
      !
      IF (.NOT. LTAMIT) RETURN
      !
      DO I=1,INUMZRN
         CLOSE(NUNDAT(I))
      ENDDO
      !
      RETURN
      !
      END SUBROUTINE CLOSEZR
