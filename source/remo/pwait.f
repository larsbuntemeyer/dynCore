      SUBROUTINE PWAIT
      !
      IMPLICIT NONE
      !
      INCLUDE "parorg.h"
      !
      ! Local Variables
      !
      INTEGER :: I
      !
      !
      DO I=2,REQUESTNR
         CALL MPI_WAIT(REQUEST(I-1),STATUS,IERROR)
      ENDDO
      !
      DO I=2,REQUESTNRR
         CALL MPI_WAIT(REQUESTR(I-1),STATUS,IERROR)
      ENDDO
      !
      REQUESTNR=1
      REQUESTNRR=1
      !
      RETURN
      !
      END SUBROUTINE PWAIT
