      SUBROUTINE PSENDRS(BUF)
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
C
      REAL, INTENT(IN) :: BUF
C
      CALL MPI_SEND(BUF,COUNT,P_REAL,DEST,TYPE,
     &                 MPI_COMM_WORLD,IERROR)
C
      RETURN
      END