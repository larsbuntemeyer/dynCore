      SUBROUTINE PSENDALLR(BUF)
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
C
      REAL, INTENT(IN) :: BUF(*)
      INTEGER :: ROOT
C
      ROOT=0
      CALL MPI_BCAST(BUF,COUNT,P_REAL,ROOT,MPI_COMM_WORLD,
     &      IERROR)
      RETURN
      END SUBROUTINE PSENDALLR