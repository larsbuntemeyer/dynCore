      SUBROUTINE PRECVRA(BUF)
      IMPLICIT NONE
      INCLUDE "parorg.h"
      REAL, INTENT(IN) :: BUF
C
      CALL MPI_IRECV(BUF,COUNT,P_REAL,SOURCE,TYPE,
     &     MPI_COMM_WORLD,REQUESTR(REQUESTNRR),IERROR)
C
      REQUESTNRR=REQUESTNRR+1
      RETURN
      END SUBROUTINE PRECVRA
