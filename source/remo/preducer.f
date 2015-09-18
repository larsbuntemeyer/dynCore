      SUBROUTINE PREDUCER(BUF,OP)
      IMPLICIT NONE
      INCLUDE "parorg.h"
      INTEGER, INTENT(IN)    :: OP
      REAL,    INTENT(INOUT) :: BUF(COUNT)
      REAL    :: BUF2(COUNT)
      INTEGER :: I,ROOT
C
      ROOT=0
      CALL MPI_REDUCE(BUF,BUF2,COUNT,P_REAL,
     &     OP,ROOT,MPI_COMM_WORLD,IERROR)
C
      IF (MYID .EQ. 0) THEN
         DO I =1,COUNT
            BUF(I)=BUF2(I)
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE PREDUCER