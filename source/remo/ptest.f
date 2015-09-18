      SUBROUTINE PTEST
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
C
      LOGICAL :: FLAG
      INTEGER :: I
C
      IF (TAGCOUNT.EQ.0) THEN
         CALL MPI_IPROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,
     &        FLAG,STATUS,IERROR)
         IF (FLAG .EQV. .TRUE.) THEN
            SOURCE = STATUS(MPI_SOURCE)
            TYPE   = STATUS(MPI_TAG)

            IF (TYPE .EQ. 2) THEN
               CALL MPI_GET_COUNT(STATUS,MPI_INTEGER,COUNT,IERROR)
            ELSEIF ((TYPE .GT. 500) .AND. (TYPE .LT. 510)) THEN
               CALL MPI_GET_COUNT(STATUS,MPI_LOGICAL,COUNT,IERROR)
            ELSE
               CALL MPI_GET_COUNT(STATUS,P_REAL,COUNT,IERROR)
            ENDIF

         ELSE
            TYPE   = -1
         ENDIF
      ELSE
         TYPE = -1
         DO WHILE(TYPE.EQ.-1)
            DO I=1,TAGCOUNT
               CALL MPI_IPROBE(MPI_ANY_SOURCE,TAGTABLE(I),
     &              MPI_COMM_WORLD,FLAG,STATUS,IERROR)
               IF (FLAG .EQV. .TRUE.) THEN
                  SOURCE = STATUS(MPI_SOURCE)
                  TYPE   = STATUS(MPI_TAG)

                  IF (TYPE .EQ. 2) THEN
                     CALL MPI_GET_COUNT(STATUS,MPI_INTEGER,COUNT,IERROR)
                  ELSEIF ((TYPE .GT. 500) .AND. (TYPE .LT. 510)) THEN
                     CALL MPI_GET_COUNT(STATUS,MPI_LOGICAL,COUNT,IERROR)
                  ELSE
                     CALL MPI_GET_COUNT(STATUS,P_REAL,COUNT,IERROR)
                  ENDIF

                  EXIT
               ELSE
                  TYPE   = -1
               ENDIF
            ENDDO
         ENDDO
      ENDIF
C     
      RETURN
      END SUBROUTINE PTEST
