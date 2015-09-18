      SUBROUTINE SENDSUBGRID(IVLOC,IPDB,UFELD,VFELD)
      IMPLICIT NONE
      
      INCLUDE "parorg.h"

      INTEGER, INTENT(IN)    :: IVLOC
      INTEGER, INTENT(IN)    :: IPDB(37)
      REAL,    INTENT(IN)    :: UFELD(MOIE,MOJE)
      REAL,    INTENT(INOUT) :: VFELD(IE*JE)

      INTEGER, PARAMETER :: ZBW=2     
      REAL    :: SENDFELD(IALLOGI*IALLOGJ+38)
      INTEGER :: IZNEWI,IZNEWJ,I,J,K,IZREBI

C     IZNEWI ,IZNEWJ   = CALCULATED GLOBAL INDICES
C     I,J,K	       = LOOP VARIABLES
C     ZBW	       = BOUNDARY WIDTH
C     IZREBI	       = LOCAL INDEX

      
      DO K=1,NPROC
         IZREBI = 0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ES WERDEN AUCH DIE RANDWERTE GESETZT.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

         DO J=ISUBPOS(K,2)-ZBW,ISUBPOS(K,4)+ZBW
            DO I=ISUBPOS(K,1)-ZBW,ISUBPOS(K,3)+ZBW
               IZREBI        = IZREBI + 1
               IF (I.LE.1) THEN
                  IZNEWI = 1
               ELSE IF (I.GE.MOIE) THEN
                  IZNEWI = MOIE
               ELSE
                  IZNEWI = I
               ENDIF
               IF (J.LE.1) THEN
                  IZNEWJ = 1
               ELSE IF (J.GE.MOJE) THEN
                  IZNEWJ = MOJE
               ELSE
                  IZNEWJ = J
               ENDIF
               IF (K-1.NE.MYID) THEN
                  SENDFELD(IZREBI) = UFELD(IZNEWI,IZNEWJ)
               ELSE
                  VFELD(IZREBI) = UFELD(IZNEWI,IZNEWJ)
               ENDIF
            ENDDO
         ENDDO
         IF (K-1.NE.MYID) THEN
            IZREBI = IZREBI + 1
            DO I=1,37
               SENDFELD(IZREBI)=IPDB(I)
               IZREBI=IZREBI+1
            END DO
            SENDFELD(IZREBI)=IVLOC
            COUNT=IZREBI*1
            DEST=K-1
            TYPE=1
            CALL PSENDRS(SENDFELD(1))
         ENDIF
      ENDDO
      
      RETURN
      END SUBROUTINE SENDSUBGRID
