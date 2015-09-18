      SUBROUTINE COLLECTDATA( GUFELD, UFELD, IVLOC, K )
C
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
C
C     Dummy Arguments
C
      REAL,    INTENT(IN)    :: UFELD(IE,JE)
      INTEGER, INTENT(IN)    :: IVLOC,K
      REAL,    INTENT(INOUT) :: GUFELD(MOIE,MOJE)
C
C     Local Variables
C
      REAL    :: REBUF(MOIE*MOJE)
      INTEGER :: ILO, IUP, JLO, JUP,I,J,I1D,IX
C
C
      TAGCOUNT=1
      TAGTABLE(1)=IVLOC+10000+200*K

      IF (NEIGHBOR(1) .EQ. -1) THEN
         ILO =  2
      ELSE
         ILO = 3 
      ENDIF

      IF (NEIGHBOR(2) .EQ. -1) THEN
         JUP = JE-1
      ELSE
         JUP = JE-2
      ENDIF

      IF (NEIGHBOR(3) .EQ. -1) THEN
         IUP = IE-1
      ELSE
         IUP = IE-2
      ENDIF

      IF (NEIGHBOR(4) .EQ. -1) THEN
         JLO = 2
      ELSE
         JLO = 3
      ENDIF

      I1D = 0
      DO J = JLO,JUP
         DO I = ILO,IUP
            I1D = I1D + 1
            REBUF(I1D) = UFELD(I,J)
         ENDDO
      ENDDO

      IF (MYID /= 0) THEN
         TAGTABLE(1) = IVLOC + 10000 + 200*K
         TYPE  = TAGTABLE(1)
         COUNT = I1D*1
         DEST  = 0
         CALL PSENDR(REBUF(1))
      ELSE
         IX     = 0
         SOURCE = 0
         DO WHILE(IX.LT.NPROC)
            IF (IX .NE. 0) THEN

C               DO I = 1,MOIE*MOJE
C                  REBUF(I) = IX
C               ENDDO
               CALL PTEST
               CALL PRECVR(REBUF(1))
            ENDIF

            IF (ISUBNEIGH(SOURCE+1,1) .EQ. -1) THEN
               ILO = ISUBPOS(SOURCE+1,1) - 1
            ELSE
               ILO = ISUBPOS(SOURCE+1,1)
            ENDIF
            
            IF (ISUBNEIGH(SOURCE+1,2) .EQ. -1) THEN
               JUP = ISUBPOS(SOURCE+1,4) + 1
            ELSE
               JUP = ISUBPOS(SOURCE+1,4)
            ENDIF
            
            IF (ISUBNEIGH(SOURCE+1,3) .EQ. -1) THEN
               IUP = ISUBPOS(SOURCE+1,3) + 1
            ELSE
               IUP = ISUBPOS(SOURCE+1,3)
            ENDIF
            
            IF (ISUBNEIGH(SOURCE+1,4) .EQ. -1) THEN
               JLO = ISUBPOS(SOURCE+1,2) - 1
            ELSE
               JLO = ISUBPOS(SOURCE+1,2)
            ENDIF

            I1D = 0
            DO J = JLO ,JUP
               DO I = ILO ,IUP
                  I1D = I1D + 1
                  GUFELD(I,J) = REBUF(I1D)
               ENDDO
            ENDDO
            IX = IX+1
         ENDDO
      ENDIF

      CALL PWAIT

      CALL PSTOP

      END SUBROUTINE COLLECTDATA
