      SUBROUTINE READBOD(ZALB, ZVGR, ZVLT)
C
C**** READBOD  -   UP:EINLESEN DER MONATSWEISE VARIABLEN BODENFELDER
C**
C**   AUFRUF   :   CALL READBOD IN UP *EC4ORG*
C**
C**   ENTRIES  :   KEINE
C**
C**   ZWECK    :   EINLESEN DER MONATSWEISE VARIABLEN BODENFELDER
C**                ALBECH, VGRAT UND VLT
C**
C**   VERSIONS-
C**   DATUM    :   22.02.02
C**                2007
C**
C**   EXTERNALS:   KEINE
C**
C**   EINGABE-
C**   PARAMETER:   KEINE
C**
C**   AUSGABE-
C**   PARAMETER:   ZALB, ZVGR, ZVLT
C**
C**   COMMON-
C**   BLOECKE  :   KEINE
C**
C**   METHODE  :   -
C**
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   R.PODZUN
C
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
      INCLUDE "corg.h"
C
C     Dummy Arguments
C
      DOUBLE PRECISION, INTENT(OUT) :: ZALB(IEJE,12), ZVGR(IEJE,12),
     &                                 ZVLT(IEJE,12)
C
C     Local Declarations
C
      INTEGER, PARAMETER :: KNT=8
      DOUBLE PRECISION   :: UFELD(MOIEJE), VFELD(IEJE)
      INTEGER            :: IH(KNT), NUN(NZMXVB)   
      REAL               :: ZZALB(MOIEJE,12), ZZVGR(MOIEJE,12),
     &                      ZZVLT(MOIEJE,12)
      INTEGER            :: I, IJ, M
C
C     KONSTANTEN BELEGEN
C
      TAGCOUNT=1
      TAGTABLE(1)=1
C
C
C     FILES OEFFNEN
C
      IF (MYID.EQ.0) THEN
         DO I=1,NZMXVB
            NUN(I)=I+12
            CALL GETD(NUN(I),YBDNAM(I),YBDCAT)
         ENDDO
      ENDIF

C
C     FELDER EINLESEN
C
      DO I=1,NZMXVB
         DO M=1,12
            IF (MYID.EQ.0) THEN
               READ(NUN(I)) IH
            ENDIF

            COUNT=KNT
            CALL PSENDALLI(IH(1))
C
C     DIMENSIONEN PRUEFEN
C

            IF (IH(5).NE.MOIE.AND.IH(6).NE.MOJE) THEN
               IF (MYID.NE.0) THEN
                  PRINT *,'READBOD: DIMENSION ERROR FELD NR:',IH(1)
                  PRINT *,'REMO: IE=',MOIE,' INPUT: IE=',IH(5)
                  PRINT *,'REMO: JE=',MOJE,' INPUT: JE=',IH(6)
               ENDIF
               STOP 'ERROR'
            ENDIF

C
            IF (IH(1).EQ.174) THEN
               IF (MYID.EQ.0) THEN
                  READ(NUN(I)) (ZZALB(IJ,M),IJ=1,MOIEJE)

                  DO IJ=1,MOIEJE
                     UFELD(IJ)=ZZALB(IJ,M)
                  ENDDO

                  CALL SSGB(UFELD,VFELD)
               ENDIF

               IF (MYID.NE.0) THEN
                  CALL PTEST
                  CALL PRECVR(VFELD(1))
               ENDIF
               CALL PSTOP

               DO IJ=1,IEJE
                  ZALB(IJ,M)=VFELD(IJ)
               ENDDO
            ENDIF

            IF (IH(1).EQ.198) THEN
               IF (MYID.EQ.0) THEN
                  READ(NUN(I)) (ZZVGR(IJ,M),IJ=1,MOIEJE)

                  DO IJ=1,MOIEJE
                     UFELD(IJ)=ZZVGR(IJ,M)
                  ENDDO

                  CALL SSGB(UFELD,VFELD)
               ENDIF

               IF (MYID.NE.0) THEN
                  CALL PTEST
                  CALL PRECVR(VFELD(1))
               ENDIF
               CALL PSTOP

               DO IJ=1,IEJE
                  ZVGR(IJ,M)=VFELD(IJ)
               ENDDO
            ENDIF

            IF (IH(1).EQ.200) THEN
               IF (MYID.EQ.0) THEN
                  READ(NUN(I)) (ZZVLT(IJ,M),IJ=1,MOIEJE)

                  DO IJ=1,MOIEJE
                     UFELD(IJ)=ZZVLT(IJ,M)
                  ENDDO

                  CALL SSGB(UFELD,VFELD)
               ENDIF

               IF (MYID.NE.0) THEN
                  CALL PTEST
                  CALL PRECVR(VFELD(1))
               ENDIF
               CALL PSTOP

               DO IJ=1,IEJE
                  ZVLT(IJ,M)=VFELD(IJ)
               ENDDO
            ENDIF

         ENDDO
      ENDDO

      IF (MYID.EQ.0) THEN
         DO I=1,NZMXVB
            NUN(I)=I+12
            CLOSE(NUN(I))
         ENDDO
      ENDIF
C
      RETURN
      END SUBROUTINE READBOD
