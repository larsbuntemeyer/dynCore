      SUBROUTINE READO3(ZOZACT, NOZ, YJDAT)
C
C**** READO3  -   UP:EINLESEN DER MONATSWEISE VARIABLEN OZONFELDER
C**
C**   AUFRUF   :   CALL READO3 IN UP *EC4ORG*
C**
C**   ENTRIES  :   KEINE
C**
C**   ZWECK    :   EINLESEN DER OZONFELDER ZOZACT
C**
C**   VERSIONS-
C**   DATUM    :   27.04.06
C**                2007
C**
C**   EXTERNALS:   KEINE
C**
C**   EINGABE-
C**   PARAMETER:   KEINE
C**
C**   AUSGABE-
C**   PARAMETER:   ZOZACT
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
      INTEGER,          INTENT(IN)  :: NOZ
      CHARACTER,        INTENT(IN)  :: YJDAT*10
      DOUBLE PRECISION, INTENT(OUT) :: ZOZACT(IEJE,NOZ,14)
C
C     Local Declarations
C      
      INTEGER, PARAMETER :: KNT=8
      REAL               :: UFELD(MOIEJE), VFELD(IEJE)
      INTEGER            :: IH(KNT)
      REAL               :: ZZOZACT(MOIEJE,NOZ,14)    
      REAL               :: AMD, AMO3, PPM2VGG
      INTEGER            :: IDA, IHA, IJ, IMA, IOST, IYA, 
     &                      JM, JY, K, M, NOZACT
C
C     KONSTANTEN BELEGEN
C
      TAGCOUNT=1
      TAGTABLE(1)=1
C
C     AKTUELLES DATUM EINLESEN
C
      READ(YJDAT,'(I4,3I2)') IYA,IMA,IDA,IHA
C
C     UNIT NUMMERN SETZEN
C
      NOZACT=25
C
C     FILE OEFFNEN
C
      IF (MYID.EQ.0) THEN
         CALL GETD(NOZACT,YO3DNAM,YO3DCAT)
      ENDIF
C
C     FELDER EINLESEN
C
C     LESEN VON DEZEMBER DES VORJAHRES BIS INKL. JANUAR DES FOLGEJAHRES
C
      DO
C
         IF (MYID.EQ.0) THEN
            DO M=1,14
               DO K=1,NOZ
                  READ(NOZACT,IOSTAT=IOST) IH
                  IF (IOST.EQ.0) THEN
                     READ(NOZACT) (ZZOZACT(IJ,K,M),IJ=1,MOIEJE)
                  ELSEIF (IOST.GT.0) THEN
                     PRINT *,'READO3: FEHLER BEIM LESEN VON OZON'
                     STOP 'ERROR'
                  ELSE
                     PRINT *,'READO3: DATEIENDE VON OZON ERREICHT'
                     STOP 'ERROR'
                  ENDIF
               ENDDO
            ENDDO
         ENDIF

         COUNT=KNT
         CALL PSENDALLI(IH)

         JY=IH(3)/10000
         JM=(IH(3)/100)-(JY*100)
         IF (JY.EQ.(IYA+1).AND.JM.EQ.1) EXIT
C
      ENDDO
      IF (MYID.EQ.0) THEN
         PRINT *,'READO3: DATUM=',JM,' ',JY
      ENDIF
C
C     DIMENSIONEN PRUEFEN
C
      IF (IH(5).NE.MOIE.AND.IH(6).NE.MOJE) THEN
         IF (MYID.EQ.0) THEN
            PRINT *,'READO3: DIMENSION ERROR'
            PRINT *,'REMO: IE=',MOIE,' INPUT: IE=',IH(5)
            PRINT *,'REMO: JE=',MOJE,' INPUT: JE=',IH(6)
         ENDIF
         STOP 'OZACT'
      ENDIF
C
C     UMSPEICHERN VON 32-BIT AUF 64-BIT UND
C     UMRECHNEN VON PPMV AUF G/G
C
      AMO3=47.9982
      AMD=28.970
      PPM2VGG=1.E-6*AMO3/AMD
C
      DO M=1,14
         DO K=1,NOZ

            IF (MYID.EQ.0) THEN

               DO IJ=1,MOIEJE
                  UFELD(IJ)=ZZOZACT(IJ,K,M)
               ENDDO

               CALL SSGB(UFELD,VFELD)
            ENDIF

            IF (MYID.NE.0) THEN
               CALL PTEST
               CALL PRECVR(VFELD)
            ENDIF
            CALL PSTOP
C
C           UMSPEICHERN
C
            DO IJ=1,IEJE
               ZOZACT(IJ,K,M)=VFELD(IJ)*PPM2VGG
            ENDDO
         ENDDO
      ENDDO
C
C     UNITS SCHLIESSEN
C
      IF (MYID.NE.0) THEN
         CLOSE(NOZACT)
      ENDIF
C
      RETURN
      END SUBROUTINE READO3
