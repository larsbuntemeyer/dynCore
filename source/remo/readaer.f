      SUBROUTINE READAER(ZSO4ALL, ZSO4NAT, YJDAT)
C
C**** READAER  -   UP:EINLESEN DER MONATSWEISE VARIABLEN AEROSOLE
C**
C**   AUFRUF   :   CALL READAER IN UP *EC4ORG*
C**
C**   ENTRIES  :   KEINE
C**
C**   ZWECK    :   EINLESEN DER AEROSOLE SO4ALL, SO4NAT
C**
C**   VERSIONS-
C**   DATUM    :   04.04.06
C**                2007
C**
C**   EXTERNALS:   KEINE
C**
C**   EINGABE-
C**   PARAMETER:   KEINE
C**
C**   AUSGABE-
C**   PARAMETER:   ZSO4ALL, ZSO4NAT
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
      CHARACTER,        INTENT(IN)  :: YJDAT*10
      DOUBLE PRECISION, INTENT(OUT) :: ZSO4ALL(IEJE,KE,14)
      DOUBLE PRECISION, INTENT(OUT) :: ZSO4NAT(IEJE,KE,14)
C
C     Local Declarations
C
      INTEGER, PARAMETER :: KNT=8
      REAL    :: UFELD(MOIEJE), VFELD(IEJE)
      INTEGER :: IH(KNT)
      REAL    :: ZZSO4ALL(MOIEJE,KE,14), ZZSO4NAT(MOIEJE,KE,14)
      INTEGER :: IDA, IHA, IJ, IMA, IOST, IYA, JM, JY, K, M, NSO4ALL, 
     &           NSO4NAT
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
      NSO4ALL=26
      NSO4NAT=27
C
C     FILES OEFFNEN
C
      IF (MYID.EQ.0) THEN
         CALL GETD(NSO4ALL,YSADNAM,YSADCAT)
         CALL GETD(NSO4NAT,YSNDNAM,YSNDCAT)
      ENDIF
C
      DO M=2,13
         DO K=1,KE
            IF (MYID.EQ.0) THEN
               READ(NSO4NAT) IH
               READ(NSO4NAT) (ZZSO4NAT(IJ,K,M),IJ=1,MOIEJE)
            ENDIF

            COUNT=KNT
            CALL PSENDALLI(IH)
C      CALL PSTOP
C
C           DIMENSIONEN PRUEFEN
C
            IF (IH(5).NE.MOIE.AND.IH(6).NE.MOJE) THEN
               IF (MYID.EQ.0) THEN
                  PRINT *,'READAER: DIMENSION ERROR'
                  PRINT *,'REMO: IE=',MOIE,' INPUT: IE=',IH(5)
                  PRINT *,'REMO: JE=',MOJE,' INPUT: JE=',IH(6)
               ENDIF
               STOP 'SO4NAT'
            ENDIF

         ENDDO
         IF (MYID.EQ.0) THEN
            PRINT*,'READAER: SO4NAT',IH
         ENDIF
      ENDDO

C
C     DEZ IN JAN KOPIEREN BZW. UMGEKEHRT
C
      IF (MYID.EQ.0) THEN
         DO K=1,KE
            DO IJ=1,IEJE
               ZZSO4NAT(IJ,K,1)=ZZSO4NAT(IJ,K,13)
               ZZSO4NAT(IJ,K,14)=ZZSO4NAT(IJ,K,2)
            ENDDO
         ENDDO
      ENDIF
C
C     EINTEILEN
C
      DO M=1,14
         DO K=1,KE

            IF (MYID.EQ.0) THEN

               DO IJ=1,MOIEJE
                  UFELD(IJ)=ZZSO4NAT(IJ,K,M)
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
               ZSO4NAT(IJ,K,M)=VFELD(IJ)
            ENDDO
         ENDDO
      ENDDO
C
C     LESEN VON DEZEMBER DES VORJAHRES BIS INKL. JANUAR DES FOLGEJAHRES
C
      DO
C
         IF (MYID.EQ.0) THEN
            DO M=1,14
               DO K=1,KE
                  READ(NSO4ALL,IOSTAT=IOST) IH
                  IF (IOST.EQ.0) THEN
                     READ(NSO4ALL) (ZZSO4ALL(IJ,K,M),IJ=1,MOIEJE)
                  ELSEIF (IOST.GT.0) THEN
                     PRINT *,'READAER: FEHLER BEIM LESEN VON SO4ALL'
                     STOP 'ERROR'
                  ELSE
                     PRINT *,'READAER: DATEIENDE VON SO4ALL ERREICHT'
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
         PRINT *,'READAER: DATUM=',JM,' ',JY
      ENDIF
C
C     DIMENSIONEN PRUEFEN
C
      IF (IH(5).NE.MOIE.AND.IH(6).NE.MOJE) THEN
         IF (MYID.EQ.0) THEN
            PRINT *,'READAER: DIMENSION ERROR'
            PRINT *,'REMO: IE=',MOIE,' INPUT: IE=',IH(5)
            PRINT *,'REMO: JE=',MOJE,' INPUT: JE=',IH(6)
         ENDIF
         STOP 'SO4ALL'
      ENDIF
C
C     EINTEILEN
C
      DO M=1,14
         DO K=1,KE

            IF (MYID.EQ.0) THEN

               DO IJ=1,MOIEJE
                  UFELD(IJ)=ZZSO4ALL(IJ,K,M)
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
               ZSO4ALL(IJ,K,M)=VFELD(IJ)
            ENDDO
         ENDDO
      ENDDO
C
C     UNITS SCHLIESSEN
C
      IF (MYID.NE.0) THEN
         CLOSE(NSO4ALL)
         CLOSE(NSO4NAT)
      ENDIF
C
      RETURN
      END SUBROUTINE READAER
