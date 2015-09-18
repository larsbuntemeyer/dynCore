      SUBROUTINE ECKONT(
     &    PVERV, PDQDT,
     &    AK , BK , AKH   , BKH   , DAK   ,
     &    DBK  , ALPHABOUND, A1T, A2T, ACPHIR, CPHI  , PS    ,
     &    QDB  , U    , V  , T  , QD    , TMKVMH, TMCH)
      !
      IMPLICIT NONE
      !
C
C**** PROGORG  -   UP: STEUERUNG DER KONVEKTIONSPARAMETRISIERUNG
C**   AUFRUF   :   CALL KONTORG IN PROGORG
C**   ENTRIES  :      ---
C**   ZWECK    :   SCHEIBENWEISE BERECHNUNG DER KONVEKTION
C**   VERSIONS-
C**   DATUM    :   12.05.1992
C**                2007
C**
C**   EXTERNALS:   INITIE, EMTIED
C**
C**   EINGABE-
C**   PARAMETER:      ---
C**
C**   AUSGABE-
C**   PARAMETER:      ---
C**
C**   COMMON-
C**   BLOECKE  :   PARAM  , ORG, COMDYN, COMPHY, PHYKON, HIGKON, PARKON,
C**                COMPCST
C**
C**   METHODE  :   DER BEREICH VON JAH BIS JEH FUER DIE KONVEKTIONSPARA-
C**                METRISIERUNG WIRD IN NTASKS UNTERBEREICHE UNTERTEILT,
C**                DIE PARALLEL ABGEARBEITET WERDEN. INNERHALB EINES UN-
C**                TERBEREICHS WIRD DIE KONVEKTIONSRECHNUNG "SCHEIBEN-
C**                WEISE" DURCHGEFUEHRT. IN UP *KONTORG* WERDEN DIE
C**                SPEICHERBEREICHE FUER ALLE NTASKS UNTERBEREICHE
C**                ANGEFORDERT.
C**   FEHLERBE-
C**   HANDLUNG :      ---
C**   VERFASSER:   D.MAJEWSKI
C
      INCLUDE "parorg.h"
      INCLUDE "org.h"
C
C     VERTIKAL-KOORDINATEN-PARAMETER
      REAL, INTENT(IN) ::   AK(KE1), BK(KE1), AKH(KE), 
     &                      BKH(KE), DAK(KE), DBK(KE)
C
C     EXTERNE PARAMETER
      REAL, INTENT(IN) ::   ALPHABOUND(IE,JE,3)
C
C     VERTIKAL VARIIERENDER IMPLIZITHEITSGRAD DER V-DIFFUSION
      REAL, INTENT(IN) ::   A1T(KE1), A2T(KE1)
C
C     EXTERNE PARAMETER
      REAL, INTENT(IN) ::   ACPHIR(JE,2), CPHI(JE,2)
C
C     PROGNOSTISCHE BODENFELDER
      REAL, INTENT(IN) ::   PS(IE,JE,3), QDB(IE,JE,3)
C
C     ATMOSPHAEREN-FELDER
      REAL, INTENT(IN) ::   U(IE,JE,KE,3), V(IE,JE,KE,3),
     &                      T(IE,JE,KE,3), QD(IE,JE,KE,3)
C
C     FELDER DER KOEFFIZIENTEN UND FLUESSE (PHYSIKALISCHE UP'S)
      REAL, INTENT(IN) :: TMKVMH(IE*(KE-1),JE,2)
      REAL, INTENT(IN) :: TMCH(IE,JE)
C
C     UEBERGABEFELDER
      REAL, INTENT(INOUT) :: PVERV(IE,JE,KE),PDQDT(IE,JE,KE)
C
      !
      ! Local Variables
      !
      INTEGER :: I,J,K
      INTEGER :: IZJEA,IZJAA
      INTEGER :: JATKON,JETKON
      !
      JATKON = JAH
      JETKON = JEH
C
      CALL ECTIED (
     &    JATKON, JETKON, PVERV, PDQDT,
     &    AK , BK , AKH   ,
     &    BKH   , DAK   , DBK  , ALPHABOUND, A1T, A2T, ACPHIR,
     &    CPHI  , PS    , QDB  , U    , V  , T  , QD    ,
     &    TMKVMH, TMCH)
C
C
C     INNEREN RAND IN DEN AEUSSEREN SPIEGELN
C     --------------------------------------
C
C     1. UNTEREN RAND BELEGEN
C
      IF (NEIGHBOR(4) .EQ. -1) THEN
         DO K=1,KE
            DO I=IAH,IEH
               PVERV(I,JAA,K)=PVERV(I,JAH,K)
               PDQDT(I,JAA,K)=PDQDT(I,JAH,K)
               PVERV(I,JAA+1,K)=PVERV(I,JAH,K)
               PDQDT(I,JAA+1,K)=PDQDT(I,JAH,K)
            ENDDO
         ENDDO
      ENDIF
C
C     2. OBEREN RAND BELEGEN
C
      IF (NEIGHBOR(2) .EQ. -1) THEN
         DO K=1,KE
            DO I=IAH,IEH
               PVERV(I,JEA,K)=PVERV(I,JEH,K)
               PDQDT(I,JEA,K)=PDQDT(I,JEH,K)
               PVERV(I,JEA-1,K)=PVERV(I,JEH,K)
               PDQDT(I,JEA-1,K)=PDQDT(I,JEH,K)
            ENDDO
         ENDDO
      ENDIF
C
C     3. LINKEN RAND BELEGEN
C
      IF (NEIGHBOR(1) .EQ. -1) THEN
         DO K=1,KE
            IZJAA=JAH
            IZJEA=JEH
            IF (NEIGHBOR(4) .EQ. -1) IZJAA=JAA
            IF (NEIGHBOR(2) .EQ. -1) IZJEA=JEA
            DO J=IZJAA,IZJEA
               PVERV(IAA,J,K)=PVERV(IAH,J,K)
               PDQDT(IAA,J,K)=PDQDT(IAH,J,K)
               PVERV(IAA+1,J,K)=PVERV(IAH,J,K)
               PDQDT(IAA+1,J,K)=PDQDT(IAH,J,K)
            ENDDO
         ENDDO
      ENDIF
C
C     4. RECHTEN RAND BELEGEN
C
      IF (NEIGHBOR(3) .EQ. -1) THEN
         DO K=1,KE
            IZJAA=JAH
            IZJEA=JEH
            IF (NEIGHBOR(4) .EQ. -1) IZJAA=JAA
            IF (NEIGHBOR(2) .EQ. -1) IZJEA=JEA
            DO J=IZJAA,IZJEA
               PVERV(IEA,J,K)=PVERV(IEH,J,K)
               PDQDT(IEA,J,K)=PDQDT(IEH,J,K)
               PVERV(IEA-1,J,K)=PVERV(IEH,J,K)
               PDQDT(IEA-1,J,K)=PDQDT(IEH,J,K)
            ENDDO
         ENDDO
      ENDIF
C
      RETURN
      END SUBROUTINE ECKONT
