C
C
C**** MAKEPDB  -   UP:ERSTELLUNG DES PDB FUER EM-ERGEBNISDATEIEN
C**   AUFRUF   :   CALL MAKEPDB (YTYP, NZTX, NTRI)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   PRODUCT DEFINITION BLOCK (PDB) BIS AUF ELEMENTNUMMER
C**                UND LEVEL-ANGABEN BESETZEN FUER EM-ERGEBNIS-DATEIEN.
C**   VERSIONS-
C**   DATUM    :   16.03.89
C**
C**   EXTERNALS:   DATE, TIME
C**
C**   EINGABE-
C**   PARAMETER:   YTYP : TYP DER ERGEBNISDATEN:
C**                      'E': NORMALE ERGEBNISDATEN FUER GESAMTGEBIET
C**                      'T': TRAJEKTORIEN-DATEN    FUER GESAMTGEBIET
C**                      'D':         ERGEBNISDATEN FUER TEILGEBIET
C**                      'F': FORTSETZUNGSDATEN     FUER GESAMTGEBIET
C**                NZTX : ZEITSCHRITT
C**                NTRI : TIME RANGE INDICATOR IM GRIB-CODE (TAB. 5)
C**   AUSGABE-
C**   PARAMETER:   KEINE
C**
C**   COMMON-
C**   BLOECKE  :   ORG, CORG, COMDYN, GRIB
C**
C**   METHODE  :   BESETZEN ALLER PDB-WERTE BIS AUF ELEMENT-NUMMER UND
C**                LEVEL-ANGABEN
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   D.MAJEWSKI
C
      SUBROUTINE MAKEPDB (YTYP, NZTX, NTRI)
C
      IMPLICIT NONE
C
      INCLUDE "org.h"
      INCLUDE "corg.h"
      INCLUDE "comdia.h"
      INCLUDE "comdyn.h"
      INCLUDE "grib.h"
C
C     Dummy Arguments
C
      INTEGER,   INTENT(IN) :: NZTX, NTRI
      CHARACTER, INTENT(IN) :: YTYP*(*)
C
C     Local Variables
C
      CHARACTER :: YDATE*8, YTIME*10
      INTEGER   :: MONAT(12)
      DATA         MONAT  / 31 ,  28 ,  31 ,  30 ,  31 ,  30 ,
     &                      31 ,  31 ,  30 ,  31 ,  30 ,  31 /
C
      INTEGER :: I100, I4, I400, ISCHLT, JJERST, JJGB, M, MMERST, MMGB, 
     &           MNERST, NDERST, NDGB, NHERST, NHGB, NMLDAY, NVGB, NVGBA
C
C     PDB MIT -1 VORBESETZEN
      DO M  = 1,37
         IPDB (M) = -1
      ENDDO

C     BLOCK-FLAG SETZEN
      IPDB(6) = 128

C     FUER FORTSETZUNGSDATEIEN (YTYP='F'): VORHERSAGEZEIT IN MINUTEN
C     SONST IN STUNDEN
      IF(YTYP.EQ.'F') THEN
         IPDB(16) = 0
         IF(NZTX.GT.NDMXN) THEN
            NVGBA = NINT((((NZTX-1)/NDMXN)*NDMXN)*DT/60.0)
         ELSE
            NVGBA = NINT( NANF*DT/60.0)
         ENDIF
         NVGB     = NINT (NZTX*DT/60.0)
      ELSE IF(YTYP.EQ.'I') THEN
         IPDB(16) = 1
         NVGB     = 0
      ELSE
         IPDB(16) = 1
         IF(NZTX.GT.NDMXN) THEN
            NVGBA  = NINT((((NZTX-1)/NDMXN)*NDMXN)*DT/3600.0)
         ELSE
            NVGBA  = NINT (NANF*DT/3600.0)
         ENDIF
         NVGB      = NINT (NZTX*DT/3600.0)
      ENDIF
C     AKTUELLES DATUM IN IPDB SCHREIBEN

C     AKTUELLES DATUM LESEN
      READ(YAKDAT1,'(I4,3I2)') JJGB, MMGB, NDGB, NHGB

      NHGB=NHGB+1
      IF (NHGB.EQ.24) THEN
         NHGB=0
         NDGB=NDGB+1
      ENDIF
      IF (LMOMON) THEN
         IF (NDGB.EQ.31) THEN
            NDGB=1
            MMGB=MMGB+1
            IF (MMGB.EQ.13) THEN
               MMGB=1
               JJGB=JJGB+1
            ENDIF
         ENDIF
      ELSE
         NMLDAY=MONAT(MMGB)
C
C        SCHALTJAHR BESTIMMEN
         ISCHLT=0
         I400=MOD(JJGB,400)
         I100=MOD(JJGB,100)
         I4=MOD(JJGB,4)
         IF (I4.EQ.0) ISCHLT=1
         IF (I100.EQ.0) ISCHLT=0
         IF (I400.EQ.0) ISCHLT=1
         IF (MMGB.EQ.2) NMLDAY=MONAT(2)+ISCHLT
         IF (NDGB.EQ.NMLDAY+1) THEN
            NDGB=1
            MMGB=MMGB+1
            IF (MMGB.EQ.13) THEN
               MMGB=1
               JJGB=JJGB+1
            ENDIF
         ENDIF
      ENDIF
      IPDB(11) = JJGB
      IPDB(12) = MMGB
      IPDB(13) = NDGB
      IPDB(14) = NHGB

C     VORHERSAGEZEITPUNKT/ZEITRAUM
      IF(NTRI.EQ.1 .OR. NTRI.EQ.3 .OR. NTRI.EQ.4 .OR. NTRI.EQ.10) THEN
         IPDB(17) = 0
         IPDB(18) = NVGB
      ELSE IF(NTRI.EQ.2) THEN
         IPDB(17) = NVGBA
         IPDB(18) = NVGB
      ENDIF

      IPDB(19) = NTRI

C     ERSTELLUNGSDATUM UND VERSIONSNUMMER
C      IDATE    = DATE ( )
C      ITIME    = CLOCK( )
C      WRITE(YDATE,'(A8)') IDATE
C      WRITE(YTIME,'(A8)') ITIME
C      READ(YDATE(7:8),'(I2)') JJERST
C      READ(YDATE(1:2),'(I2)') MMERST
C      READ(YDATE(4:5),'(I2)') NDERST
C      READ(YTIME(1:2),'(I2)') NHERST
C      READ(YTIME(4:5),'(I2)') MNERST
C
      CALL DATE_AND_TIME(YDATE, YTIME)
      READ(YDATE(1:4),'(I4)') JJERST
      READ(YDATE(5:6),'(I2)') MMERST
      READ(YDATE(7:8),'(I2.2)') NDERST
      READ(YTIME(1:2),'(I2.2)') NHERST
      READ(YTIME(3:4),'(I2.2)') MNERST
C
      IPDB(32) = JJERST
      IPDB(33) = MMERST
      IPDB(34) = NDERST
      IPDB(35) = NHERST
      IPDB(36) = MNERST
      IPDB(37) = 1

      RETURN
      END SUBROUTINE MAKEPDB
