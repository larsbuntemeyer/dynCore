      SUBROUTINE HORINT(FPFL, FIB, FIBTEST)
C
C**** HORINT   -   UP:KORREKTUR DER AUF P-FLAECHEN INTERPOLIERTEN FELDER
C****                 MITTELS HORIZONTAL INTERPOLIERTER WERTE
C**   AUFRUF   :   CALL HORINT (FPFL, FIB, IE, JE, FIBTEST)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   KORREKTUR DER AUF P-FLAECHEN INTERPOLIERTEN FELDER
C**                MITTELS HORIZONTAL INTERPOLIERTER WERTE
C**   VERSIONS-
C**   DATUM    :   20.01.91
C**   EXTERNALS:   KEINE
C**   EINGABE-
C**   PARAMETER:   FPFL     FELD AUF P-FLAECHE
C**                FIB      FELD DER MODELLOROGRAPHIE
C**                IE,JE    DIMENSIONEN VON FPFL UND FIB
C**                FIBTEST  WENN FIB(I,J)>FIBTEST, SO WIRD DER WERT
C**                         FPFL(I,J) DURCH HORIZONTALE INTERPOLATION
C**                         (CRESSMAN-SCHEMA) ERSETZT
C**   AUSGABE-
C**   PARAMETER:   FPFL     KORRIGIERTES FELD AUF P-FLAECHE
C**   COMMON-
C**   BLOECKE  :   KEINE
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   D.MAJEWSKI
C
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
C
C     Dummy Arguments
C
      REAL, INTENT(IN)    :: FIBTEST
      REAL, INTENT(IN)    :: FIB(MOIE,MOJE)
      REAL, INTENT(INOUT) :: FPFL(MOIE,MOJE)
C
C     Local Variables
C
      REAL, DIMENSION(25) :: WGT, DIF, FIN
      REAL    :: FPFLINT,C
      INTEGER :: J,I,NUMI,N,M,L,KINCR,IDIF
C
      DO J = 1,MOJE
         INNERMO: DO I = 1,MOIE
C
            IF (FIB(I,J).LE.FIBTEST) CYCLE INNERMO
C
            KINCR = 1
C
            ENDLESS: DO
C
C     HERAUSSUCHEN BENACHBARTER PUNKTE
C
               NUMI = 0
C
               DO L = I-KINCR,I+KINCR
                  DO M = J-KINCR,J+KINCR
                     IF(L.LT.1.OR.L.GT.MOIE) CYCLE
                     IF(M.LT.1.OR.M.GT.MOJE) CYCLE
                     IF(FIB(L,M).GT.FIBTEST) CYCLE
                     IDIF = (L-I)**2 + (M-J)**2
                     IF(IDIF .LT. KINCR**2) CYCLE
                     NUMI = NUMI + 1
                     FIN(NUMI) = FPFL(L,M)
                     DIF(NUMI) = FLOAT(IDIF)
                     IF(KINCR.EQ.1 .AND. NUMI.GE. 6) EXIT ENDLESS
                     IF(KINCR.EQ.2 .AND. NUMI.GE.15) EXIT ENDLESS
                     IF(NUMI.GE.25) EXIT ENDLESS
                  ENDDO
               ENDDO
               KINCR = KINCR + 1
               IF(KINCR.GT.10) CYCLE INNERMO
C
            ENDDO ENDLESS
C
C     BERECHNUNG DER GEWICHTSFUNKTION
C
            C = 0.0
            DO N = 1,NUMI
               C = C + 1.0/DIF(N)
            ENDDO
            C = 1.0/C
            DO N = 1,NUMI
               WGT(N) = C/DIF(N)
            ENDDO
C
C     INTERPOLATION DES FEHLENDEN GITTERPUNKTES
C
            FPFLINT = 0.0
            DO N = 1,NUMI
               FPFLINT = FPFLINT + WGT(N)*FIN(N)
            ENDDO
            FPFL(I,J) = 0.25*FPFL(I,J) + 0.75*FPFLINT
C
         ENDDO INNERMO
      ENDDO
C
      RETURN
      END SUBROUTINE HORINT
