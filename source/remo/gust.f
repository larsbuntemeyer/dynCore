C
C     SUBROUTINE GUST 
C
C**** NEARSFC  -   UP:BERECHNUNG DER MAXIMALEN WINDBOE
C**   AUFRUF   :   CALL GUST IN UP *EC4ORG*
C**   ENTRIES  :      ---
C**   ZWECK    :   BERECHNUNG DER MAXIMALEN WINDBOE
C**   VERSIONS-
C**   DATUM    :   16.02.05
C**
C**   EXTERNALS:   KEINE
C**
C**   EINGABE-
C**   PARAMETER:   NX     (INT)   ZEITEBENE
C**   AUSGABE-
C**   PARAMETER:   VBM10M
C**
C**   COMMON-
C**   BLOECKE  :   ORG, PHYKON, HIGKON, COMPHY
C**
C**   METHODE  :   PROFILE IN DER PRANDTLSCHICHT VERWENDEN
C**   FEHLERBE-
C**   HANDLUNG :      ---
C**   VERFASSER:   R.PODZUN (NACH D.MAJEWSKI)
C
C
      SUBROUTINE GUST(NX    ,
     &   T     , QD    , U     , V     , PS    , TB    , TG    ,
     &   QDB   , TMCM  , FIB   , AKH   , BKH   , VBM10M)
C
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
      INCLUDE "org.h"
      INCLUDE "phykon.h"
      INCLUDE "higkon.h"
      INCLUDE "comphy.h"
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN) :: NX
C
C     EINGABEFELDER VON *GUST*
      REAL,    INTENT(IN) ::
     &          T (IEJE ,KE,3), QD (IEJE ,KE,3),
     &          U (IE,JE,KE,3), V  (IE,JE,KE,3),
     &          PS(IEJE    ,3), TB (IEJE    ,3),
     &          TG(IEJE    ,3), QDB(IEJE    ,3)
C
      REAL, INTENT(IN)    ::   TMCM(IEJE), FIB(IEJE)
C
      REAL, INTENT(IN)    ::   AKH(KE), BKH(KE)
C
C     AUSGABEFELDER VON *GUST*
C
      REAL, INTENT(INOUT) ::   VBM10M(IEJE)
C
C     Local Variables
C
C     LOKALE FELDER VON *GUST*
C
      REAL ::
     &          ZCM  (IEJE), ZSQCM(IEJE), ZVPB(IEJE),
     &          ZUMA (IEJE), ZVMA (IEJE), ZHKE(IEJE),
     &          ZFIHF(IEJE)
C
      INTEGER :: I,J,IJ
      REAL :: ZVMIN,ZTVKE,ZTVG,ZTVB,ZRIS
C
C     MINIMALE WINDGESCHWINDIGKEIT
C
      ZVMIN=0.01
C
C     WINDGESCHWINDIGKEIT IN DER PRANDTL-SCHICHT (K = KE) FUER
C     ZEITPUNKTE NA UND NX
C
      DO J  = 2, JE
         DO I  = 2, IE
            IJ       = I + (J - 1)*IE
            ZUMA(IJ) = 0.5*(U(I,J,KE,NA) + U(I-1,J  ,KE,NA))
            ZVMA(IJ) = 0.5*(V(I,J,KE,NA) + V(I  ,J-1,KE,NA))
            ZVPB(IJ) = MAX(SQRT(ZUMA(IJ)**2 + ZVMA(IJ)**2) , ZVMIN)
         ENDDO
      ENDDO
C
C     LINKEN UND UNTEREN RAND AUFFUELLEN
C
      DO I  = 2, IE
         ZUMA(I) = ZUMA(I + IE)
         ZVMA(I) = ZVMA(I + IE)
         ZVPB(I) = ZVPB(I + IE)
      ENDDO
      DO J  = 1, JE
         ZUMA(1 + (J - 1)*IE) = ZUMA(2 + (J - 1)*IE)
         ZVMA(1 + (J - 1)*IE) = ZVMA(2 + (J - 1)*IE)
         ZVPB(1 + (J - 1)*IE) = ZVPB(2 + (J - 1)*IE)
      ENDDO
C
C     VBM10M BERECHNEN
C     VBM10M IST AM T, PS-GITTERPUNKT DEFINIERT.
C     VBM10M ENTHAELT DIE MAXIMALE ERWARTETE WINDBOE, DIE BERECHNET
C     WIRD ALS 3.0*(WURZEL AUS DER DOPPELTEN TURBULENTEN KINETISCHEN
C     ENERGIE IN DER PRANDTL-SCHICHT)
C
C     CM UND CH BERECHNEN
C
      DO IJ  = 1,IEJE
        IF (LPHY) THEN
          ZTVB       = TB(IJ,NA)*(1.0 + RDDRM1*QDB(IJ,NA))
          ZCM   (IJ) = TMCM(IJ)/(G*ZVPB(IJ)*PS(IJ,NA)/(R*ZTVB))
          ZCM   (IJ) = MAX ( ZCM(IJ), 5.0 E -4 )
C
C     CM VORGEBEN, WENN KEINE BERECHNUNG DER GROESSEN
C     GEWUENSCHT WAR
C
        ELSE
          ZCM   (IJ) = 1.0 E-4
        ENDIF
C
C     GEOPOTENTIAL FUER DIE HAUPTFLAECHE K=KE
C
         ZFIHF (IJ) = FIB(IJ) + R*T(IJ,KE,NX)*
     &                (1.0 + RDDRM1*QD(IJ,KE,NX))*
     &                ALOG( PS(IJ,NX)/(AKH(KE) + BKH(KE)*PS(IJ,NX)) )
C
C     TROCKENSTATISCHE ENERGIE AM BODEN (S) UND IN DER PRANTL-SCHICHT
C
         ZSQCM (IJ) = SQRT(ZCM(IJ))
         ZHKE  (IJ) = ZFIHF(IJ) - FIB(IJ)
C
C     VIRTUELLE TEMPERATUR AM BODEN UND IN DER PRANDTL-SCHICHT
C
         ZTVG       = TG(IJ,   NX)*( 1.0 + RDDRM1*QDB(IJ,   NX) )
         ZTVKE      = T (IJ,KE,NX)*( 1.0 + RDDRM1*QD (IJ,KE,NX) )
         ZRIS       = ZTVKE - ZTVG + WCPR*ZHKE(IJ)
C
C     WIND IN 10 M
C
         IF( ZRIS.LE.0.0 ) THEN
C
C     INSTABILER FALL
C
            ZSQCM(IJ) = MAX(ZSQCM (IJ), 0.01)
         ENDIF
C
         VBM10M(IJ) = MAX(VBM10M(IJ), (1. + 3.0*2.4*ZSQCM(IJ))*ZVPB(IJ))
C
      ENDDO
C
      RETURN
      END SUBROUTINE GUST
