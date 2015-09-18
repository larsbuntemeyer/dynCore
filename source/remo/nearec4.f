C
C
C**** NEARSFC  -   UP:BERECHNUNG DER BODENNAHEN FELDER
C**   AUFRUF   :   CALL NEARSFC  IN HP *DWDORG*
C**   ENTRIES  :      ---
C**   ZWECK    :   BERECHNUNG DER BODENNAHEN FELDER
C**   VERSIONS-
C**   DATUM    :   17.04.89
C**
C**   EXTERNALS:   KEINE
C**
C**   EINGABE-
C**   PARAMETER:   NX     (INT)   ZEITEBENE
C**   AUSGABE-
C**   PARAMETER:   KEINE
C**
C**   COMMON-
C**   BLOECKE  :   ORG, PHYKON, HIGKON, PARKON, COMDIA, COMPHY,
C**                COMDYN
C**
C**   METHODE  :   PROFILE IN DER PRANDTLSCHICHT VERWENDEN
C**   FEHLERBE-
C**   HANDLUNG :      ---
C**   VERFASSER:   D.MAJEWSKI
C
      SUBROUTINE NEAREC4(NX    ,
     &   T     , QD    , U     , V     , PS    , TB    , TG    ,
     &   QDB   , TMCM  , TMCH  , AZ0   , FIB   ,
     &   AKH   , BKH   , U10M  , V10M  , VBM10M, T2M   , TD2M  ,
     &   TMIN2M, TMAX2M)
C
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
      INCLUDE "org.h"
      INCLUDE "phykon.h"
      INCLUDE "higkon.h"
      INCLUDE "parkon.h"
      INCLUDE "comphy.h"
      INCLUDE "comdyn.h"
C
C     EINGABEFELDER VON *NEARSFC*
      INTEGER, INTENT(IN)    :: NX
      REAL,    INTENT(IN)    ::
     &          T (IEJE ,KE,3), QD (IEJE ,KE,3),
     &          U (IE,JE,KE,3), V  (IE,JE,KE,3),
     &          PS(IEJE    ,3), TB (IEJE    ,3),
     &          TG(IEJE    ,3), QDB(IEJE    ,3)
      REAL,    INTENT(IN)    ::
     &          TMCM(IEJE), TMCH(IEJE), FIB (IEJE),
     &          AZ0 (IEJE)
      REAL,    INTENT(IN)    :: AKH(KE), BKH(KE)

C     AUSGABEFELDER VON *NEARSFC*
      REAL,    INTENT(INOUT) ::
     &          U10M (IEJE), V10M  (IEJE), VBM10M(IEJE), T2M   (IEJE),
     &          TD2M (IEJE), TMIN2M(IEJE), TMAX2M(IEJE)

C     LOKALE FELDER VON *NEARSFC*
      REAL ::
     &          ZCM  (IEJE), ZCH  (IEJE), ZH  (IEJE), ZHS   (IEJE),
     &          ZSQCM(IEJE), ZH2M (IEJE), ZVPB(IEJE), ZCHDCM(IEJE),
     &          ZUMA (IEJE), ZVMA (IEJE), ZHKE(IEJE), ZUMX  (IEJE),
     &          ZVMX (IEJE), ZFIHF(IEJE), GZ0 (IEJE)
C
      INTEGER :: I, IJ, J
      REAL    :: Z10G, Z2G, ZEW2M, ZEWKE, ZFRAC, ZGQD2M, ZGQDKE, ZH05M, 
     &           ZLNZ3, ZP2M, ZPKE, ZQD2M, ZRIS, ZTVB, ZTVG, ZTVKE, 
     &           ZLNZ1, ZVMIN
CC
C     MINIMALE WINDGESCHWINDIGKEIT
      ZVMIN=  0.01

      Z10G = 10.0*G
      Z2G  =  2.0*G

      DO IJ=1,IEJE
         GZ0(IJ)=AZ0(IJ)*G
      ENDDO

C     WINDGESCHWINDIGKEIT IN DER PRANDTL-SCHICHT (K = KE) FUER
C     ZEITPUNKTE NA UND NX
      DO J  = 2, JE
         DO I  = 2, IE
            IJ       = I + (J - 1)*IE
            ZUMA(IJ) = 0.5*(U(I,J,KE,NA) + U(I-1,J  ,KE,NA))
            ZVMA(IJ) = 0.5*(V(I,J,KE,NA) + V(I  ,J-1,KE,NA))
            ZVPB(IJ) = MAX(SQRT(ZUMA(IJ)**2 + ZVMA(IJ)**2) , ZVMIN)
            ZUMX(IJ) = 0.5*(U(I,J,KE,NX) + U(I-1,J  ,KE,NX))
            ZVMX(IJ) = 0.5*(V(I,J,KE,NX) + V(I  ,J-1,KE,NX))
         ENDDO
      ENDDO

C     LINKEN UND UNTEREN RAND AUFFUELLEN
      DO I  = 2, IE
         ZUMA(I) = ZUMA(I + IE)
         ZVMA(I) = ZVMA(I + IE)
         ZVPB(I) = ZVPB(I + IE)
         ZUMX (I) = ZUMX (I + IE)
         ZVMX (I) = ZVMX (I + IE)
      ENDDO
      DO J  = 1, JE
         ZUMA(1 + (J - 1)*IE) = ZUMA(2 + (J - 1)*IE)
         ZVMA(1 + (J - 1)*IE) = ZVMA(2 + (J - 1)*IE)
         ZVPB(1 + (J - 1)*IE) = ZVPB(2 + (J - 1)*IE)
         ZUMX(1 + (J - 1)*IE) = ZUMX(2 + (J - 1)*IE)
         ZVMX(1 + (J - 1)*IE) = ZVMX(2 + (J - 1)*IE)
      ENDDO

C     U10M, V10M, VBM10M, T2M, TD2M, TMIN2M, TMAX2M BERECHNEN
C     U10M, V10M UND VBM10M SIND AM T, PS-GITTERPUNKT DEFINIERT.
C     VBM10M ENTHAELT DIE MAXIMALE ERWARTETE WINDBOE, DIE BERECHNET
C     WIRD ALS 3.0*(WURZEL AUS DER DOPPELTEN TURBULENTEN KINETISCHEN
C     ENERGIE IN DER PRANDTL-SCHICHT)

C     CM UND CH BERECHNEN
      DO IJ  = 1,IEJE
         ZTVB       = TB(IJ,NA)*(1.0 + RDDRM1*QDB(IJ,NA))
         ZCM   (IJ) = TMCM(IJ)/(G*ZVPB(IJ)*PS(IJ,NA)/(R*ZTVB))
         ZCH   (IJ) = TMCH(IJ)/(G*ZVPB(IJ)*PS(IJ,NA)/(R*ZTVB))
         ZCM   (IJ) = MAX ( ZCM(IJ), 5.0 E -4 )
         ZCH   (IJ) = MAX ( ZCH(IJ), 4.0 E -5 )

C        CM, CH UND GZ0 VORGEBEN, WENN KEINE BERECHNUNG DER GROESSEN
C        GEWUENSCHT WAR
         IF(.NOT.LPHY) THEN
            ZCM   (IJ) = 1.0 E-4
            ZCH   (IJ) = 1.0 E-4
            GZ0   (IJ) = MAX(GZ0(IJ), 1.0 E-3)
         ENDIF

C        GEOPOTENTIAL FUER DIE HAUPTFLAECHE K=KE
         ZFIHF (IJ) = FIB(IJ) + R*T(IJ,KE,NX)*
     &        (1.0 + RDDRM1*QD(IJ,KE,NX))*
     &        ALOG( PS(IJ,NX)/(AKH(KE) + BKH(KE)*PS(IJ,NX)) )

C        TROCKENSTATISCHE ENERGIE AM BODEN (S) UND IN DER PRANTL-SCHICHT
         ZH    (IJ) = WCP*T (IJ,KE,NX) + ZFIHF(IJ)
         ZHS   (IJ) = WCP*TG(IJ   ,NX) + FIB  (IJ)
         ZSQCM (IJ) = SQRT(ZCM(IJ))
         ZCHDCM(IJ) = ZCH(IJ)/(AKT*ZSQCM(IJ))
         ZHKE  (IJ) = ZFIHF(IJ) - FIB(IJ)

C     VIRTUELLE TEMPERATUR AM BODEN UND IN DER PRANDTL-SCHICHT
         ZTVG       = TG(IJ,   NX)*( 1.0 + RDDRM1*QDB(IJ,   NX) )
         ZTVKE      = T (IJ,KE,NX)*( 1.0 + RDDRM1*QD (IJ,KE,NX) )
         ZRIS       = ZTVKE - ZTVG + WCPR*ZHKE(IJ)

C     WIND IN 10 M UND TEMPERATUR IN 2 M
         IF( ZRIS.GT.0.0 ) THEN

C     STABILER FALL
            ZLNZ1     = ALOG  ((Z10G+GZ0(IJ))/GZ0(IJ)) - Z10G/ZHKE(IJ)*
     &           ALOG  ((ZHKE(IJ)+GZ0(IJ))/GZ0(IJ))
            U10M (IJ) = ZUMX(IJ)*(Z10G/ZHKE(IJ) + ZSQCM(IJ)/AKT*ZLNZ1)
            V10M (IJ) = ZVMX(IJ)*(Z10G/ZHKE(IJ) + ZSQCM(IJ)/AKT*ZLNZ1)

            ZH05M     = ZHS(IJ) + 0.25*( ZH(IJ) - ZHS(IJ) )
            ZH2M (IJ) = ZH05M  + 1.5*G*( ZH(IJ) - ZH05M )/
     &                                 ( ZHKE(IJ) - 0.5*G )
         ELSE

C     INSTABILER FALL
            ZSQCM(IJ) = MAX(ZSQCM (IJ), 0.01)
            ZCHDCM(IJ)= MAX(ZCHDCM(IJ), 0.01)
            ZLNZ3     = 1.0-ZSQCM(IJ)/AKT*
     &           ALOG(1.0+(EXP(AKT/ZSQCM(IJ))-1.0)
     &                *GZ0(IJ)*(ZHKE(IJ) - Z10G)/
     &                (ZHKE(IJ)*(Z10G + GZ0(IJ))))
            U10M (IJ) = ZUMX(IJ)*ZLNZ3
            V10M (IJ) = ZVMX(IJ)*ZLNZ3

            ZH2M (IJ) = ZHS(IJ) + (ZH(IJ) - ZHS(IJ))*
     &           (1.0-ZCHDCM(IJ)*ALOG  (1.0+(EXP  (1.0/ZCHDCM(IJ))-
     &           1.0)*(GZ0(IJ)*(ZHKE(IJ) - GZ0(IJ))/
     &           (ZHKE(IJ)*(Z2G + GZ0(IJ))))))
         ENDIF

         T2M (IJ)      = (ZH2M(IJ) - Z2G - FIB(IJ))*WCPR

C     TAUPUNKT IN 2 M
         ZP2M       = PS(IJ,NX)*(1.0 - Z2G/(R*T2M(IJ)*
     &        (1.0 + RDDRM1*QD(IJ,KE,NX))))
         ZPKE       = AKH(KE) + BKH(KE)*PS(IJ,NX)
         ZEW2M      = B1*EXP  (B2W*(T2M(IJ      )-B3)/
     &        (T2M(IJ      )-B4W))
         ZEWKE      = B1*EXP  (B2W*(T  (IJ,KE,NX)-B3)/
     &        (T  (IJ,KE,NX)-B4W))
         ZGQD2M     = RDRD*ZEW2M/(ZP2M - EMRDRD*ZEW2M)
         ZGQDKE     = RDRD*ZEWKE/(ZPKE - EMRDRD*ZEWKE)
         ZQD2M      = MAX( 1.0 E-10 , QD(IJ,KE,NX)*ZGQD2M/ZGQDKE )
         ZFRAC      = ALOG  (ZP2M*ZQD2M/(B1*(RDRD + EMRDRD*ZQD2M)))

         TD2M  (IJ) = (B2W*B3 - B4W*ZFRAC)/(B2W - ZFRAC)
         TD2M  (IJ) = MIN(TD2M(IJ) , T2M(IJ))

         TMIN2M(IJ) = MIN(TMIN2M(IJ), T2M(IJ))
         TMAX2M(IJ) = MAX(TMAX2M(IJ), T2M(IJ))
         VBM10M(IJ) = MAX(VBM10M(IJ), (1. + 3.0*2.4*ZSQCM(IJ))*ZVPB(IJ))

      ENDDO
C
      RETURN
      END SUBROUTINE NEAREC4
