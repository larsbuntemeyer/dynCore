      SUBROUTINE ECRANDUP(  A1   , A2     ,
     &     U      , V     , T    , QD     , QW     , QI     ,
     &     UR     , VR    , TR   , QDR    , QWR    , PS     , TSWECH,
     &     TSIECH , SEAICE, INFRL, INFRW  , INFRI  , QDB    , QDBL  ,
     &     QDBW   , QDBI  , PSR  , TSWECHR, TSIECHR, SEAICER, QDBLR ,
     &     SICED  , SICEDR, QIR  , PINT   , DWDT   , W      , AK    ,
     &     BK)
C
C     SUBROUTINE ECRANDUP
C
C**** ECRANDUP -   UP:VORBESETZUNG MIT RANDWERTEN
C**   AUFRUF   :   CALL ECRANDUP(A1,A2) IN UP *PROGEC4*
C**   ENTRIES  :      ---
C**   ZWECK    :   VORBESETZUNG ALLER PROGNOSTISCHEN VARIABLEN MIT
C**                INTERPOLIERTEN RANDDATEIWERTEN
C**   VERSIONS-
C**   DATUM    :   21.06.01
C**                2007
C**
C**   EXTERNALS:      ---
C**
C**   EINGABE-
C**   PARAMETER:   A1,A2 : PARAMETER FUER LINEARE INTERPOLATION
C**                   A1 : FAKTOR FUER DEN ERSTEN  STUETZSTELLENWERT
C**                   A2 : FAKTOR FUER DEN ZWEITEN STUETZSTELLENWERT
C**   AUSGABE-
C**   PARAMETER:      ---
C**
C**   COMMON-
C**   BLOECKE  :   ORG, PHYKON, PARKON, COMDIA, HIGKON, COMPHY,
C**
C**   METHODE  :   LINEARE ZEITLICHE INTERPOLATION DER RANDWERTE
C**                F(T) = A1(T)*F(T1) + A2(T)*F(T2)
C**   FEHLERBE-
C**   HANDLUNG :      ---
C**   VERFASSER:   R.PODZUN
C
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
      INCLUDE "org.h"
      INCLUDE "parkon.h"
      INCLUDE "higkon.h"
      INCLUDE "faktinf.h"
C
C     Dummy Arguments
C     ERFORDERLICHE FELDER DIMENSIONIEREN
C
C
      INTEGER, INTENT(IN) ::
     &          INFRL (IEJE  ), INFRW (IEJE  ),
     &          INFRI (IEJE)
      REAL, INTENT(IN) :: A1, A2
      REAL, INTENT(IN) :: AK(KE1), BK(KE1)

      REAL, INTENT(INOUT) ::
     &          U (IEJEKE,3), V (IEJEKE,3), T (IEJEKE,3),
     &          QD(IEJEKE,3), QW(IEJEKE,3), QI(IEJEKE,3)
C
      REAL,    INTENT(IN)    ::
     &          UR (IEJEKE,2), VR (IEJEKE,2), TR (IEJEKE,2),
     &          QDR(IEJEKE,2), QWR(IEJEKE,2)
C
      REAL,    INTENT(INOUT) :: QIR(IEJEKE,2)
C
      REAL, INTENT(INOUT) ::
     &          PS  (IEJE,3), TSWECH(IEJE,3), TSIECH(IEJE,3),
     &          SEAICE(IEJE), QDB   (IEJE,3), SICED (IEJE  ),
     &          QDBL(IEJE,3), QDBW  (IEJE,3), QDBI  (IEJE,3)
C
      REAL,    INTENT(IN)    ::
     &          PSR    (IEJE,2), QDBLR  (IEJE,2), SICEDR (IEJE,2),
     &          TSWECHR(IEJE,2), TSIECHR(IEJE,2), SEAICER(IEJE,2)
C
      REAL, INTENT(INOUT) :: PINT(IEJE,KE1,3), DWDT(IEJE,KE,3), 
     &                       W(IEJE,KE1,3)
C
C     Local Varibles
C
      LOGICAL :: LQIR
C
      INTEGER :: IJ,K
C
C     STATEMENT - FUNCTIONS ZUR BERECHNUNG DER SAEETIGUNGSFEUCHTE
C     MAGNUS-FORMEL FUER WASSER
C      ZGEW(TT)        = B1 * EXP  ( B2W*(TT - B3)/(TT - B4W) )
C      ZGEE(TT)        = B1 * EXP  ( B2E*(TT - B3)/(TT - B4E) )
C     SPEZIFISCHE FEUCHTE AUS DAMPFDRUCK UND LUFTDRUCK
C      ZGQD(GE,PP)     = RDRD*GE/(PP - EMRDRD*GE)

C     BODENDRUCK VORBESETZEN
      DO IJ = 1 , IEJE
         PS (IJ,NE) = A1* PSR(IJ,NRD1) + A2* PSR(IJ,NRD2)
      ENDDO
C
C     ATMOSPHAERENFELDER VORBESETZEN; WENN QW ALS RANDFELD GEGEBEN,
C     AUCH QW BESETZEN, SONST IST QW=0.
      DO IJ  = 1 , IEJEKE
         U (IJ,NE) = A1* UR (IJ,NRD1) + A2* UR (IJ,NRD2)
         V (IJ,NE) = A1* VR (IJ,NRD1) + A2* VR (IJ,NRD2)
         T (IJ,NE) = A1* TR (IJ,NRD1) + A2* TR (IJ,NRD2)
         QD(IJ,NE) = A1* QDR(IJ,NRD1) + A2* QDR(IJ,NRD2)
      ENDDO
      IF (LQWR) THEN
         DO IJ  = 1 , IEJEKE
            QW(IJ,NE) = A1* QWR(IJ,NRD1) + A2* QWR(IJ,NRD2)
         ENDDO
      ELSE
         DO IJ  = 1 , IEJEKE
            QW(IJ,NE) = 0.0
         ENDDO
      ENDIF
C
CKS   QI IS NOT A BOUNDARY FIELD YET. IT NEEDS TO BE SET ANYWAY, BECAUSE
C     OF THE OUTFLOW CONDITION FORMULATION. OTHERWISE IT CAN HAPPEN THAT
C     SPURIOUS CLOUD ICE IS FORMED (SEE ECRANDAS).
C     THE FORMULATION IS ALREADY PREPARED FOR A VERSION WHERE CLOUD ICE
C     CAN OPTIONALLY BE READ FROM A FORCING FIELD.
C
      LQIR = .FALSE.
      IF (LQIR) THEN
        QI(:,NE) = A1* QIR(:,NRD1) + A2* QIR(:,NRD2)
      ELSE
        QI(:,NE) = 0.0
        QIR(:,:) = 0.0
      ENDIF
C
C     BODENFELDER VORBESETZEN
C
      IF (LSICED) THEN
         DO IJ = 1 , IEJE
            SICED (IJ   )= A1* SICEDR (IJ,NRD1) + A2* SICEDR (IJ,NRD2)
         ENDDO
      ENDIF

      DO IJ = 1 , IEJE
        TSWECH(IJ,NE) = A1* TSWECHR(IJ,NRD1) + A2* TSWECHR(IJ,NRD2)
        TSIECH(IJ,NE) = A1* TSIECHR(IJ,NRD1) + A2* TSIECHR(IJ,NRD2)
        SEAICE(IJ   ) = A1* SEAICER(IJ,NRD1) + A2* SEAICER(IJ,NRD2)
        QDBL  (IJ,NE) = A1* QDBLR  (IJ,NRD1) + A2* QDBLR  (IJ,NRD2)
        QDBW  (IJ,NE) = ZGQD( ZGEW( TSWECH(IJ,NE) ), PS(IJ,NE) )
        QDBI  (IJ,NE) = ZGQD( ZGEE( TSIECH(IJ,NE) ), PS(IJ,NE) )
        QDB   (IJ,NE) = (REAL(INFRL(IJ))*QDBL(IJ,NE)
     &                +  REAL(INFRW(IJ))*QDBW(IJ,NE)
     &                +  REAL(INFRI(IJ))*QDBI(IJ,NE))*EDFAKINF
      ENDDO
C
C     INITIALIZE NON-HYDROSTATIC VARIABLES
C
      DO K  = 1 , KE1
        W(1:IEJE,K,NE) = 0.
        PINT(1:IEJE,K,NE) = AK(K) + BK(K)*PS(1:IEJE,NE)
      ENDDO

      DO K  = 1  , KE
        DWDT(1:IEJE,K,NE) = 0.
      ENDDO
         !
         !
         !
         CONTAINS
         !
         ! STATEMENT - FUNCTIONS ZUR BERECHNUNG DER SAEETIGUNGSFEUCHTE
         ! MAGNUS-FORMEL FUER WASSER
         !
         REAL FUNCTION ZGEW(TT)
            IMPLICIT NONE
            REAL :: TT
            ZGEW = B1 * EXP  ( B2W*(TT - B3)/(TT - B4W) )
         END FUNCTION ZGEW
         !
         REAL FUNCTION ZGEE(TT)
            IMPLICIT NONE
            REAL :: TT
            ZGEE = B1 * EXP  ( B2E*(TT - B3)/(TT - B4E) )
         END FUNCTION ZGEE
         !
         ! SPEZIFISCHE FEUCHTE AUS DAMPFDRUCK UND LUFTDRUCK
         !  
         REAL FUNCTION ZGQD(GE,PP)
            IMPLICIT NONE
            REAL :: GE,PP
            ZGQD = RDRD*GE/(PP - EMRDRD*GE)
         END FUNCTION ZGQD
         !
         !
         ! 
      END SUBROUTINE ECRANDUP
