      SUBROUTINE ECRANDAS
     &    (AKH , BKH , ALPHABOUND, PS  ,
     &     U   , V   , T   , QD  , QW    , UR    ,
     &     VR  , TR  , QDR , QWR , PSR   , QI    ,
     &     QIR , PINT, DWDT, W   , AK    , BK    , DAK   , DBK)
C
C     SUBROUTINE ECRANDAS
C
C**** RANDASS  -   UP: RANDRELAXATION UND ASSELIN-FILTERUNG
C**   AUFRUF   :   CALL RANDASS IN UP *PROGEC4*
C**   ENTRIES  :      ---
C**   ZWECK    :   RANDRELAXATION NACH DAVIES (1976) UND
C**                ASSELIN-ZEITFILTER
CKS                UPDATE OF HALOS
C**   VERSIONS-
C**   DATUM    :   21.11.03
C**                2007
C**
C**   EXTERNALS:   KEINE
C**
C**   EINGABE-
C**   PARAMETER:   KEINE
C**
C**   AUSGABE-
C**   PARAMETER:      ---
C**
C**   COMMON-
C**   BLOECKE  :   ORG, COMDYN, PHYKON, HIGKON, PARKON, COMDIA,
C**                PROGCHK
C**
C**   METHODE  :   NACH DAVIES (1976) IN DER VEREINFACHTEN FORM VON
C**                KALLBERG (1979)
C**   FEHLERBE-
C**   HANDLUNG :      ---
C**   VERFASSER:   G. DOMS / D.MAJEWSKI / T.SEMMLER / R.PODZUN / K.SIECK
C
CKS   REPLACED JAHCOMP, JEHCOMP ETC... BY 1, JE FOR FULL FIELD TREATMENT
CKS   IN ASSELIN FILTER. SAVES COMMUNICATION
C
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
      INCLUDE "org.h"
      INCLUDE "comdyn.h"
      INCLUDE "phykon.h"
      INCLUDE "higkon.h"
      INCLUDE "parkon.h"
      INCLUDE "progchk.h"
      INCLUDE "comdia.h"
C
      INCLUDE "haloexchorg"
C
C     Dummy Arguments
C-----------------------------------------------------------------------
C     ERFORDERLICHE EM-FELDER AUS DEM LANGZEITSPEICHER DIMENSIONIEREN
C     ---------------------------------------------------------------
C     VERTIKAL-KOORDINATEN-PARAMETER
      REAL, INTENT(IN)    :: AKH(KE), BKH(KE)

C     EXTERNE PARAMETER
      REAL, INTENT(IN)    :: ALPHABOUND(IEJE,3)

C     PROGNOSTISCHE BODENFELDER
      REAL, INTENT(INOUT) :: PS(IEJE,3)

      REAL, INTENT(INOUT) ::
     1          PINT (IEJE,KE1,3), DWDT(IEJE,KE,3),
     2          W    (IEJE,KE1,3), T(IEJE,KE,3)
      REAL, INTENT(IN) ::
     3          AK   (KE1)       , BK  (KE1)       ,
     4          DAK  (KE)        , DBK (KE)

C     ATMOSPHAEREN-FELDER
      REAL, INTENT(INOUT) :: U(IEJE,KE,3)
      REAL, INTENT(INOUT) :: V(IEJE,KE,3)
      REAL, INTENT(INOUT) :: QD(IEJE,KE,3)
      REAL, INTENT(INOUT) :: QW(IEJE,KE,3)
      REAL, INTENT(INOUT) :: QI(IEJE,KE,3)

C     RANDWERT-FELDER
      REAL, INTENT(IN)    ::  UR (IEJE,KE,2),
     &                        VR (IEJE,KE,2),
     &                        TR (IEJE,KE,2),
     &                        QDR(IEJE,KE,2),
     &                        QWR(IEJE,KE,2),
     &                        PSR(IEJE,2   ),
     &                        QIR(IEJE,KE,2)
C
C     LOKALE FELDER DIMENSIONIEREN
C     ----------------------------
      REAL    :: ZMYA (IEJE   ), ZMYB (IEJE   ),
     &           ZMYPS(IEJE   ), ZZMY (IE,JE,KE),
     &           ZPSRE(IEJE   )
      REAL :: ZDP, ZDPI

      REAL    :: ZWGTVD(KE)
      REAL    :: ZMYQDW
      REAL    :: ZQDR, ZQWR, ZQIR, ZTR

C     FELDER FUER SUBROUTINEN EXRANDASS1 UND EXRANDASS2
      LOGICAL :: LSCHALT1(JE*KE),LSCHALT2(JE*KE),
     &           LSCHALT3(IE*KE),LSCHALT4(IE*KE)

      LOGICAL :: LSCHALT5(KE),LSCHALT6(KE),LSCHALT7(KE),LSCHALT8(KE)
C
      REAL    :: ZVRE,ZVFK,ZURE,ZEMZDR,ZDX,ZDTIRD
      INTEGER :: I,J,K,JA,JI,IA,II,IJ,IJA,IJI
C
      INTEGER :: KSIZE(10), ITYPE
C  
C     CHECK FOR OUTFLOW CONDITIONS AT THE BOUNDARIES
C     TRACER LIKE VARIABLES WILL NOT BE RELAXED TO THE BOUNDARY VALUE
C     IF OUTFLOW IS PRESENT
C
      CALL EXRANDASS1(U,V,NE,
     &                LSCHALT1,LSCHALT2,LSCHALT3,LSCHALT4,
     &                LSCHALT5,LSCHALT6,LSCHALT7,LSCHALT8)

C-----------------------------------------------------------------------
C     VORBEREITENDE MASSNAHMEN
C     KOEFFIZIENTEN FUER LINEARE ZEITLICHE RANDWERTINTERPOLATION
      ZDTIRD = REAL( NZT + 1 - NZTRD ) / REAL( NDR )
      ZEMZDR = 1. - ZDTIRD

C     ERFORDERLICHE LOKALE FELDER ADRESSIEREN
C     ---------------------------------------

C     VERSTAERKUNGSFAKTOR FUER DIE VERTIKALE NESTUNG
      IF ( LVNEST ) THEN
         ZDX = RERD/MAX ( EDDLAM, EDDPHI )
         DO K = 1,KE
            IF ( K.LT.KFL850 ) THEN
               ZVFK = MAX ( 0.0, DT*SQRT(2.0)*VBMXV/ZDX - 1.0 )
               ZVFK = MIN ( ZVFK, 1.0 )
               ZWGTVD(K) =
     &              ( 1.0 - MIN ( (AKH(K)+BKH(K)*1.0 E5), 0.85 E5 )/
     &              0.85 E5 )**2*ZVFK*0.5
            ELSE
               ZWGTVD(K) = 0.0
            ENDIF
         ENDDO
      ELSE
         DO K = 1,KE
            ZWGTVD(K) = 0.0
         ENDDO
      ENDIF
C
C     COPY BOUNDARY RELAXATION FACTORS TO LOCAL ARRAYS
C
      ZMYPS(:) = ALPHABOUND(:,1) ! SURFACE PRESSURE
      ZMYA(:)  = ALPHABOUND(:,2) ! U VELOCITY
      ZMYB(:)  = ALPHABOUND(:,3) ! V VELOCITY
C
C     COPY TO ARRAY THAT TAKES IN-/OUTFLOW CONDITIONS INTO ACCOUNT
C     (SEE EXRANDASS2)
C
      DO J = 1, JE
         DO I = 1, IE
            IJ = I + (J-1)*IE
            ZZMY(I,J,:) = ALPHABOUND(IJ,1)
         ENDDO
      ENDDO
C
C     COPY INNER (THIRD ROW/LINE) TO BOUNDARIES FOR OUTFLOW CONDITION
C     TREATMENT OF TRACER LIKE VARIABLES
C
CRP
C
CKS      DO K=1,KE ! REPLACED BY F90 NOTATION
C
C     3.REIHE IN 2. REIHE SPIEGELN
C     --------------------------------------
C
C     1. LINKEN RAND BELEGEN
C
      IF (NEIGHBOR(1) .EQ. -1) THEN
        IA=2
        II=3
        DO J=2,JE-1
          IJA=IA+(J-1)*IE
          IJI=II+(J-1)*IE
          QD(IJA,:,NE)=QD(IJI,:,NE)
          QW(IJA,:,NE)=QW(IJI,:,NE)
          QI(IJA,:,NE)=QI(IJI,:,NE)
          T (IJA,:,NE)=T (IJI,:,NE)
        END DO
      ENDIF
C
C     2. OBEREN RAND BELEGEN
C
      IF (NEIGHBOR(2) .EQ. -1) THEN
        JA=JE-1
        JI=JE-2
        DO I=2,IE-1
          IJA=I+(JA-1)*IE
          IJI=I+(JI-1)*IE
          QD(IJA,:,NE)=QD(IJI,:,NE)
          QW(IJA,:,NE)=QW(IJI,:,NE)
          QI(IJA,:,NE)=QI(IJI,:,NE)
          T (IJA,:,NE)=T (IJI,:,NE)
        END DO
      ENDIF
C
C     3. RECHTEN RAND BELEGEN
C
      IF (NEIGHBOR(3) .EQ. -1) THEN
        IA=IE-1
        II=IE-2
        DO J=2,JE-1
          IJA=IA+(J-1)*IE
          IJI=II+(J-1)*IE
          QD(IJA,:,NE)=QD(IJI,:,NE)
          QW(IJA,:,NE)=QW(IJI,:,NE)
          QI(IJA,:,NE)=QI(IJI,:,NE)
          T (IJA,:,NE)=T (IJI,:,NE)
        END DO
      ENDIF
C
C     4. UNTEREN RAND BELEGEN
C
      IF (NEIGHBOR(4) .EQ. -1) THEN
        JA=2
        JI=3
        DO I=2,IE-1
          IJA=I+(JA-1)*IE
          IJI=I+(JI-1)*IE
          QD(IJA,:,NE)=QD(IJI,:,NE)
          QW(IJA,:,NE)=QW(IJI,:,NE)
          QI(IJA,:,NE)=QI(IJI,:,NE)
          T (IJA,:,NE)=T (IJI,:,NE)
        END DO
      ENDIF
C
C     ECKPUNKT LINKS UNTEN BELEGEN
C
      IF ( (NEIGHBOR(1) .EQ. -1) .AND. (NEIGHBOR(4) .EQ. -1) ) THEN
        I=2
        J=2
        IJA=I+(J-1)*IE
        I=3
        J=3
        IJI=I+(J-1)*IE
        QD(IJA,:,NE)=QD(IJI,:,NE)
        QW(IJA,:,NE)=QW(IJI,:,NE)
        QI(IJA,:,NE)=QI(IJI,:,NE)
        T (IJA,:,NE)=T (IJI,:,NE)
      ENDIF
C
C     ECKPUNKT LINKS OBEN BELEGEN
C
      IF ( (NEIGHBOR(1) .EQ. -1) .AND. (NEIGHBOR(2) .EQ. -1) ) THEN
        I=2
        J=JE-1
        IJA=I+(J-1)*IE
        I=3
        J=JE-2
        IJI=I+(J-1)*IE
        QD(IJA,:,NE)=QD(IJI,:,NE)
        QW(IJA,:,NE)=QW(IJI,:,NE)
        QI(IJA,:,NE)=QI(IJI,:,NE)
        T (IJA,:,NE)=T (IJI,:,NE)
      ENDIF
C
C     ECKPUNKT RECHTS OBEN BELEGEN
C
      IF ( (NEIGHBOR(2) .EQ. -1) .AND. (NEIGHBOR(3) .EQ. -1) ) THEN
        I=IE-1
        J=JE-1
        IJA=I+(J-1)*IE
        I=IE-2
        J=JE-2
        IJI=I+(J-1)*IE
        QD(IJA,:,NE)=QD(IJI,:,NE)
        QW(IJA,:,NE)=QW(IJI,:,NE)
        QI(IJA,:,NE)=QI(IJI,:,NE)
        T (IJA,:,NE)=T (IJI,:,NE)
      ENDIF
C
C     ECKPUNKT RECHTS UNTEN BELEGEN
C
      IF ( (NEIGHBOR(3) .EQ. -1) .AND. (NEIGHBOR(4) .EQ. -1) ) THEN
        I=IE-1
        J=2
        IJA=I+(J-1)*IE
        I=IE-2
        J=3
        IJI=I+(J-1)*IE
        QD(IJA,:,NE)=QD(IJI,:,NE)
        QW(IJA,:,NE)=QW(IJI,:,NE)
        QI(IJA,:,NE)=QI(IJI,:,NE)
        T (IJA,:,NE)=T (IJI,:,NE)
      ENDIF
C
C     2.REIHE IN 1. REIHE SPIEGELN
C     --------------------------------------
C
C     1. LINKEN RAND BELEGEN
C
      IF (NEIGHBOR(1) .EQ. -1) THEN
        IA=1
        II=2
        DO J=1,JE
          IJA=IA+(J-1)*IE
          IJI=II+(J-1)*IE
          QD(IJA,:,NE)=QD(IJI,:,NE)
          QW(IJA,:,NE)=QW(IJI,:,NE)
          QI(IJA,:,NE)=QI(IJI,:,NE)
          T (IJA,:,NE)=T (IJI,:,NE)
        END DO
      ENDIF
C
C     2. OBEREN RAND BELEGEN
C
      IF (NEIGHBOR(2) .EQ. -1) THEN
        JA=JE
        JI=JE-1
        DO I=1,IE
          IJA=I+(JA-1)*IE
          IJI=I+(JI-1)*IE
          QD(IJA,:,NE)=QD(IJI,:,NE)
          QW(IJA,:,NE)=QW(IJI,:,NE)
          QI(IJA,:,NE)=QI(IJI,:,NE)
          T (IJA,:,NE)=T (IJI,:,NE)
        END DO
      ENDIF
C
C     3. RECHTEN RAND BELEGEN
C
      IF (NEIGHBOR(3) .EQ. -1) THEN
        IA=IE
        II=IE-1
        DO J=1,JE
          IJA=IA+(J-1)*IE
          IJI=II+(J-1)*IE
          QD(IJA,:,NE)=QD(IJI,:,NE)
          QW(IJA,:,NE)=QW(IJI,:,NE)
          QI(IJA,:,NE)=QI(IJI,:,NE)
          T (IJA,:,NE)=T (IJI,:,NE)
        END DO
      ENDIF
C
C     4. UNTEREN RAND BELEGEN
C
      IF (NEIGHBOR(4) .EQ. -1) THEN
        J=2
        DO I=1,IE
          IJ=I+(J-1)*IE
          QD(I,:,NE)=QD(IJ,:,NE)
          QW(I,:,NE)=QW(IJ,:,NE)
          QI(I,:,NE)=QI(IJ,:,NE)
          T (I,:,NE)=T (IJ,:,NE)
        END DO
      ENDIF
C
C     ECKPUNKT LINKS UNTEN BELEGEN
C
      IF ( (NEIGHBOR(1) .EQ. -1) .AND. (NEIGHBOR(4) .EQ. -1) ) THEN
        I=1
        J=1
        IJA=I+(J-1)*IE
        I=2
        J=2
        IJI=I+(J-1)*IE
        QD(IJA,:,NE)=QD(IJI,:,NE)
        QW(IJA,:,NE)=QW(IJI,:,NE)
        QI(IJA,:,NE)=QI(IJI,:,NE)
        T (IJA,:,NE)=T (IJI,:,NE)
      ENDIF
C
C     ECKPUNKT LINKS OBEN BELEGEN
C
      IF ( (NEIGHBOR(1) .EQ. -1) .AND. (NEIGHBOR(2) .EQ. -1) ) THEN
        I=1
        J=JE
        IJA=I+(J-1)*IE
        I=2
        J=JE-1
        IJI=I+(J-1)*IE
        QD(IJA,:,NE)=QD(IJI,:,NE)
        QW(IJA,:,NE)=QW(IJI,:,NE)
        QI(IJA,:,NE)=QI(IJI,:,NE)
        T (IJA,:,NE)=T (IJI,:,NE)
      ENDIF
C
C     ECKPUNKT RECHTS OBEN BELEGEN
C
      IF ( (NEIGHBOR(2) .EQ. -1) .AND. (NEIGHBOR(3) .EQ. -1) ) THEN
        I=IE
        J=JE
        IJA=I+(J-1)*IE
        I=IE-1
        J=JE-1
        IJI=I+(J-1)*IE
        QD(IJA,:,NE)=QD(IJI,:,NE)
        QW(IJA,:,NE)=QW(IJI,:,NE)
        QI(IJA,:,NE)=QI(IJI,:,NE)
        T (IJA,:,NE)=T (IJI,:,NE)
      ENDIF
C
C     ECKPUNKT RECHTS UNTEN BELEGEN
C
      IF ( (NEIGHBOR(3) .EQ. -1) .AND. (NEIGHBOR(4) .EQ. -1) ) THEN
        I=IE
        J=1
        IJA=I+(J-1)*IE
        I=IE-1
        J=2
        IJI=I+(J-1)*IE
        QD(IJA,:,NE)=QD(IJI,:,NE)
        QW(IJA,:,NE)=QW(IJI,:,NE)
        QI(IJA,:,NE)=QI(IJI,:,NE)
        T (IJA,:,NE)=T (IJI,:,NE)
      ENDIF
C
CKS      END DO
C
CRP
C
C     COMPUTE MODIFIED RELAXATION FACTORS FOR OUTFLOW CONDITIONS
C     OF TRACER LIKE VARIABLES
C
      CALL EXRANDASS2(ZZMY,
     &                LSCHALT1,LSCHALT2,LSCHALT3,LSCHALT4,
     &                LSCHALT5,LSCHALT6,LSCHALT7,LSCHALT8)
C     1. RANDRELAXATION
C
C     U- UND V-GLEICHUNG
C     ------------------
C
      DO K = 1,KE

         DO J = JAHCOMP, JEHCOMP
            DO I = IAHCOMP, IEHCOMP
               IJ = I + (J-1)*IE
               ZURE        = ZEMZDR*UR(IJ,K,NRD1) + ZDTIRD*UR(IJ,K,NRD2)
               U(IJ,K,NE)  = U(IJ,K,NE) -
     &              ZMYA (IJ)*( U(IJ,K,NE) - ZURE )
               U(IJ,K,NE)  = U(IJ,K,NE) -
     &              ZWGTVD(K)*( U(IJ,K,NE) - ZURE )
               ZVRE        = ZEMZDR*VR(IJ,K,NRD1) + ZDTIRD*VR(IJ,K,NRD2)
               V(IJ,K,NE)  = V(IJ,K,NE) -
     &              ZMYB (IJ)*( V(IJ,K,NE) - ZVRE )
               V(IJ,K,NE)  = V(IJ,K,NE) -
     &              ZWGTVD(K)*( V(IJ,K,NE) - ZVRE )
            ENDDO
         ENDDO

      ENDDO !K = 1,KE
C
C     PS - GLEICHUNG
C     ----------------------
C
      DO J = JAHCOMP, JEHCOMP
        DO I = IAHCOMP, IEHCOMP
          IJ = I + (J-1)*IE
          ZPSRE(IJ)   = ZEMZDR*PSR(IJ,NRD1) + ZDTIRD*PSR(IJ,NRD2)
          PS(IJ,NE) = PS(IJ,NE) - ZMYPS(IJ)*(PS(IJ,NE) - ZPSRE(IJ))
        ENDDO
      ENDDO

CHG   DRUECKE UND VERTIKALGESCHWINDIGKEIT
      DO K = 1, KE1
        DO J = JAHCOMP, JEHCOMP
          DO I = IAHCOMP, IEHCOMP
            IJ = I + (J-1)*IE
            PINT(IJ,K,NE) = PINT(IJ,K,NE)
     &        - ZMYPS(IJ)*(PINT(IJ,K,NE)
     &                    - (AK(K)+BK(K)*ZPSRE(IJ)))
            W(IJ,K,NE) = (1. - ZMYPS(IJ))*W(IJ,K,NE)
          ENDDO
        ENDDO
      ENDDO
C
C     TEMPERATURE AND HUMIDITY
C     ------------------------
C
      DO K = 1,KE

         DO J = JAHCOMP, JEHCOMP
            DO I = IAHCOMP, IEHCOMP
               IJ = I + (J-1)*IE
               ZMYQDW      = ZZMY(I,J,K)
               ZQDR        = ZEMZDR*QDR(IJ,K,NRD1) +
     &                       ZDTIRD*QDR(IJ,K,NRD2)
               ZQWR        = ZEMZDR*QWR(IJ,K,NRD1) +
     &                       ZDTIRD*QWR(IJ,K,NRD2)
               ZQIR        = ZEMZDR*QIR(IJ,K,NRD1) +
     &                       ZDTIRD*QIR(IJ,K,NRD2)
               ZTR         = ZEMZDR*TR (IJ,K,NRD1) +
     &                       ZDTIRD*TR (IJ,K,NRD2)
               QD(IJ,K,NE) = QD(IJ,K,NE) -
     &                       ZMYQDW * (QD(IJ,K,NE) - ZQDR)
               QD(IJ,K,NE) = QD(IJ,K,NE) -
     &                       ZWGTVD(   K)*(QD(IJ,K,NE) - ZQDR)
               QW(IJ,K,NE) = QW(IJ,K,NE) -
     &                       ZMYQDW * (QW(IJ,K,NE) - ZQWR)
               QW(IJ,K,NE) = QW(IJ,K,NE) -
     &                       ZWGTVD(   K)*(QW(IJ,K,NE) - ZQWR)
               QI(IJ,K,NE) = QI(IJ,K,NE) -
     &                       ZMYQDW * (QI(IJ,K,NE) - ZQIR)
               QI(IJ,K,NE) = QI(IJ,K,NE) -
     &                       ZWGTVD(   K)*(QI(IJ,K,NE) - ZQIR)
               T (IJ,K,NE) = T (IJ,K,NE) -
     &                       ZMYQDW * (T (IJ,K,NE) - ZTR)
               T (IJ,K,NE) = T (IJ,K,NE) -
     &                       ZWGTVD(   K)*(T (IJ,K,NE) - ZTR)
            ENDDO
         ENDDO

      ENDDO
C
C     2. HALO UPDATE
C
      ITYPE = 90
      KSIZE(1:10) = (/KE,KE,KE,1,KE,KE,KE,KE1,KE,KE1/)
      CALL HALOEXCH(KSIZE, ITYPE, U(:,:,NE), V(:,:,NE), T(:,:,NE),
     &     PS(:,NE), QD(:,:,NE), QW(:,:,NE), QI(:,:,NE),
     &     PINT(:,:,NE), DWDT(:,:,NE), W(:,:,NE))
C
C     3. ASSELIN - ZEITFILTERUNG
C-----------------------------------------------------------------------
      DO IJ = 1, IEJE
        PS(IJ,NJ) = PS(IJ,NJ)
     &            + EPSASS*(PS(IJ,NE)- 2.*PS(IJ,NJ) + PS(IJ,NA))
      ENDDO

      DO K = 1, KE
        DO IJ = 1, IEJE
          U(IJ,K,NJ)  = U(IJ,K,NJ)
     &                + EPSASS*( U(IJ,K,NE)-2.* U(IJ,K,NJ)+ U(IJ,K,NA))
          V(IJ,K,NJ)  = V(IJ,K,NJ)
     &                + EPSASS*( V(IJ,K,NE)-2.* V(IJ,K,NJ)+ V(IJ,K,NA))
          T (IJ,K,NJ) = T (IJ,K,NJ)
     &                + EPSASS*( T(IJ,K,NE)-2.* T(IJ,K,NJ)+ T(IJ,K,NA))
          QD(IJ,K,NJ) = QD(IJ,K,NJ)
     &                + EPSASS*(QD(IJ,K,NE)-2.*QD(IJ,K,NJ)+QD(IJ,K,NA))
          QW(IJ,K,NJ) = QW(IJ,K,NJ)
     &                + EPSASS*(QW(IJ,K,NE)-2.*QW(IJ,K,NJ)+QW(IJ,K,NA))
          QI(IJ,K,NJ) = QI(IJ,K,NJ)
     &                + EPSASS*(QI(IJ,K,NE)-2.*QI(IJ,K,NJ)+QI(IJ,K,NA))
          W(IJ,K,NJ)  = W(IJ,K,NJ) 
     &                + EPSASS*( W(IJ,K,NE)-2.* W(IJ,K,NJ)+ W(IJ,K,NA))
        ENDDO
      ENDDO

      DO K = 1, KE1
        DO IJ = 1, IEJE
          PINT(IJ,K,NJ) = PINT(IJ,K,NJ) + EPSASS
     &         * (PINT(IJ,K,NE) - 2.*PINT(IJ,K,NJ) + PINT(IJ,K,NA))
        ENDDO
      ENDDO

CHG   DWDT ENTSPRECHEND DEN NEUEN DRUECKEN ANPASSEN
      DO K  = 1, KE
        DO IJ = 1, IEJE
          ZDPI          = DAK(K) + DBK(K)*PS(IJ,NE)
          ZDP           = PINT(IJ,K+1,NE) - PINT(IJ,K,NE)
          DWDT(IJ,K,NE) = ZDP/ZDPI
          ZDPI          = DAK(K) + DBK(K)*PS(IJ,NJ)
          ZDP           = PINT(IJ,K+1,NJ) - PINT(IJ,K,NJ)
          DWDT(IJ,K,NJ) = ZDP/ZDPI
        ENDDO
      ENDDO

      END SUBROUTINE ECRANDAS
