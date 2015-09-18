      SUBROUTINE FINAL(
     &    AK     , BK    , DAK    , DBK  ,    
     &    DWDT   , PINT  , W      , T      , QD     , QW   , PS    ,
     &    DIFW   , DIFD  , QI)

      INCLUDE "phykon.h"
      INCLUDE "higkon.h"
      INCLUDE "org.h"
      INCLUDE "parorg.h"

      REAL, INTENT(IN)    ::
     &         QD(IE,JE,KE,3) ,
     &         QW(IE,JE,KE,3) ,
     &         QI(IE,JE,KE,3) ,
     &         PS(IE,JE,3), 
     &         DIFD(IE,JE,KE),
     &         DIFW(IE,JE,KE)

      REAL, INTENT(IN)    ::
     &         AK  (KE1)  , BK  (KE1)  ,
     &         DAK (KE)   , DBK (KE) 

      REAL, INTENT(INOUT) :: DWDT(IE,JE,KE,3)
      REAL, INTENT(INOUT) :: T(IE,JE,KE,3)
      REAL, INTENT(INOUT) :: PINT(IE,JE,KE1,3)
      REAL, INTENT(INOUT) :: W(IE,JE,KE1,3)

      REAL, PARAMETER :: EPS   = 1.
      REAL, PARAMETER :: WP    = 0.05
      REAL, PARAMETER :: WP1   = 1. - WP
      REAL, PARAMETER :: WA    = 0.05

      INTEGER :: I, J, K
      INTEGER :: KM1, KP1
      REAL    :: CAPPA, GDT, GDT2, FCC, WGHT
      REAL    :: DPTL, DPTU, DELP, DPSTR
      REAL    :: PP1, RDPLDN, RDPLUP
      REAL    :: TTFC
      REAL    :: RTOP
      REAL ::
     &         B1(KE)   , B2(KE)   , B3(KE)   , C0(KE)    ,
     &         PONE(KE1), PSTR(KE1), PNP1(KE1), COFF(KE1) ,
     &         CHI(KE1)

      REAL ::  P1(IE,KE1)

C-----------------------------------------------------------------------

      CAPPA = R/WCP
      GDT   = G*DT2
      GDT2  = GDT*GDT
      FCC   = -R/GDT2
      WGHT  = CAPPA


CHG   SICHERHEITSABFANG
!
!     HERE DWDT BECOMES EPSILON!
!
      DO K = 1 , KE
         DO J = JAH , JEH
            DO I = IAH, IEH
               DWDT(I,J,K,NE) = DWDT(I,J,K,NE) 
     &              + WA*(DIFD(I,J,K) + DIFW(I,J,K))
               DWDT(I,J,K,NE) = MAX( DWDT(I,J,K,NE),-EPS )
               DWDT(I,J,K,NE) = MIN( DWDT(I,J,K,NE), EPS )
               DWDT(I,J,K,NE) = (1. + DWDT(I,J,K,NE)/G)*WP1
     &              + DWDT(I,J,K,NA)*WP
            ENDDO
         ENDDO
      ENDDO

      DO J = JAH , JEH !J-LOOP
        DO I = IAH , IEH
          PDP = PS(I,J,NE)
          CHI(1) = 0.
C
          PONE(1) = AK(1) + BK(1)*PDP
          PSTR(1) = PONE(1)
          PNP1(1) = PONE(1)
          P1(I,1) = PONE(1)
C
          DO K = 2 , KE1
            P1(I,K)=PINT(I,J,K,NE)
            PONE(K)=PINT(I,J,K,NE)
            KM1 = K - 1
            DPSTR   = DWDT(I,J,KM1,NE)*(DAK(KM1)+DBK(KM1)*PDP)
            PSTR(K) = PSTR(KM1) + DPSTR
            PP1     = PNP1(KM1) + DPSTR
            DP      = (PP1 - PONE(K))*WGHT
            PNP1(K) = PONE(K) + DP
          ENDDO
C
          DO K = 2 , KE1
            KM1 = K - 1
            TFAC    = 1. + 0.608 * QD(I,J,KM1,NE) -
     &           (QW(I,J,KM1,NE) + QI(I,J,KM1,NE))
            TTFC    = 1. - CAPPA * TFAC
            COFF(KM1) = TFAC * T(I,J,KM1,NE) * (DAK(KM1)+DBK(KM1) * PDP)
     &           * FCC * TTFC / ((PNP1(KM1)+PNP1(K)))**2
          ENDDO

          DO K = 2 , KE
            KM1 = K - 1
            KP1 = K + 1
            DFRC=((PSTR(KM1)+PSTR(K)-PONE(KM1)-PONE(K))*COFF(KM1)
     &           +(PSTR(K)+PSTR(KP1)-PONE(K)-PONE(KP1))*COFF(K))*0.5
            RDPLDN = 1./(DAK(K)+DBK(K)*PDP)
            RDPLUP = 1./(DAK(KM1)+DBK(KM1)*PDP)

            B1(K) = COFF(KM1)*0.5 + RDPLUP
            B2(K) = (COFF(KM1) + COFF(K))*0.5 - (RDPLUP + RDPLDN)
            B3(K) = COFF(K)*0.5 + RDPLDN

            C0(K) = -DFRC
          ENDDO

          B2(KE) = B2(KE) + B3(KE)

          DO K = 3 , KE
            KM1 = K - 1
            TMP   = -B1(K)/B2(KM1)
            B2(K) = B3(KM1)*TMP + B2(K)
            C0(K) = C0(KM1)*TMP + C0(K)
          ENDDO

          CHI(KE)  = C0(KE)/B2(KE)
          CHI(KE1) = CHI(KE)

          DO K = KE-1 , 2 , -1
            CHI(K) = (-B3(K)*CHI(K+1) + C0(K))/B2(K)
          ENDDO

          PNP1(1:KE1)   = CHI(1:KE1) + PSTR(1:KE1)
          PINT(I,J,1:KE1,NE) = PNP1(1:KE1)
        ENDDO

C----------BACKSUBSTITUTION---------------------------------------------
        DO I = IAH , IEH
          DPTU = 0.
          DO K = 1 , KE
            KP1  = K + 1
            PDP=PS(I,J,NE)
C
            DPTL = PINT(I,J,KP1,NE) - P1(I,KP1)
            TFAC = 1. + 0.608 * QD(I,J,K,NA) -
     &           (QW(I,J,K,NA) + QI(I,J,K,NA))
            RTOP = T(I,J,K,NA) * TFAC * R /
     &           ((PINT(I,J,K,NA) + PINT(I,J,KP1,NA)) * 0.5)
            T(I,J,K,NE) = (DPTU + DPTL) * RTOP * 0.5 / WCP + T(I,J,K,NE)
            DELP      = (PINT(I,J,KP1,NE) - PINT(I,J,K,NE))/
     &           (DAK(K) + DBK(K) * PDP)
            DWDTT = DWDT(I,J,K,NE)
            DWDT(I,J,K,NE) = DELP 
            WIL = W(I,J,K,NE) + (DELP-DWDTT) * GDT
            W(I,J,K,NE) = WIL
C
            DPTU = DPTL
          ENDDO
        ENDDO
       
      ENDDO ! J-LOOP

      END SUBROUTINE FINAL
