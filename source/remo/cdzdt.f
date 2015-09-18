      SUBROUTINE CDZDT(
     &    DAK    , DBK   , U      , V      ,
     &    ACPHIR , DT    , FIB    , QI     ,
     &    DWDT   , PINT  , W      , T      , QD     , QW   ,
     &    PS     , ETAS  , ZPW    )
!
!     COMPUTES THE VERTICAL VELOCITY DIAGNOSTICALLY (W=DZ/DT)
!
      INCLUDE "phykon.h"
      INCLUDE "parorg.h"
      INCLUDE "org.h"

      REAL, INTENT(IN) ::
     &         DAK (KE)   , DBK (KE)   , ACPHIR(JE,2)

      REAL, INTENT(IN) ::
     &         FIB (IE,JE),
     &         DWDT(IE,JE,KE ,3), PINT(IE,JE,KE1,3),
     &         T(IE,JE,KE ,3)  , QI(IE,JE,KE,3),
     &         QD(IE,JE,KE ,3) , QW(IE,JE,KE,3),
     &         PS(IE,JE    ,3) , ETAS(IE,JE,KE1),
     &         U (IE,JE,KE ,3) , V (IE,JE,KE ,3)

      REAL, INTENT(INOUT) :: ZPW(IE,JE,KE)
      REAL, INTENT(INOUT) :: W(IE,JE,KE1,3)

      REAL ::   ADVZ(IE,JE,KE)

      REAL, PARAMETER :: WA    = 0.0         ! ZEITFILTER
      REAL, PARAMETER :: WAM1  = 1. - WA     ! ZEITFILTER
      LOGICAL, PARAMETER ::  LLOG = .FALSE.

      INTEGER :: I, J, K, KP1
      REAL :: TWODT
      REAL ::
     &         Z   (IE,JE,KE1  ) ,
     &         TTB (IE,JE      )

      TWODT = 2.*DT

      ADVZ(:,:,:)   = 0.0

      DO J = JAA , JEA
        DO I = IAA, IEA
          Z(I,J,KE1)    = FIB(I,J)/G
          W(I,J,KE1,NE) = 0.
        END DO
      END DO

      IF (LLOG) THEN
        DO K = KE , 1  , -1
          DO J = JAA , JEA
            DO I = IAA , IEA
              ALP1  = ALOG(PINT(I,J,K+1,NE))
              ALP1O = ALOG(PINT(I,J,K+1,NA))
              IF (K==1) THEN
                ALP1U  = ALOG(2.0)
                ALP1UO = ALOG(2.0)
              ELSE
                ALP1U  = ALOG(PINT(I,J,K,NE))
                ALP1UO = ALOG(PINT(I,J,K,NA))
              ENDIF
              TFAC = 1. + 0.608*QD(I,J,K,NE) - QW(I,J,K,NE)
              DZ = TFAC*T(I,J,K,NE)*R*(ALP1 - ALP1U)/(DWDT(I,J,K,NJ)*G)
              TFAC = 1. + 0.608*QD(I,J,K,NA) - QW(I,J,K,NA)
              DZ2 = TFAC*T(I,J,K,NA)*R*(ALP1O - ALP1UO)
     &                   /(DWDT(I,J,K,NA)*G)
              Z(I,J,K) = DZ + Z(I,J,K+1)
              W(I,J,K,NE) = W(I,J,K+1,NE) + DZ - DZ2
            END DO
          END DO
        END DO
      ELSE
        DO K = KE , 1  , -1
          KP1 = K + 1
          DO J = JAA , JEA
            DO I = IAA , IEA
              TFAC = 1. + 0.608*QD(I,J,K,NE) 
     &                  -(QW(I,J,K,NE) + QI(I,J,K,NE))
              DZ = TFAC*T(I,J,K,NE)*R*2.
     &             /((PINT(I,J,KP1,NE) + PINT(I,J,K,NE))*G)
     &             * (DAK(K) + DBK(K)*PS(I,J,NE))
              TFAC = 1. + 0.608*QD(I,J,K,NA) 
     &                  - (QW(I,J,K,NA) + QI(I,J,K,NA))
              DZ2 = TFAC*T(I,J,K,NA)*R*2.
     &              /((PINT(I,J,KP1,NA) + PINT(I,J,K,NA))*G)
     &              * (DAK(K) + DBK(K)*PS(I,J,NA))
              Z(I,J,K) = DZ + Z(I,J,KP1)
              W(I,J,K,NE) = W(I,J,KP1,NE) + DZ - DZ2
            END DO
          END DO
        END DO
      ENDIF

C-----------------------------------------------------------------------
      DO K = 1 , KE
        KP1 = K + 1
        DO J = JAA , JEA
          DO I = IAA , IEA
            ZPW(I,J,K)  = (Z(I,J,K) + Z(I,J,KP1))*0.5
            W(I,J,K,NE) = (W(I,J,K,NE) + W(I,J,KP1,NE))*0.5/TWODT
          END DO
        END DO
      END DO
C--------------VERTICAL ADVECTION---------------------------------------
      DO J = JAA , JEA
        DO I = IAA , IEA
          TTB(I,J) = 0.
        END DO
      END DO

      DO K = 1 , KE-1
        DO J = JAA , JEA
          DO I = IAA , IEA
            TTAL = (ZPW(I,J,K+1) - ZPW(I,J,K))*ETAS(I,J,K)*0.5
            W(I,J,K,NE) = W(I,J,K,NE)
     &       + (TTAL + TTB(I,J))/(DAK(K) + DBK(K)*PS(I,J,NJ))
            TTB(I,J)    = TTAL
          END DO
        END DO
      END DO

      DO J = JAA , JEA
        DO I = IAA , IEA
          W(I,J,KE,NE) = 
     &    TTB(I,J)/(DAK(KE) + DBK(KE)*PS(I,J,NJ)) + W(I,J,KE,NE)
        END DO
      END DO

      IF( WA /= 0. ) THEN
        DO K = 1 , KE
          DO J = JAA , JEA
            W(IAA:IEA,J,K,NE) = 
     &      W(IAA:IEA,J,K,NA)*WA + W(IAA:IEA,J,K,NE)*WAM1
          END DO
        END DO
      ENDIF

!
!     CALCULATES DIAGNOSTICALLY THE HORIZONTAL ADVECTION OF HEIGHT
!
      CALL ADV_X(  NE    , DAK    , DBK    ,
     &    ACPHIR , U  , V   , PS  , ADVZ   , ZPW   )

      DO K = 1, KE
        DO J = JAH, JEH
          DO I = IAH, IEH
            W(I,J,K,NE) = W(I,J,K,NE) + ADVZ(I,J,K)
          END DO
        END DO
      END DO

      END SUBROUTINE CDZDT
