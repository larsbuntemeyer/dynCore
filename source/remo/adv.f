      SUBROUTINE ADV_X(    NT     , DAK    , DBK  ,
     &    ACPHIR , U     , V      , PS     , ADVT , ADV   )
C
C     COMPUTES THE HORIZONTAL ADVECTION OF ADVT
C
      INCLUDE "higkon.h"
      INCLUDE "org.h"
      INCLUDE "parorg.h"
C
      REAL, INTENT(INOUT) :: ADVT(IE,JE,KE)
      REAL, INTENT(IN)    ::
     &         U   (IE,JE,KE ,3) , V   (IE,JE,KE ,3) ,
     &         PS  (IE,JE    ,3) ,
     &         ADV (IE,JE,KE)
C
      REAL, INTENT(IN)    ::
     &         DAK (KE)   , DBK (KE)   , ACPHIR(JE,2)
C
      DIMENSION
     &         ZGU (IE, KE, 2) , ZGV (IE, KE, 3) ,
     &         ZDP (IE, KE, 5) , ZEW (IE, KE)    ,
     &         ZSN (IE, KE, 5) , ZEWF(IE, KE)    ,
     &         ZSNF(IE, KE, 3) , RDX (5)
C
      REAL, PARAMETER :: WCA = 1.
      REAL, PARAMETER :: WCB = 1. - WCA
      REAL, PARAMETER :: WCC = .25*WCB

      JU2 = 1
      JM2 = 2

      JU3 = 1
      JM3 = 2
      JO3 = 3

      JS5 = 1
      JU5 = 2
      JM5 = 3
      JO5 = 4
      JN5 = 5

      JZS = JAH - 2
      JZU = JAH - 1
      JZM = JAH
      JZO = JAH + 1
      JZN = JAH + 2

      RDX(JS5) = ACPHIR(JZS,1)*EDDLAM
      RDX(JU5) = ACPHIR(JZU,1)*EDDLAM
      RDX(JM5) = ACPHIR(JZM,1)*EDDLAM
      RDX(JO5) = ACPHIR(JZO,1)*EDDLAM
      RDY      = EDADPHI 

      DO K = 1 , KE
        DO I = IAA , IEA
          ZSNF(I,K,JU3) = 0.
          ZSNF(I,K,JM3) = 0.
          ZSNF(I,K,JO3) = 0.
          ZEWF(I,K)     = 0.

          ZDP(I,K,JS5) = DAK(K) + DBK(K)*PS(I,JZS,NT)
          ZDP(I,K,JU5) = DAK(K) + DBK(K)*PS(I,JZU,NT)
          ZDP(I,K,JM5) = DAK(K) + DBK(K)*PS(I,JZM,NT)
          ZDP(I,K,JO5) = DAK(K) + DBK(K)*PS(I,JZO,NT)
        END DO
      END DO

      DO K = 1 , KE
        DO I = IAA , IEH+1
          ZGU(I,K,JU2) =
     &       0.5 * (ZDP(I+1,K,JU5)+ZDP(I,K,JU5)) * U(I,JZU,K,NT)
          ZGV(I,K,JU3) =
     &       0.5 * (ZDP(I,K,JU5)+ZDP(I,K,JM5))   * V(I,JZU,K,NT)
          ZGV(I,K,JM3) =
     &       0.5 * (ZDP(I,K,JM5)+ZDP(I,K,JO5))   * V(I,JZM,K,NT)

          ZSN(I,K,JS5) =
     &       0.5 * (ZDP(I,K,JS5) + ZDP(I,K,JU5)) * V(I,JZS,K,NT)
     &           * (ADV(I,JZU,K) - ADV(I,JZS,K)) *.5*RDY
          ZSN(I,K,JU5) =
     &      ZGV(I,K,JU3) * (ADV(I,JZM,K) - ADV(I,JZU,K)) *.5*RDY
          ZSN(I,K,JM5) =
     &      ZGV(I,K,JM3) * (ADV(I,JZO,K) - ADV(I,JZM,K)) *.5*RDY

          ZSNF(I,K,JU2)= ZSN(I,K,JU5)*WCA
     &       + WCC*(ZSN(I,K,JS5)-2.*ZSN(I,K,JU5)+ZSN(I,K,JM5))
        END DO
      END DO

      J_LOOP :
     & DO J = JAH, JEH
        JZN = J + 2
        RDX(JN5) = ACPHIR(JZN,1)*EDDLAM

        DO K = 1 , KE
          DO I = IAA , IEA
            ZDP(I,K,JN5) = DAK(K) + DBK(K)*PS(I,JZN,NT)
          END DO
        END DO

        DO K = 1 , KE
          DO I = IAA , IEH+1
            ZGU(I,K,JM2) =
     &         0.5 * (ZDP(I+1,K,JM5) + ZDP(I,K,JM5)) * U(I,J,K,NT)
          END DO
CKS          DO I = IAH , IEH
          DO I = IAH - 1, IEH + 1
            ZGV(I,K,JO3) =
     &         0.5 * (ZDP(I,K,JO5) + ZDP(I,K,JN5)) * V(I,J+1,K,NT)
          END DO
        END DO

        DO K = 1 , KE
          DO I = IAA , IEH+1
            ZEW(I,K) = 
     &        ZGU(I,K,JM2)*(ADV(I+1,J,K) - ADV(I,J,K))  *.5*RDX(JM5)
          END DO

CKS          DO I = IAH , IEH Eventuell ohne +1
          DO I = IAH - 1, IEH + 1
            ZSN(I,K,JO5) =
     &         ZGV(I,K,JO3)*(ADV(I,JZN,K) - ADV(I,J+1,K))*.5*RDY
          END DO
        END DO

        IF( WCC == 0. ) THEN
          DO K = 1 , KE
            ZEWF(IAH-1:IEH,K) = ZEW(IAH-1:IEH,K)
            ZSNF(IAH-1:IEH,K,JM2) = ZSN(IAH-1:IEH,K,JM5)
          END DO
        ELSE
          DO K = 1 , KE
            DO I = IAH-1 , IEH
            ZEWF(I,K) = ZEW(I,K)*WCA
     &                + WCC*(ZEW(I-1,K)-2.*ZEW(I,K)+ZEW(I+1,K))
            END DO

CKS            DO I = IAH , IEH
            DO I = IAH - 1, IEH
            ZSNF(I,K,JM2) = ZSN(I,K,JM5)*WCA
     &        + WCC*(ZSN(I,K,JU5)-2.*ZSN(I,K,JM5)+ZSN(I,K,JO5))
            END DO
          END DO
        ENDIF

        DO K = 1   , KE
          DO I = IAH , IEH
            ADVT(I,J,K) = .5*(ZEWF(I-1,K) + ZEWF(I,K)
     &                        + ZSNF(I,K,JU2) + ZSNF(I,K,JM2) )
     &                    /ZDP(I,K,JM5)
          END DO
        END DO

        JA  = JU2
        JU2 = JM2
        JM2 = JA

        JA  = JU3
        JU3 = JM3
        JM3 = JA
      
        JA  = JS5
        JS5 = JU5
        JU5 = JM5
        JM5 = JO5
        JO5 = JN5
        JN5 = JA

      END DO J_LOOP

      END SUBROUTINE ADV_X
