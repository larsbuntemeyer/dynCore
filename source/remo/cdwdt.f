      SUBROUTINE CDWDT(
     &    DAK    , DBK   , 
     &    ACPHIR , DT    , ETAS ,
     &    DWDT   , W     , U      , V      , PS   )

      INCLUDE "org.h"
      INCLUDE "parorg.h"
      INCLUDE "higkon.h"
      INCLUDE "parkon.h"

      REAL, INTENT(IN) :: DT
      REAL, INTENT(IN) ::
     &         U   (IE, JE, KE, 3) , V   (IE, JE, KE, 3) ,
     &         W   (IE, JE,KE1, 3) , PS  (IE, JE,     3) ,
     &         ETAS(IE, JE, KE1  )

      REAL, INTENT(IN) ::
     &         DAK (KE)   , DBK (KE)   , ACPHIR(JE,2) 

      REAL, INTENT(INOUT) :: DWDT(IE,JE,KE,3)

      REAL ::
     &         TDWDT(IE,JE,KE)     , ADVW(IE, JE, KE),
     &         WHELP(IE, JE, KE)

      REAL ::  TWB(IE,JE)

      INTEGER :: I, J, K, KP1
      REAL    :: TWODTI

      TWODTI = .5 / DT        ! 1/(2*DT)
C--------------LOCAL DERIVATIVE OF W------------------------------------
      DO K = 1 , KE
        DO J = JAH, JEH
          DO I = IAH,IEH
            DWDT(I,J,K,NE) = 0.
            TDWDT(I,J,K)   = (W(I,J,K,NE) - W(I,J,K,NA))*TWODTI
          END DO
        END DO
      END DO

C--------------VERTICAL ADVECTION OF W----------------------------------
      TWB(:,:) = 0.

      DO K = 1 , KE-1
        KP1 = K + 1
        DO J = JAH, JEH
          DO I = IAH,IEH
            TWAL = (W(I,J,KP1,NE) - W(I,J,K,NE))*ETAS(I,J,KP1)*0.5
            DWDT(I,J,K,NE) = DWDT(I,J,K,NE)
     &         + (TWAL + TWB(I,J))/(DAK(K) + DBK(K)*PS(I,J,NJ))
            TWB(I,J) = TWAL
          END DO
        END DO
      END DO

      DO J = JAH, JEH
        DO I = IAH,IEH
          DWDT(I,J,KE,NE) = TWB(I,J)/(DAK(KE)+DBK(KE)*PS(I,J,NJ))
     &                    + DWDT(I,J,KE,NE)
        END DO
      END DO

C-----------------------------------------------------------------------
      WHELP(:,:,1:KE) = W(:,:,1:KE,NE)

      CALL ADV_X(  NE    , DAK    , DBK    ,
     &    ACPHIR , U     , V      , PS     , ADVW   , WHELP  )

      DO K = 1 , KE
        DO J = JAH, JEH
          DO I = IAH , IEH
            DWDT(I,J,K,NE) = DWDT(I,J,K,NE) + ADVW(I,J,K) + TDWDT(I,J,K)
          END DO
        END DO
      END DO

      END SUBROUTINE CDWDT
