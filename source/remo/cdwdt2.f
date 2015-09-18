      SUBROUTINE CDWDT2( DIFW  , DIFD   , DWDT   , W )
C
      INCLUDE "org.h"
      INCLUDE "parorg.h"
      INCLUDE "higkon.h"
      INCLUDE "parkon.h"
      INCLUDE "comdyn.h"
C
      REAL, INTENT(IN)    :: W(IE,JE,KE1,3) 
      REAL, INTENT(IN)    :: DWDT(IE,JE,KE,3) 
      REAL, INTENT(INOUT) :: DIFW(IE,JE,KE)
      REAL, INTENT(INOUT) :: DIFD(IE,JE,KE)

      INTEGER :: KM1
      REAL    :: AK0
      
      AK0 = .5 * AKS2/DT

      DO K = 1 , KE
        DO J = JAH, JEH
          DO I = IAH , IEH
            DIFW(I,J,K) =
     &       (W(I-1,J,K,NE) - 2.*W(I,J,K,NE) + W(I+1,J,K,NE))
     &                    * AK0
            DIFW(I,J,K) =  DIFW(I,J,K)
     &      + (W(I,J-1,K,NE) - 2.*W(I,J,K,NE) + W(I,J+1,K,NE))
     &                    * AK0
          END DO
        END DO
      END DO

      DO K = 2 , KE-1
        DO J = JAH, JEH
          DO I = IAH , IEH
            DIFW(I,J,K) =  DIFW(I,J,K)
     &      + (W(I,J,K-1,NE) - 2.*W(I,J,K,NE) + W(I,J,K+1,NE))
     &                    * AK0
          END DO
        END DO
      END DO

      DO K = 1 , KE
        DO J = JAH, JEH
          DO I = IAH , IEH
            DIFD(I,J,K) = DWDT(I-1,J,K,NE) - 2.*DWDT(I,J,K,NE)
     &                  + DWDT(I+1,J,K,NE)
            DIFD(I,J,K) = DWDT(I  ,J-1,K,NE) - 2.*DWDT(I,J,K,NE)
     &                  + DWDT(I  ,J+1,K,NE)
     &                  + DIFD(I,J,K)
          END DO
        END DO
      END DO

      DO K = 2 , KE-1
        KM1 = K - 1
        DO J = JAH, JEH
          DO I = IAH , IEH
            DIFD(I,J,K) = DWDT(I,J,KM1,NE)
     &            -2.*DWDT(I,J,K,NE) + DWDT(I,J,K+1,NE)
     &                  + DIFD(I,J,K)
          END DO
        END DO
      END DO

      END SUBROUTINE CDWDT2
