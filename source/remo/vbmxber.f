C
C     SUBROUTINE VBMXBER
C
C**** PROGAUS  -   UP:BERECHNUNG DER MAXIMALEN WINDGESCHWINDIGKEIT
C**   AUFRUF   :   CALL VBMXBER(U,V)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   BERECHNUNG DER MAXIMALEN WINDGESCHWINDIGKEIT
C**
C**   VERSIONS-
C**   DATUM    :   15.01.01
C**                2007
C**
C**   EXTERNALS:   KEINE
C**
C**   EINGABE-
C**   PARAMETER:   U, V
C**   AUSGABE-
C**   PARAMETER:   KEINE  (VBMXV IN COMDYN)
C**
C**   COMMON-
C**   BLOECKE  :   ORG, PARORG, COMDYN
C**
C**   METHODE  :   BERECHNUNG DER MAXIMALEN WINDGESCHWINDIGKEIT
C**                ZU JEDEM ZEITSCHRITT
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   R.PODZUN/T.DIEHL
C
      SUBROUTINE VBMXBER(U, V)
C
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
      INCLUDE "org.h"
      INCLUDE "comdyn.h"
C
C     Dummy Arguments
C
      REAL, DIMENSION(IE,JE,KE,3), INTENT(IN) :: U,V
C
C     Local Variables
C
      REAL    :: V2
      REAL    :: ZVB,ZVBMAX,ZU,ZV
      INTEGER :: NT,K,J,I
C
      V2     = 0.
      ZVBMAX = 0.
      NT = NE
C
      DO K = 1  ,KE
         DO J = JAH,JEH
            DO I = IAH,IEH
               ZU=0.5*(U(I,J,K,NT) + U(I-1,J  ,K,NT))
               ZV=0.5*(V(I,J,K,NT) + V(I  ,J-1,K,NT))
               ZVB      = SQRT( ZU*ZU + ZV*ZV )
               ZVBMAX   = MAX(ZVBMAX, ZVB)
            ENDDO
         ENDDO
      ENDDO

      COUNT = 1
      CALL PALLREDUCER(ZVBMAX,MPI_MAX)

      COUNT = 1
      CALL PALLREDUCER(V2,MPI_SUM)

      VBMXV = ZVBMAX

      END SUBROUTINE VBMXBER
