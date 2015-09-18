      SUBROUTINE GEOPOT(NT, NT2, ILO, IUP, JLO, JUP,
     &    T  , QD, QW  , FI    , FIB ,
     &    QI , PINT    , DWDT)
      !
      IMPLICIT NONE
      !
C
C**** GEOPOT   -   UP:BERECHNUNG DES GEOPOTENTIALS FI AN DEN NEBENFL.
C**   AUFRUF   :   CALL GEOPOT(NT, NT2)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   BERECHNUNG DES GEOPOTENTIALS FI AN DEN NEBENFLAECHEN.
C**   VERSIONS-
C**   DATUM    :   14.02.89
C**                2007
C**
C**   EXTERNALS:   COPYRE
C**
C**   EINGABE-
C**   PARAMETER:   NT:  ZEITINDEX DER PROGNOSTISCHEN VARIABLEN
C**                NT2: ZEITINDEX FUER FI
C**   AUSGABE-
C**   PARAMETER:   KEINE
C**
C**   COMMON-
C**   BLOECKE  :   PHYKON, HIGKON, ORG
C**                COMDYN
C**
C**   METHODE  :   INTEGRATION DER STATISCHEN GRUNDGLEICHUNG
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   D.MAJEWSKI

      INCLUDE "parorg.h"
      INCLUDE "phykon.h"
      INCLUDE "org.h"
      INCLUDE "higkon.h"
      INCLUDE "comdyn.h"
      !
      ! Dummy Arguments
      !
      INTEGER, INTENT(IN) :: NT, NT2, ILO, IUP, JLO, JUP
      !
      REAL, INTENT(IN) :: DWDT(IEJE,KE,3)
      REAL, INTENT(IN) :: FIB(IEJE)
      REAL, INTENT(IN) :: QI(IEJE,KE,3)
      REAL, INTENT(IN) :: QD(IEJE,KE,3)
      REAL, INTENT(IN) :: QW(IEJE,KE,3)
      REAL, INTENT(IN) :: PINT(IEJE,KE1,3)
      REAL, INTENT(IN) :: T(IEJE,KE,3)
      REAL, INTENT(INOUT) :: FI(IEJE,KE,2)
      !
      ! Local Variables
      !
      REAL :: ZPU (IEJE), ZPO (IEJE), ZTV(IEJE)
      REAL :: ZFIU(IEJE), ZFIO(IEJE)
      INTEGER  :: I,J,K,IJ
      REAL :: ALOG2

      ALOG2 = ALOG( 2.0 )
C     BERECHNUNG DES GEOPOTENTIALS FI DURCH INTEGRATION DER HYDROSTATI-
C     SCHEN APPROXIMATION
      DO J = JLO , JUP
        DO I = ILO , IUP
          IJ = I + (J - 1)*IE
          ZFIU(IJ) = FIB(IJ)
          ZPU (IJ) = PINT (IJ,KE1,NT)
        ENDDO
      ENDDO
      DO K  = KE,2, - 1
        DO J = JLO , JUP
          DO I = ILO , IUP
            IJ = I + (J - 1)*IE
            ZPO(IJ)   = PINT(IJ,K,NT)
            ZTV(IJ)   = T(IJ,K,NT)*(1.0 + RDDRM1*QD(IJ,K,NT)
     &                  - (QW(IJ,K,NT) + QI(IJ,K,NT)))
            ZFIO (IJ) = ZFIU(IJ) + R*ZTV(IJ)*ALOG(ZPU(IJ)/ZPO(IJ))
     &                             /DWDT(IJ,K,NT)
            ZPU(IJ)  = ZPO (IJ)
            ZFIU(IJ) = ZFIO(IJ)
            FI(IJ,K,NT2) = ZFIO(IJ)
          ENDDO
        ENDDO
      ENDDO
      DO J = JLO , JUP
        DO I = ILO , IUP
          IJ = I + (J - 1)*IE
          ZTV(IJ) = T(IJ,1,NT)*(1.0 + RDDRM1*QD(IJ,1,NT)
     &              - (QW(IJ,1,NT) + QI(IJ,1,NT)))
          FI(IJ,1,NT2) = FI(IJ,2,NT2) + R*ZTV(IJ)*ALOG2
        ENDDO
      ENDDO

      END SUBROUTINE GEOPOT
