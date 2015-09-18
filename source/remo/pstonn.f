      SUBROUTINE PSTONN(AK,BK,PS,T,FIB,PSRED)
C
C     BERECHNUNG DES AUF NN REDUZIERTEN BODENDRUCKES
C
      IMPLICIT NONE
C
      INCLUDE "parorg.h"

      REAL, INTENT(IN)    :: AK(KE+1), BK(KE+1)
      REAL, INTENT(IN)    :: T(IEJE), PS(IEJE), FIB(IEJE)
      REAL, INTENT(INOUT) :: PSRED(IEJE)
C
C     LOKALE FELDER
C
      REAL    :: ZTSTAR(IEJE), ZALPHA(IEJE), ZPHFKE(IEJE)
      REAL    :: G, R, ZT0
      INTEGER :: IJ
C
      R=287.05
      G=9.80665
C
C     DRUCK IM NIVEAU K = KE (UNTERSTE HAUPTFLAECHE) BERECHNEN
C     TEMPERATUR AM BODEN EXTRAPOLIEREN AUS UNTERSTER MODELLFLAECHE
C
      DO IJ  = 1,IEJE
         ZPHFKE(IJ) = 0.5*(AK(KE)+ AK(KE+1)+
     &        (BK(KE)+BK(KE+1))*PS(IJ))
         ZTSTAR(IJ) = T(IJ)+0.0065*R/G*T(IJ)*(PS(IJ)/ZPHFKE(IJ)-1.0)
      ENDDO
C
C     ZALPHA (LAPSE RATE) BERECHNEN
C
      DO IJ  = 1,IEJE
         ZALPHA(IJ) = 0.0065*R/G
         ZT0        = ZTSTAR(IJ) + 0.0065*FIB(IJ)/G
         IF (ZTSTAR(IJ).LE.290.5 .AND. ZT0.GT.290.5) THEN
            ZALPHA(IJ) = R*(290.5 - ZTSTAR(IJ))/FIB(IJ)
         ELSE IF (ZTSTAR(IJ).GT.290.5 .AND. ZT0.GT.290.5) THEN
            ZTSTAR(IJ) = 0.5*(290.5 + ZTSTAR(IJ))
            ZALPHA(IJ) = 0.0
         ELSE IF (ZTSTAR(IJ).LT.255.0) THEN
            ZTSTAR(IJ) = 0.5*(255.0 + ZTSTAR(IJ))
         ENDIF
      ENDDO
C
C     AUF NN REDUZIERTEN BODENDRUCK BERECHNEN
C
      DO IJ  = 1,IEJE
         IF (ABS(FIB(IJ)).GT.0.1) THEN
            PSRED(IJ) = PS(IJ)*       EXP(FIB(IJ)/(R*ZTSTAR(IJ))*
     &             (1.0 - 0.5*(ZALPHA(IJ)*FIB(IJ)/(R*ZTSTAR(IJ))) +
     &                  0.333*(ZALPHA(IJ)*FIB(IJ)/(R*ZTSTAR(IJ)))**2))
         ELSE
            PSRED(IJ) = PS(IJ)
         ENDIF
      ENDDO
C
      RETURN
      END SUBROUTINE PSTONN
