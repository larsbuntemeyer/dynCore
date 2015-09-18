      INTEGER FUNCTION MDAT2C(KD,KM,KY)
      !
      IMPLICIT NONE
      !
C
C             THIS FUNCTION RETURNS FORECAST DAY
C             STARTING AT 1ST OF JAN  OF FIRST FORECAST-YEAR
C             KD        DAY
C             KM        MONTH
C             KY        YEAR (19XX OR XX)
C
C  *****  VERSION FOR 30-DAY MONTHS , 360-DAYS YEAR  (U. SCHLESE MAR-89)
C----------------------------------------------------------------------
C
      !
      ! Dummy arguments
      !
      INTEGER :: KD,KM,KY
      !
      ! Local Variables
      !
      INTEGER :: IY,IYDAY
      !
      IY=MOD(KY,100)
      IYDAY=(KM-1)*30+KD
      MDAT2C=(IY-1)*360+IYDAY
      !
      RETURN
      !
      END FUNCTION MDAT2C
