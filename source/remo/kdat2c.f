      INTEGER FUNCTION KDAT2C(KD,KM,KY)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: KD, KM, KY
C
C             THIS FUNCTION RETURNS CENTURY DAY GIVEN
C             KD        DAY
C             KM        MONTH
C             KY        YEAR (19XX OR XX)
C
C----------------------------------------------------------------------
C
      INTEGER :: IDAYS(12),JDAYS(12)
      INTEGER :: IY,ILY,IYDAY
C
      DATA IDAYS/0,31,59,90,120,151,181,212,243,273,304,334/
      DATA JDAYS/0,31,60,91,121,152,182,213,244,274,305,335/
      IY=MOD(KY,100)
      IYDAY=IDAYS(KM)+KD
      IF (MOD(IY,4).EQ.0) IYDAY=JDAYS(KM)+KD
      ILY=(IY-1)/4
      KDAT2C=IY*365+ILY+IYDAY
      RETURN
      END FUNCTION KDAT2C
