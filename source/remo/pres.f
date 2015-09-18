C
C
C**** *PRES* - *CALCULATES HALF-LEVEL PRESSURES.
C
C     A.J.SIMMONS     E.C.M.W.F.    16/11/81.
C     MODIFIED BY:
C     R.PODZUN        D.K.R.Z.      25/08/94.
C
C     PURPOSE
C     -------
C
C             *TO CALCULATE HALF-LEVEL PRESSURES AT ALL MODEL LEVELS
C     FOR A GIVEN SURFACE PRESSURE.
C
C**   INTERFACE
C     ---------
C
C             *THIS SUBROUTINE IS CALLED FROM MANY POINTS WITHIN THE
C     FORECASTING SYSTEM. PARAMETERS ARE:
C
C             *PH*        *COMPUTED HALF-LEVEL PRESSURES.
C             *KLEN*     *FIRST DIMENSION OF 2-D ARRAY *PH.*
C             *PS*        *SURFACE PRESSURE.
C             *KLEN*      *NUMBER OF POINTS FOR WHICH CALCULATION IS
C                         PERFORMED.
C
C     RESULTS
C     -------
C
C             *RESULTS ARE COMPUTED FOR *KLEN* CONSECUTIVE POINTS AT
C     EACH MODEL HALF LEVEL.
C
C     METHOD
C     ------
C
C             *CALCULATIONS ARE PERFORMED FOR ALL LEVELS
C
C     REFERENCE
C     ---------
C
C             *EXTERNAL DOCUMENTATION OF THE MODEL EQUATIONS AND THE
C     ORGANIZATION OF THE VERTICAL CALCULATION.
C
C
      SUBROUTINE PRES(KLEN,KSTART,KSTOP,NLEVP1,VCT,PH,PS)
C
      IMPLICIT NONE
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: KLEN,KSTART,KSTOP,NLEVP1
      REAL,    INTENT(IN)    :: PS(KLEN),VCT(2*NLEVP1)
      REAL,    INTENT(INOUT) :: PH(KLEN,NLEVP1)
C
C     Local Variables
C
      INTEGER :: JK,JL
      REAL    :: ZP,ZB
C
C     ------------------------------------------------------------------
!DIR$ NOBOUNDS
C
C
      DO JK=1,NLEVP1
         ZP=VCT(JK)
         ZB=VCT(JK+NLEVP1)
C
         DO JL=KSTART,KSTOP
            PH(JL,JK)=ZP+ZB*PS(JL)
         ENDDO
C
      ENDDO
C
C     ------------------------------------------------------------------
C
      RETURN
      END SUBROUTINE PRES
