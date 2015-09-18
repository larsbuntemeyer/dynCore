C
C
C**** *PRESF* *COMPUTES FULL-LEVEL PRESSURES.
C
C          M.JARRAUD     E.C.M.W.F.     7/5/1982.
C
C     PURPOSE
C     -------
C
C          *TO COMPUTE FULL-LEVEL PRESSURES FROM HALF-LEVEL VALUES.
C
C**   INTERFACE
C     ---------
C
C          *PRESF* IS CALLED FROM *PHECHAM*. PARAMETERS ARE:
C             *PF*        *COMPUTED FULL-LEVEL PRESSURES.
C             *PH*        *HALF-LEVEL PRESSURES.
C             *KDIMP*     *FIRST DIMENSION OF 2-D ARRAYS *PF* AND *PH*
C
C
C     METHOD
C     ------
C
C          *FULL-LEVEL PRESSURES ARE DEFINED AS THE ARITHMETIC
C     AVERAGE OF THE TWO ADJOINING HALF-LEVEL PRESSURES.
C
C     REFERENCE
C     ---------
C
C
      SUBROUTINE PRESF(KDIMP,KSTART,KSTOP,NLEV,PF,PH)
C
      IMPLICIT NONE
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: KDIMP,KSTART,KSTOP,NLEV
      REAL,    INTENT(IN)    :: PH(KDIMP,NLEV+1)
      REAL,    INTENT(INOUT) :: PF(KDIMP,NLEV)
C
C     Local Variables
C
      INTEGER :: JLEV,JL
C
!DIR$ NOBOUNDS
C     ------------------------------------------------------------
C
C*        1.    COMPUTE FULL-LEVEL PRESSURE VALUES.
C               ------- ---- ----- -------- -------
C
      DO JLEV=1,NLEV
C
         DO JL=KSTART,KSTOP
            PF(JL,JLEV)=(PH(JL,JLEV)+PH(JL,JLEV+1))*.5
         ENDDO
C
      ENDDO
C
C     ------------------------------------------------------------
C
      RETURN
      END SUBROUTINE PRESF
