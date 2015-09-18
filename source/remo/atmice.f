      SUBROUTINE ATMICE(KLON, KSTART, KSTOP,  SICED, INFRI)
CTS 250100
      IMPLICIT NONE
C
C**** *ATMICE*  - COMPUTES SEAICE COVER AND DEPTH FOR UNCOUPLED RUNS
C
C      U.SCHLESE       UNI HAMBURG     12.07.89
CDJ    D. JACOB        MPI HAMBURG      1.06.94
C
      INCLUDE "COMPH2"
C
C     Dummy Arguments
C
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
      INTEGER,                  INTENT(IN)    :: KLON, KSTART, KSTOP
      INTEGER, DIMENSION(KLON), INTENT(IN)    :: INFRI
      REAL,    DIMENSION(KLON), INTENT(INOUT) :: SICED
C
C     Local Variables
C
      INTEGER :: JL
C
CTS 250100
C
C*       2.   COMPUTE SEAICE.
C             ------- -------
C
      DO JL=KSTART,KSTOP
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
         IF (INFRI(JL).GT.0) THEN
C
CDJ  SOUTHPOL SICED=1. ----- NORTHPOL SICED=2.
C
C            SICED(JL)=1.
            SICED(JL)=2.
         ELSE
            SICED(JL)=0.
         ENDIF
CTS 250100
      ENDDO
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
CTS SKIN TEMPERATURE NOT NEEDED
CTS 250100
C
      RETURN
      END SUBROUTINE ATMICE
