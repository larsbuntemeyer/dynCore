C
C**** *SWTT* - COMPUTES THE SHORTWAVE TRANSMISSION FUNCTIONS
C
C     PURPOSE.

C     --------
C           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
C     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN THE TWO SPECTRAL
C     INTERVALS.
C
C**   INTERFACE.
C     ----------
C          *SWTT* IS CALLED FROM *SW*.
C
C     *CALL*     SWTT (KLON,KNU,KA,PU,PTR)
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C KNU    :                     ; INDEX OF THE SPECTRAL INTERVAL
C KA     :                     ; INDEX OF THE ABSORBER
C PU     : (KLON)             ; ABSORBER AMOUNT
C     ==== OUTPUTS ===
C PTR    : (KLON)             ; TRANSMISSION FUNCTION
C
C        IMPLICIT ARGUMENTS :   NONE
C        --------------------
C
C     METHOD.
C     -------
C
C          TRANSMISSION FUNCTION ARE COMPUTED USING PADE APPROXIMANTS
C     AND HORNER'S ALGORITHM.
C
C     EXTERNALS.
C     ----------
C
C          NONE
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "IN CORE MODEL"
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 88-12-15
C-----------------------------------------------------------------------
C
      SUBROUTINE SWTT (KLON,KNU,KA,PU,PTR)
C
      IMPLICIT NONE
C
      INCLUDE "YOMSW"
C-----------------------------------------------------------------------
C
C*       0.1   DUMMY ARGUMENTS
C              ---------------
C
      INTEGER, INTENT(IN)    :: KLON,KNU,KA
      REAL,    INTENT(IN)    :: PU(KLON)
      REAL,    INTENT(INOUT) :: PTR(KLON)
C
C     Local Variables
C
      REAL    :: ZR1,ZR2
      INTEGER :: JL
C
C-----------------------------------------------------------------------
C
C*         1.      HORNER'S ALGORITHM TO COMPUTE TRANSMISSION FUNCTION
C
      DO JL = 1 , KLON
         ZR1     = APAD(KNU,KA,1) + PU(JL) * (APAD(KNU,KA,2) + PU(JL)
     &         * ( APAD(KNU,KA,3) + PU(JL) * (APAD(KNU,KA,4) + PU(JL)
     &         * ( APAD(KNU,KA,5) + PU(JL) * (APAD(KNU,KA,6) + PU(JL)
     &         * ( APAD(KNU,KA,7) ))))))
C
         ZR2     = BPAD(KNU,KA,1) + PU(JL) * (BPAD(KNU,KA,2) + PU(JL)
     &         * ( BPAD(KNU,KA,3) + PU(JL) * (BPAD(KNU,KA,4) + PU(JL)
     &         * ( BPAD(KNU,KA,5) + PU(JL) * (BPAD(KNU,KA,6) + PU(JL)
     &         * ( BPAD(KNU,KA,7) ))))))
C
C*         2.      ADD THE BACKGROUND TRANSMISSION
C
         PTR(JL) = (ZR1 / ZR2) * (1. - D(KNU,KA)) + D(KNU,KA)
      ENDDO
C
      RETURN
      END SUBROUTINE SWTT
