      SUBROUTINE SUSWX
      IMPLICIT NONE
C
C**** *SUSW*   - INITIALIZE COMMON YOMSW
C
C     PURPOSE.
C     --------
C           INITIALIZE YOMSW, THE COMMON THAT CONTAINS COEFFICIENTS
C           NEEDED TO RUN THE SHORTWAVE RADIATION SUBROUTINES
C
C**   INTERFACE.
C     ----------
C        *CALL* *SUSW
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C        NONE
C
C        IMPLICIT ARGUMENTS :
C        --------------------
C        COMMON YOMSW
C
C     METHOD.
C     -------
C        SEE DOCUMENTATION
C
C     EXTERNALS.
C     ----------
C
C     REFERENCE.
C     ----------
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE
C     "IN CORE MODEL"
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 88-12-15
C     ------------------------------------------------------------------
      INCLUDE "COMCON"
      INCLUDE "YOMSW"
      REAL :: ZH2O, ZPDH2O, ZPDUMG, ZPRH2O, ZPRUMG, ZUMG
CRP
      DATA ZPDH2O,ZPDUMG / 0.8   ,  0.75 /
      DATA ZPRH2O,ZPRUMG / 30000., 30000./
C      ----------------------------------------------------------------
C
C*       1.    SET VALUES.
C              -----------
C
C
      RPDH1=ZPDH2O+1.
      RPDU1=ZPDUMG+1.
      ZH2O=1./( 10.* G * RPDH1 )
      ZUMG=1./( 10.* G * RPDU1 )
      RPNU = ZUMG/(ZPRUMG**ZPDUMG)
      RPNH = ZH2O/(ZPRH2O**ZPDH2O)
C
      RSWCP=0.002*RSWCE
C
      RETURN
      END
