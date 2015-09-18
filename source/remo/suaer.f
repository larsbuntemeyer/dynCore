      BLOCK DATA SUAER
      IMPLICIT NONE
C
C**** *SUAER*   - INITIALIZE COMMON YOMAER
C
C     PURPOSE.
C     --------
C           INITIALIZE YOMAER, THE COMMON THAT CONTAINS THE
C           RADIATIVE CHARACTERISTICS OF THE AEROSOLS
C
C**   INTERFACE.
C     ----------
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C        NONE
C
C        IMPLICIT ARGUMENTS :
C        --------------------
C        COMMON YOMAER
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
C        ORIGINAL : 88-02-15
C     ------------------------------------------------------------------
      INCLUDE "YOMAER"
C      ----------------------------------------------------------------
      INTEGER :: JA,IN
C
C*       1.    SHORTWAVE COEFFICIENTS
C              ----------------------
C
      DATA ((TAUA(IN,JA),JA=1,5),IN=1,2) /
     & .730719, .912819, .725059, .745405, .682188 ,
     & .730719, .912819, .725059, .745405, .682188 /
      DATA ((PIZA(IN,JA),JA=1,5),IN=1,2) /
     & .872212, .982545, .623143, .944887, .997975 ,
     & .872212, .982545, .623143, .944887, .997975 /
      DATA ((CGA (IN,JA),JA=1,5),IN=1,2) /
     & .647596, .739002, .580845, .662657, .624246 ,
     & .647596, .739002, .580845, .662657, .624246 /
C      ----------------------------------------------------------------
C
C*       2.    LONGWAVE COEFFICIENTS
C              ---------------------
C
      DATA CAER / .038520, .037196, .040532, .054934, .038520
     &          , .12613 , .18313 , .10357 , .064106, .126130
     &          , .012579, .013649, .018652, .025181, .012579
     &          , .011890, .016142, .021105, .028908, .011890
     &          , .013792, .026810, .052203, .066338, .013792 /
C
      END BLOCK DATA SUAER
