      SUBROUTINE SULWX
      IMPLICIT NONE
C
C**** *SULW*   - INITIALIZE COMMON YOMLW
C
C     PURPOSE.
C     --------
C           INITIALIZE YOMLW, THE COMMON THAT CONTAINS COEFFICIENTS
C           NEEDED TO RUN THE LONGWAVE RADIATION ROUTINES
C
C**   INTERFACE.
C     ----------
C        *CALL* *SULW
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C        NONE
C
C        IMPLICIT ARGUMENTS :
C        --------------------
C        COMMON YOMLW
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
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ROB VAN DORLAND *KNMI*
C        MODIFIED : 92-01-27  N2O+CH4, GEISA1984
C
C        MARCO GIORGETTA *MPI*
C        MODIFIED : 93-06-22  H2O,N2O,CH4:HITRAN91 ; CO2:AFGL80
C     ------------------------------------------------------------------
      INCLUDE "YOMLW"
C      ----------------------------------------------------------------
C
C*       1.    SET VALUES.
C              -----------
C
      NG1=2
      NG1P1=NG1+1
      NINT=6
      NIPD=8
      NIPD2=2*NIPD
      NTR=11
      NTRA=15
      NUA=29
      PVGCO2=60.
      PVGH2O=30.
      PVGO3 =400.
      RETURN
      END SUBROUTINE SULWX
