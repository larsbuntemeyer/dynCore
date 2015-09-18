      BLOCK DATA SUSW
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
      INTEGER :: I,J,K
C      ----------------------------------------------------------------
C
C*       1.    SET VALUES.
C              -----------
C
      DATA RTDH2O,RTDUMG / 0.40  , 0.375 /
      DATA RTH2O ,RTUMG  / 240.  , 240.  /
      DATA RSWCE ,RSWCP                / 0.000 , 0.000 /
C
      DATA SUN(1) / 0.441676 /
      DATA (D(1,K),K = 1,3) / 0.00, 0.00, 0.00 /
C* DERIVED FROM HITRAN APRIL 1991
C       H2O:  PREF=300 HPA, TREF=240K, PDEP=0.8
C       O3 :  UNCHANGED
C
      DATA ((APAD(1,I,J),I=1,3),J=1,7) /
     & 0.912418292E+05, 0.000000000E-00, 0.925887084E-04,
     & 0.723613782E+05, 0.000000000E-00, 0.129353723E-01,
     & 0.596037057E+04, 0.000000000E-00, 0.800821928E+00,
     & 0.000000000E-00, 0.000000000E-00, 0.242715973E+02,
     & 0.000000000E-00, 0.000000000E-00, 0.878331486E+02,
     & 0.000000000E-00, 0.000000000E-00, 0.191559725E+02,
     & 0.000000000E-00, 0.000000000E-00, 0.000000000E+00 /
C
      DATA ((BPAD(1,I,J),I=1,3),J=1,7) /
     & 0.912418292E+05, 0.000000000E-00, 0.925887084E-04,
     & 0.724555318E+05, 0.000000000E-00, 0.131812683E-01,
     & 0.602593328E+04, 0.000000000E-00, 0.812706117E+00,
     & 0.100000000E+01, 0.000000000E-00, 0.249863591E+02,
     & 0.000000000E-00, 0.000000000E-00, 0.931071925E+02,
     & 0.000000000E-00, 0.000000000E-00, 0.252233437E+02,
     & 0.000000000E-00, 0.000000000E-00, 0.100000000E+01 /
C
C
      DATA (CRAY(1,K),K=1,6) /
     & .428937E-01, .890743E+00,-.288555E+01,
     & .522744E+01,-.469173E+01, .161645E+01/
C
      DATA SUN(2) / 0.558324 /
C
      DATA (D(2,K),K=1,3) / 0.000000000, 0.000000000, 0.800000000 /
C
C* INTERVAL 2:  0.68 - 4.00 MICRONS
C* DERIVED FROM HITRAN APRIL 1991
C       H2O:  PREF=300 HPA, TREF=240K, PDEP=0.80
C       UMG:  PREF=300 HPA, TREF=240K, PDEP=0.75 (CO2+O2+CH4+N2O+CO)
C       O3 :  UNCHANGED
C
      DATA ((APAD(2,I,J),I=1,3),J=1,7) /
     & 0.376655383E-08, 0.739646016E-08, 0.410177786E+03,
     & 0.978576773E-04, 0.131849595E-03, 0.672595424E+02,
     & 0.387714006E+00, 0.437772681E+00, 0.000000000E-00,
     & 0.118461660E+03, 0.151345118E+03, 0.000000000E-00,
     & 0.119079797E+04, 0.233628890E+04, 0.000000000E-00,
     & 0.293353397E+03, 0.797219934E+03, 0.000000000E-00,
     & 0.000000000E+00, 0.000000000E+00, 0.000000000E+00 /
C
      DATA ((BPAD(2,I,J),I=1,3),J=1,7) /
     & 0.376655383E-08, 0.739646016E-08, 0.410177786E+03,
     & 0.979023421E-04, 0.131861712E-03, 0.731185438E+02,
     & 0.388611139E+00, 0.437949001E+00, 0.100000000E+01,
     & 0.120291383E+03, 0.151692730E+03, 0.000000000E+00,
     & 0.130531005E+04, 0.237071130E+04, 0.000000000E+00,
     & 0.415049409E+03, 0.867914360E+03, 0.000000000E+00,
     & 0.100000000E+01, 0.100000000E+01, 0.000000000E+00 /
C
C
      DATA (CRAY(2,K),K=1,6) /
     & .697200E-02, .173297E-01,-.850903E-01,
     & .248261E+00,-.302031E+00, .129662E+00/
CRP   NOW IN SUSWX
C     RPDH1=ZPDH2O+1.
C     RPDU1=ZPDUMG+1.
C     ZH2O=1./( 10.* G * RPDH1 )
C     ZUMG=1./( 10.* G * RPDU1 )
C     RPNU = ZUMG/(ZPRUMG**ZPDUMG)
C     RPNH = ZH2O/(ZPRH2O**ZPDH2O)
C
C     RSWCP=0.002*RSWCE
C
      END BLOCK DATA SUSW
