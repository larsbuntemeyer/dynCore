C      SUBROUTINE SURADI(PCO2FAC)
C
C**** *SURADI*   - INITIALIZE COMMON YOMRDI CONTROLLING RADINT
C
C     PURPOSE.
C     --------
C           INITIALIZE YOMRDI, THE COMMON THAT CONTROLS THE
C           RADIATION INTERFACE
C
C**   INTERFACE.
C     ----------
C        *CALL* *SURADI
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C        NONE
C
C        IMPLICIT ARGUMENTS :
C        --------------------
C        COMMON YOMRDI
C
C     METHOD.
C     -------
C        SEE DOCUMENTATION
C
C     EXTERNALS.
C     ----------
C        NONE
C
C     REFERENCE.
C     ----------
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE
C     "IN CORE MODEL"
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 88-12-15
C     ------------------------------------------------------------------
      SUBROUTINE SURADI(PCO2FAC)
      IMPLICIT NONE
      INCLUDE "YOMRDI"
C
      REAL, INTENT(IN) :: PCO2FAC
      REAL :: ZAIRMWG, ZCH4MWG, ZCO2MWG, ZN2OMWG
C      ----------------------------------------------------------------
C
C*       1.    SET DEFAULT VALUES.
C              -------------------
C
      CSDTSN  = 5.0
C
      ZALBICE = 0.55
      ZALBSEA = 0.07
      ZALBSNO = 0.80
      ZALBSNM = 0.40
      ZSNOWAL = 0.01
      ZEMISS  = 0.996
      ZVLBDC  = 0.5
C
C* CONCENTRATION OF VARIOUS TRACE GASES (IPCC/SACC VALUES FOR 1990)
C    CO2     CH4       N2O     CFC11    CFC12    HCFC113  CFC114
C  353PPMV 1.72PPMV  310PPBV  280PPTV  484PPTV   60PPTV   15PPTV
C  CFC115   HCFC22   HCFC123  HCFC124   HFC125   HFC134A HCFC141B
C    5PPTV  122PPTV     ?        ?        ?         ?       ?
C HCFC142B  HFC143A  HFC152A   CCL4     H3CCL3
C    ?        ?         ?     146PPTV     ?
C
      ZAIRMWG = 28.970
      ZCO2MWG = 44.011
      ZCH4MWG = 16.043
      ZN2OMWG = 44.013
C
      ZCARDI  = 353.E-06 * ZCO2MWG/ZAIRMWG * PCO2FAC
      ZMETHA  = 1.72E-06 * ZCH4MWG/ZAIRMWG
      ZNITOX  = 310.E-09 * ZN2OMWG/ZAIRMWG
      ZCFC( 1)= 280.E-12
      ZCFC( 2)= 484.E-12
      ZCFC( 3)=  60.E-12
      ZCFC( 4)=  15.E-12
      ZCFC( 5)=   5.E-12
      ZCFC( 6)= 122.E-12
      ZCFC( 7)=   0.
      ZCFC( 8)=   0.
      ZCFC( 9)=   0.
      ZCFC(10)=   0.
      ZCFC(11)=   0.
      ZCFC(12)=   0.
      ZCFC(13)=   0.
      ZCFC(14)=   0.
      ZCFC(15)= 146.E-12
      ZCFC(16)=   0.
C
      ZEPSEC=1.E-12
      ZEPCLC=1.E-12
      ZEPH2O=1.E-12
      ZEPALB=1.E-12
C
C     -----------------------------------------------------------------
C
      RETURN
      END SUBROUTINE SURADI
