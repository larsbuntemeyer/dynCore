C
C     SUBROUTINE SWDE
C
C**** *SWDE* - DELTA-EDDINGTON IN A CLOUDY LAYER
C
C     PURPOSE.
C     --------
C           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY OF A CLOUDY
C     LAYER USING THE DELTA-EDDINGTON'S APPROXIMATION.
C
C**   INTERFACE.
C     ----------
C          *SWDE* IS CALLED BY *SWR*, *SW2S*
C
C     *CALL*     SWDE (KLON,PGG,PREF,PRMUZ,PTO1,PW
C    S                ,      PRE1,PRE2,PTR1,PTR2         )
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C PGG    : (NDLON)             ; ASSYMETRY FACTOR
C PREF   : (NDLON)             ; REFLECTIVITY OF THE UNDERLYING LAYER
C PRMUZ  : (NDLON)             ; COSINE OF SOLAR ZENITH ANGLE
C PTO1   : (NDLON)             ; OPTICAL THICKNESS
C PW     : (NDLON)             ; SINGLE SCATTERING ALBEDO
C     ==== OUTPUTS ===
C PRE1   : (NDLON)             ; LAYER REFLECTIVITY ASSUMING NO
C                              ; REFLECTION FROM UNDERLYING LAYER
C PTR1   : (NDLON)             ; LAYER TRANSMISSIVITY ASSUMING NO
C                              ; REFLECTION FROM UNDERLYING LAYER
C PRE2   : (NDLON)             ; LAYER REFLECTIVITY ASSUMING
C                              ; REFLECTION FROM UNDERLYING LAYER
C PTR2   : (NDLON)             ; LAYER TRANSMISSIVITY ASSUMING
C                              ; REFLECTION FROM UNDERLYING LAYER
C
C        IMPLICIT ARGUMENTS :   NONE
C        --------------------
C
C     METHOD.
C     -------
C
C          STANDARD DELTA-EDDINGTON LAYER CALCULATIONS.
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
C     ------------------------------------------------------------------
C
      SUBROUTINE SWDE (KLON,PGG,PREF,PRMUZ,PTO1,PW
     &                ,      PRE1,PRE2,PTR1,PTR2         )
c
      IMPLICIT NONE
C
C*       0.1   ARGUMENTS
C              ---------
C
      INTEGER, INTENT(IN)    ::  KLON
      REAL,    INTENT(IN)    ::  PGG(KLON),PREF(KLON),PRMUZ(KLON),
     &                           PTO1(KLON),PW(KLON)
      REAL,    INTENT(INOUT) ::  PRE1(KLON),PRE2(KLON),
     &                           PTR1(KLON),PTR2(KLON)
C
C     LOCAL VARIABLES
C
      INTEGER :: JL
      REAL    :: ZFF,ZGP,ZTOP,ZEXMU0,ZEXKP,ZDT,ZEXKM,ZDENB,ZDENA,ZC2B,
     &           ZC2A,ZC1B,ZC1A,ZRM2,ZRK,ZRI1C,ZRI1B,ZRI1A,ZRI0D,ZRI0C,
     &           ZRI0B,ZRI0A,ZXP2P,ZXM2P,ZX1,ZX2,ZWM,ZWCP,ZRP,ZRI1D,
     &           ZBETA,ZB23,ZB22,ZB21,ZARG2,ZARG,ZAP2B,ZALPHA,ZA11,ZA12,
     &           ZA23,ZA22,ZAM2B,ZA21,ZA13
C
C     ------------------------------------------------------------------
C
C*       0.2   FUNCTIONS
C              ---------
C
C
C     ------------------------------------------------------------------
C
C*         1.      DELTA-EDDINGTON CALCULATIONS
C
      DO JL   =   1 , KLON
C
C*         1.1     SET UP THE DELTA-MODIFIED PARAMETERS
C
         ZFF = PGG(JL)*PGG(JL)
         ZGP = PGG(JL)/(1.+PGG(JL))
         ZTOP = (1.- PW(JL) * ZFF) * PTO1(JL)
         ZWCP = (1-ZFF)* PW(JL) /(1.- PW(JL) * ZFF)
         ZDT = 2./3.
         ZX1 = 1.-ZWCP*ZGP
         ZWM = 1.-ZWCP
         ZRM2 =  PRMUZ(JL) * PRMUZ(JL)
         ZRK = SQRT(3.*ZWM*ZX1)
         ZX2 = 4.*(1.-ZRK*ZRK*ZRM2)
         ZRP=ZRK/ZX1
         ZALPHA = 3.*ZWCP*ZRM2*(1.+ZGP*ZWM)/ZX2
         ZBETA = 3.*ZWCP* PRMUZ(JL) *(1.+3.*ZGP*ZRM2*ZWM)/ZX2
         ZARG=MIN(ZTOP/PRMUZ(JL),200.)
         ZEXMU0=EXP  (-ZARG)
         ZARG2=MIN(ZRK*ZTOP,200.)
         ZEXKP=EXP  (ZARG2)
         ZEXKM = 1./ZEXKP
         ZXP2P = 1.+ZDT*ZRP
         ZXM2P = 1.-ZDT*ZRP
         ZAP2B = ZALPHA+ZDT*ZBETA
         ZAM2B = ZALPHA-ZDT*ZBETA
C
C*         1.2     WITHOUT REFLECTION FROM THE UNDERLYING LAYER
C
         ZA11 = ZXP2P
         ZA12 = ZXM2P
         ZA13 = ZAP2B
         ZA22 = ZXP2P*ZEXKP
         ZA21 = ZXM2P*ZEXKM
         ZA23 = ZAM2B*ZEXMU0
         ZDENA = ZA11 * ZA22 - ZA21 * ZA12
         ZC1A = (ZA22*ZA13-ZA12*ZA23)/ZDENA
         ZC2A = (ZA11*ZA23-ZA21*ZA13)/ZDENA
         ZRI0A = ZC1A+ZC2A-ZALPHA
         ZRI1A = ZRP*(ZC1A-ZC2A)-ZBETA
         PRE1(JL) = (ZRI0A-ZDT*ZRI1A)/ PRMUZ(JL)
         ZRI0B = ZC1A*ZEXKM+ZC2A*ZEXKP-ZALPHA*ZEXMU0
         ZRI1B = ZRP*(ZC1A*ZEXKM-ZC2A*ZEXKP)-ZBETA*ZEXMU0
         PTR1(JL) = ZEXMU0+(ZRI0B+ZDT*ZRI1B)/ PRMUZ(JL)
C
C*         1.3     WITH REFLECTION FROM THE UNDERLYING LAYER
C
         ZB21 = ZA21- PREF(JL) *ZXP2P*ZEXKM
         ZB22 = ZA22- PREF(JL) *ZXM2P*ZEXKP
         ZB23 = ZA23- PREF(JL) *ZEXMU0*(ZAP2B - PRMUZ(JL) )
         ZDENB = ZA11 * ZB22 - ZB21 * ZA12
         ZC1B = (ZB22*ZA13-ZA12*ZB23)/ZDENB
         ZC2B = (ZA11*ZB23-ZB21*ZA13)/ZDENB
         ZRI0C = ZC1B+ZC2B-ZALPHA
         ZRI1C = ZRP*(ZC1B-ZC2B)-ZBETA
         PRE2(JL) = (ZRI0C-ZDT*ZRI1C) / PRMUZ(JL)
         ZRI0D = ZC1B*ZEXKM + ZC2B*ZEXKP - ZALPHA*ZEXMU0
         ZRI1D = ZRP * (ZC1B*ZEXKM - ZC2B*ZEXKP) - ZBETA*ZEXMU0
         PTR2(JL) = ZEXMU0 + (ZRI0D + ZDT*ZRI1D) / PRMUZ(JL)
C
      ENDDO
C
      RETURN
      END SUBROUTINE SWDE
