C
C**** *LWTT* - LONGWAVE TRANSMISSION FUNCTIONS
C
C     PURPOSE.
C     --------
C           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
C     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN ALL SIX SPECTRAL
C     INTERVALS.
C
C**   INTERFACE.
C     ----------
C          *LWTT* IS CALLED FROM *LWVN*, *LWVD*, *LWVB*
C
C     *CALL*     LWTT (KLON,KAER,KCFC,PGA,PGB,PGC,PGD,PUU,PTT)
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C     ==== INPUTS ===
C KND    :                     ; WEIGHTING INDEX
C PUU    : (KLON,NUA)         ; ABSORBER AMOUNTS
C     ==== OUTPUTS ===
C PTT    : (KLON,NTRA)        ; TRANSMISSION FUNCTIONS
C
C        IMPLICIT ARGUMENTS :   NONE
C        --------------------
C
C     METHOD.
C     -------
C
C          1. TRANSMISSION FUNCTION BY H2O AND UNIFORMLY MIXED GASES ARE
C     COMPUTED USING PADE APPROXIMANTS AND HORNER'S ALGORITHM.
C          2. TRANSMISSION BY O3 IS EVALUATED WITH MALKMUS'S BAND MODEL.
C        3. TRANSMISSION BY H2O CONTINUUM, CFC'S AND AEROSOLS FOLLOW AN
C           A SIMPLE EXPONENTIAL DECREASE WITH ABSORBER AMOUNT.
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
C        MODIFIED: 92-06-12  ROB VAN DORLAND *KNMI*
C        REWRITTEN: U.SCHLESE  MAY-93
C-----------------------------------------------------------------------
      SUBROUTINE LWTT (KLON,KAER,KCFC,PGA,PGB,PGC,PGD,PUU,PTT)
C
      IMPLICIT NONE
C
      INCLUDE "YOMLW"
C     ------------------------------------------------------------------
C
C*        0.1   ARGUMENTS
C               ---------
C
      INTEGER, INTENT(IN)    :: KLON,KAER,KCFC
      REAL,    INTENT(IN)    :: PUU(KLON,NUA)
     &                       ,  PGA(KLON,8,2), PGB(KLON,8,2)
     &                       ,  PGC(KLON,5,2),PGD(KLON,5,2)
      REAL,    INTENT(INOUT) :: PTT(KLON,NTRA)
C
C     ------------------------------------------------------------------
C
C*        0.2   LOCAL ARRAYS
C               ------------
C
      REAL    :: ZTTNC(KLON,5)
C
      INTEGER :: ITMP1, ITMP2, JL
      REAL    :: 
     &     ZCOAC, ZEU, ZPEU10, ZPEU11, ZPEU12, ZPEU13, ZPU, ZSOZ, ZSQ1, 
     &     ZSQ2, ZTO1, ZTO2, ZTOZ, ZUXY, ZUXZ, ZVXY, ZVXZ, ZX, ZXD, ZXI2
      REAL    :: ZXN, ZY, ZYI2, ZZ
C
C     ------------------------------------------------------------------
C
C*        0.3   FUNCTIONS
C               ---------
C
C
!DIR$ NOBOUNDS

C     ------------------------------------------------------------------
C
C*         1.1    HORNER'S ALGORITHM FOR H2O AND CO2 TRANSMISSION
C                 -----------------------------------------------
C
      DO itmp2 = 1, 8
         DO itmp1 = 1, KLON
            ZZ=SQRT  (PUU(itmp1,itmp2))
            ZXD=PGB(itmp1,itmp2,1)+ZZ*(PGB(itmp1,itmp2,2)+ZZ)
            ZXN=PGA(itmp1,itmp2,1)+ZZ*(PGA(itmp1,itmp2,2))
            PTT(itmp1,itmp2)=ZXN/ZXD
         ENDDO
      ENDDO

C
C
C*         1.2    HORNER'S ALGORITHM FOR N2O AND CH4 TRANSMISSION
C                 -----------------------------------------------
C
      DO itmp2 = 1, 5
         DO itmp1 = 1, KLON
            ZZ=SQRT  (PUU(itmp1,14+itmp2-1) )
            ZXN=PGC(itmp1,itmp2,1)+ZZ*(PGC(itmp1,itmp2,2))
            ZXD=PGD(itmp1,itmp2,1)+ZZ*(PGD(itmp1,itmp2,2)+ZZ)
            ZTTNC(itmp1,itmp2)=ZXN/ZXD
         ENDDO
      ENDDO
C     ------------------------------------------------------------------
C
C*      2.     CONTINUUM, OZONE, AEROSOLS AND TRACEGASES
C              -----------------------------------------
C
C    CODE FOR KAER=1 AND KCFC=1
C
      IF(KAER.EQ.1.AND.KCFC.EQ.1) THEN
         DO JL = 1 , KLON
            PTT(JL, 9) = PTT(JL, 8)
C
C-  CONTINUUM ABSORPTION: E- AND P-TYPE
C
            ZPU=PUU(JL,10)
            ZEU=PUU(JL,11)
C
C-  OZONE ABSORPTION
C
            ZX = PUU(JL,12)
            ZY = PUU(JL,13)
            ZUXY = 4. * ZX * ZX / (PIALF0 * ZY)
            ZVXY=(PIALF0*ZY)/(ZX+ZX)
            ZSQ1 = SQRT  (1. + O1H * ZUXY ) - 1.
            ZSQ2 = SQRT  (1. + O2H * ZUXY ) - 1.
C
            ZXI2 = PUU(JL,28)
            ZYI2 = PUU(JL,29)
            ZUXZ = 4. * ZXI2 * ZXI2 / (PIAOD2 * ZYI2)
            ZSOZ = SQRT  (1. + SAVOD2 * ZUXZ ) - 1.
            ZVXZ = PIAOD2 * ZYI2 / (2. * ZXI2)
C
C*      INTERVAL 0-350 CM-1 + 1440-1880 CM-1
C
            PTT(JL,10)=EXP  (-PUU(JL,23))
C
C*      INTERVAL 500-800 CM-1
C
            ZCOAC=47.7*(0.017*ZPU+ZEU)+PUU(JL,19)+PUU(JL,24)
            ZTOZ=EXP  (-ZVXZ*ZSOZ-ZCOAC)
            PTT(JL,11)=ZTTNC(JL,1)*ZTOZ
C
C*      INTERVAL 800-970 CM-1 + 1110-1250 CM-1
C
            ZCOAC=8.31*(0.0025*ZPU+ZEU)+PUU(JL,20)+PUU(JL,25)
            PTT(JL,12)=ZTTNC(JL,2)*ZTTNC(JL,4)*EXP  (-ZCOAC)
C
C*      INTERVAL 970-1110 CM-1
C
            ZCOAC=5.87*(0.0018*ZPU+ZEU)+PUU(JL,21)+PUU(JL,26)
            ZTO1=EXP  (-ZVXY*ZSQ1-ZCOAC)
            ZTO2=EXP  (-ZVXY*ZSQ2-ZCOAC)
            PTT(JL,13)=0.7554*ZTO1+0.2446*ZTO2
C
C*      INTERVAL 350-500 CM-1
C
            ZCOAC=209.*(0.059*ZPU+ZEU)+PUU(JL,27)
            PTT(JL,14)=EXP  (-ZCOAC)
C
C*      INTERVAL 1250-1440 CM-1 + 1880-2820 CM-1
C
            ZCOAC=PUU(JL,22)+PUU(JL,23)
            PTT(JL,15)=ZTTNC(JL,3)*ZTTNC(JL,5)*EXP  (-ZCOAC)
C
         ENDDO
C
C----------------
      ELSEIF (KAER.EQ.0.AND.KCFC.EQ.0) THEN
C ---------------

         DO JL = 1 , KLON
            PTT(JL, 9) = PTT(JL, 8)
C
C-  CONTINUUM ABSORPTION: E- AND P-TYPE
C
            ZPU=PUU(JL,10)
            ZEU=PUU(JL,11)
            ZPEU10=47.7*(0.017*ZPU+ZEU)
            ZPEU11=8.31*(0.0025*ZPU+ZEU)
            ZPEU12=5.87*(0.0018*ZPU+ZEU)
            ZPEU13=209.*(0.059*ZPU+ZEU)
C
C-  OZONE ABSORPTION
C
            ZX = PUU(JL,12)
            ZY = PUU(JL,13)
            ZUXY = 4. * ZX * ZX / (PIALF0 * ZY)
            ZVXY=(PIALF0*ZY)/(ZX+ZX)
            ZSQ1 = SQRT  (1. + O1H * ZUXY ) - 1.
            ZSQ2 = SQRT  (1. + O2H * ZUXY ) - 1.
C
            ZXI2 = PUU(JL,28)
            ZYI2 = PUU(JL,29)
            ZUXZ = 4. * ZXI2 * ZXI2 / (PIAOD2 * ZYI2)
            ZSOZ = SQRT  (1. + SAVOD2 * ZUXZ ) - 1.
            ZVXZ = PIAOD2 * ZYI2 / (2. * ZXI2)
C
C*      INTERVAL 0-350 CM-1 + 1440-1880 CM-1
C
            PTT(JL,10)=1.
C
C*      INTERVAL 500-800 CM-1
C
            ZCOAC=ZPEU10
            ZTOZ=EXP  (-ZVXZ*ZSOZ-ZCOAC)
            PTT(JL,11)=ZTTNC(JL,1)*ZTOZ
C
C*      INTERVAL 800-970 CM-1 + 1110-1250 CM-1
C
            ZCOAC=ZPEU11
            PTT(JL,12)=ZTTNC(JL,2)*ZTTNC(JL,4)*EXP  (-ZCOAC)
C
C*      INTERVAL 970-1110 CM-1
C
            ZCOAC=ZPEU12
            ZTO1=EXP  (-ZVXY*ZSQ1-ZCOAC)
            ZTO2=EXP  (-ZVXY*ZSQ2-ZCOAC)
            PTT(JL,13)=0.7554*ZTO1+0.2446*ZTO2
C
C*      INTERVAL 350-500 CM-1
C
            ZCOAC=ZPEU13
            PTT(JL,14)=EXP  (-ZCOAC)
C
C*      INTERVAL 1250-1440 CM-1 + 1880-2820 CM-1
C
            PTT(JL,15)=ZTTNC(JL,3)*ZTTNC(JL,5)
C
         ENDDO
C------------
      ELSE
C -----------

         DO JL = 1 , KLON
            PTT(JL, 9) = PTT(JL, 8)
C
C-  CONTINUUM ABSORPTION: E- AND P-TYPE
C

            ZPU=PUU(JL,10)
            ZEU=PUU(JL,11)
            ZPEU10=47.7*(0.017*ZPU+ZEU)
            ZPEU11=8.31*(0.0025*ZPU+ZEU)
            ZPEU12=5.87*(0.0018*ZPU+ZEU)
            ZPEU13=209.*(0.059*ZPU+ZEU)
C
C-  OZONE ABSORPTION
C
            ZX = PUU(JL,12)
            ZY = PUU(JL,13)
            ZUXY = 4. * ZX * ZX / (PIALF0 * ZY)
            ZVXY=(PIALF0*ZY)/(ZX+ZX)
            ZSQ1 = SQRT  (1. + O1H * ZUXY ) - 1.
            ZSQ2 = SQRT  (1. + O2H * ZUXY ) - 1.
C
            ZXI2 = PUU(JL,28)
            ZYI2 = PUU(JL,29)
            ZUXZ = 4. * ZXI2 * ZXI2 / (PIAOD2 * ZYI2)
            ZSOZ = SQRT  (1. + SAVOD2 * ZUXZ ) - 1.
            ZVXZ = PIAOD2 * ZYI2 / (2. * ZXI2)
C
C*      INTERVAL 0-350 CM-1 + 1440-1880 CM-1
C
            ZCOAC=KAER*PUU(JL,23)
            PTT(JL,10)=EXP  (-ZCOAC)
C
C*      INTERVAL 500-800 CM-1
C
            ZCOAC=ZPEU10+KCFC*PUU(JL,19)+KAER*PUU(JL,24)
            ZTOZ=EXP  (-ZVXZ*ZSOZ-ZCOAC)
            PTT(JL,11)=ZTTNC(JL,1)*ZTOZ
C
C*      INTERVAL 800-970 CM-1 + 1110-1250 CM-1
C
            ZCOAC=ZPEU11+KCFC*PUU(JL,20)+KAER*PUU(JL,25)
            PTT(JL,12)=ZTTNC(JL,2)*ZTTNC(JL,4)*EXP  (-ZCOAC)
C
C*      INTERVAL 970-1110 CM-1
C
            ZCOAC=ZPEU12+KCFC*PUU(JL,21)+KAER*PUU(JL,26)
            ZTO1=EXP  (-ZVXY*ZSQ1-ZCOAC)
            ZTO2=EXP  (-ZVXY*ZSQ2-ZCOAC)
            PTT(JL,13)=0.7554*ZTO1+0.2446*ZTO2
C
C*      INTERVAL 350-500 CM-1
C
            ZCOAC=ZPEU13+KAER*PUU(JL,27)
            PTT(JL,14)=EXP  (-ZCOAC)
C
C*      INTERVAL 1250-1440 CM-1 + 1880-2820 CM-1
C
            ZCOAC=KCFC*PUU(JL,22)+KAER*PUU(JL,23)
            PTT(JL,15)=ZTTNC(JL,3)*ZTTNC(JL,5)*EXP  (-ZCOAC)
C
         ENDDO
      ENDIF
C
      RETURN
      END SUBROUTINE LWTT
