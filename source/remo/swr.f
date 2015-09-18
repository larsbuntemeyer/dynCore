C
C     SUBROUTINE SWR
C
C**** *SWR* - CONTINUUM SCATTERING COMPUTATIONS
C
C     PURPOSE.
C     --------
C           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
C     CONTINUUM SCATTERING
C
C**   INTERFACE.
C     ----------
C     *CALL*     SWR ( KLON, KLEV, KAER, KNU
C    S  , PAER, PALBS, PCG, PCLDSW, PDSIG, POMEGA, PRAYL, PSEC, PTAU
C    S  , PCGAZ, PPIZAZ, PRAY1, PRAY2, PREFZ, PRJ, PRK, PRMUE
C    S  , PTAUAZ, PTRA1, PTRA2 )
C
C          *SWR* IS CALLED EITHER FROM *SW1S*
C                              OR FROM *SW2S*
C
C        IMPLICIT ARGUMENTS :
C        --------------------
C
C     ==== INPUTS ===
C     ==== OUTPUTS ===
C
C     METHOD.
C     -------
C
C          1. COMPUTES CONTINUUM FLUXES CORRESPONDING TO AEROSOL
C     OR/AND RAYLEIGH SCATTERING (NO MOLECULAR GAS ABSORPTION)
C
C     EXTERNALS.
C     ----------
C
C          *SWDE*
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
C        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C       MODIFIED   U. SCHLESE   DKRZ-HAMBURG  JUL-93
C       MODIFIED   T. DIEHL     DKRZ-HAMBURG  SEP-98
C     ------------------------------------------------------------------
C
      SUBROUTINE SWR ( KLON, KLEV, KNU
     &  , PAER, PALBS, PCG, PCLDSW, PDSIG, POMEGA, PRAYL, PSEC, PTAU
     &  , PCGAZ, PPIZAZ, PRAY1, PRAY2, PREFZ, PRJ, PRK, PRMUE
     &  , PTAUAZ, PTRA1, PTRA2 )
C
      IMPLICIT NONE
C
      INCLUDE "YOMAER"
      INCLUDE "YOMRDU"
      INCLUDE "YOMRDI"
C     ------------------------------------------------------------------
C
C*       0.1   DUMMY ARGUMENTS
C              ---------------
C
      INTEGER, INTENT(IN)    :: KLON, KLEV, KNU
C
      REAL,    INTENT(IN)    ::
     &     PAER(KLON,KLEV,5), PALBS(KLON,2), PCG(KLON,2,KLEV)
     &  ,  PCLDSW(KLON,KLEV), PDSIG(KLON,KLEV)
     &  ,  POMEGA(KLON,2,KLEV), PRAYL(KLON)
     &  ,  PSEC(KLON), PTAU(KLON,2,KLEV)
C
      REAL,    INTENT(INOUT) :: 
     &     PRAY1(KLON,KLEV+1), PRAY2(KLON,KLEV+1)
     &  ,  PREFZ(KLON,2,KLEV+1), PRJ(KLON,6,KLEV+1)
     &  ,  PRK(KLON,6,KLEV+1), PRMUE(KLON,KLEV+1)
     &  ,  PCGAZ(KLON,KLEV),PPIZAZ(KLON,KLEV),PTAUAZ(KLON,KLEV)
     &  ,  PTRA1(KLON,KLEV+1), PTRA2(KLON,KLEV+1)
C
C     ------------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS
C              ------------
C
      REAL ZC1I(KLON,KLEV+1), ZGG(KLON), ZREF(KLON)
     &  ,  ZRE1(KLON), ZRE2(KLON)
     &  ,  ZRMUZ(KLON), ZRNEB(KLON)
     &  ,  ZTO1(KLON), ZTR(KLON,2,KLEV+1)
     &  ,  ZTR1(KLON), ZTR2(KLON), ZW(KLON)
C
      INTEGER :: JK,JA,JL,JAE,JAJ,JKL,JKLP1,JKM1
      REAL    :: ZBMU0,ZBMU1,ZWW,ZTRAY,ZTO,ZSS1,ZRE11,ZRATIO,ZR22,
     &           ZR21,ZMUE,ZMU1,ZGAR,ZGAP,ZFF,ZFACOC,ZFACOA,ZDEN1,
     &           ZDEN,ZCORCD,ZCORAE,ZCLTEST
C
C
C     ------------------------------------------------------------------
C
C*         1.    OPTICAL PARAMETERS FOR AEROSOLS AND RAYLEIGH
C                --------------------------------------------
C
      DO JK = 1 , KLEV+1
         DO JA = 1 , 6
            DO JL = 1 , KLON
               PRJ(JL,JA,JK) = 0.
               PRK(JL,JA,JK) = 0.
            ENDDO
         ENDDO
      ENDDO
C
      DO JK = 1 , KLEV
         DO JL = 1 , KLON
            PCGAZ(JL,JK) = 0.
            PPIZAZ(JL,JK) =  0.
            PTAUAZ(JL,JK) = 0.
         ENDDO
         DO JAE=1,5
            DO JL = 1 , KLON
               PTAUAZ(JL,JK)=PTAUAZ(JL,JK)
     &              +PAER(JL,JK,JAE)*TAUA(KNU,JAE)
               PPIZAZ(JL,JK)=PPIZAZ(JL,JK)+PAER(JL,JK,JAE)
     &              * TAUA(KNU,JAE)*PIZA(KNU,JAE)
               PCGAZ(JL,JK) =  PCGAZ(JL,JK) +PAER(JL,JK,JAE)
     &              * TAUA(KNU,JAE)*PIZA(KNU,JAE)*CGA(KNU,JAE)
            ENDDO
         ENDDO
C
         DO JL = 1 , KLON
            PCGAZ (JL,JK)=PCGAZ (JL,JK)/PPIZAZ(JL,JK)
            PPIZAZ(JL,JK)=PPIZAZ(JL,JK)/PTAUAZ(JL,JK)
            ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
            ZRATIO = ZTRAY / (ZTRAY + PTAUAZ(JL,JK))
            ZGAR = PCGAZ(JL,JK)
            ZFF = ZGAR * ZGAR
            PTAUAZ(JL,JK)=ZTRAY+PTAUAZ(JL,JK)*(1.-PPIZAZ(JL,JK)*ZFF)
            PCGAZ(JL,JK) = ZGAR * (1. - ZRATIO) / (1. + ZGAR)
            PPIZAZ(JL,JK) =ZRATIO+(1.-ZRATIO)*PPIZAZ(JL,JK)*(1.-ZFF)
     &           / (1. - PPIZAZ(JL,JK) * ZFF)
         ENDDO
C
      ENDDO
C
C
C     ------------------------------------------------------------------
C
C*         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
C                ----------------------------------------------
C
      DO JL = 1 , KLON
         ZC1I(JL,KLEV+1) = 0.
      ENDDO
C
      DO JK = 1 , KLEV
         JKL = KLEV+1 - JK
         JKLP1 = JKL + 1
         DO JL = 1 , KLON
            ZFACOA = 1.-PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
            ZFACOC = 1. - POMEGA(JL,KNU,JKL) * PCG(JL,KNU,JKL)
     &                                       * PCG(JL,KNU,JKL)
            ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
            ZCORCD = ZFACOC * PTAU(JL,KNU,JKL) * PSEC(JL)
            ZR21 = EXP  (-ZCORAE   )
            ZR22 = EXP  (-ZCORCD   )
            ZSS1 = PCLDSW(JL,JKL)*(1.0-ZR21*ZR22)
     &           + (1.0-PCLDSW(JL,JKL))*(1.0-ZR21)
CRAN  ZC1I(JL,JKL) = 1.0-(1.0-ZSS1)*(1.0-ZC1I(JL,JKLP1))
CMAX  ZC1I(JL,JKL) = AMAX1( ZSS1 , ZC1I(JL,JKLP1) )
            ZC1I(JL,JKL) = 1.0-(1.0-ZSS1)
     &           *(1.0-AMAX1(ZSS1,ZC1I(JL,JKLP1)))
     &           /(1.0-AMIN1(ZC1I(JL,JKLP1),1.-ZEPSEC))
         ENDDO
C
      ENDDO
C
C
C     ------------------------------------------------------------------
C
C*         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
C                -----------------------------------------------
C
      DO JL = 1 , KLON
         PRAY1(JL,KLEV+1) = 0.
         PRAY2(JL,KLEV+1) = 0.
         PREFZ(JL,2,1) = PALBS(JL,KNU)
         PREFZ(JL,1,1) = PALBS(JL,KNU)
         PTRA1(JL,KLEV+1) = 1.
         PTRA2(JL,KLEV+1) = 1.
      ENDDO
C
      DO JK = 2 , KLEV+1
         JKM1 = JK-1
         DO JL = 1 , KLON
            ZRNEB(JL)= PCLDSW(JL,JKM1)
            ZRE1(JL)=0.
            ZTR1(JL)=0.
            ZRE2(JL)=0.
            ZTR2(JL)=0.
C
C
C     ------------------------------------------------------------------
C
C*         3.1  EQUIVALENT ZENITH ANGLE
C               -----------------------
C
            ZMUE = (1.-ZC1I(JL,JK)) * PSEC(JL)
     &           + ZC1I(JL,JK) * 1.66
            PRMUE(JL,JK) = 1./ZMUE
C
C
C     ------------------------------------------------------------------
C
C*         3.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS
C               ----------------------------------------------------
C
            ZGAP = PCGAZ(JL,JKM1)
            ZBMU0 = 0.5 - 0.75 * ZGAP / ZMUE
            ZWW = PPIZAZ(JL,JKM1)
            ZTO = PTAUAZ(JL,JKM1)
            ZDEN = 1. + (1. - ZWW + ZBMU0 * ZWW) * ZTO * ZMUE
     &           + (1-ZWW) * (1. - ZWW +2.*ZBMU0*ZWW)*ZTO*ZTO*ZMUE*ZMUE
            PRAY1(JL,JKM1) = ZBMU0 * ZWW * ZTO * ZMUE / ZDEN
            PTRA1(JL,JKM1) = 1. / ZDEN
C
            ZMU1 = 0.5
            ZBMU1 = 0.5 - 0.75 * ZGAP * ZMU1
            ZDEN1= 1. + (1. - ZWW + ZBMU1 * ZWW) * ZTO / ZMU1
     &           + (1-ZWW) * (1. - ZWW +2.*ZBMU1*ZWW)*ZTO*ZTO/ZMU1/ZMU1
            PRAY2(JL,JKM1) = ZBMU1 * ZWW * ZTO / ZMU1 / ZDEN1
            PTRA2(JL,JKM1) = 1. / ZDEN1
C
C
C     ------------------------------------------------------------------
C
C*         3.3  EFFECT OF CLOUD LAYER
C               ---------------------
C
            ZW(JL) = POMEGA(JL,KNU,JKM1)
            ZTO1(JL) = PTAU(JL,KNU,JKM1)/ZW(JL)
     &           + PTAUAZ(JL,JKM1)/PPIZAZ(JL,JKM1)
            ZR21 = PTAU(JL,KNU,JKM1) + PTAUAZ(JL,JKM1)
            ZR22 = PTAU(JL,KNU,JKM1) / ZR21
            ZGG(JL) = ZR22 * PCG(JL,KNU,JKM1)
     &           + (1. - ZR22) * PCGAZ(JL,JKM1)
            ZW(JL) = ZR21 / ZTO1(JL)
            ZREF(JL) = PREFZ(JL,1,JKM1)
            ZRMUZ(JL) = PRMUE(JL,JK)
         ENDDO
C
C      ZCLTEST=SSUM(KLON,ZRNEB,1)
         ZCLTEST=SUM(ZRNEB)
         IF(ZCLTEST.GT.ZEPSC*1.E3) THEN
C
            CALL SWDE ( KLON
     &           , ZGG,ZREF,ZRMUZ,ZTO1,ZW
     &           , ZRE1,ZRE2,ZTR1,ZTR2     )
C
         ENDIF
C
         DO JL = 1 , KLON
C
            PREFZ(JL,1,JK) = (1.-ZRNEB(JL)) * (PRAY1(JL,JKM1)
     &           + PREFZ(JL,1,JKM1) * PTRA1(JL,JKM1)
     &           * PTRA2(JL,JKM1)
     &           / (1.-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))
     &           + ZRNEB(JL) * ZRE2(JL)
C
            ZTR(JL,1,JKM1) = ZRNEB(JL) * ZTR2(JL) + (PTRA1(JL,JKM1)
     &           / (1.-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))
     &           * (1.-ZRNEB(JL))
C
            PREFZ(JL,2,JK) = (1.-ZRNEB(JL)) * (PRAY1(JL,JKM1)
     &           + PREFZ(JL,2,JKM1) * PTRA1(JL,JKM1)
     &           * PTRA2(JL,JKM1) )
     &           + ZRNEB(JL) * ZRE1(JL)
C
            ZTR(JL,2,JKM1) = ZRNEB(JL) * ZTR1(JL)
     &           + PTRA1(JL,JKM1) * (1.-ZRNEB(JL))
C
         ENDDO
      ENDDO
      DO JL = 1 , KLON
         ZMUE = (1.-ZC1I(JL,1))*PSEC(JL)+ZC1I(JL,1)*1.66
         PRMUE(JL,1)=1./ZMUE
      ENDDO


C
C
C     ------------------------------------------------------------------
C
C*         3.5    REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
C                 -------------------------------------------------
C
      IF (KNU.EQ.1) THEN
         JAJ = 2
         DO JL = 1 , KLON
            PRJ(JL,JAJ,KLEV+1) = 1.
            PRK(JL,JAJ,KLEV+1) = PREFZ(JL, 1,KLEV+1)
         ENDDO
C
         DO JK = 1 , KLEV
            JKL = KLEV+1 - JK
            JKLP1 = JKL + 1
            DO JL = 1 , KLON
               ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,  1,JKL)
               PRJ(JL,JAJ,JKL) = ZRE11
               PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,  1,JKL)
            ENDDO
         ENDDO
C
      ELSE
C
         DO JAJ = 1 , 2
            DO JL = 1 , KLON
               PRJ(JL,JAJ,KLEV+1) = 1.
               PRK(JL,JAJ,KLEV+1) = PREFZ(JL,JAJ,KLEV+1)
            ENDDO
C
            DO JK = 1 , KLEV
               JKL = KLEV+1 - JK
               JKLP1 = JKL + 1
               DO JL = 1 , KLON
                  ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,JAJ,JKL)
                  PRJ(JL,JAJ,JKL) = ZRE11
                  PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,JAJ,JKL)
               ENDDO
            ENDDO
         ENDDO
C
      ENDIF
C
C     ------------------------------------------------------------------
C
      RETURN
      END
