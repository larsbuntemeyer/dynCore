      SUBROUTINE AOSW2S ( KLON, KLEV, KEWAER, KAERH, KNU
     &  ,  PPMB
     &  ,  PAER,PAKI,PALBS,PCG,PCLDSW,PDSIG,POMEGA,PRMU,PSEC,PTAU
     &  ,  PUD,PUM,PWV,PQS
     &  ,  PFDOWN,PFUP                                                )
C
      IMPLICIT NONE
C
C**** *SW2* - SHORTWAVE RADIATION, 2ND SPECTRAL INTERVAL
C
C     PURPOSE.
C     --------
C
C          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN THE
C     SECOND SPECTRAL INTERVAL FOLLOWING FOUQUART AND BONNEL (1980).
C
C**   INTERFACE.
C     ----------
C     *CALL*     SW2S ( KLON, KLEV, KAER, KNU
C    1  ,  PAER,PAKI,PALBS,PCG,PCLDSW,PDSIG,POMEGA,PRMU,PSEC,PTAU
C    2  ,  PUD,PUM
C    3  ,  PFDOWN,PFUP                                                )
C
C          *SW2S* IS CALLED FROM *SW*.
C
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
C          1. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING TO
C     CONTINUUM SCATTERING
C          2. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING FOR
C     A GREY MOLECULAR ABSORPTION
C          3. LAPLACE TRANSFORM ON THE PREVIOUS TO GET EFFECTIVE AMOUNTS
C     OF ABSORBERS
C          4. APPLY H2O AND U.M.G. TRANSMISSION FUNCTIONS
C          5. MULTIPLY BY OZONE TRANSMISSION FUNCTION
C
C     EXTERNALS.
C     ----------
C
C          *SWR*, *SWDE*, *SWTT*
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
C         MODIFIED  U. SCHLESE   DKRZ HAMBURG  JUL-93
C       MODIFIED: ROB VAN DORLAND, KNMI, 95-05-10
C     ------------------------------------------------------------------
C
C
      INCLUDE "YOMSW"
      INCLUDE "YOMRDU"
C     ------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
      INTEGER, INTENT(IN) :: KEWAER, KLEV, KLON, KNU
C     
      INTEGER, INTENT(IN) :: KAERH(KLON,KLEV)
C    
      REAL,    INTENT(IN) :: 
     &     PAER(KLON,KLEV,5+KEWAER)
     &  ,  PAKI(KLON,2), PALBS(KLON,2)
     &  ,  PPMB (KLON,KLEV+1)
     &  ,  PCG(KLON,2,KLEV), PCLDSW(KLON,KLEV), PDSIG(KLON,KLEV)
     &  ,  POMEGA(KLON,2,KLEV), PRMU(KLON), PSEC(KLON)
     &  ,  PTAU(KLON,2,KLEV), PUD(KLON,3,KLEV+1)
     &  ,  PUM(KLON,KLEV+1), PWV(KLON,KLEV), PQS(KLON,KLEV)
C
      REAL, INTENT(INOUT) :: PFDOWN(KLON,KLEV+1),PFUP(KLON,KLEV+1)    
C
C     ------------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS
C              ------------
C     Local Variables
C
      INTEGER :: IKL, JABS, JAJ, JAJP, JK, JKKI, JKKP4, JKL, JKLP1,
     &           JL, JN, JN2J, KREF, JKM1
      REAL    :: ZAA, ZBB, ZCLTEST, ZCNEB, ZR21, ZR22, ZRE11, ZRKI, 
     &           ZWH2O, ZRMUM1
C
      REAL ZCGAZ(KLON,KLEV), ZG(KLON), ZGG(KLON)
     &  ,  ZPIZAZ(KLON,KLEV), ZRAYL(KLON), ZRAY1(KLON,KLEV+1)
     &  ,  ZRAY2(KLON,KLEV+1), ZREF(KLON), ZREFZ(KLON,2,KLEV+1)
     &  ,  ZRE1(KLON), ZRE2(KLON), ZRJ(KLON,6,KLEV+1)
     &  ,  ZRK(KLON,6,KLEV+1), ZRL(KLON,8),  ZRMUE(KLON,KLEV+1)
     &  ,  ZRMUZ(KLON), ZRNEB(KLON),  ZR1(KLON)
     &  ,  ZR2(KLON),   ZS(KLON)
     &  ,  ZTAUAZ(KLON,KLEV), ZTO1(KLON), ZTR(KLON,2,KLEV+1)
     &  ,  ZTRA1(KLON,KLEV+1), ZTRA2(KLON,KLEV+1)
     &  ,  ZTR1(KLON), ZTR2(KLON), ZW(KLON), ZW1(KLON), ZW2(KLON)
C
C
!DIR$ NOBOUNDS
C     ------------------------------------------------------------------
C
C*         1.     SECOND SPECTRAL INTERVAL (0.68-4.00 MICRON)
C                 -------------------------------------------
C
C
C*         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
C                 -----------------------------------------
C
      DO JL = 1 , KLON
         ZRMUM1 = 1. - PRMU(JL)
         ZRAYL(JL) = CRAY(KNU,1) + ZRMUM1   * (CRAY(KNU,2) + ZRMUM1
     &            * (CRAY(KNU,3) + ZRMUM1   * (CRAY(KNU,4) + ZRMUM1
     &            * (CRAY(KNU,5) + ZRMUM1   *  CRAY(KNU,6)       ))))
      ENDDO
C
C
C     ------------------------------------------------------------------
C
C*         2.    CONTINUUM SCATTERING CALCULATIONS
C                ---------------------------------
C
      CALL AOSWR ( KLON, KLEV, KEWAER, KAERH, KNU
     &  ,  PPMB
     &  , PAER, PALBS, PCG, PCLDSW, PDSIG, POMEGA, ZRAYL, PSEC, PTAU
     &  , ZCGAZ, ZPIZAZ, ZRAY1, ZRAY2, ZREFZ, ZRJ, ZRK, ZRMUE
     &  , ZTAUAZ, ZTRA1, ZTRA2 )
C
C
C     ------------------------------------------------------------------
C
C*         3.    SCATTERING CALCULATIONS WITH GREY MOLECULAR ABSORPTION
C                ------------------------------------------------------
C
      JN = 2
C
      DO JABS=1,2
C
C
C*         3.1  SURFACE CONDITIONS
C               ------------------
C
         DO JL = 1 , KLON
            ZREFZ(JL,2,1) = PALBS(JL,KNU)
            ZREFZ(JL,1,1) = PALBS(JL,KNU)
         ENDDO
C
C
C*         3.2  INTRODUCING CLOUD EFFECTS
C               -------------------------
C
         DO JK = 2 , KLEV+1
            JKM1 = JK - 1
            IKL=KLEV+1-JKM1
            DO JL = 1 , KLON
               ZRNEB(JL) = PCLDSW(JL,JKM1)
               IF (JABS.EQ.1 .AND. ZRNEB(JL).GT.2.*ZEELOG) THEN
                  ZWH2O=AMAX1(PWV(JL,IKL),ZEELOG)
                  ZCNEB=AMAX1(ZEELOG,AMIN1(ZRNEB(JL),1.-ZEELOG))
                  ZBB=PUD(JL,JABS,JKM1)*PQS(JL,IKL)/ZWH2O
                  ZAA=AMAX1((PUD(JL,JABS,JKM1)-ZCNEB*ZBB)/(1.-ZCNEB)
     &                      ,ZEELOG)
               ELSE
                  ZAA=PUD(JL,JABS,JKM1)
                  ZBB=ZAA
               END IF
               ZRKI = PAKI(JL,JABS)
               ZS(JL) = EXP(-ZRKI * ZAA * 1.66)
               ZG(JL) = EXP(-ZRKI * ZAA / ZRMUE(JL,JK))
               ZTR1(JL) = 0.
               ZRE1(JL) = 0.
               ZTR2(JL) = 0.
               ZRE2(JL) = 0.
C
               ZW(JL)= POMEGA(JL,KNU,JKM1)
               ZTO1(JL) = PTAU(JL,KNU,JKM1) / ZW(JL)
     &              + ZTAUAZ(JL,JKM1) / ZPIZAZ(JL,JKM1)
     &              + ZBB * ZRKI
               ZR21 = PTAU(JL,KNU,JKM1) + ZTAUAZ(JL,JKM1)
               ZR22  = PTAU(JL,KNU,JKM1) / ZR21
               ZGG(JL) = ZR22 * PCG(JL,KNU,JKM1)
     &              + (1. - ZR22) * ZCGAZ(JL,JKM1)
               ZW(JL) = ZR21 / ZTO1(JL)
               ZREF(JL) = ZREFZ(JL,1,JKM1)
               ZRMUZ(JL) = ZRMUE(JL,JK)
            ENDDO
C
            ZCLTEST=SUM(ZRNEB)
            IF(ZCLTEST.GT.ZEPSC*1.E3) THEN
C
               CALL SWDE ( KLON
     &              , ZGG,ZREF,ZRMUZ,ZTO1,ZW
     &              , ZRE1,ZRE2,ZTR1,ZTR2     )
C
            ENDIF
C
            DO JL = 1 , KLON
C
               ZREFZ(JL,2,JK) = (1.-ZRNEB(JL)) * (ZRAY1(JL,JKM1)
     &              + ZREFZ(JL,2,JKM1) * ZTRA1(JL,JKM1)
     &              * ZTRA2(JL,JKM1) ) * ZG(JL) * ZS(JL)
     &              + ZRNEB(JL) * ZRE1(JL)
C
               ZTR(JL,2,JKM1)=ZRNEB(JL)*ZTR1(JL)
     &              + ZTRA1(JL,JKM1) * ZG(JL) * (1.-ZRNEB(JL))
C
               ZREFZ(JL,1,JK)=(1.-ZRNEB(JL))*(ZRAY1(JL,JKM1)
     &              +ZREFZ(JL,1,JKM1)*ZTRA1(JL,JKM1)*ZTRA2(JL,JKM1)
     &              /(1.-ZRAY2(JL,JKM1)*ZREFZ(JL,1,JKM1)))*ZG(JL)*ZS(JL)
     &              + ZRNEB(JL) * ZRE2(JL)
C
               ZTR(JL,1,JKM1)= ZRNEB(JL) * ZTR2(JL)
     &              + (ZTRA1(JL,JKM1)/(1.-ZRAY2(JL,JKM1)
     &              * ZREFZ(JL,1,JKM1)))
     &              * ZG(JL) * (1. -ZRNEB(JL))
C
            ENDDO
         ENDDO
C
C*         3.3  REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
C               -------------------------------------------------
C
         DO KREF=1,2
C
            JN = JN + 1
C
            DO JL = 1 , KLON
               ZRJ(JL,JN,KLEV+1) = 1.
               ZRK(JL,JN,KLEV+1) = ZREFZ(JL,KREF,KLEV+1)
            ENDDO
C
            DO JK = 1 , KLEV
               JKL = KLEV+1 - JK
               JKLP1 = JKL + 1
               DO JL = 1 , KLON
                  ZRE11 = ZRJ(JL,JN,JKLP1) * ZTR(JL,KREF,JKL)
                  ZRJ(JL,JN,JKL) = ZRE11
                  ZRK(JL,JN,JKL) = ZRE11 * ZREFZ(JL,KREF,JKL)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C
C     ------------------------------------------------------------------
C
C*         4.    INVERT GREY AND CONTINUUM FLUXES
C                --------------------------------
C
C
C*         4.1   UPWARD (ZRK) AND DOWNWARD (ZRJ) PSEUDO-FLUXES
C                ---------------------------------------------
C
      DO JK = 1 , KLEV+1
         DO JAJ = 1 , 5 , 2
            JAJP = JAJ + 1
            DO JL = 1 , KLON
               ZRJ(JL,JAJ,JK)=  ZRJ(JL,JAJ,JK) - ZRJ(JL,JAJP,JK)
               ZRK(JL,JAJ,JK)=  ZRK(JL,JAJ,JK) - ZRK(JL,JAJP,JK)
            ENDDO
         ENDDO
      ENDDO
C
      DO JL = 1 , KLON*6*(KLEV+1)
         ZRJ(JL,1,1)= AMAX1( ZRJ(JL,1,1) , ZEELOG )
         ZRK(JL,1,1)= AMAX1( ZRK(JL,1,1) , ZEELOG )
      ENDDO
C
C*         4.2    EFFECTIVE ABSORBER AMOUNTS BY INVERSE LAPLACE
C                 ---------------------------------------------
C
      DO JK = 1 , KLEV+1
         JKKI = 1
         DO JAJ = 1 , 2
            DO JN = 1 , 2
               JN2J = JN + 2 * JAJ
               JKKP4 = JKKI + 4
C
C*         4.2.1  EFFECTIVE ABSORBER AMOUNTS
C                 --------------------------
C
               DO JL = 1 , KLON
                  ZW1(JL) = ALOG( ZRJ(JL,JN,JK) / ZRJ(JL,JN2J,JK))
     &                 / PAKI(JL,JAJ)
               ENDDO
C
C*         4.2.2  TRANSMISSION FUNCTION
C                 ---------------------
C
               CALL SWTT ( KLON, KNU, JAJ, ZW1, ZR1 )
C
               DO JL = 1 , KLON
                  ZRL(JL,JKKI) = ZR1(JL)
                  ZW2(JL) = ALOG( ZRK(JL,JN,JK) / ZRK(JL,JN2J,JK))
     &                 / PAKI(JL,JAJ)
               ENDDO
C
               CALL SWTT ( KLON, KNU, JAJ, ZW2, ZR2 )
C
               DO JL = 1 , KLON
                  ZRL(JL,JKKP4) = ZR2(JL)
               ENDDO
C
               JKKI=JKKI+1
            ENDDO
         ENDDO
C
C*         4.3    UPWARD AND DOWNWARD FLUXES WITH H2O AND UMG ABSORPTION
C                 ------------------------------------------------------
C
         DO JL = 1 , KLON
            PFDOWN(JL,JK) = ZRJ(JL,1,JK) * ZRL(JL,1) * ZRL(JL,3)
     &                    + ZRJ(JL,2,JK) * ZRL(JL,2) * ZRL(JL,4)
            PFUP(JL,JK)   = ZRK(JL,1,JK) * ZRL(JL,5) * ZRL(JL,7)
     &                    + ZRK(JL,2,JK) * ZRL(JL,6) * ZRL(JL,8)
         ENDDO
      ENDDO
C
C
C     ------------------------------------------------------------------
C
C*         5.     INTRODUCTION OF OZONE ABSORPTION
C                 --------------------------------
C
      JABS=3
      DO JK = 1 , KLEV+1
         DO JL = 1 , KLON
            ZW1(JL) = PUD(JL,JABS,JK)
         ENDDO
C
         CALL SWTT ( KLON, KNU, JABS, ZW1, ZR1 )
C
         DO JL = 1 , KLON
            PFDOWN(JL,JK) = ZR1(JL) * PFDOWN(JL,JK) * SUN(KNU)
            ZW2(JL) = PUM(JL,JK)
         ENDDO
C
         CALL SWTT ( KLON, KNU, JABS, ZW2, ZR2 )
C
         DO JL = 1 , KLON
            PFUP(JL,JK) = ZR2(JL) * PFUP(JL,JK) * SUN(KNU)
         ENDDO
      ENDDO
C
C     ------------------------------------------------------------------
C
      RETURN
      END SUBROUTINE AOSW2S
