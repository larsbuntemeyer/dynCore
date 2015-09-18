      SUBROUTINE AOSW1S ( KLON, KLEV, KEWAER, KAERH, KNU
     &  ,  PPMB
     &  ,  PAER,PALBS,PCG,PCLDSW,PDSIG,POMEGA,POZ,PRMU,PSEC,PTAU
     &  ,  PUD,PUM
     &  ,  PFD,PFU                                                 )
C
      IMPLICIT NONE
C
C**** *SW1S* - SHORTWAVE RADIATION, FIRST SPECTRAL INTERVAL
C
C     PURPOSE.
C     --------
C
C          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
C     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).
C
C**   INTERFACE.
C     ----------
C     *CALL*     SW1S ( KLON, KLEV, KAER, KNU
C    &  ,  PAER,PALBS,PCG,PCLDSW,PDSIG,POMEGA,PRMU,PSEC,PTAU
C    &  ,  PUD,PUM
C    &  ,  PFD,PFU                                                 )
C
C          *SW1S* IS CALLED FROM *SW*.
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
C          1. COMPUTES UPWARD AND DOWNWARD FLUXES CORRESPONDING TO
C     CONTINUUM SCATTERING
C          2. MULTIPLY BY OZONE TRANSMISSION FUNCTION
C
C     EXTERNALS.
C     ----------
C
C          *SWR*, *SWTT*
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
C       MODIFIED: ROB VAN DORLAND, KNMI, 95-05-10
C     ------------------------------------------------------------------
C
C      IMPLICIT LOGICAL (L)
C
C
      INCLUDE "YOMSW"
      INCLUDE "COMCON"
C     ------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN) :: KEWAER, KLEV, KLON, KNU, KAERH(KLON,KLEV)
C
      REAL,    INTENT(IN) :: 
     &     PAER(KLON,KLEV,5+KEWAER)
     &  ,  PALBS(KLON,2), PCG(KLON,2,KLEV)
     &  ,  PPMB (KLON,KLEV+1)
     &  ,  PCLDSW(KLON,KLEV), PDSIG(KLON,KLEV)
     &  ,  POMEGA(KLON,2,KLEV), POZ(KLON,KLEV)
     &  ,  PRMU(KLON)        , PSEC(KLON)
     &  ,  PTAU(KLON,2,KLEV)
      REAL, INTENT(INOUT) :: PUD(KLON,3,KLEV+1), PUM(KLON,KLEV+1),
     &                       PFD(KLON,KLEV+1),PFU(KLON,KLEV+1)
C
C     ------------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS
C 
C      Local Variables
C 
      INTEGER :: IKL, IKM1, JAJ, JK, JL
      REAL    :: ZOZFAC
C
      REAL ZCGAZ(KLON,KLEV), ZPIZAZ(KLON,KLEV)
     &  ,  ZRAYL(KLON), ZRAY1(KLON,KLEV+1), ZRAY2(KLON,KLEV+1)
     &  ,  ZREFZ(KLON,2,KLEV+1), ZRJ(KLON,6,KLEV+1)
     &  ,  ZRK(KLON,6,KLEV+1), ZRMUE(KLON,KLEV+1)
     &  ,  ZR1(KLON), ZR2(KLON), ZTAUAZ(KLON,KLEV)
     &  ,  ZTRA1(KLON,KLEV+1), ZTRA2(KLON,KLEV+1)
     &  ,  ZW1(KLON), ZW2(KLON)
C
C     ------------------------------------------------------------------
C
C*         1.     FIRST SPECTRAL INTERVAL (0.25-0.68 MICRON)
C                 ----------------------- ------------------
C
C
C*         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
C                 -----------------------------------------
C
      DO JL = 1 , KLON
         ZRAYL(JL) = CRAY(KNU,1) + PRMU(JL) * (CRAY(KNU,2) + PRMU(JL)
     &            * (CRAY(KNU,3) + PRMU(JL) * (CRAY(KNU,4) + PRMU(JL)
     &            * (CRAY(KNU,5) + PRMU(JL) *  CRAY(KNU,6)       ))))
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
     &  , PAER,PALBS,PCG,PCLDSW,PDSIG,POMEGA,ZRAYL,PSEC,PTAU
     &  , ZCGAZ,ZPIZAZ,ZRAY1,ZRAY2,ZREFZ,ZRJ,ZRK,ZRMUE
     &  , ZTAUAZ,ZTRA1,ZTRA2                                  )
C
C
C     ------------------------------------------------------------------
C
C*         3.    OZONE ABSORPTION
C                ----------------
C
C
C*         3.1   DOWNWARD FLUXES
C                ---------------
C
      JAJ = 2
      ZOZFAC=46.6968/G
C
      DO JL = 1 , KLON
         ZW1(JL)=0.
         ZW2(JL)=0.
         ZR1(JL)=1.
         PFD(JL,KLEV+1)=ZRJ(JL,JAJ,KLEV+1) * SUN(KNU)
         PUD(JL,3,KLEV+1)=0.
      ENDDO
      DO JK = 1 , KLEV
         IKL = KLEV+1-JK
         DO JL = 1 , KLON
            ZW1(JL)=ZW1(JL)+PUD(JL,1,IKL)/ZRMUE(JL,IKL)
            ZW2(JL)=ZW2(JL)+POZ(JL,JK)*ZOZFAC/ZRMUE(JL,IKL)
         ENDDO
C
         CALL SWTT ( KLON, KNU, 1, ZW1, ZR1 )
         CALL SWTT ( KLON, KNU, 3, ZW2, ZR2 )
C
         DO JL = 1 , KLON
            PFD(JL,IKL) = ZR1(JL)*ZR2(JL) * ZRJ(JL,JAJ,IKL) * SUN(KNU)
            PUD(JL,3,IKL) = ZW2(JL)
         ENDDO
      ENDDO
C
C
C*         3.2   UPWARD FLUXES
C                -------------
C
      DO JL = 1 , KLON
         PFU(JL,1) = PALBS(JL,KNU) * PFD(JL,1)
         PUM(JL,1)=PUD(JL,3,1)
      ENDDO
C
      DO JK = 2 , KLEV+1
         IKM1=JK-1
         IKL = KLEV+2-JK
         DO JL = 1 , KLON
            ZW1(JL)=ZW1(JL)+PUD(JL,1,IKM1)*1.66
            ZW2(JL)=ZW2(JL)+POZ(JL,IKL)*ZOZFAC*1.66
         ENDDO
C
         CALL SWTT ( KLON, KNU, 1, ZW1, ZR1 )
         CALL SWTT ( KLON, KNU, 3, ZW2, ZR2 )
C
         DO JL = 1 , KLON
            PFU(JL,JK) = ZR1(JL)*ZR2(JL) * ZRK(JL,JAJ,JK) * SUN(KNU)
            PUM(JL,JK) = ZW2(JL)
         ENDDO
      ENDDO
C
C     ------------------------------------------------------------------
C
      RETURN
      END SUBROUTINE AOSW1S
