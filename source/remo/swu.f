C
C     SUBROUTINE SWU
C
C**** *SWU* - SHORTWAVE RADIATION, ABSORBER AMOUNTS
C
C     PURPOSE.
C     --------
C           COMPUTES THE ABSORBER AMOUNTS USED IN SHORTWAVE RADIATION
C     CALCULATIONS
C
C**   INTERFACE.
C     ----------
C          *SWU* IS CALLED BY *SW*
C
C     *CALL*     SWU ( KLON, KLEV, KAER
C    S              , PSCT,PCARDI,POZ,PPMB,PPSOL,PRMU0,PTAVE,PWV
C    S              , PAKI,PDSIG,PFACT,PRMU,PSEC,PUD,PUM              )
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
C          1. COMPUTES ABSORBER AMOUNTS WITH TEMPERATURE AND PRESSURE
C     SCALING.
C
C     EXTERNALS.
C     ----------
C
C          *SWTT*
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
C     ------------------------------------------------------------------
C
      SUBROUTINE SWU ( KLON, KLEV
     &              , PSCT,PCARDI,PPMB,PPSOL,PRMU0,PTAVE,PWV
     &              , PAKI,PDSIG,PFACT,PRMU,PSEC,PUD )
C
      IMPLICIT NONE
C
      INCLUDE "YOMSW"
      INCLUDE "YOMRDU"
C     ------------------------------------------------------------------
C
C*       0.1   DUMMY ARGUMENTS
C              ---------
C
      INTEGER, INTENT(IN)    :: KLON, KLEV
C
      REAL,    INTENT(IN)    :: PSCT, PCARDI
C
      REAL,    INTENT(IN)    :: PPMB(KLON,KLEV+1),PPSOL(KLON),
     &                          PRMU0(KLON),PTAVE(KLON,KLEV),
     &                          PWV(KLON,KLEV)
C
      REAL,    INTENT(INOUT) :: PAKI(KLON,2),PDSIG(KLON,KLEV),
     &                          PFACT(KLON),PRMU(KLON),
     &                          PSEC(KLON), PUD(KLON,3,KLEV+1)
C
C     ------------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS
C              ------------
C
      REAL    :: ZR1(KLON), ZR2(KLON), ZU1D(KLON)
     &  ,        ZU2D(KLON), ZO175(KLON)
     &  ,        ZO190(KLON)
     &  ,        ZSIGO(KLON)
      REAL    :: ZWH2O,ZSIGN,ZRTU,ZRTH,ZN190,ZN175,
     &           ZDSH2O,ZDSCO2
      INTEGER :: JL,JKP1,JK,JKL
C
C     ------------------------------------------------------------------
C
C*       0.3   FUNCTIONS
C              ---------
C
C     ------------------------------------------------------------------
C
C*         1.     COMPUTES AMOUNTS OF ABSORBERS
C                 -----------------------------
C
C
C*         1.1    INITIALIZES QUANTITIES
C                 ----------------------
C
      DO JL = 1 , KLON
         PUD(JL,1,KLEV+1)=0.
         PUD(JL,2,KLEV+1)=0.
         PUD(JL,3,KLEV+1)=0.
         PFACT(JL)= PRMU0(JL) * PSCT
         PRMU(JL)=SQRT(1224.* PRMU0(JL) * PRMU0(JL) + 1.) / 35.
         PSEC(JL)=1./PRMU(JL)
      ENDDO
C
C*          1.3    AMOUNTS OF ABSORBERS
C                  --------------------
C
      DO JL= 1 , KLON
         ZU1D(JL) = 0.
         ZU2D(JL) = 0.
         ZO175(JL) = PPSOL(JL)** RPDU1
         ZO190(JL) = PPSOL(JL)** RPDH1
         ZSIGO(JL) = PPSOL(JL)
      ENDDO
C
      DO JK = 1 , KLEV

         JKP1 = JK + 1
         JKL = KLEV+1 - JK
         DO JL = 1 , KLON
            ZRTH=(RTH2O/PTAVE(JL,JK))**RTDH2O
            ZRTU=(RTUMG/PTAVE(JL,JK))**RTDUMG
            ZWH2O = AMAX1 (PWV(JL,JKL) , ZEPSCQ )
            ZSIGN = 100. * PPMB(JL,JKP1)
            PDSIG(JL,JK) = (ZSIGO(JL) - ZSIGN)/PPSOL(JL)
            ZN175 = ZSIGN ** RPDU1
            ZN190 = ZSIGN ** RPDH1
            ZDSCO2 = ZO175(JL) - ZN175
            ZDSH2O = ZO190(JL) - ZN190
            PUD(JL,1,JK) = RPNH * ZDSH2O * ZWH2O  * ZRTH
            PUD(JL,2,JK) = RPNU * ZDSCO2 * PCARDI * ZRTU
            ZU1D(JL) = ZU1D(JL) + PUD(JL,1,JK)
            ZU2D(JL) = ZU2D(JL) + PUD(JL,2,JK)
            ZSIGO(JL) = ZSIGN
            ZO175(JL) = ZN175
            ZO190(JL) = ZN190
         ENDDO
      ENDDO
C
C*         1.4    COMPUTES CLEAR-SKY GREY ABSORPTION COEFFICIENTS
C                 -----------------------------------------------
C
      DO JL = 1 , KLON
         ZU1D(JL) = ZU1D(JL) * PSEC(JL)
         ZU2D(JL) = ZU2D(JL) * PSEC(JL)
      ENDDO
C
      CALL SWTT ( KLON, 2, 1, ZU1D, ZR1 )
C
      DO JL = 1 , KLON
         PAKI(JL,1) = -ALOG  ( ZR1 (JL)) / ZU1D(JL)
      ENDDO
C
      CALL SWTT ( KLON, 2, 2, ZU2D, ZR2 )
C
      DO JL = 1 , KLON
         PAKI(JL,2) = -ALOG  ( ZR2 (JL)) / ZU2D(JL)
      ENDDO
C
C     ------------------------------------------------------------------
C
      RETURN
      END SUBROUTINE SWU
