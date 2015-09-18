C
C     SUBROUTINE SW
C
C**** *SW* - COMPUTES THE SHORTWAVE RADIATION FLUXES.
C
C     PURPOSE.
C     --------
C
C          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
C     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).
C
C**   INTERFACE.
C     ----------
C     *CALL*     SW ( KLON, KLEV, KAER
C    S              , PSCT, PCARDI, PPSOL, PALBS, PWV, PRMU0, PCG
C    S              , PCLDSW, PDP, POMEGA, POZ, PPMB, PTAU, PTAVE, PAER
C    S              , PHEAT, PFDOWN, PFUP, PFDNN, PFDNV, PFUPN, PFUPV  )
C
C          *SW* IS CALLED FROM *RADITE*
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
C          1. COMPUTES ABSORBER AMOUNTS
C          2. COMPUTES FLUXES IN 1ST SPECTRAL INTERVAL
C          3. COMPUTES FLUXES IN 2ND SPECTRAL INTERVAL
C
C     EXTERNALS.
C     ----------
C
C          *SWU*, *SW1S*, *SW2S*
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
      SUBROUTINE SW ( KLON, KLEV
     &              , PSCT, PCARDI, PPSOL, PALBS, PWV, PQS, PRMU0, PCG
     &              , PCLDSW, POMEGA, POZ, PPMB, PTAU, PTAVE, PAER
     &              , PFDOWN, PFUP  )
C
      IMPLICIT NONE
C
C     ------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
      INTEGER, INTENT(IN)    :: KLON, KLEV
      REAL,    INTENT(IN)    :: PSCT, PCARDI
C
      REAL,    INTENT(IN)    :: 
     &     PPSOL(KLON), PAER(KLON,KLEV,5), PRMU0(KLON)
     &  ,  PWV(KLON,KLEV),PQS(KLON,KLEV)
C
      REAL,    INTENT(IN)    :: 
     &     PALBS(KLON,2), PCG(KLON,2,KLEV), PCLDSW(KLON,KLEV),
     &     POMEGA(KLON,2,KLEV), POZ(KLON,KLEV),
     &     PPMB(KLON,KLEV+1), PTAU(KLON,2,KLEV), PTAVE(KLON,KLEV)
C
      REAL,    INTENT(INOUT) :: PFDOWN(KLON,KLEV+1),PFUP(KLON,KLEV+1)
C
C     ------------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS
C              ------------
C
      REAL :: ZAKI(KLON,2), ZDSIG(KLON,KLEV), ZFACT(KLON)
     &  ,     ZFD(KLON,KLEV+1), ZFDOWN(KLON,KLEV+1)
     &  ,     ZFU(KLON,KLEV+1), ZFUP(KLON,KLEV+1)
     &  ,     ZRMU(KLON), ZSEC(KLON)
     &  ,     ZUD(KLON,3,KLEV+1), ZUM(KLON,KLEV+1)
C
      INTEGER :: INU,JK,JL
C
C     ------------------------------------------------------------------
C
C*         1.     ABSORBER AMOUNTS AND OTHER USEFUL QUANTITIES
C                 --------------------------------------------
C
      CALL SWU (      KLON, KLEV
     &              , PSCT,PCARDI,PPMB,PPSOL,PRMU0,PTAVE,PWV
     &              , ZAKI,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD )
C
C
C     ------------------------------------------------------------------

C
C*         2.     FIRST SPECTRAL INTERVAL (0.25-0.68 MICRON)
C                 ----------------------- ------------------
C
      INU = 1
C
      CALL SW1S ( KLON, KLEV, INU
     &  ,  PAER,PALBS,PCG,PCLDSW,ZDSIG,POMEGA,POZ,ZRMU,ZSEC,PTAU
     &  ,  ZUD,ZUM
     &  ,  ZFD,ZFU                                                 )
C
C
C     ------------------------------------------------------------------
C
C*         3.     SECOND SPECTRAL INTERVAL (0.68-4.00 MICRON)

C                 -------------------------------------------
C
      INU = 2
C
      CALL SW2S ( KLON, KLEV, INU
     &  ,  PAER,ZAKI,PALBS,PCG,PCLDSW,ZDSIG,POMEGA,ZRMU,ZSEC,PTAU
     &  ,  ZUD,ZUM,PWV,PQS
     &  ,  ZFDOWN,ZFUP                                                )
C
C
C     ------------------------------------------------------------------
C
C*         4.     FILL THE DIAGNOSTIC ARRAYS
C                 --------------------------
C
      DO JK = 1 , KLEV+1
         DO JL = 1 , KLON
            PFUP(JL,JK)   = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
            PFDOWN(JL,JK) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
         ENDDO
      ENDDO
C
C     ------------------------------------------------------------------
C
      RETURN
      END SUBROUTINE SW
