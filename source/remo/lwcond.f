      SUBROUTINE LWCOND (KIDIA,KFDIA,KLON,KTDIA,KLEV,KLEVP1
     &                 , TWODT, NSTART, NSTEP
C----------------------------------------------------------------------
C - INPUT  2D .
     &                 , KLAB,     PAPHM1, PAPHP1
     &                 , PAPM1,    PAPP1,    PGEOM1, PQM1,   PTM1
     &                 , PXM1,     PXTEC
C - INPUT  1D .
     &                 , LALAND,   KTYPE
C - OUTPUT 2D .
     &                 , PACLC,    PACLCAC
C - OUTPUT 1D .
     &                 , PACLCOV,  PALWCVI,  PAPRL,  PQVI, PSSFL
C - INPUT/OUTPUT 2D .
     &                 , PTTE,     PQTE,     PXTE
C - INPUT/OUTPUT 1D .
     &                 , PAPRS,    PRSFL
     &                 , PXIM1  ,  PQIVI,  PXITE,  PACDNC, PVERVEL, VCT
C - DIAGNOSTIC FIELD
     &                 , RPRAC)
C
C
C**** *COND* - COMPUTES LARGE-SCALE WATER PHASE CHANGES AND CLOUD COVER.
C
C
C     SUBJECT.
C     --------
C
C          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE THREE
C     PROGNOSTIC VARIABLES T, Q, AND X (CLOUD WATER) DUE TO WATER PHASE
C     CHANGES (CONDENSATION, EVAPORATION OF FALLING PRECIPITATION IN
C     UNSATURATED LAYERS AND MELTING OR FREEZING OF THE FALLING WATER)
C     AND PRECIPITAION FORMATION (COALESCENSE, SEDIMENTATION).
C
C     RAIN, SNOWFALL, SURFACE FLUXES, CLOUD WATER AND CLOUD COVER
C     LATER TO BE USED FOR SOIL PROCESSES AND RADIATION ARE STORED.
C
C
C**   INTERFACE.
C     ----------
C
C     *CALL* *COND*
C
C     INPUT ARGUMENTS.
C     ----- ----------
C  - 2D
C  KLAB     : CONVECTION FLAG (0: NO CONVECTION, 1:  ....
C  PAPHM1   : PRESSURE AT HALF LEVELS (T-DT)
C  PAPHP1   : PRESSURE AT HALF LEVELS (T+DT)
C  PAPM1    : PRESSURE AT FULL LEVELS (T-DT)
C  PAPP1    : PRESSURE AT FULL LEVELS (T+DT)
C  PGEOM1   : GEOPOTENTIAL AT FULL LEVELS (T-DT)
C  PQM1     : SPECIFIC HUMIDITY (T-DT)
C  PTM1     : TEMPERATURE (T-DT)
C  PXM1     : CLOUD WATER (T-DT)
C  PXTEC    : TENDENCY OF DETRAINED CONVECTIVE CLOUD WATER
C  - 1D
C  KTYPE    : TYPE OF CONVECTION
C  LALAND   : LOGICAL MASK FOR LAND
C
C     OUTPUT ARGUMENTS.
C     ------ ----------
C  - 2D
C  PACLC    : CLOUD COVER
C  PACLCAC  : CLOUD COVER, ACCUMULATED
C  - 1D
C  PACLCOV  : TOTAL CLOUD COVER
C  PALWCVI  : VERTICALLY INTEGRATED CLOUD WATER, ACCUMULATED
C  PAPRL    : STRATIFORM PRECIPITATION, ACCUMULATED
C  PQVI     : VERTICALLY INTEGRATED SPEC. HUMIDITY, ACCUMULATED
C  PSSFL    : SURFACE SNOW FLUX
C
C     INPUT/OUTPUT ARGUMENTS.
C     ------------ ----------
C  - 2D
C  PTTE     : TENDENCY OF TEMPERATURE
C  PQTE     : TENDENCY OF SPECIFIC HUMIDITY
C  PXTE     : TENDENCY OF CLOUD WATER
C  - 1D
C  PAPRS    : SNOW FALL, ACCUMULATED
C  PRSFL    : SURFACE RAIN FLUX
C
C
C     EXTERNALS.
C     ----------
C
C     METHOD.
C     -------
C
C          SEE REFERENCES
C
C     REFERENCES.
C     ----------
C
C     1. LARGE SCALE PHASE CHANGES' PART OF THE MODEL'S DOCUMENTATION
C     2. ROECKNER, E. AND U. SCHLESE (1985), ECMWF-WORKSHOP ON
C        "CLOUD COVER PARAMETERIZATION IN NUMERICAL MODELS", PP87-108.
C     3. ROECKNER ET AL. (1991) ECMWF/WRCP WORKSHOP ON
C        "CLOUDS, RADIATIVE TRANSFER AND THE HYDROLOGICAL CYCLE",
C         199-222, ECMWF, READING, U.K.
C     4. SUNDQVIST, H. (1978), QJRMS, 104, 677-690.
C
C     AUTHOR.
C     -------
C     U.SCHLESE     DKRZ-HAMBURG   NOV-92
C
C     MODIFICATIONS.
C     --------------
C
      IMPLICIT NONE
C
      INCLUDE "COMCON"
      INCLUDE "YOTLUC"
C---------------------------------------------
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: KIDIA,KFDIA,KLON,KTDIA,KLEV,KLEVP1,
     &                          NSTART,NSTEP
      REAL,    INTENT(IN)    :: TWODT
C
      INTEGER, INTENT(IN)    :: KLAB(KLON,KLEV), KTYPE(KLON)
      REAL,    INTENT(IN)    :: PAPHM1(KLON,KLEVP1), PAPHP1(KLON,KLEVP1)
     &    ,PAPM1(KLON,KLEV),    PAPP1(KLON,KLEV), PQM1(KLON,KLEV)
     &    ,PGEOM1(KLON,KLEV)
     &    ,PTM1(KLON,KLEV),     PXM1(KLON,KLEV),PXTEC(KLON,KLEV)
      REAL,    INTENT(INOUT) :: PACLC(KLON,KLEV), PACLCAC(KLON,KLEV)
      REAL,    INTENT(INOUT) :: PACLCOV(KLON), PALWCVI(KLON), 
     &                          PQVI(KLON), PSSFL(KLON), PAPRL(KLON)
      REAL,    INTENT(INOUT) :: PTTE(KLON,KLEV), PQTE(KLON,KLEV),
     &                          PXTE(KLON,KLEV)
      REAL,    INTENT(INOUT) :: PAPRS(KLON), PRSFL(KLON)
      REAL,    INTENT(IN)    :: PXIM1(KLON,KLEV), PVERVEL(KLON,KLEV)
      REAL,    INTENT(INOUT) :: PXITE(KLON,KLEV), PQIVI(KLON)
      REAL,    INTENT(IN)    :: PACDNC(KLON,KLEV)
      REAL,    INTENT(IN)    :: VCT(2*KLEVP1)
      LOGICAL, INTENT(IN)    :: LALAND(KLON)
      REAL,    INTENT(INOUT) :: RPRAC(KLON,KLEV)
C---------------------------------------------
C     Local Variables
C
C     TEMPORARY ARRAYS
C
      INTEGER :: 
     &     INVB(KLON)
      REAL    ::
     &     ZCLCP1(KLON,KLEV),
     &     ZSAT(KLON,KLEV)  , ZACLCOV(KLON)    , ZQVI(KLON)        ,
     &     ZDP(KLON)        ,
     &     ZGEOH(KLON,KLEVP1),
     &     ZLSDCP(KLON)     , ZLVDCP(KLON)     ,
     &     ZQP1(KLON,KLEV)  ,
     &     ZQSP1(KLON,KLEV) ,
     &     ZRFL(KLON)       , ZRHC(KLON,KLEV)  , ZSFL(KLON)        ,
     &     ZTP1(KLON,KLEV)  ,
     &     ZDTHMIN(KLON)    , ZTHETA(KLON,KLEV)
      REAL    :: 
     &     ZCND(KLON)       ,
     &     ZDEP(KLON)       , ZRHO(KLON,KLEV)  , ZFRL(KLON)        ,
     &     ZXLEVAP(KLON)     ,
     &     ZXLTE(KLON)      ,  ZCLCAUX(KLON)   , ZXLB(KLON)        ,
     &     ZTP1TMP(KLON)    ,  ZQP1TMP(KLON)   ,
     &     ZCLCPRE(KLON)    ,  ZSACL(KLON)     , ZXIMLT(KLON)      ,
     &     ZXISUB(KLON)     ,  ZDZ(KLON)

      REAL    :: 
     &     ZXIEVAP(KLON)    , ZIMLT(KLON)      , ZXITE(KLON)       ,
     &     ZXIFLUX(KLON)    ,
     &     ZXIB(KLON)       , ZQIVI(KLON)      , ZRPR(KLON)        ,
     &     ZQLVI(KLON)      , ZEVP(KLON)       , ZSUB(KLON)        ,
     &     ZSPR(KLON)       , ZSMLT(KLON)
C
      REAL :: ZA, ZB, ZPH(KLEVP1), ZP(KLEV), ZH(KLEV)      
      REAL :: ZQCDIF, ZQRHO, ZDEPCOR
C
CJ-PP
      LOGICAL :: LO1,LO
      LOGICAL :: LO2,LOCC
CSP
      REAL :: CAULOC, CCRAUT, CCSACL, CCSAUT, CCWMIN, CEFFMAX, CEFFMIN, 
     &        CLMAX, CLMIN, CN0S, CQTMIN, CRHOI, CRHOSNO, CSECFRL,
     &        CVTFALL, RDCPD, ZAL1, ZAL2, ZAST, CTHOMI
      REAL :: ZAULOC, ZB1, ZB2, ZBST, ZC1, ZCFAC4C, ZCLAMBS, ZCLCSTAR, 
     &        ZCLOIA, ZCLP1, ZCNDCOR, ZCOEFF, ZCOLLEFFI, ZCOND, ZCONS, 
     &        ZCONS1, ZCONS2, ZCOR, ZDEPOS, ZDIAGT
      REAL :: ZDIAGW, ZDPG, ZDQSAT, ZDQSDT, ZDT2, ZDTDT, ZDTDTSTAR,  
     &        ZDV, ZDXICOR, ZDXLCOR, ZEPQP1, ZEPSEC, ZEPZWP, ZES, ZESAT, 
     &        ZESI, ZESW, ZEXM1, ZEXP, ZDTHDP
      REAL :: ZF1, ZIFRAC, ZLAMSM, ZLC, ZLCDQSDT, ZOVERSAT, ZPREDEL, 
     &        ZPRESUM, ZPRETOT, ZQCON, ZQR, ZQSEC, ZQSI, ZQSM1, 
     &        ZQST1, ZQSW, ZQVDT, ZRAC1, ZRAC2, ZQSP1TMP
      REAL :: ZRADL, ZRAUT, ZRCP, ZRELHUM, ZRHSC, ZRHTEST, ZRI, ZRIEFF, 
     &        ZRIH, ZRS, ZRTC, ZRTL, ZSACI1, ZSACI2, ZSACL1, ZSACL2, 
     &        ZSATSC, ZSAUT, ZSNMLT, ZSUBI
      REAL :: ZSUSATI, ZSUSATW, ZTDIF, ZTMST, ZXIDT, ZXIDTSTAR, ZXIFALL, 
     &        ZXILB, ZXIM1EVP, ZXIMELT, ZXIOLD, ZXIP1, ZXISED, ZXLDT, 
     &        ZXLDTSTAR, ZXLM1EVP, ZXLOLD, ZXLP1, ZXRP1, ZXSEC
      REAL :: ZXSP1, ZXSP2, ZZDRR, ZZDRS, ZZEPR, ZZEPS
C
      INTEGER :: IT, IT1, JB, JBMAX, JBMIN, JK, JL, KLEV2, KLEV2P1,
     &           NEXC, NEXL, KLEVM1
C
C     TODO : Make the following variables save, so we do not have to initialize
C            them at each call!
C
      CQTMIN=1.E-12
      CRHOSNO=100.
      CN0S=3.E6
      CSECFRL=5.E-7
      CCWMIN=5.E-7
      CTHOMI=TMELT-35.
      CCSACL=0.1
      CCSAUT=95.
      CRHOI=500.
CBE   CAULOC=2.0 5 ori
      CAULOC=5.0     !RESOLUTION DEPENDENT IN ECHAM5
CBE CCRAUT=7.  15 ori
      CCRAUT=15.     !TUNING FACTOR FOR AUTOCONVERSION RATE
      CEFFMIN=10.    !MINIMUM EFFECTIVE RADIUS
      CEFFMAX=150.   !MAXIMUM EFFECTIVE RADIUS
      CLMAX=0.5
      CLMIN=0.0
CSP
C
C
C
C
C*    SECURITY PARAMETERS.
C     --------------------
C
C     ZEPFLM IS A MINIMUM FLUX TO AVOID DIVIDING BY ZERO IN THE ICE
C     PROPORTION CALCULATIONS.
C
      ZEPSEC=1.E-12
      ZEPQP1=0.
      ZEPZWP=1.E-20
CSP
      ZXSEC=1.0-ZEPSEC
      ZQSEC=1.0-CQTMIN
      CVTFALL = 3.29
CSP
C
C*    COMPUTATIONAL CONSTANTS.
C     ------------- ----------
C
      ZTMST=TWODT
      IF (NSTEP.EQ.NSTART) ZTMST=0.5*TWODT
      ZDIAGT=0.5*TWODT
      ZDIAGW=ZDIAGT/RHOH2O
C
      ZCONS1=CPD*VTMPC2
      ZCONS2=1./(ZTMST*G)
C
C ATTENTION: ZSATSC>ZRHSC
C
      ZRS=0.99
C
CBE    von 0.7 zu 0.98
      ZRHSC=0.7
      ZSATSC=0.99
CBE   von 0.7 zu 0.98
      ZRTC=0.7
      ZRTL=0.7
      NEXC=4
      NEXL=4
C
      ZC1=2.*ZTMST
      KLEV2=KLEV/2
      KLEVM1=KLEV-1
      KLEV2P1=KLEV2+1
      RDCPD=RD/CPD
      ZCLOIA=1.0E+02
C
C
C
COLS      IROW=NROW(1)
C
C     ------------------------------------------------------------------
C
C*         2.  TOP BOUNDARY CONDITIONS AND QUANTITIES NEEDED FOR
C*             --- -------- ---------- --- ---------- ------ ---
C*                CONDENSATION AND PRECIPITATION CALCULATIONS.
C*                ------------ --- ------------- ------------
C
C*         2.1     SET TO ZERO PRECIPITATION FLUXES AT THE TOP.
C
      DO JL=KIDIA,KFDIA
         ZCLCPRE(JL)=0.
         ZRFL(JL)=0.
         ZSFL(JL)=0.
         ZDTHMIN(JL)=0.0
         INVB(JL)=1
         ZXIFLUX(JL)=0.
      ENDDO
C
C     TODO: 2.2 and 2.3 might be merged.
C
C*        2.2      CALCULATE POTENTIAL TEMPERATURES
C
C
      DO JK=KLEV,KLEV2,-1
         DO JL=KIDIA,KFDIA
            ZTHETA(JL,JK)=PTM1(JL,JK)*(1.0E5/PAPM1(JL,JK))**RDCPD
         ENDDO
      ENDDO
C
C*         2.3    CHECK FOR PRESENCE OF LOW-LEVEL INVERSION
C                  ( SEA POINTS ONLY )
C
      KLEV2P1=KLEV2+1
      DO JK=KLEV,KLEV2P1,-1
         DO JL=KIDIA,KFDIA
            IF((.NOT.LALAND(JL)) .AND. KTYPE(JL).EQ.0) THEN
               ZDTHDP=(ZTHETA(JL,JK)-
     &              ZTHETA(JL,JK-1))*ZCLOIA/(PAPM1(JL,JK)
     &              -PAPM1(JL,JK-1))
               LO=ZDTHDP.LT.ZDTHMIN(JL)
               IF (LO) THEN
                  ZDTHMIN(JL)=ZDTHDP
                  INVB(JL)=JK
               ENDIF
            ENDIF
         ENDDO
      ENDDO
C

C Todo: Merge these K loops
CSP HALF LEVEL PRESSURE VALUES, ASSUMING 101320 PA SURFACE PRESSURE

      DO JK=1,KLEVP1
         ZA=VCT(JK)
         ZB=VCT(JK+KLEVP1)
         ZPH(JK)=ZA+ZB*101320.
      ENDDO

CSP   FULL LEVEL PRESSURE

      DO JK = 1, KLEV
         ZP(JK)=(ZPH(JK)+ZPH(JK+1))*0.5
      ENDDO

      DO JK = 1, KLEV
         ZH(JK)=(ZPH(KLEVP1)-ZP(JK))/(G*1.25)
      ENDDO

CSP SEARCH FOR HIGHEST INVERSION LEVEL (FIRST FULL LEVEL BELOW 1000 M)

      DO JK = 1, KLEV
         JBMIN=JK
         IF(ZH(JK).LT.1000.) EXIT
      ENDDO

CSP -- SEARCH FOR LOWEST INVERSION LEVEL (FIRST FULL LEVEL BELOW 500 M)

      DO JK = 1, KLEV
         JBMAX=JK
         IF(ZH(JK).LT.500.) EXIT
      ENDDO


C
      DO JK=KTDIA,KLEV
         DO JL=KIDIA,KFDIA

C
C*     2.4   T, Q AND QS PROVISIONAL VALUES AT T+DT
C*            EFFECTIVES L FOR VAPORISATION AND SUBLIMATION
C

            ZTP1(JL,JK)=PTM1(JL,JK)+ZTMST*PTTE(JL,JK)
            ZQP1(JL,JK)=PQM1(JL,JK)+ZTMST*PQTE(JL,JK)
            ZQP1(JL,JK)=AMAX1(ZQP1(JL,JK),ZEPQP1)
C
            ZRCP=1./(CPD+ZCONS1*ZQP1(JL,JK))
            ZLVDCP(JL)=ALV*ZRCP
            ZLSDCP(JL)=ALS*ZRCP
C
            IT=INT(ZTP1(JL,JK)*1000.)
            ZQSP1(JL,JK)=TLUCUA(IT)/PAPP1(JL,JK)
            ZQSP1(JL,JK)=AMIN1(ZQSP1(JL,JK),0.5)
            ZCOR=1./(1.-VTMPC1*ZQSP1(JL,JK))
            ZQSP1(JL,JK)=ZQSP1(JL,JK)*ZCOR
C
C*     2.5    THRESHOLD RELATIVE HUMIDITY
C
            JB=INVB(JL)
            LO=(JB.GT.KLEV-6.AND. JB.LT.KLEVM1)
            IF(KLAB(JL,JK).EQ.2) THEN
               ZRHC(JL,JK)=ZRTC+(ZRS-ZRTC)*EXP(1.-
     &              (PAPHP1(JL,KLEVP1)/PAPP1(JL,JK))**NEXC)
               ZSAT(JL,JK)=1.
            ELSE
CCC  ACHTUNG AENDERUNG
               IF(JK.GE.JB+1 .AND. LO) THEN
                  ZRHC(JL,JK)=ZRHSC
                  ZSAT(JL,JK)=ZSATSC
               ELSE
                  ZRHC(JL,JK)=ZRTL+(ZRS-ZRTL)*EXP(1.-
     &                 (PAPHP1(JL,KLEVP1)/PAPP1(JL,JK))**NEXL)
                  ZSAT(JL,JK)=1.
               ENDIF
            ENDIF
C
C*    2.6  CLOUD COVER AT T+DT
C
            ZQR=ZQSP1(JL,JK)*ZSAT(JL,JK)*ZRHC(JL,JK)
            ZCLCP1(JL,JK)=(ZQP1(JL,JK)-ZQR)/
     &                   (ZQSP1(JL,JK)*ZSAT(JL,JK)-ZQR)
            ZCLCP1(JL,JK)=AMAX1(ZCLCP1(JL,JK),0.)
            ZCLCP1(JL,JK)=AMIN1(ZCLCP1(JL,JK),1.)
            ZCLCP1(JL,JK)=1.-SQRT(1.-ZCLCP1(JL,JK))
C
C
         ENDDO
      ENDDO
C
C
C*      2.8 CLOUD TOPS/BASES AND THICKNESS.
C
C      2.81 GEOPOTENTIAL AT HALF LEVELS
C
      DO JK=KTDIA,KLEV-1
         DO JL=KIDIA,KFDIA
            ZGEOH(JL,JK+1)=(PGEOM1(JL,JK+1)+PGEOM1(JL,JK))*0.5
         ENDDO
      ENDDO

      DO JL=KIDIA,KFDIA
         ZGEOH(JL,1)=2.*PGEOM1(JL,1)-ZGEOH(JL,2)
         ZGEOH(JL,KLEVP1)=0.
      ENDDO
C
C***
      DO JK=KTDIA,KLEV

C***
CSP     2.9. SET TO ZERO LOCAL TENDENCIES
         DO JL=KIDIA,KFDIA
            ZCND(JL)=0.0
            ZDEP(JL)=0.0
            ZFRL(JL)=0.0
            ZIMLT(JL)=0.0
            ZXIEVAP(JL)=0.0
            ZXLEVAP(JL)=0.0
            ZRPR(JL)=0.0
            ZSPR(JL)=0.0
            ZSMLT(JL)=0.0
            ZEVP(JL)=0.0
            ZSUB(JL)=0.0
            ZXIMLT(JL)=0.0
            ZXISUB(JL)=0.0
            PACLC(JL,JK)=ZCLCP1(JL,JK)
            ZSACL(JL)=0.0
            ZDP(JL)=PAPHM1(JL,JK+1)-PAPHM1(JL,JK)
            ZRHO(JL,JK)=PAPP1(JL,JK)/(RD*ZTP1(JL,JK))
            ZDZ(JL)=(ZGEOH(JL,JK)-ZGEOH(JL,JK+1))/G
         ENDDO

CSP
CBE  CHANGE THE LOCATION OF THE ADV FROM 2.9 TO 7.0
C
C     ------------------------------------------------------------------
C
C*         3.   MODIFICATION OF INCOMING PRECIPITATION FLUXES BY MELTING,
C               SUBLIMATION AND EVAPORATION
C               ---- ---- --- ---------- -------- ---------- -- ----
C
!DIR$ IVDEP
         IF (JK .GT. 1) THEN
            DO JL=KIDIA,KFDIA
C
C*         3.1   MELTING OF SNOW AND ICE
C

               ZCONS=ZCONS2*ZDP(JL)/(ZLSDCP(JL)-ZLVDCP(JL))
               ZTDIF=MAX(0.0,PTM1(JL,JK)-TMELT)
               ZSNMLT=MIN(ZXSEC*ZSFL(JL),ZCONS*ZTDIF)
               ZSFL(JL)=ZSFL(JL)-ZSNMLT
               ZRFL(JL)=ZRFL(JL)+ZSNMLT
               ZSMLT(JL)=ZSNMLT/(ZCONS2*ZDP(JL))
               ZXIMELT=MIN(ZXSEC*ZXIFLUX(JL),ZCONS*ZTDIF)
               ZXIFLUX(JL)=ZXIFLUX(JL)-ZXIMELT
               ZXIMLT(JL)=ZXIMELT/(ZCONS2*ZDP(JL))
               IF (ZTDIF .GT. 0.0) THEN
                  ZIMLT(JL)=MAX(0.0,PXIM1(JL,JK)+PXITE(JL,JK)*ZTMST)
               ELSE
                  ZIMLT(JL)=0.0
               END IF


C        3.2 SUBLIMATION OF SNOW AND ICE (LIN ET. AL, 1983)

               IF (ZCLCPRE(JL) .GT. 0.0) THEN
                  ZCLCSTAR=ZCLCPRE(JL)
                  ZDPG=ZDP(JL)/G
                  ZQRHO=1.3/ZRHO(JL,JK)
                  IT=NINT(PTM1(JL,JK)*1000.)
CSP HIER IST IN ECHAM EINE ABFRAGE AUF LOOKUPOVERFLOW, WENN IT
CSP JPTLUCU1 ODER > LPTLUCU2
                  IT=MAX(MIN(IT, JPTLUCU2),JPTLUCU1)
                  ZESI=TLUCUA(IT)/PAPM1(JL,JK)
                  ZESI=MIN(ZESI, 0.5)
                  ZQSI=ZESI/(1.-VTMPC1*ZESI) !*VTMPC1=RV/RD-1.
                  ZSUSATI=MIN(PQM1(JL,JK)/ZQSI-1.0,0.0)
                  ZB1=ZLSDCP(JL)**2/(2.43E-2*RV*(PTM1(JL,JK)**2))
                  ZB2=1./(ZRHO(JL,JK)*ZQSI*0.211E-4)
                  ZCOEFF=3.E6*2.*API*ZSUSATI/(ZRHO(JL,JK)*(ZB1+ZB2))
                  IF (ZSFL(JL) .GT. CQTMIN) THEN
                     ZXSP1=(ZSFL(JL)/(ZCLCPRE(JL)*CVTFALL))**(1./1.16)
                     ZCLAMBS=(ZXSP1/(API*CRHOSNO*CN0S))**0.25
                     ZCFAC4C=0.78*ZCLAMBS**2+232.19*ZQRHO**0.25
     &                    *ZCLAMBS**2.625
                     ZZEPS=MAX(-ZXSEC*ZSFL(JL)/ZCLCPRE(JL),
     &                    ZCOEFF*ZCFAC4C*ZDPG)
                     ZSUB(JL)=-ZZEPS/ZDPG*ZTMST*ZCLCSTAR
                     ZSUB(JL)=MIN(ZSUB(JL),
     &                    MAX(ZXSEC*(ZQSI-PQM1(JL,JK)),0.0))
                     ZSUB(JL)=MAX(ZSUB(JL),0.0)
                  ENDIF
C
                  IF (ZXIFLUX(JL) .GT. CQTMIN) THEN
                     ZXSP1=(ZXIFLUX(JL)/
     &                    (ZCLCPRE(JL)*CVTFALL))**(1./1.16)
                     ZCLAMBS=(ZXSP1/(API*CRHOSNO*CN0S))**0.25
                     ZCFAC4C=0.78*ZCLAMBS**2+232.19*ZQRHO**0.25
     &                    *ZCLAMBS**2.625
                     ZZEPS=MAX(-ZXSEC*ZXIFLUX(JL)/ZCLCPRE(JL),
     &                    ZCOEFF*ZCFAC4C*ZDPG)
                     ZSUBI=-ZZEPS/ZDPG*ZTMST*ZCLCSTAR
                     ZSUBI=MIN(ZSUBI,MAX(ZXSEC*(ZQSI-PQM1(JL,JK)),0.0))
                     ZSUBI=MAX(ZSUBI,0.0)
                     ZXIFLUX(JL)=ZXIFLUX(JL)-ZSUBI*ZCONS2*ZDP(JL)
                     ZXISUB(JL)=ZSUBI
                  ENDIF
               ENDIF
C
C        3.3  EVAPORATION OF RAIN  (ROTSTAYN, 1997)
C
               IF (ZCLCPRE(JL) .GT. 0.0 .AND. ZRFL(JL) .GT. CQTMIN) THEN
                  ZCLCSTAR=ZCLCPRE(JL)
                  ZDPG=ZDP(JL)/G
                  ZQRHO=1.3/ZRHO(JL,JK)
                  ZXRP1=(ZRFL(JL)/(ZCLCPRE(JL)*12.45*SQRT(ZQRHO)))
     &                 **(8./9.)
                  IT=NINT(PTM1(JL,JK)*1000.)
CSP HIER IST IN ECHAM EINE ABFRAGE AUF LOOKUPOVERFLOW, WENN
CSP   IT < JPTLUCU1 ODER > LPTLUCU2
                  IT=MAX(MIN(IT,JPTLUCU2),JPTLUCU1)
                  ZESW=TLUCUAW(IT)/PAPM1(JL,JK)
                  ZESAT=ZESW*PAPM1(JL,JK)*RV/RD
                  ZESW=MIN(ZESW,0.5)
                  ZQSW=ZESW/(1.-VTMPC1*ZESW)
                  ZSUSATW=MIN(PQM1(JL,JK)/ZQSW-1.0,0.0)
                  ZDV=2.21/PAPM1(JL,JK)
                  ZAST=ALV*(ALV/(RV*PTM1(JL,JK))-1.0)/
     &                 (0.024*PTM1(JL,JK))
                  ZBST=RV*PTM1(JL,JK)/(ZDV*ZESAT)
                  ZZEPR=870.*ZSUSATW*(ZRFL(JL)/ZCLCPRE(JL))**0.61
     &                 /(SQRT(ZRHO(JL,JK))*(ZAST+ZBST))
                  ZZEPR=MAX(-ZXSEC*ZRFL(JL)/ZCLCPRE(JL),ZZEPR*ZDPG)
                  ZEVP(JL)=-ZZEPR/ZDPG*ZTMST*ZCLCSTAR
                  ZEVP(JL)=MIN(ZEVP(JL),
     &                 MAX(ZXSEC*(ZQSW-PQM1(JL,JK)),0.0))
                  ZEVP(JL)=MAX(ZEVP(JL),0.0)
               ENDIF
            ENDDO
         ENDIF

         DO JL=KIDIA,KFDIA

C      4.SEDIMENTATION OF CLOUD ICE FROM GRID-MEAN VALUES.
C        UPDATING THE TENDENCY 'PXITE' TO INCLUDE SEDIMENTATION
C        AT JK=KLEV THE SEDIMENTATION SINK IS BALANCED BY
C        PRECIPITATION AT THE SURFACE (THROUGH ZZDRS, SEE 7.3).
C        FINALLY: IN-CLOUD CLOUD WATER/ICE.
C
C TODO: Save inverse of ZCLCAUX, ZTMST, PAPP1
C

            ZXIP1=PXIM1(JL,JK)+PXITE(JL,JK)*ZTMST-ZIMLT(JL)
            ZXIP1=MAX(ZXIP1,EPSILON(1.))
            ZXIFALL=CVTFALL*(ZRHO(JL,JK)*ZXIP1)**0.16
            ZAL1=ZXIFALL*G*ZRHO(JL,JK)*ZTMST/ZDP(JL)
            ZAL2=ZXIFLUX(JL)/(ZRHO(JL,JK)*ZXIFALL)
            ZXISED=ZXIP1*EXP(-ZAL1)+ZAL2*(1.-EXP(-ZAL1))
            ZXIFLUX(JL)=ZXIFLUX(JL)+(ZXIP1-ZXISED)*ZCONS2*ZDP(JL)
            PXITE(JL,JK)=(ZXISED-PXIM1(JL,JK))/ZTMST


C      IN-CLOUD WATER/ICE CALCULATED FROM RESPECTIVE GRID-MEANS
C      PARTIAL CLOUD COVER, ADVECTIVE/DIFFUSIVE TENDENCIES, DETRAINED
C      CLOUD WATER/ICE AND ICE SEDIMENTATION.
C      IN-CLOUD VALUES ARE REQUIRED FOR CLOUD MICROPHYSICS
C
            ZCLCAUX(JL)=PACLC(JL,JK)
            LOCC=ZCLCAUX(JL) .GT. 0.0
            LO2=(PTM1(JL,JK) .LT. CTHOMI) .OR.
     &           (PTM1(JL,JK) .LT. TMELT .AND. ZXISED .GT. CSECFRL)
            IF (LO2) THEN
               ZXITE(JL)=PXTEC(JL,JK)
               ZXLTE(JL)=0.0
               IF (LOCC) THEN
                  ZXIB(JL)=PXIM1(JL,JK)/ZCLCAUX(JL)
                  ZXLB(JL)=PXM1(JL,JK)/ZCLCAUX(JL)
                  ZXIM1EVP=0.0
                  ZXLM1EVP=0.0
                  ZXIDT=(PXITE(JL,JK)+ZXITE(JL))*ZTMST
                  ZXLDT=PXTE(JL,JK)*ZTMST+ZXIMLT(JL)+ZIMLT(JL)
                  IF(ZXIDT .GT. 0.0) THEN
                     ZXIDTSTAR=ZXIDT
                     ZXIB(JL)=ZXIB(JL)+ZXIDT
                  ELSE
                     ZXIDTSTAR=0.0
                     ZXIB(JL)=ZXIB(JL)+MAX(ZXIDT/ZCLCAUX(JL),-ZXIB(JL))
                     PXITE(JL,JK)=MAX(PXITE(JL,JK),
     &                    -(PXIM1(JL,JK)/ZTMST+ZXITE(JL)))
                  ENDIF
                  IF(ZXLDT .GT. 0.0) THEN
                     ZXLDTSTAR=ZXLDT
                     ZXLB(JL)=ZXLB(JL)+ZXLDT
                  ELSE
                     ZXLDTSTAR=0.0
                     ZXLB(JL)=ZXLB(JL)+MAX(ZXLDT/ZCLCAUX(JL),-ZXLB(JL))
                     PXTE(JL,JK)=MAX(PXTE(JL,JK),-PXM1(JL,JK)/ZTMST)
                  ENDIF
               ELSE
                  ZXIB(JL)=0.0
                  ZXLB(JL)=0.0
                  ZXIDT=(PXITE(JL,JK)+ZXITE(JL))*ZTMST
                  ZXLDT=PXTE(JL,JK)*ZTMST+ZXIMLT(JL)+ZIMLT(JL)
                  IF (ZXIDT .GT. 0.0) THEN
                     ZXIDTSTAR=ZXIDT
                     ZXIM1EVP=PXIM1(JL,JK)
                  ELSE
                     ZXIDTSTAR=0.0
                     PXITE(JL,JK)=MAX(PXITE(JL,JK),
     &                    -(PXIM1(JL,JK)/ZTMST+ZXITE(JL)))
                     ZXIM1EVP=PXIM1(JL,JK)+
     &                    (PXITE(JL,JK)+ZXITE(JL))*ZTMST
                  ENDIF
                  IF (ZXLDT .GT. 0.0) THEN
                     ZXLDTSTAR=ZXLDT
                     ZXLM1EVP=PXM1(JL,JK)
                  ELSE
                     ZXLDTSTAR=0.0
                     PXTE(JL,JK)=MAX(PXTE(JL,JK),-PXM1(JL,JK)/ZTMST)
                     ZXLM1EVP=PXM1(JL,JK)+PXTE(JL,JK)*ZTMST
                  ENDIF
               ENDIF

            ELSE
               ZXLTE(JL)=PXTEC(JL,JK)
               ZXITE(JL)=0.0
               IF (LOCC) THEN
                  ZXLB(JL)=PXM1(JL,JK)/ZCLCAUX(JL)
                  ZXIB(JL)=PXIM1(JL,JK)/ZCLCAUX(JL)
                  ZXLM1EVP=0.0
                  ZXIM1EVP=0.0
                  ZXLDT=(PXTE(JL,JK)+ZXLTE(JL))*ZTMST+ZXIMLT(JL)
     &                 +ZIMLT(JL)
                  ZXIDT=PXITE(JL,JK)*ZTMST
                  IF (ZXLDT .GT. 0.0) THEN
                     ZXLDTSTAR=ZXLDT
                     ZXLB(JL)=ZXLB(JL)+ZXLDT
                  ELSE
                     ZXLDTSTAR=0.0
                     ZXLB(JL)=ZXLB(JL)+MAX(ZXLDT/ZCLCAUX(JL),
     &                    -ZXLB(JL))
                     PXTE(JL,JK)=MAX(PXTE(JL,JK),
     &                    -(PXM1(JL,JK)/ZTMST+ZXLTE(JL)))
                  ENDIF
                  IF (ZXIDT .GT. 0.0) THEN
                     ZXIDTSTAR=ZXIDT
                     ZXIB(JL)=ZXIB(JL)+ZXIDT
                  ELSE
                     ZXIDTSTAR=0.0
                     ZXIB(JL)=ZXIB(JL)+MAX(ZXIDT/ZCLCAUX(JL),
     &                    -ZXIB(JL))
                     PXITE(JL,JK)=MAX(PXITE(JL,JK),
     &                    -PXIM1(JL,JK)/ZTMST)
                  ENDIF
               ELSE
                  ZXLB(JL)=0.0
                  ZXIB(JL)=0.0
                  ZXLDT=(PXTE(JL,JK)+ZXLTE(JL))*ZTMST+ZXIMLT(JL)
     &                 +ZIMLT(JL)
                  ZXIDT=PXITE(JL,JK)*ZTMST
                  IF (ZXLDT .GT. 0.0) THEN
                     ZXLDTSTAR=ZXLDT
                     ZXLM1EVP=PXM1(JL,JK)
                  ELSE
                     ZXLDTSTAR=0.0
                     PXTE(JL,JK)=MAX(PXTE(JL,JK),
     &                    -(PXM1(JL,JK)/ZTMST+ZXLTE(JL)))
                     ZXLM1EVP=PXM1(JL,JK)+
     &                    (PXTE(JL,JK)+ZXLTE(JL))*ZTMST
                  ENDIF
                  IF (ZXIDT .GT. 0.0) THEN
                     ZXIDTSTAR=ZXIDT
                     ZXIM1EVP=PXIM1(JL,JK)
                  ELSE
                     ZXIDTSTAR=0.0
                     PXITE(JL,JK)=MAX(PXITE(JL,JK),
     &                    -PXIM1(JL,JK)/ZTMST)
                     ZXIM1EVP=PXIM1(JL,JK)+PXITE(JL,JK)*ZTMST
                  ENDIF
               ENDIF
            ENDIF


C     5. CONDENSATION/DEPOSITION AND EVAPORATION/SUBLIMATION
C
C     ZLC = L_{V/S} / C _ P
C     ZLCDQSDT = L DQ_{SAT} / C_P DT
C     ZDQSDT = DQ_{SAT}/DT
C
            ZRCP=1./(CPD+ZCONS1*MAX(PQM1(JL,JK),0.0))
            ZLVDCP(JL)=ALV*ZRCP
            ZLSDCP(JL)=ALS*ZRCP

            IF (LO2) THEN
               ZLC=ZLSDCP(JL)
            ELSE
               ZLC=ZLVDCP(JL)
            ENDIF
            IT=NINT(PTM1(JL,JK)*1000.)
CSP         HIER IST IN ECHAM EINE ABFRAGE AUF LOOKUPOVERFLOW, WENN
CSP         IT < JPTLUCU1 ODER > LPTLUCU2
            IT=MAX(MIN(IT,JPTLUCU2),JPTLUCU1)
            IF(LO2) THEN
               ZQSM1=TLUCUA(IT)/PAPM1(JL,JK)
            ELSE
               ZQSM1=TLUCUAW(IT)/PAPM1(JL,JK)
            ENDIF
            ZQSM1=MIN(ZQSM1,0.5)
            ZQSM1=ZQSM1/(1.-VTMPC1*ZQSM1)
            IT1=IT+1
            IT1=MAX(MIN(IT1,JPTLUCU2),JPTLUCU1)
            IF(LO2) THEN
               ZQST1=TLUCUA(IT1)/PAPM1(JL,JK)
            ELSE
               ZQST1=TLUCUAW(IT1)/PAPM1(JL,JK)
            ENDIF
            ZQST1=MIN(ZQST1,0.5)
            ZQST1=ZQST1/(1.-VTMPC1*ZQST1)
            ZDQSDT=(ZQST1-ZQSM1)*1000.
            ZLCDQSDT=ZLC*ZDQSDT
            ZXIEVAP(JL)=(1.0-ZCLCAUX(JL))*ZXIDTSTAR+ZXIM1EVP
            ZXLEVAP(JL)=(1.0-ZCLCAUX(JL))*ZXLDTSTAR+ZXLM1EVP
            ZQVDT=PQTE(JL,JK)*ZTMST+ZEVP(JL)+ZSUB(JL)
     &           +ZXIEVAP(JL)+ZXLEVAP(JL)+ZXISUB(JL)
            ZDTDT=PTTE(JL,JK)*ZTMST-ZLVDCP(JL)*(ZEVP(JL)+ZXLEVAP(JL))
     &           -ZLSDCP(JL)*(ZSUB(JL)+ZXIEVAP(JL)+ZXISUB(JL))
     &           -(ZLSDCP(JL)-ZLVDCP(JL))
     &           *(ZSMLT(JL)+ZXIMLT(JL)+ZIMLT(JL))
            ZQP1(JL,JK)=MAX(PQM1(JL,JK)+ZQVDT,0.0)
            ZTP1(JL,JK)=PTM1(JL,JK)+ZDTDT
            ZDTDTSTAR=ZDTDT+ZCLCAUX(JL)*(ZLC*PQTE(JL,JK)*ZTMST
     &           +ZLVDCP(JL)*(ZEVP(JL)+ZXLEVAP(JL))
     &           +ZLSDCP(JL)*(ZSUB(JL)+ZXIEVAP(JL)+ZXISUB(JL)))
            ZDQSAT=ZDTDTSTAR*ZDQSDT/(1.+ZCLCAUX(JL)*ZLCDQSDT)
            ZXIB(JL)=MAX(ZXIB(JL),0.0)
            ZXLB(JL)=MAX(ZXLB(JL),0.0)
            ZXILB=ZXIB(JL)+ZXLB(JL)

            ZQCDIF=(ZQVDT-ZDQSAT)*ZCLCAUX(JL)
            ZQCDIF=MAX(ZQCDIF,-ZXILB*ZCLCAUX(JL))
            ZQCDIF=MIN(ZQCDIF,ZQSEC*ZQP1(JL,JK))
            IF (ZQCDIF .LT. 0.0) THEN
               ZIFRAC=ZXIB(JL)/MAX(ZEPSEC,ZXILB)
               ZIFRAC=MAX(MIN(ZIFRAC,1.0),0.0)
               ZDEP(JL)=ZQCDIF*ZIFRAC
               ZCND(JL)=ZQCDIF*(1.0-ZIFRAC)
            ELSE
               IF (LO2) THEN
                  ZDEP(JL)=ZQCDIF
                  ZCND(JL)=0.0
               ELSE
                  ZCND(JL)=ZQCDIF
                  ZDEP(JL)=0.0
               ENDIF
            ENDIF

C     5.4 ACCOUNTING FOR CLOUD EVAPORATION IN CLEAR AIR AND
C         CHECKING FOR SUPERSATURATION
C

            ZTP1TMP(JL)=ZTP1(JL,JK)+ZLVDCP(JL)
     &           *ZCND(JL)+ZLSDCP(JL)*ZDEP(JL)
            ZQP1TMP(JL)=ZQP1(JL,JK)-ZCND(JL)-ZDEP(JL)
            ZXIP1=MAX(ZXISED+ZXITE(JL)*ZTMST-ZXIEVAP(JL)
     &           +ZDEP(JL),0.0)
            LO2=(ZTP1TMP(JL) .LT. CTHOMI) .OR.
     &           (ZTP1TMP(JL) .LT. TMELT .AND. ZXIP1 .GT. CSECFRL)
            IT=NINT(ZTP1TMP(JL)*1000.)
CSP HIER IST IN ECHAM EINE ABFRAGE AUF LOOKUPOVERFLOW, WENN
CSP   IT < JPTLUCU1 ODER > LPTLUCU2
            IT=MAX(MIN(IT,JPTLUCU2),JPTLUCU1)
            IF(LO2) THEN
               ZES=TLUCUA(IT)/PAPP1(JL,JK)
            ELSE
               ZES=TLUCUAW(IT)/PAPP1(JL,JK)
            ENDIF
            ZES=MIN(ZES,0.5)
            LO=ZES .LT. 0.4
            ZCOR=1./(1.-VTMPC1*ZES)
            ZQSP1TMP=ZES*ZCOR
            ZOVERSAT=ZQSP1TMP*0.01
            ZRHTEST=MIN(PQM1(JL,JK)/ZQSM1,1.)*ZQSP1TMP
            IT1=IT+1
            IT1=MAX(MIN(IT1,JPTLUCU2),JPTLUCU1)
            IF(LO2) THEN
               ZQST1=TLUCUA(IT1)/PAPP1(JL,JK)
            ELSE
               ZQST1=TLUCUAW(IT1)/PAPP1(JL,JK)
            ENDIF
            ZQST1=MIN(ZQST1,0.5)
            ZQST1=ZQST1/(1.-VTMPC1*ZQST1)
            ZDQSDT=(ZQST1-ZQSP1TMP)*1000.
            IF(LO2) THEN
               ZLC=ZLSDCP(JL)
            ELSE
               ZLC=ZLVDCP(JL)
            ENDIF
            IF(LO) THEN
               ZLCDQSDT=ZLC*ZDQSDT
            ELSE
               ZLCDQSDT=ZQSP1TMP*ZCOR*TLUCUB(IT)
            ENDIF
            ZQCON=1./(1.+ZLCDQSDT)
            IF (LO2) THEN       ! ICE CLOUD
               IF (ZQP1TMP(JL) .GT. ZQSP1TMP+ZOVERSAT) THEN
                  ZDEPCOR = (ZQP1TMP(JL)-ZQSP1TMP-ZOVERSAT)*ZQCON
                  ZDEP(JL)=ZDEP(JL)+ZDEPCOR
               ENDIF
               IF (ZDEP(JL) .GT. 0.0 .AND. ZQP1TMP(JL) .LT. ZRHTEST
     &              .AND. ZQSP1TMP .LE. ZQSM1) THEN
                  ZDEP(JL)=ZQP1(JL,JK)-ZRHTEST
                  ZDEP(JL)=MAX(ZDEP(JL),0.0)
               ENDIF
            ELSE                ! WATER CLOUD
               IF (ZQP1TMP(JL) .GT. ZQSP1TMP+ZOVERSAT) THEN
                  ZCNDCOR=(ZQP1TMP(JL)-ZQSP1TMP-ZOVERSAT)*ZQCON
                  ZCND(JL)=ZCND(JL)+ZCNDCOR
               ENDIF
               IF (ZCND(JL) .GT. 0.0 .AND. ZQP1TMP(JL) .LT. ZRHTEST
     &              .AND. ZQSP1TMP .LE. ZQSM1) THEN
                  ZCND(JL)=ZQP1(JL,JK)-ZRHTEST
                  ZCND(JL)=MAX(ZCND(JL),0.0)
               ENDIF
            ENDIF

C     5.5 CHANGE OF IN-CLOUD WATER DUE TO DEPOSITION/SUBLIMATION
C         AND CONDENSATION/EVAPORATION (INPUT FOR CLOUD MICROPHYSICS)
C
            ZRELHUM=ZQP1TMP(JL)/ZQSP1TMP
            ZDEPOS=MAX(ZDEP(JL),0.0)
            ZCOND=MAX(ZCND(JL),0.0)
            IF (LOCC) THEN
               ZXIB(JL)=MAX(ZXIB(JL)+ZDEP(JL)/ZCLCAUX(JL),0.0)
               ZXLB(JL)=MAX(ZXLB(JL)+ZCND(JL)/ZCLCAUX(JL),0.0)
            ELSE IF (ZDEPOS .GT. 0.0 .OR. ZCOND .GT. 0.0) THEN
               ZCLCAUX(JL)=MAX(MIN(ZRELHUM,1.0),0.01)
               ZXIB(JL)=ZDEPOS/ZCLCAUX(JL)
               ZXLB(JL)=ZCOND/ZCLCAUX(JL)
            ENDIF

            ZTP1TMP(JL)=ZTP1(JL,JK)+ZLVDCP(JL)*
     &           ZCND(JL)+ZLSDCP(JL)*ZDEP(JL)


C     6.0 FREEZING OD CLOUD WATER
C
C     6.2 FREEZING OF CLOUD WATER FOR T < 238 K

            IF(ZTP1TMP(JL) .LE. CTHOMI) THEN
               ZFRL(JL)=ZXLB(JL)*ZCLCAUX(JL)
               ZXIB(JL)=ZXIB(JL)+ZXLB(JL)
               ZXLB(JL)=0.0
            ENDIF

C     6.3 FREEZING OF CLOUD WATER BETWEEN 238 AND 273 K

            LO = ZXLB(JL).GT. 0.0 .AND. ZTP1TMP(JL) .LT. TMELT
     &           .AND. ZTP1TMP(JL) .GT. CTHOMI
            IF (LO) THEN
               ZFRL(JL)=100.*(EXP(0.66*(TMELT-ZTP1TMP(JL)))-1.)
     &              *ZRHO(JL,JK)/(RHOH2O*PACDNC(JL,JK))
               ZFRL(JL)=ZXLB(JL)*(1.-1./(1.+ZFRL(JL)
     &              *ZTMST*ZXLB(JL)))
               ZRADL=(0.75*ZXLB(JL)*ZRHO(JL,JK)
     &              /(API*RHOH2O*PACDNC(JL,JK)))**(1./3.)
               ZF1=4.*API*ZRADL*PACDNC(JL,JK)*2.E5
     &              *(TMELT-3.-ZTP1TMP(JL))/ZRHO(JL,JK)
               ZF1=MAX(0.0,ZF1)
               ZFRL(JL)=ZFRL(JL)+ZTMST*1.4E-20*ZF1
               ZFRL(JL)=MAX(0.0,MIN(ZFRL(JL),ZXLB(JL)))
               ZXLB(JL)=ZXLB(JL)-ZFRL(JL)
               ZXIB(JL)=ZXIB(JL)+ZFRL(JL)
               ZFRL(JL)=ZFRL(JL)*ZCLCAUX(JL)
            ENDIF

         ENDDO

CBE  CHANGE THE LOCATION OF THE ADV FROM 2.9 TO 7.0
C
CJ-PP   2.91 ADD THE ADVECTED PRECIPITATION FROM PREVIOUS TIME STEP IF
C            ADDVECTION IS ON AND IF (WHEN PRECIPIATION PRESENT) THE
C            CURRENT OR NEXT LEVEL IS IN CLOUD (CLOUD COVER OVER ZEPSEC)
CKS         IF(LADVV .AND. JK .LT. KLEV) THEN
C           DO JL=KIDIA,KFDIA
C             IF((PADVE(JL,JK,3)+PADVE(JL,JK,4)) .NE. 0.0) THEN
C               IF(ZCLCP1(JL,JK) .GE. ZEPSEC .OR.
C     &              ZCLCP1(JL,JK+1) .GE. ZEPSEC) THEN
C                 ZRFL(JL)       = ZRFL(JL) + PADVE(JL,JK,3)
C                 PADVE(JL,JK,3) = 0.0
C                 ZSFL(JL)       = ZSFL(JL) + PADVE(JL,JK,4)
C                 PADVE(JL,JK,4) = 0.0
C               ENDIF
C             ENDIF
C           ENDDO
C         ENDIF
CJ-PP
C
C     7. CLOUD PHYSICS AND PRECIPITATION FLUXES AT THE SURFACE

         DO JL=KIDIA,KFDIA
            LOCC=(ZCLCAUX(JL) .GT. 0.0)
            ZCLCSTAR=MIN(ZCLCAUX(JL),ZCLCPRE(JL))
            ZAULOC=CAULOC*ZDZ(JL)/5000.
            ZAULOC=MAX(MIN(ZAULOC,CLMAX),CLMIN)
            JB=INVB(JL)
            LO=(JB .GE. JBMIN .AND. JB .LE. JBMAX
     &          .AND.PVERVEL(JL,JK).GT.0.)
            LO1=(JK .EQ. JB .OR. JK .EQ. JB+1)
            IF (LO .AND. LO1) THEN
               ZAULOC=0.0
            ENDIF
            ZQRHO=1.3/ZRHO(JL,JK)
            ZXLB(JL)=MAX(ZXLB(JL),1.E-20)
            ZXIB(JL)=MAX(ZXIB(JL),1.E-20)
            IF (ZCLCPRE(JL) .GT. 0.0) THEN
               ZXRP1=(ZRFL(JL)/(ZCLCPRE(JL)*12.45
     &              *SQRT(ZQRHO)))**(8./9.)
               ZXSP1=(ZSFL(JL)/(ZCLCPRE(JL)*CVTFALL))**(1./1.16)
            ELSE
               ZXRP1=0.0
               ZXSP1=0.0
            ENDIF



C       7.1 WARM CLOUDS: COALESCENCE PROCESS AFTER BEHENG(1994).
C         AUTOCONVERSION OF CLOUD DROPLETS AND COLLECTION OF CLOUD DROPLETS BY
C         FALLING RAIN. ACRETION OF CLOUD DROPLETS BY FALLING SNOW (ZSACL) IS
C         CALCULATED UNDER 7.2

            IF (LOCC .AND. ((ZXLB(JL) .GT. CQTMIN) .OR.
     &           (ZXIB(JL) .GT. CQTMIN))) THEN
               ZRAUT=CCRAUT*1.2E27/ZRHO(JL,JK)*(PACDNC(JL,JK)*1.E-6)
     &              **(-3.3)*(ZRHO(JL,JK)*1.E-3)**4.7
               ZEXM1=4.7-1.0
               ZEXP=-1./ZEXM1
               ZRAUT=ZXLB(JL)*(1.-(1.+ZRAUT*ZTMST*ZEXM1*ZXLB(JL)
     &              **ZEXM1)**ZEXP)
               ZXLB(JL)=ZXLB(JL)-ZRAUT
               ZRAC1=6.*ZXRP1*ZTMST
               ZRAC1=ZXLB(JL)*(1.-EXP(-ZRAC1))
               ZXLB(JL)=ZXLB(JL)-ZRAC1
               ZRAC2=6.*ZAULOC*ZRHO(JL,JK)*ZRAUT*ZTMST
               ZRAC2=ZXLB(JL)*(1.-EXP(-ZRAC2))
               ZXLB(JL)=ZXLB(JL)-ZRAC2
               ZRPR(JL)=ZRPR(JL)+ZCLCAUX(JL)*(ZRAUT+ZRAC2)+
     &              ZCLCSTAR*ZRAC1
C
C       7.2 COLD CLOUDS: CONVERSION OF CLOUD ICE TO SNOW AFTER
C          LEVKOV ET AL, 1992: AGGREGATION OF ICE CRYSTALS TO SNOW
C          AND ACCRETION OF ICE BY FALLING SNOW. ACCRETION OF CLOUD DROPLETS
C
               ZRIEFF=83.8*(ZXIB(JL)*ZRHO(JL,JK)*1000.)**0.216
               ZRIEFF=MIN(MAX(ZRIEFF,CEFFMIN),CEFFMAX)
               ZRIH=-2261.+SQRT(5113188.+2809.*ZRIEFF*ZRIEFF*ZRIEFF)
               ZRI=1.E-6*ZRIH**(1./3.)
               ZCOLLEFFI=EXP(0.025*(ZTP1TMP(JL)-TMELT))
               ZC1=17.5*ZRHO(JL,JK)/CRHOI*ZQRHO**0.33
               ZDT2=-6./ZC1*LOG10(ZRI*1.E4)
               ZSAUT=CCSAUT/ZDT2
               ZSAUT=ZXIB(JL)*(1.-1./(1.+ZSAUT*ZTMST*ZXIB(JL)))
               ZXIB(JL)=ZXIB(JL)-ZSAUT
               ZSACI1=0.0
               ZSACI2=0.0
               ZSACL1=0.0
               ZSACL2=0.0
               IF (ZXSP1 .GT. CQTMIN) THEN
                  ZLAMSM=(ZXSP1/(API*CRHOSNO*CN0S))**0.8125
                  ZSACI1=API*CN0S*3.078*ZLAMSM*ZQRHO**0.5
                  ZSACL1=ZXLB(JL)*(1.-EXP(-ZSACI1*CCSACL*ZTMST))
                  ZXLB(JL)=ZXLB(JL)-ZSACL1
                  ZSACL1=ZCLCSTAR*ZSACL1
                  ZSACI1=ZSACI1*ZCOLLEFFI*ZTMST
                  ZSACI1=ZXIB(JL)*(1.-EXP(-ZSACI1))
                  ZXIB(JL)=ZXIB(JL)-ZSACI1
               ENDIF
               ZXSP2=ZAULOC*ZRHO(JL,JK)*ZSAUT
               IF (ZXSP2 .GT. CQTMIN) THEN
                  ZLAMSM=(ZXSP2/(API*CRHOSNO*CN0S))**0.8125
                  ZSACI2=API*CN0S*3.078*ZLAMSM*ZQRHO**0.5
                  ZSACL2=ZXLB(JL)*(1.-EXP(-ZSACI2*CCSACL*ZTMST))
                  ZXLB(JL)=ZXLB(JL)-ZSACL2
                  ZSACL2=ZCLCAUX(JL)*ZSACL2
                  ZSACI2=ZSACI2*ZCOLLEFFI*ZTMST
                  ZSACI2=ZXIB(JL)*(1.-EXP(-ZSACI2))
                  ZXIB(JL)=ZXIB(JL)-ZSACI2
               ENDIF
               ZSACL(JL)=ZSACL1+ZSACL2
               ZSPR(JL)=ZSPR(JL)+ZCLCAUX(JL)*(ZSAUT+ZSACI2)
     &              +ZCLCSTAR*ZSACI1

            ENDIF
C
C       7.3 UPDATING PRECIPITATION FLUXES. IN THE LOWEST LAYER (KLEV)
C        THE SEDIMENTATION SINK OF CLOUD ICE IS BALANCED
C        BY PRECIPITATION AT THE SURFACE (THROUGH ZZDRS).
C        FRACTION OF PRECIPITATION CLOUDS (ZCLCPRE) USED FOR THE
C        CALCULATION OF EVAPORATION/SUBLIMATION OF RAIN/SNOW IN THE NEXT LAYER
C
            ZZDRR=ZCONS2*ZDP(JL)*ZRPR(JL)
            ZZDRS=ZCONS2*ZDP(JL)*(ZSPR(JL)+ZSACL(JL))
            IF (JK .EQ. KLEV) THEN
               ZZDRS=ZZDRS+ZXIFLUX(JL)
               ZCONS=ZCONS2*ZDP(JL)/(ZLSDCP(JL)-ZLVDCP(JL))
               ZSNMLT=MIN(ZXSEC*ZZDRS,ZCONS*MAX(0.,(ZTP1TMP(JL)-TMELT)))
               ZZDRR=ZZDRR+ZSNMLT
               ZZDRS=ZZDRS-ZSNMLT
               ZSMLT(JL)=ZSMLT(JL)+ZSNMLT/(ZCONS2*ZDP(JL))
            ENDIF
            ZPRETOT=ZRFL(JL)+ZSFL(JL)
            ZPREDEL=ZZDRR+ZZDRS
            LO=(ZPRETOT .GT. ZPREDEL)
            ZCLCPRE(JL)=MERGE(ZCLCPRE(JL),ZCLCAUX(JL),LO)
            ZPRESUM=ZPRETOT+ZPREDEL
            IF (ZPRESUM .GT. CQTMIN) THEN
               ZCLCPRE(JL)=MAX(ZCLCPRE(JL),(ZCLCAUX(JL)*ZPREDEL
     &              +ZCLCPRE(JL)*ZPRETOT)/ZPRESUM)
               ZCLCPRE(JL)=MIN(ZCLCPRE(JL),1.0)
               ZCLCPRE(JL)=MAX(ZCLCPRE(JL),0.0)
            ELSE
               ZCLCPRE(JL)=0.0
            ENDIF

            ZSFL(JL)=ZSFL(JL)+ZZDRS-ZCONS2*ZDP(JL)*ZSUB(JL)
            ZRFL(JL)=ZRFL(JL)+ZZDRR-ZCONS2*ZDP(JL)*ZEVP(JL)

         ENDDO
C
C    8.   UPDATING TENDENCIES OF T, Q, XL, XI

         DO JL=KIDIA,KFDIA

C    8.1. TENDENCIES OF THERMODYNAMIC VARIABLES
C             ATTN: THE TERMS ZXISUB AND ZXIMLT DO NOT APPEAR IN
C                   PXITE BECAUSE THESE PROCESSES HAVE ALREADY BEEN
C                   INCLUDED IN PXITE VIA CHANGES IN CLOUD ICE
C                   SEDIMENTATION (SEE 3.1, 3.2 AND 4)
C

            PQTE(JL,JK)=PQTE(JL,JK)
     &           +(-ZCND(JL)+ZEVP(JL)+ZXLEVAP(JL)
     &             -ZDEP(JL)+ZSUB(JL)+ZXIEVAP(JL)
     &           + ZXISUB(JL))       /ZTMST
            PTTE(JL,JK)=PTTE(JL,JK)+(ZLVDCP(JL)
     &          *(ZCND(JL)-ZEVP(JL)-ZXLEVAP(JL))
     &                              +ZLSDCP(JL)
     &          *(ZDEP(JL)-ZSUB(JL)-ZXIEVAP(JL)-ZXISUB(JL))
     &                             +(ZLSDCP(JL)-ZLVDCP(JL))
     &           *(-ZSMLT(JL)-ZIMLT(JL)-ZXIMLT(JL)+ZFRL(JL)
     &           +ZSACL(JL)))       /ZTMST
            PXTE(JL,JK)=PXTE(JL,JK)+ZXLTE(JL)
     &           +(ZIMLT(JL)+ZXIMLT(JL)-ZFRL(JL)-ZRPR(JL)
     &           -ZSACL(JL)+ZCND(JL)-ZXLEVAP(JL))/ZTMST
            PXITE(JL,JK)=PXITE(JL,JK)+ZXITE(JL)
     &           +(ZFRL(JL)-ZSPR(JL)
     &           +ZDEP(JL)-ZXIEVAP(JL))/ZTMST
            ZTP1(JL,JK)=PTM1(JL,JK)+PTTE(JL,JK)*ZTMST
            ZQP1(JL,JK)=PQM1(JL,JK)+PQTE(JL,JK)*ZTMST
            ZXLP1=PXM1(JL,JK)+PXTE(JL,JK)*ZTMST
            ZXIP1=PXIM1(JL,JK)+PXITE(JL,JK)*ZTMST

C  ***** WIEDER ORIGINAL ECHAM4/REMO ZUR BERECHNUNG DER WOLKENBEDECKUNG:

            IT=INT(ZTP1(JL,JK)*1000.)
            ZQSP1(JL,JK)=TLUCUA(IT)/PAPP1(JL,JK)
            ZQSP1(JL,JK)=MIN(ZQSP1(JL,JK),0.5)
            ZCOR=1./(1.-VTMPC1*ZQSP1(JL,JK))
            ZQSP1(JL,JK)=ZQSP1(JL,JK)*ZCOR

            ZQR=ZQSP1(JL,JK)*ZSAT(JL,JK)*ZRHC(JL,JK)
            ZCLP1=(ZQP1(JL,JK)-ZQR)/(ZQSP1(JL,JK)*ZSAT(JL,JK)-ZQR)
            ZCLP1=MAX(ZCLP1,0.)
            ZCLP1=MIN(ZCLP1,1.)
            ZCLP1=1.-SQRT(1.-ZCLP1)
            IF (ZXIP1 <=CCWMIN .AND. ZXLP1 <=CCWMIN) THEN
               PACLC(JL,JK)=0.
            ELSE
               PACLC(JL,JK)=ZCLP1
            ENDIF

C***************  BIS HIER



CSP     8.2 CORRECTIONS: AVOID NEGATIVE CLOUD WATER/ICE
CSP         DIE ABFRAGE ERFOLGT JETZT EINZELN FUER WASSER UND EIS,
CSP         ANSONSTEN WIE REMO_ORIGINAL

            ZXLOLD=ZXLP1
            LO=(ZXLP1 .LT. ZEPZWP)
            IF(LO) THEN
               ZXLP1=0.0
            ENDIF
            ZDXLCOR=(ZXLP1-ZXLOLD)/ZTMST
            ZXIOLD=ZXIP1
            LO1=(ZXIP1 .LT. ZEPZWP)
            IF(LO1) THEN
               ZXIP1=0.0
            ENDIF
            ZDXICOR=(ZXIP1-ZXIOLD)/ZTMST
            IF(LO .AND. LO1) THEN
               PACLC(JL,JK)=0.0
            ELSE
               PACLC(JL,JK)=PACLC(JL,JK)
            ENDIF
            PACLCAC(JL,JK)=PACLCAC(JL,JK)+PACLC(JL,JK)*ZDIAGT
            PXTE(JL,JK)=PXTE(JL,JK)+ZDXLCOR
            PXITE(JL,JK)=PXITE(JL,JK)+ZDXICOR
            PQTE(JL,JK)=PQTE(JL,JK)-ZDXLCOR-ZDXICOR
            PTTE(JL,JK)=PTTE(JL,JK)+ZLVDCP(JL)*ZDXLCOR
     &           +ZLSDCP(JL)*ZDXICOR
            RPRAC(JL,JK)=RPRAC(JL,JK)+ZRPR(JL)*ZDIAGT
         ENDDO
      ENDDO

CSP    9. DIAGNOSTICS ORIGINAL REMO

C
C*         8.5     SURFACE FLUXES.
C
      DO JL=KIDIA,KFDIA
         PRSFL(JL)=PRSFL(JL)+ZRFL(JL)
         PSSFL(JL)=ZSFL(JL)
         PAPRL(JL)=PAPRL(JL)+ZDIAGW*(PRSFL(JL)+PSSFL(JL))
         PAPRS(JL)=PAPRS(JL)+ZDIAGW*PSSFL(JL)
      ENDDO
C
C*          8.6    ACCUMULATED TOTAL CLOUDCOVER
C
      DO JL=KIDIA,KFDIA
         ZACLCOV(JL)=1.-PACLC(JL,1)
      ENDDO
      DO JK=2,KLEV
         DO JL=KIDIA,KFDIA
            ZACLCOV(JL)=
     &           ZACLCOV(JL)*(1.-AMAX1(PACLC(JL,JK),PACLC(JL,JK-1)))
     &           /(1.-AMIN1(PACLC(JL,JK-1),1.-ZEPSEC))
         ENDDO
      ENDDO
      DO JL=KIDIA,KFDIA
         ZACLCOV(JL)=1.-ZACLCOV(JL)
         PACLCOV(JL)=PACLCOV(JL)+ZDIAGT*ZACLCOV(JL)
      ENDDO

C
C*          8.7     VERTICALLY INTEGRATED HUMIDITY AND CLOUD WATER
C
      DO JL=KIDIA,KFDIA
         ZQVI(JL)=0.
         ZQLVI(JL)=0.
         ZQIVI(JL)=0.
      ENDDO
C
      DO JK=KTDIA,KLEV
         DO JL=KIDIA,KFDIA
            ZDPG=(PAPHM1(JL,JK+1)-PAPHM1(JL,JK))/G
            ZQVI(JL)=ZQVI(JL)+PQM1(JL,JK)*ZDPG
            ZQLVI(JL)=ZQLVI(JL)+PXM1(JL,JK)*ZDPG
            ZQIVI(JL)=ZQIVI(JL)+PXIM1(JL,JK)*ZDPG
         ENDDO
      ENDDO
C
      DO JL=KIDIA,KFDIA
         PQVI(JL)=PQVI(JL)+ZDIAGT*ZQVI(JL)
         PALWCVI(JL)=PALWCVI(JL)+ZDIAGT*ZQLVI(JL)
         PQIVI(JL)=PQIVI(JL)+ZDIAGT*ZQIVI(JL)
      ENDDO
C
      RETURN
      END SUBROUTINE LWCOND
