      SUBROUTINE PHYORG
         !
         ! WARNING!! Here the actual parameter names differ from
         ! the names of the dummy arguments when called in:
         ! phyec
         !
     &  (NSTEP  , NHOR   , NLEV   , NLEVP1 , CONACC , VCT    ,
     &   NRADIA , CEVAPCU, CVDAES , CVDAEL , CVDAEU ,
     &   CVDAED , DTIME  , TWODT  , EPS    , NCBASE ,
     &   NTBASE , NTIMST , NMONTH , LVDIFF , LSURF  ,
     &   LGWDRAG, LCONV  , LCOND  , NRADFR , LRAD   ,
     &   LSOLC  , LAER   , LCFC   , LSICED , IEXC   ,
     &   NOZ    , LGADSRH, LAEROZ , L5LAY  , LWDIF  ,
C        -- PROGNOSTIC FIELDS AND TENDENCIES --
     &   UM1    , VM1    , TM1    , QM1    , QE     , XM1    ,
     &   TRAD   , TGWDR  , UGWDR  , VGWDR  , TCOND  , QCOND  ,
     &   XCOND  , TCONV  , UCONV  , VCONV  , QCONV  ,
C        -- PRESSURE FIELDS GALCIERMASK AND COORDINATES --
     &   ZPS    , SINLAT , COSLAT , SINLON , COSLON ,
C        -- SKIN-TEMPERATURE RELATED FIELDS --
     &   SICED  , TEFF   ,
C        -- SURF TEMP AND MOIST FIELDS --
     &   TS     , TSM    , TSM1M  , TSL    , TSLM   , TSLM1M ,
     &   TSW    , TSWM   , TSWM1M , TSI    , TSIM   , TSIM1M ,
     &   WS     , WSM    , WSM1M  , SN     , SNM    , SNM1M  ,
     &   WL     , WLM    , WLM1M  , TD     , TDM    , TDM1M  ,
C        -- SURFACE FLUXES AND FLUX DERIVATIVES --
     &   SRFL   , THFL   , QHFL   , XHFL   , RSFC   , SSFC   ,
     &   RSFL   , SSFL   , AHFS   , AHFSL  , AHFSW  , AHFSI  ,
     &   AHFL   , AHFICE , ZQDB   ,
     &   ZQDBL  , ZQDBW  , ZQDBI  , DHFT   , DHFQW  , DHFQS  ,
C        -- FIELDS USED IN DEEP CONVECTION --
     &   TOPMAX , VERVEL , AZ0    , AZ0L   , AZ0W   , AZ0I   ,
C        -- RAIN AND SNOW FALL AND EVAPORATION --
     &   APRC   , APRL   , APRS   , EVAP   , EVAPM  , EVAPL  ,
     &   EVAPW  , EVAPI  , ACLC   , ACLCAC , VGRAT  , FOREST ,
     &   ALB    , ALBEDO , ALSOL  , ALSOW  , ALSOI  ,
C        -- REMAINING ELEMENTS IN *VDIFF* --
     &   DEW2   , TKE    , TKEM   , TKEM1M , WSMX   , VLT    ,
     &   RGCGN  , TEMP2  , T2MAX  , T2MIN  , USTAR3 ,
     &   USTR   , U10    , VDIS   , VSTR   , V10    , WIND10 ,
     &   WIMAX  , USTRL  , USTRW  , USTRI  , VSTRL  , VSTRW  ,
     &   VSTRI  , INFRL  , INFRW  , INFRI  , Z2ZM   ,
     &   Z2ZH   , Z2ZHL  , Z2ZHW  , Z2ZHI  ,
C        -- REMAINING ELEMENTS --
     &   ACLCOV , ALWCVI , QVI    , EMTER  , TRSOL  , EMTEF  ,
     &   TRSOF  , SCLF0  , SCLFS  , SRAF0  , SRAFS  , TCLF0  ,
     &   TCLFS  , TRAF0  , TRAFS  , ACLCV  , SRADS  , SRADSU ,
     &   SRAD0  , SRAD0U , TRADS  , TRADSU , TRAD0  , DSNAC  ,
     &   RUNOFF , SNMEL  , TDCL   , TDCLM  , TDCLM1M, TD3    ,
     &   TD3M   , TD3M1M , TD4    , TD4M   , TD4M1M , TD5    ,
     &   TD5M   , TD5M1M , QRES   , TSLIN  , TSMAX  , TSMIN  ,
     &   TSN    , TSNM   , TSNM1M , TSURF  , VAROR  , DRAIN  ,
C        -- REMAINING ELEMENTS IN *GWDRAG* --
     &   USTRGW , VAR    , VDISGW , VSTRGW ,
C        -- SOIL PARAMETER --
     &   TLAMBDA, DLAMBDA, PORVOL , FCAP   , WI3    , WI3M   ,
     &   WI3M1M , WI4    , WI4M   , WI4M1M , WI5    , WI5M   ,
     &   WI5M1M , WI     , WIM    , WIM1M  , WICL   , WICLM  ,
     &   WICLM1M, LOGLAC , LOLAND , LOSEA  , LOICE  , LALAND ,
     &   GHPBL  , BETA   , WMINLOK, WMAXLOK, CAPE   , OZONPL ,
     &   SO4ALL , SO4NAT ,
C        -- 5 LAYER SOIL SCHEME: MOISTURE LAYERS AND PARAMETERS --
     &   WS1    , WS2    , WS3    , WS4    , WS5   , DZR     ,
     &   DZS    , FKSAT  , FMPOT  , BCLAPP , VPOR  ,
C        -- EVAPORATION FLUXES --
     &   ETRANS , EBSOIL , ESNOW  , ESKIN  , ERES,
C        -- CLOUD ICE
     &   XIM1   , XICOND , XIVI   , RPRAC)
C
C     AUTHOR: R. PODZUN 07/02/2007
C     --------
C
C     ******* NOV. 2005/JAN. 2007 -- STEFAN HAGEMANN
C     ***           CHANGES FOR 5 LAYER SOIL SCHEME
C     ***     NDEEP   = NUMBER OF SOIL LAYERS (CURRENTLY HARD CODED 5)
C     ***      WSI(I) = NEW SOIL MOISTURE CONTENT IN I. SOIL LAYER [MM]
C     ***   WSIM1M(I) = OLD SOIL MOISTURE CONTENT IN I. SOIL LAYER [MM]
C
C     *** SEPT. 2007 -- STEFAN HAGEMANN
C     ***
C     *** TRANSPIRATION AND BARE SOIL EVAPORATION MAY BE REDUCED BY CALCUCLATIONS
C     ***   IN ROUTINE SURF->SOILCHANGE BY AMOUNT REDEVAP. CONSEQUENTLY EVAP AND EVAPL
C     ***   MUST BE CHANGED IN PHYORG, TOO.
C     *** BARE SOIL EVAPORATION IN VDIFF IS ONLY TAKEN FROM THE MOST UPPER SOIL LAYER,
C     ***   THUS THE FIELD CAPACITY OF THIS LAYER IS CALCULATED AND SUBMITTED TO VDIFF
C     ***   THEREFORE THE SOIL LAYER THICKNESS DEFINITION ZDEL HAS MOVED FROM SURF
C     ***   TO PHYORG.
C
C     ******* JANUARY 2008 - S. HAGEMANN
C     *** BARE SOIL EVAPORATION CAN ONLY BE TAKEN FROM THE MOST UPPER LAYER
C         IF 5 LAYER SOIL SCHEME IS USED. THUS, IF THE BUCKET IS USED, THE
C         THE OLD FORMULATION HAS TO BE USED. THEREFORE THE I5LAYER SWITCH
C         HAS TO BE INCLUDED IN VDIFF.
C
C     ******* MARCH 2008 - S. HAGEMANN
C     *** REDUCTION OF EVAPORATION IS TAKEN OUT AGAIN AS THERE WAS NO FEEDBACK
C         TO THE ATMOSPHERE. AN IMPLEMENTATION OF THIS FEEDBACK WOULD REQUIRE
C         MORE COMPLEX CHANGES TO VDIFF AND COMPREHENSIVE TESTING. IN
C         MY OPINIOPN A REDUCTION OF TRANSPIRATION IN VDIFF WOULD BE DESIRABLE.
C
C     *** UM TRANSPIRATION IN VDIFF REDUZIEREN ZU KOENNEN, MUSS 5 SOIL LAYER INFO
C         SCHON IN PHYORG VORHANDEN SEIN. DAHER WIRD CALL DER ROUTINE SOILDEF
C         NUN NACH PHYORG VERLEGT, UND SOMIT DIE FELDER NACH SURF UEBERGEBEN.

C
C     PURPOSE.
C     --------
C
C            THIS SUBROUTINE CONTROLS THE CALLS TO THE VARIOUS
C     PHYSICAL SUBROUTINES.
C
C**   INTERFACE.
C     ----------
C
C            *PHYORG* IS CALLED FROM *PHYEC*.
C
C     RESULTS.
C     --------
C
C     EXTERNALS.
C     ----------
C
C            *COPYR*    COPY REAL VALUES FROM A REAL ARRAY INTO ANOTHER.
C            *SETRA*   RESET REAL ARRAY TO GIVEN VALUES.
C            *SIGMA*   COMPUTES SCALAR PRODUCT.
C            *GEOPEC*   COMPUTES FULL LEVEL GEOPOTENTIALS.
C            *PRES*     COMPUTES HALF LEVEL PRESSURES.
C            *PRESF*    COMPUTES FULL LEVEL PRESSURES.
C            *VDIFF*    VERTICAL EXCHANGE OF U,V,T,Q BY TURBULENCE.
C            *COND*     LARGE SCALE WATER PHASE CHANGES.
C            *RADMOD*  CONTROLS RADIATION COMPUTATIONS.
C            *SURF*    COMPUTE NEW SURFACE VALUES.
C            *GWDRAG *    GRAVITY WAVE DRAG SCHEME
C
C     REFERENCE.
C     ----------
C
      IMPLICIT NONE
C
      INCLUDE "COMCON"
      INCLUDE "faktinf.h"
C
      INTEGER, INTENT(IN) :: NSTEP  , NHOR   , NLEV   , NLEVP1
      INTEGER, INTENT(IN) :: NCBASE, NMONTH, NOZ, NRADFR,
     &                       NTBASE, NTIMST, NRADIA, IEXC
      REAL,    INTENT(IN) :: CONACC, DTIME, EPS, TWODT
C
C     -------------------------------------------------------
C     DECLARATION OF REAL ARRAYS
C     TODO: These arrays are duplicated (memcpy, "pass by copy").
C           This should be avoided sinc memcpy is quite costly.
C
      REAL, INTENT(IN) ::
     &     VCT    (2*NLEVP1),
     &     CEVAPCU(NLEV)    , CVDAES(NLEVP1)   ,
     &     CVDAEL (NLEVP1)  , CVDAEU(NLEVP1)   ,
     &     CVDAED (NLEVP1)  , OZONPL(NHOR,NOZ) ,
     &     SO4ALL(NHOR,NLEV), SO4NAT(NHOR,NLEV)
      REAL, INTENT(INOUT) ::
     &     UM1   (NHOR,NLEV)  , VM1   (NHOR,NLEV)  ,
     &     TM1   (NHOR,NLEV)  , QM1   (NHOR,NLEV)  ,
     &     QE    (NHOR,NLEV)  , XM1   (NHOR,NLEV)  ,
     &     TRAD  (NHOR,NLEV)  , TGWDR (NHOR,NLEV)  ,
     &     UGWDR (NHOR,NLEV)  , VGWDR (NHOR,NLEV)  ,
     &     TCOND (NHOR,NLEV)  , QCOND (NHOR,NLEV)  ,
     &     XCOND (NHOR,NLEV)  , TCONV (NHOR,NLEV)  ,
     &     UCONV (NHOR,NLEV)  , QCONV (NHOR,NLEV)  ,
     &     VCONV (NHOR,NLEV)  , VAR   (NHOR,4)     ,
     &     ACLC  (NHOR,NLEV)  , ACLCAC(NHOR,NLEV)  ,
     &     TKE   (NHOR,NLEV)  , TKEM  (NHOR,NLEV)  ,
     &     TKEM1M(NHOR,NLEV)  , VERVEL(NHOR,NLEV)  ,
     &     Z2ZM  (NHOR,NLEV)  , Z2ZH  (NHOR,NLEV)  ,
     &     EMTER (NHOR,NLEVP1), TRSOL (NHOR,NLEVP1),
     &     EMTEF (NHOR,2)     , TRSOF (NHOR,2)     ,
     &     SINLAT(NHOR)       , COSLAT(NHOR)       ,
     &     SINLON(NHOR)       , COSLON(NHOR)
      REAL :: ALNPR (NHOR,NLEV), ALPHA (NHOR,NLEV)
      REAL, INTENT(INOUT) ::
     &     SICED (NHOR), TEFF  (NHOR),
     &     TS    (NHOR), TSM   (NHOR), TSM1M  (NHOR),
     &     TSL   (NHOR), TSLM  (NHOR), TSLM1M (NHOR),
     &     TSW   (NHOR), TSWM  (NHOR), TSWM1M (NHOR),
     &     TSI   (NHOR), TSIM  (NHOR), TSIM1M (NHOR),
     &     WS    (NHOR), WSM   (NHOR), WSM1M  (NHOR),
     &     SN    (NHOR), SNM   (NHOR), SNM1M  (NHOR),
     &     WL    (NHOR), WLM   (NHOR), WLM1M  (NHOR),
     &     TD    (NHOR), TDM   (NHOR), TDM1M  (NHOR),
     &     SRFL  (NHOR),
     &     THFL  (NHOR), QHFL  (NHOR), XHFL   (NHOR),
     &     RSFC  (NHOR), SSFC  (NHOR), RSFL   (NHOR),
     &     SSFL  (NHOR), AHFL  (NHOR), AHFS   (NHOR),
     &     ZQDB  (NHOR), ZQDBL (NHOR), ZQDBW  (NHOR),
     &     ZQDBI (NHOR), AHFSL (NHOR), AHFSW  (NHOR),
     &     AHFSI (NHOR), AHFICE(NHOR), DHFT   (NHOR),
     &     DHFQW (NHOR), DHFQS (NHOR), TOPMAX (NHOR),
     &     AZ0   (NHOR), AZ0L  (NHOR), AZ0W   (NHOR),
     &     AZ0I  (NHOR), APRC  (NHOR), APRL   (NHOR),
     &     APRS  (NHOR), EVAP  (NHOR), EVAPM  (NHOR),
     &     EVAPL (NHOR), EVAPW (NHOR), EVAPI  (NHOR),
     &     VGRAT (NHOR), FOREST(NHOR), ALB    (NHOR),
     &     ALBEDO(NHOR), ALSOL (NHOR), ALSOW  (NHOR),
     &     ALSOI (NHOR), DEW2  (NHOR), WSMX   (NHOR),
     &     VLT   (NHOR), RGCGN (NHOR), CAPE   (NHOR),
     &     TEMP2 (NHOR), T2MAX (NHOR), T2MIN  (NHOR),
     &     USTAR3(NHOR), USTR  (NHOR), U10    (NHOR),
     &     VDIS  (NHOR), VSTR  (NHOR), V10    (NHOR),
     &     WIND10(NHOR), WIMAX (NHOR), USTRL  (NHOR),
     &     USTRW (NHOR), USTRI (NHOR), VSTRL  (NHOR),
     &     VSTRW (NHOR), VSTRI (NHOR), Z2ZHL  (NHOR),
     &     Z2ZHW (NHOR), Z2ZHI (NHOR), ACLCOV (NHOR),
     &     ALWCVI(NHOR), QVI   (NHOR), SCLF0  (NHOR),
     &     SCLFS (NHOR), SRAF0 (NHOR), SRAFS  (NHOR),
     &     TCLF0 (NHOR), TCLFS (NHOR), TRAF0  (NHOR),
     &     TRAFS (NHOR), ACLCV (NHOR), SRADS  (NHOR),
     &     SRADSU(NHOR), SRAD0 (NHOR), SRAD0U (NHOR),
     &     TRADS (NHOR), TRADSU(NHOR), TRAD0  (NHOR),
     &     DSNAC (NHOR), RUNOFF(NHOR), SNMEL  (NHOR),
     &     TDCL  (NHOR), TDCLM (NHOR), TDCLM1M(NHOR),
     &     TD3   (NHOR), TD3M  (NHOR), TD3M1M (NHOR),
     &     TD4   (NHOR), TD4M  (NHOR), TD4M1M (NHOR),
     &     TD5   (NHOR), TD5M  (NHOR), TD5M1M (NHOR),
     &     QRES  (NHOR), TSLIN (NHOR), TSMAX  (NHOR),
     &     TSMIN (NHOR), TSN   (NHOR), TSNM   (NHOR),
     &     TSNM1M(NHOR), TSURF (NHOR), VAROR  (NHOR),
     &     DRAIN (NHOR), USTRGW(NHOR), VDISGW (NHOR),
     &     VSTRGW(NHOR), ZPS   (NHOR), TLAMBDA(NHOR),
     &     DLAMBDA(NHOR),PORVOL(NHOR), FCAP   (NHOR),
     &     WI3   (NHOR), WI3M  (NHOR), WI3M1M (NHOR),
     &     WI4   (NHOR), WI4M  (NHOR), WI4M1M (NHOR),
     &     WI5   (NHOR), WI5M  (NHOR), WI5M1M (NHOR),
     &     WI    (NHOR), WIM   (NHOR), WIM1M  (NHOR),
     &     WICL  (NHOR), WICLM (NHOR), WICLM1M(NHOR),
     &     GHPBL (NHOR), BETA  (NHOR), WMINLOK(NHOR),
     &     WMAXLOK(NHOR)
      REAL :: Z2ZML (NHOR), Z2ZMW (NHOR), Z2ZMI (NHOR)
      REAL, INTENT(INOUT) ::
     &     WS1   (NHOR), WS2   (NHOR), WS3    (NHOR),
     &     WS4   (NHOR), WS5   (NHOR), DZR    (NHOR),
     &     DZS   (NHOR), FKSAT (NHOR), FMPOT  (NHOR),
     &     BCLAPP(NHOR), VPOR  (NHOR), ETRANS (NHOR),
     &     EBSOIL(NHOR), ESNOW (NHOR), ESKIN  (NHOR),
     &     ERES  (NHOR)
C
      REAL, INTENT(INOUT) :: RPRAC(NHOR,NLEV)
C
C     -------------------------------------------------------
C     DECLARATION OF LOGICAL ARRAYS
C
C     -- LAND AND ICE MASK --
      LOGICAL, INTENT(IN) :: LOLAND(NHOR), LOGLAC(NHOR),
     &        LOSEA (NHOR), LOICE (NHOR), LALAND(NHOR),
     &        LVDIFF, LGWDRAG, LSURF , LCONV, LCOND,
     &        LRAD  , LSOLC , LAER , LCFC ,
     &        LSICED, LGADSRH, LAEROZ, L5LAY, LWDIF
C
C     -------------------------------------------------------
C     DECLARATION OF INTEGER ARRAYS
C
      INTEGER, INTENT(IN) :: INFRL(NHOR), INFRW(NHOR), INFRI(NHOR)
C
C     -------------------------------------------------------
C     DECLARATION OF INTEGER WORKING ARRAYS
C
      INTEGER :: ITYPE(NHOR), ILAB(NHOR,NLEV)
C
C     -------------------------------------------------------
C     DECLARATION OF REAL WORKING ARRAYS
C
      REAL ::
     &     ZGEO(NHOR)     , CVS  (NHOR)     , CVW  (NHOR)     ,
     &     WLMX(NHOR)     , ZTVM1(NHOR,NLEV), ZXTEC(NHOR,NLEV),
     &     ZUE (NHOR,NLEV), ZVE  (NHOR,NLEV), ZTE  (NHOR,NLEV),
     &     ZXE (NHOR,NLEV), ZTM1 (NHOR,NLEV), ZQM1 (NHOR,NLEV),
     &     ZXM1(NHOR,NLEV), ZQHFLA(NHOR)    , ZDHFTI(NHOR)    ,
     &     ZTHFLI(NHOR)   , APM1 (NHOR,NLEV), APHM1 (NHOR,NLEVP1),
     &     GEOM1(NHOR,NLEV)

      REAL, INTENT (IN) :: XIM1(NHOR,NLEV)
      REAL, INTENT (INOUT) :: XICOND (NHOR,NLEV), XIVI(NHOR)
      REAL :: PACDNC(NHOR,NLEV), ZXIE(NHOR,NLEV), ZXIM1(NHOR,NLEV)
CSH
C
      INTEGER, PARAMETER :: NDEEP=5
C
      REAL ::
     &     WSI(NHOR,NDEEP), WSIM1M(NHOR,NDEEP), ZDEL(NDEEP)
C
C     *** ACTUAL EVAPORATION ARRAYS FOR THE TIME STEP
      REAL ::
     &     AETRANS(NHOR), AEBSOIL(NHOR), AESNOW(NHOR),
     &     AESKIN (NHOR), AERES (NHOR)
C     *** 5 LAYER VALUES FOR EACH LAYER: ROOTED DEPTH, WATERED DEPTH
C     *** SATURATED WATER CONTENT AND FIELD CAPACITY
      REAL ::
     &     DZRSI(NHOR, NDEEP),  DZSI(NHOR, NDEEP),
     &     ZWSAT(NHOR, NDEEP), ZWSFC(NHOR, NDEEP)
CSH
      REAL :: RALPHA, RLNPR, ZN1, ZN2, ZPRAT, ZSDISS,
     &        ZSEVAP, ZSHEAT, ZSMELT, 
     &        ZSRAIN, ZDIAGS
C
      INTEGER :: I5LAYER, IK, JHOR, JK, JL, JLEV, KSTART, KSTOP, KTDIA, 
     &           NEXP, NLEVM1, 
     &           NSTART
C
C     END DECLARATIONS
C     -------------------------------------------------------
C
C
C*    1.   RESET WORK SPACE TO ZERO AND SET PHYSICAL
C          ----- ---- ----- -- ---- --- --- --------
C          TENDENCIES TO ZERO
C          ---------- -- ----
C
      KSTART = 1
      KSTOP  = NHOR
      NSTART = 0
      NLEVM1 = NLEV-1
      KTDIA  = 1
CSH
C     SWITCH FOR 5 LAYER CALCULATION OFF=0/ON=1
      IF (L5LAY) THEN
        I5LAYER=1
      ELSE
        I5LAYER=0
      ENDIF
CSH
      DO JL=1,NHOR
        ITYPE(JL)=0
      END DO
C
      DO JK=1,NLEV
        DO JL=1,NHOR
          ILAB(JL,JK)=0
        END DO
      END DO
C
      DO JL=1,NHOR
        ZGEO   (JL)=0.
        CVS    (JL)=0.
        CVW    (JL)=0.
        WLMX   (JL)=0.
        Z2ZML  (JL)=0.
        Z2ZMW  (JL)=0.
        Z2ZMI  (JL)=0.
        ZDHFTI (JL)=0.
        ZQHFLA (JL)=0.
        ZTHFLI (JL)=0.
CSH
        AEBSOIL(JL)=0.
        AETRANS(JL)=0.
        AESNOW (JL)=0.
        AESKIN (JL)=0.
        AERES  (JL)=0.
CSH
      END DO
C
      DO JK=1,NLEV
        DO JL=1,NHOR
          ZTVM1 (JL,JK)=0.
          ZXTEC (JL,JK)=0.
          TRAD  (JL,JK)=0.
          TGWDR (JL,JK)=0.
          UGWDR (JL,JK)=0.
          VGWDR (JL,JK)=0.
          ZUE   (JL,JK)=0.
          ZVE   (JL,JK)=0.
          ZTE   (JL,JK)=0.
          ZXE   (JL,JK)=0.
          TCOND (JL,JK)=0.
          QCOND (JL,JK)=0.
          XCOND (JL,JK)=0.
          TCONV (JL,JK)=0.
          UCONV (JL,JK)=0.
          VCONV (JL,JK)=0.
          QCONV (JL,JK)=0.
          ZXIE  (JL,JK)=0.
          XICOND(JL,JK)=0.
        END DO
      END DO
C
C
C*    2.    UPDATE SOME FIELDS NEEDED BY THE PHYSICS ROUITNES
C           ------ ---- ------ ------ -- --- ------- ---------
C
C     2.1   GET SEAICE COVER AND DEPTH
C           --- ------ ----- --- -----
C
      IF (.NOT.LSICED) THEN
        CALL ATMICE (NHOR, KSTART, KSTOP, SICED, INFRI)
      ENDIF
C
C
C*    3.    UPDATE SOME AUXILIARY FIELDS
C           ------ ---- --------- ------
C
C     3.1   UPDATE ALPHA AT TIME STEP T - DT.
C           ------ ----- -- ---- ---- ------
C
      CALL PRES(NHOR,KSTART,KSTOP,NLEVP1,VCT,APHM1,ZPS)
C
CDJ 30.5.94
      CALL PRESF(NHOR,KSTART,KSTOP,NLEV,APM1,APHM1)
C
      RALPHA=RD*LOG(2.)
      RLNPR=2.*RALPHA
C
C     3.2   SET TOP-LEVEL VALUES OF ALPHA AND ALNPR
C           --- --------- ------ -- ----- --- -----
C
      DO JL=1,NHOR
        ALPHA(JL,1)=RALPHA
        ALNPR(JL,1)=RLNPR
      END DO
C
C     3.3   SET REMAINING VALUES OF ALPHA AND ALNPR
C           --- --------- ------ -- ----- --- -----
C
      DO JK=2,NLEV
        DO JL=KSTART,KSTOP
          ALNPR(JL,JK)=RD*ALOG(APHM1(JL,JK+1)/APHM1(JL,JK))
          ALPHA(JL,JK)=RD-APHM1(JL,JK)*ALNPR(JL,JK)
     &         /(APHM1(JL,JK+1)-APHM1(JL,JK))
        END DO
      END DO
C
C
C*    4.    COMPUTE SOME FIELDS NEEDED BY THE PHYSICAL ROUTINES.
C           ------- ---- ------ ------ -- --- -------- ---------
C
C*    4.1   COMPUTE VIRTUAL TEMPERATURE AT T-DT AND SET *ZGEO* TO 0.
C           ------- ------- ----------- -- ---- --- --- ------ -- --
C
      DO JLEV=1,NLEV
        DO JHOR=KSTART,KSTOP
          ZTVM1(JHOR,JLEV)=TM1(JHOR,JLEV)*(1.+VTMPC1*QM1(JHOR,JLEV)
     &                    -(XM1(JHOR,JLEV)+XIM1(JHOR,JLEV)))
        END DO
      END DO
C
C*    4.4   COMPUTE (PHI-PHIS) AT T-DT USING LN(P) AT T.
C           ------- ---------- -- ---- ----- ----- -- --
C
      CALL GEOPEC
     &     (NHOR,KSTART,KSTOP,NLEV,NLEVM1,
     &     GEOM1,ZTVM1,ALNPR,ALPHA,ZGEO)
C
C
C*    6.    VERTICAL EXCHANGE OF U,V,T,Q BY TURBULENCE.
C           -------- -------- -- - - - - -- -----------
C
      DO JL=1,NHOR
        EVAPM(JL)=EVAP(JL)
      END DO
CSH
C     *** SETZEN DER 5 SOIL MOISTURE LAYERS DES VORHERIGEN ZEITSCHRITTS
C
      IF (I5LAYER.GE.1) THEN
        DO JL=1,NHOR
          WSIM1M(JL,1) = WS1(JL)
          WSIM1M(JL,2) = WS2(JL)
          WSIM1M(JL,3) = WS3(JL)
          WSIM1M(JL,4) = WS4(JL)
          WSIM1M(JL,5) = WS5(JL)
        END DO
      END IF
C
CSH   *** SOIL LAYER THICKNESSES
      ZDEL(1) = 0.065
      ZDEL(2) = 0.254
      ZDEL(3) = 0.913
      ZDEL(4) = 2.902
      ZDEL(5) = 5.700
CSH   *** INITIALIZATION OF WORKING ARRAYS FOR 5 SOIL LAYER SCHEME
      IF (I5LAYER.GE.1) THEN
        CALL SOILDEF(KSTART, KSTOP, NHOR, NDEEP,
     &       WSMX, DZR, DZS, VPOR, LOLAND, ZDEL,
     &       DZRSI, DZSI, ZWSAT, ZWSFC)
C
C     *** NO 5 LAYERS SCHEME --> DUMMY BELEGUNG
      ELSE
        CALL SETRA(DZRSI,NHOR*NDEEP,0.)
        CALL SETRA(DZSI ,NHOR*NDEEP,0.)
        CALL SETRA(ZWSAT,NHOR*NDEEP,0.)
        CALL SETRA(ZWSFC,NHOR*NDEEP,0.)
      ENDIF
CSH
      CALL VDIFF
     &    (KSTART, KSTOP , NHOR   , KTDIA  , NLEV  , NLEVM1,
     &     NLEVP1, NSTEP , NSTART , TWODT  , EPS   , LVDIFF,
C          - INPUT 2D .
     &     ACLC  , APHM1 , APM1   , GEOM1  , QM1   , TKEM  ,
     &     TKEM1M, TM1   , UM1    , VM1    , XM1   , ZTVM1 ,
C          - INPUT 1D .
     &     INFRL , INFRW , INFRI  , AHFL   , AHFS  , AHFSL ,
     &     AHFSW , AHFSI , AZ0    , AZ0L   , AZ0W  , AZ0I  ,
     &     DEW2  , EVAPM , EVAPL  , EVAPW  , EVAPI , SNM1M ,
     &     SRFL  , TEMP2 , T2MAX  , T2MIN  , USTAR3, USTR  ,
     &     U10   , VDIS  , VSTR   , V10    , WIMAX , WIND10,
     &     USTRL , USTRW , USTRI  , VSTRL  , VSTRW , VSTRI ,
     &     WLM1M , WSM1M , WSMX   , VLT    , SINLAT,
C          - OUTPUT 2D .
     &     TKE   , Z2ZM  , Z2ZML  , Z2ZMW  , Z2ZMI , Z2ZH  ,
     &     Z2ZHL , Z2ZHW , Z2ZHI  ,
C          - OUTPUT 1D .
     &     CVS   , CVW   , DHFQS  , DHFQW  , DHFT  , ZDHFTI,
     &     EVAP  , QHFL  , ZQHFLA , RSFL   , THFL  , WLMX  ,
     &     XHFL  ,
C          - INPUT/OUTPUT 1D .
     &     VGRAT , ZQDB  , ZQDBL  , ZQDBW  , ZQDBI , ZTHFLI,
     &     TSLM1M, TSWM1M, TSIM1M , WI3M1M , LOLAND, LOSEA ,
     &     LOICE , GHPBL ,
C          - SOIL MOISTURE LAYERS / EVAP. FLUXES OVER LAND
     &     NDEEP , WSIM1M, AETRANS, AEBSOIL, AESNOW, AESKIN,
C          - CLOUD ICE
     &     XIM1  , DZR   , I5LAYER, DZRSI  , ZWSFC)
C
C
C*    6.1   GRAVITY WAVE DRAG PARAMETERISATION
C           ------- ---- ---- ----------------
C
      IF(LGWDRAG) THEN
C
C     TODO: Here is one more actual arguments than dummy arguments
C           in subroutine gwdrag.f. Something is wrong here with
C           the interface.
C
        CALL GWDRAG
     &       (KSTART, KSTOP , NHOR  , KTDIA, NLEV  , NLEVM1, NLEVP1,
     &       NSTEP , NSTART, TWODT , APHM1, APM1  , GEOM1 , TM1   ,
     &       UM1   , VM1   , VAR  , USTRGW, VDISGW, VSTRGW,
C            ----- OUTPUT 2D --------
     &       TGWDR , VGWDR , UGWDR)
C
      END IF
C
C
C*    7.    CONVECTION PARAMETERISATION.
C           ---------- -----------------
C
      DO JL=1,NHOR
        ITYPE(JL)=0
        RSFC(JL)=0.
        SSFC(JL)=0.
      END DO
C
      ZSDISS=0.
      ZSRAIN=0.
      ZSEVAP=0.
      ZSHEAT=0.
      ZSMELT=0.
C
C     7.1   STORE CURRENT TENDENCIES
C           ----- ------- ----------
C
      DO JK=1,NLEV
        DO JL=1,NHOR
          TCONV(JL,JK)=ZTE(JL,JK)
          QCONV(JL,JK)=QE (JL,JK)
          UCONV(JL,JK)=ZUE(JL,JK)
          VCONV(JL,JK)=ZVE(JL,JK)
        END DO
      END DO
C
C
C*    7.2   CALL SUBROUTINE CUCALL FOR CUMULUS PARAMETERIZATION
C           ---------------------------------------------------
C
      IF (LCONV) THEN
C
        CALL CECALL
     &    (KSTART, KSTOP , NHOR  , NLEV   , NLEVP1,
     &     NLEVM1, ILAB  , NSTEP , NSTART , TWODT ,
     &     TM1   , QM1   , UM1   , VM1    , XM1   ,
C          TE    , QE    , UE    , VE     , XE    ,
     &     TCONV , QCONV , UCONV , VCONV  , ZXE   ,
     &     VERVEL, ZQHFLA, ZXTEC , CEVAPCU, APM1  ,
     &     APHM1 , GEOM1 , LALAND, RSFC   , SSFC  ,
     &     APRC  , APRS  , ITYPE , TOPMAX , ZSRAIN,
     &     ZSEVAP, ZSHEAT, ZSDISS, ZSMELT , CAPE  ,
     &     XIM1  , ZXIE  , TSM1M)
C
      END IF
C
C
C     7.3   DETERMINE DEEP CONVECTION TENDENCIES FROM MASS-FLUX SCHEME
C           --------- ---- ---------- ---------- -------------- ------
C
      DO JK=1,NLEV
        DO IK=1,NHOR
CDJ      ----  TENDENCIES FROM TIEDKE USED IN LWCOND -----------
          TCOND(IK,JK) = TCONV(IK,JK)
          QCOND(IK,JK) = QCONV(IK,JK)
CDJ
          TCONV(IK,JK) = TCONV(IK,JK) - ZTE(IK,JK)
          QCONV(IK,JK) = QCONV(IK,JK) -  QE(IK,JK)
          UCONV(IK,JK) = UCONV(IK,JK) - ZUE(IK,JK)
          VCONV(IK,JK) = VCONV(IK,JK) - ZVE(IK,JK)
        END DO
      END DO
C
C
C*    8.    LARGE SCALE WATER PHASE CHANGES
C           ----- ----- ----- ----- -------
C
C     8.1    STORE CURRENT TENDENCIES
C            ----- ------- ----------
C
      DO JK=1,NLEV
        DO JL=1,NHOR
          ZTE(JL,JK)=TCOND(JL,JK)
          QE (JL,JK)=QCOND(JL,JK)
          ZXE(JL,JK)=XCOND(JL,JK)
          ZXIE(JL,JK)=XICOND(JL,JK)
        END DO
      END DO
C*    3.12 INITIALISATION OF CLOUD DROPLET NUMBER CONCENTRATION
C          (1/M**3) USED IN RADLSW AND LWCOND
C
      DO JK=1,NLEV
        DO JL=1,NHOR
          NEXP=2
          ZPRAT=(MIN(8.,80000./APM1(JL,JK)))**NEXP
          IF (LOLAND(JL).AND.(.NOT.LOGLAC(JL))) THEN
            ZN1= 50.
            ZN2=220.
          ELSE
            ZN1= 50.
            ZN2= 100.
          ENDIF
          IF (APM1(JL,JK).LT.80000.) THEN
            PACDNC(JL,JK)=1.e6*(ZN1+(ZN2-ZN1)*(EXP(1.-ZPRAT)))
          ELSE
            PACDNC(JL,JK)=ZN2*1.e6
          ENDIF
        END DO
      END DO
C
      IF (LCOND) THEN
C
        CALL LWCOND
     &       (KSTART, KSTOP , NHOR  , KTDIA, NLEV ,
     &       NLEVP1, TWODT , NSTART, NSTEP,
C            -- INPUT 2D:
C            -- ACHTUNG: MIT ABSICHT WURDE HIER 2 MAL DER DRUCK ZUM ZEITPUNKT
C            -- T - DT  UEBERGEBEN!!!!!!!!!!!!!!!!!!!!!!!!
     &       ILAB  , APHM1 , APHM1 , APM1 , APM1 ,
     &       GEOM1 , QM1   , TM1   , XM1  , ZXTEC,
C            -- INPUT 1D:
     &       LALAND, ITYPE ,
C            -- OUTPUT 2D:
     &       ACLC  , ACLCAC,
C            -- OUTPUT 1D:
     &       ACLCOV, ALWCVI, APRL  , QVI  , SSFL ,
C            -- INPUT 2D:
     &       TCOND , QCOND , XCOND ,
C            -- INPUT/OUTPUT 1D:
     &       APRS  , RSFL, XIM1, XIVI, XICOND, PACDNC, VERVEL, VCT,
C            -- INPUT/OUTPUT RAIN DIAGNOSTIC
     &       RPRAC)
C
      ENDIF

C
C     8.2    DETERMINE TENDENCIES FROM LWCOND
C            --------- ---------- ---- ------
C
      DO JK=1,NLEV
        DO IK=1,NHOR
          ZTM1 (IK,JK) =  TM1 (IK,JK) + TWODT*TCOND(IK,JK)
          ZQM1 (IK,JK) =  QM1 (IK,JK) + TWODT*QCOND(IK,JK)
          ZXM1 (IK,JK) =  XM1 (IK,JK) + TWODT*XCOND(IK,JK)
          ZXIM1(IK,JK) = XIM1 (IK,JK) + TWODT*XICOND(IK,JK)
          TCOND(IK,JK) = TCOND(IK,JK) - ZTE(IK,JK)
          QCOND(IK,JK) = QCOND(IK,JK) -  QE(IK,JK)
          XCOND(IK,JK) = XCOND(IK,JK) - ZXE(IK,JK)
          XICOND(IK,JK)=XICOND(IK,JK) - ZXIE(IK,JK)
        END DO
      END DO
C
C
C*    5.    RADIATION PARAMETERISATION.
C           --------- -----------------
C
      IF (NRADIA.EQ.0) THEN
C
        IF (LAEROZ) THEN
C
          CALL AORAD
     &         (NSTEP , NHOR  , NLEV  , NLEVP1, KSTART ,
     &         KSTOP , CONACC, DTIME , TWODT , NCBASE,
     &         NTBASE, NTIMST, NMONTH, NRADFR,
     &         LRAD  , LSOLC , LAER  , LCFC  , LGADSRH,
C              -- PROGNOSTIC FIELDS AND TENDENCIES --
     &         ZTM1  , ZQM1  , TRAD  , ZXM1  , APM1  , APHM1  ,
     &         LALAND, LOGLAC, SINLAT, COSLAT, SINLON, COSLON ,
C              -- SURF TEMP AND MOIST FIELDS  --
     &         TSM1M , SNM1M , SRFL  , ACLC  , TSLM1M,
     &         TSIM1M, INFRL , INFRW , INFRI ,
C              -- ALBEDOS ---
     &         FOREST, ALB   , ALBEDO, ALSOL , ALSOW , ALSOI  ,
C              -- EMISSIVITIES --
     &         EMTER , TRSOL , EMTEF , TRSOF , ACLCV , SRADS  ,
     &         SRADSU, SRAD0 , SRAD0U, TRADS , TRADSU, TRAD0  ,
     &         SCLF0 , SCLFS , SRAF0 , SRAFS , TCLF0 , TCLFS  ,
     &         TRAF0 , TRAFS , CVDAES, CVDAEL, CVDAEU, CVDAED ,
     &         OZONPL, NOZ   , SO4ALL, SO4NAT)
C
        ELSE
C
          CALL RAD
     &         (NSTEP , NHOR  , NLEV  , NLEVP1, KSTART ,
     &         KSTOP , CONACC, DTIME , TWODT , NCBASE,
     &         NTBASE, NTIMST, NMONTH, NRADFR,
     &         LRAD  , LSOLC , LAER  , LCFC  ,
C              -- PROGNOSTIC FIELDS AND TENDENCIES --
     &         ZTM1  , ZQM1  , TRAD  , ZXM1  , APM1  , APHM1  ,
     &         LALAND, LOGLAC, SINLAT, COSLAT, SINLON, COSLON ,
C              -- SURF TEMP AND MOIST FIELDS  --
     &         TSM1M , SNM1M , SRFL  , ACLC  , TSLM1M,
     &         TSIM1M, INFRL , INFRW , INFRI ,
C              -- ALBEDOS ---
     &         FOREST, ALB   , ALBEDO, ALSOL , ALSOW , ALSOI  ,
C              -- EMISSIVITIES --
     &         EMTER , TRSOL , EMTEF , TRSOF , ACLCV , SRADS  ,
     &         SRADSU, SRAD0 , SRAD0U, TRADS , TRADSU, TRAD0  ,
     &         SCLF0 , SCLFS , SRAF0 , SRAFS , TCLF0 , TCLFS  ,
     &         TRAF0 , TRAFS , CVDAES, CVDAEL, CVDAEU, CVDAED ,
     &         ZXIM1)
C
        ENDIF
C
      ENDIF
C
C
C*    9.    COMPUTATION OF NEW SURFACE VALUES.
C           ----------- -- --- ------- -------
C
      CALL SURF
     &     (KSTART , KSTOP  , NHOR  , NLEVP1 , NDEEP  ,
     &     NSTEP  , NSTART , IEXC   , WSMX   , RGCGN  ,
     &     TLAMBDA, DLAMBDA, PORVOL , FCAP   , LWDIF  ,
     &     TWODT  , LSURF  , LOGLAC , LOLAND , TSL    , TSLM1M ,
     &     WS     , WSM1M  , SN     , SNM1M  , WL     , WLM1M  ,
     &     TD     , TDM1M  , SRFL   , THFL   , QHFL   , XHFL   ,
     &     RSFC   , SSFC   , RSFL   , SSFL   , ALSOL  , ALBEDO ,
     &     DHFT   , DHFQW  , DHFQS  , EVAP   , EVAPM  , VGRAT  ,
     &     CVS    , WLMX   , EMTER  , DSNAC  , RUNOFF ,
     &     SNMEL  , TDCL   , TDCLM1M, TD3    , TD3M1M , TD4    ,
     &     TD4M1M , TD5    , TD5M1M , WI3    , WI3M1M , WI4    ,
     &     WI4M1M , WI5    , WI5M1M , WI     , WIM1M  , WICL   ,
     &     WICLM1M, TSLIN  , TSN    , TSNM1M ,
     &     TSURF  , VAROR  , DRAIN  , BETA   , WMINLOK, WMAXLOK,
     &     I5LAYER, WSI    , WSIM1M , DZR    , FKSAT  ,
     &     FMPOT  , BCLAPP , VPOR   , AETRANS, AEBSOIL, AESNOW ,
     &     AESKIN , AERES  , ZDEL   , DZRSI  , DZSI   , ZWSAT  ,
     &     ZWSFC)
C
CSH
      DO JL=KSTART,KSTOP
        IF (LOLAND(JL)) THEN
          ETRANS(JL)=ETRANS(JL)+AETRANS(JL)
          EBSOIL(JL)=EBSOIL(JL)+AEBSOIL(JL)
          ESNOW(JL)=ESNOW(JL)+AESNOW(JL)
          ESKIN(JL)=ESKIN(JL)+AESKIN(JL)
          ERES(JL)=ERES(JL)+AERES(JL)
        ENDIF
      END DO
CSH
C
C*    9.1   CALL FOR SKINTEM
C           ---- --- -------
C
      CALL SKINTEM
     &     (KSTART, KSTOP , NHOR, NLEVP1, TWODT , EMTER ,
     &     TSM1M , SICED , TEFF, ALBEDO, ALSOI , SRFL  ,
     &     ZTHFLI, ZDHFTI, TSI , TSIM1M, AHFICE, QRES  ,
     &     TSLIN , INFRL)
C
C
C*    9.2    TIME FILTER
C
C
      IF (NSTEP.NE.NSTART) THEN
C***
        DO JL=KSTART,KSTOP
          IF (LOLAND(JL)) THEN
            TSLM(JL)=TSLM(JL)+EPS*(TSLM1M(JL)-2.*TSLM(JL)+TSL(JL))
          ENDIF
          IF (LOSEA(JL)) THEN
            TSWM(JL)=TSWM(JL)+EPS*(TSWM1M(JL)-2.*TSWM(JL)+TSW(JL))
          ENDIF
          IF (LOICE(JL)) THEN
            TSIM(JL)=TSIM(JL)+EPS*(TSIM1M(JL)-2.*TSIM(JL)+TSI(JL))
          ENDIF
          TDM(JL)=TDM(JL)+EPS*(TDM1M(JL)-2.*TDM(JL)+TD(JL))
          TSNM(JL)=TSNM(JL)+EPS*(TSNM1M(JL)-2.*TSNM(JL)+TSN(JL))
          TD3M(JL)=TD3M(JL)+EPS*(TD3M1M(JL)-2.*TD3M(JL)+TD3(JL))
          TD4M(JL)=TD4M(JL)+EPS*(TD4M1M(JL)-2.*TD4M(JL)+TD4(JL))
          TD5M(JL)=TD5M(JL)+EPS*(TD5M1M(JL)-2.*TD5M(JL)+TD5(JL))
          TDCLM(JL)=TDCLM(JL)+EPS*(TDCLM1M(JL)-2.*TDCLM(JL)+TDCL(JL))
          WI3M(JL)=WI3M(JL)+EPS*(WI3M1M(JL)-2.*WI3M(JL)+WI3(JL))
          WI4M(JL)=WI4M(JL)+EPS*(WI4M1M(JL)-2.*WI4M(JL)+WI4(JL))
          WI5M(JL)=WI5M(JL)+EPS*(WI5M1M(JL)-2.*WI5M(JL)+WI5(JL))
          WIM(JL)=WIM(JL)+EPS*(WIM1M(JL)-2.*WIM(JL)+WI(JL))
          WICLM(JL)=WICLM(JL)+EPS*(WICLM1M(JL)-2.*WICLM(JL)+WICL(JL))
C
CSH       *** WS, WL UND SN OHNE ASSELIN
C         *** BEI 5 LAYER SOIL SCHEME IST ASSELINFILTER FUER WS
C         *** NICHT MEHR NOTWENDIG, DA IN SURF WS AUS WSI ERSTELLT WIRD.
          ZDIAGS=0.5
C
          IF (I5LAYER.GE.1) THEN
            WSM(JL)=WS(JL)
          ELSE
            WSM(JL)=WSM1M(JL) + (WS(JL)-WSM1M(JL)) *ZDIAGS
            WS (JL)=WSM(JL)
          ENDIF
C         ***  --> IN EC4ORG:  WSM1M=WS', WSM=WS' UND WSM WIRD NICHT
C         ***      GESONDERT GENUTZT ODER RAUSGESCHRIEBEN,
          WLM(JL)=WLM1M(JL)+(WL(JL)-WLM1M(JL))*ZDIAGS
          WL (JL)=WLM(JL)
          SNM(JL)=SNM1M(JL)+(SN(JL)-SNM1M(JL))*ZDIAGS
          SN (JL)=SNM(JL)
CSH
        END DO
C***
      ENDIF
C
C
CSH  *** STORE UPDATED 5 LAYER SOIL MOISTURE VALUES IN SINGLE ARRAYS
      IF (I5LAYER.GE.1) THEN
        DO JL=1,NHOR
          WS1(JL)=WSI(JL,1)
          WS2(JL)=WSI(JL,2)
          WS3(JL)=WSI(JL,3)
          WS4(JL)=WSI(JL,4)
          WS5(JL)=WSI(JL,5)
        END DO
      ENDIF
CSH
C
C*    9.3    NEW AVERAGED SURFACE TEMPERATURES, SURFACE MAXIMUM
C            AND MINIMUM TEMPERATURES
C
C
      DO JL=KSTART,KSTOP
        TSM(JL)= (FLOAT(INFRL(JL)) * TSLM(JL)
     &         +  FLOAT(INFRW(JL)) * TSWM(JL)
     &         +  FLOAT(INFRI(JL)) * TSIM(JL))*EDFAKINF
        TS(JL) = (FLOAT(INFRL(JL)) * TSL(JL)
     &         +  FLOAT(INFRW(JL)) * TSW(JL)
     &         +  FLOAT(INFRI(JL)) * TSI(JL))*EDFAKINF
        TSMAX(JL)=AMAX1(TS(JL),TSMAX(JL))
        TSMIN(JL)=AMIN1(TS(JL),TSMIN(JL))
      END DO
C
      RETURN
      END SUBROUTINE PHYORG
