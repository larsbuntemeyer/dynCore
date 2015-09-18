C
C**** *SURF* - UPDATES LAND VALUES OF TEMPERATURE, MOISTURE AND SNOW.
C
C     J.F.GELEYN     E.C.M.W.F.     08/06/82.
C     MODIFIED BY
C     C.A.BLONDIN    E.C.M.W.F.    18/12/86.
C     MODIFIED BY L.DUMENIL      MET.INST.HH     20/05/88
C
C     MODIFICATION
C     ------------
C
C     1. SCHEME MODIFIED FOR CLIMATE STUDIES BY USING A FIVE LAYER
C        SCHEME AS PROPOSED BY WARRILOW ET AL.(1986)
C        SOIL TEMPERATURES ARE TD3,TD4,TD5,TD,TDCL,
C     2. EXTRA TEMPERATURE VARIABLE FOR SNOW LAYER TSN IF SNOW IS
C        DEEPER THAN ZSNCRI
C     3. SOIL HYDROLOGY SCHEME BASED ON CATCHMENT CONSIDERATIONS
C
C     ****** IMPLEMENTATION OF IMPROVED ARNO SCHEME
C            STEFAN HAGEMANN, MPI-M, FEBRUARY 2005
C            HAGEMANN AND DUMENIL GATES, 2003,
C            IMPROVING A SUBGRID RUNOFF PARAMETERIZATION SCHEME FOR CLIMATE
C            MODELS BY THE USE OF HIGH RESOLUTION DATA DERIVED FROM SATELLITE
C            OBSERVATIONS
C            CLIM. DYN. 21, PP. 349-359
C
C
C     *** NOV. 2005 - APRIL 2006 -- STEFAN HAGEMANN --
C     *** PREPARAPTORY CHANGES FOR 5 LAYER SOIL SCHEME
C     *** PLUS A VERY FEW CHANGES TO AVOID UNNESSARY CONFUSING CODE
C
C
C     *** JULI 2006 -- STEFAN HAGEMANN
C     *** IMPLEMENTATION OF 5 LAYER SOIL SCHEME
C     *** SOIL MOISTURE CHANGES BELONG TO ONE TIME STEP = 0.5 * TWODT
C
C     *** NOV. 2006 -- STEFAN HAGEMANN
C     *** SMALL INCONSISTENCY CORRECTED IN THE CALCULATION OF THE
C         ROOT ZONE SOIL MOISTURE. ONLY MATTERS IN GRIDBOXES WHERE THE
C         BEDROCK IS WITHIN THE SAME SOIL LEVEL (BUT BELOW) AS THE
C         ROOTING DEPTH:
C         IN WS(JL) = WS(JL) + WSI(JL,JK) * DZRSI(JL,JK)/DZSI(JL,JK)
C         DZSI(JL,JK) HAS REPLACED ZDEL(JK)
C
C     *** SEPT. 2007 -- STEFAN HAGEMANN
C     ***
C     *** IN ROUTINE SOILCHANGE, THE FOLLOWING IS NOW HAPPENING:
C     *** 1. IF WSI < WWILT FOR ONE LAYER, TRANSPIRATION AMOUNT WILL BE REDUCED
C     ***     INSTEAD OF RE-DISTRIBUTION TO OTHER LAYERS.
C     ***     THUS, AETRANS IS CHANGED IN SOILCHANGE, BUT ALSO EVAP AND EVAPL
C     ***     MUST BE CHANGED.
C     *** BARE SOIL EVAPORATION WILL ONLY BE TAKEN FROM THE MOST UPPER LAYER.
C     *** 2.  OVERSHOOTING AMOUNTS WILL BE DISREGARDED AND AEBSOIL WILL BE
C     ***     REDUCED. THUS, AEBSOIL IS CHANGED IN SOILCHANGE, BUT ALSO EVAP
C     ***     AND EVAPL MUST BE CHANGED.
C     *** 3.  OVERSHOOTING AMOUNTS IN AESKIN WILL NOT BE TAKEN OUT OF THE
C     ***     SOIL USING AEBSOIL. INSTEAD EVAP AND EVAPL
C     ***     WILL BE REDUCED ACCORDINGLY.
C     *** 1. & 2. & 3.  --> NEW ARRAY SUBMITTED TO ROUTINE: REDEVAP
C     *** SOIL LAYER THICKNESS DEFINITION ZDEL HAS MOVED FROM SURF
C     ***   TO PHYORG.
C
C     *** MARCH 2008 - S. HAGEMANN
C     *** REDUCTION OF EVAPORATION IS TAKEN OUT AGAIN AS THERE WAS NO FEEDBACK
C         TO THE ATMOSPHERE. AN IMPLEMENTATION OF THIS FEEDBACK WOULD REQUIRE
C         MORE COMPLEX CHANGES TO VDIFF AND COMPREHENSIVE TESTING. IN
C         MY OPINION A REDUCTION OF TRANSPIRATION IN VDIFF WOULD BE DESIRABLE.
C         FELD REDEVAP BELIBT ERHALTEN, UM ALS DIAGNOSTIC ZU FUNKTIONIEREN,
C         WIRD ABER NICHT MEHR AUS PHYORG UEBERGEBEN.
C     *** UM TRANSPIRATION IN VDIFF REDUZIEREN ZU KNNEN, MUSS 5 SOIL LAYER INFO
C         SCHON IN PHYORG VORHANDEN SEIN. DAHER WIRD CALL DER ROUTINE SOILDEF
C         NUN NACH PHYORG VERLEGT, UND SOMIT DIE FELDER NACH SURF BERGEBEN.
C
C     PURPOSE.
C     --------
C
C          THIS ROUTINE UPDATES THE LAND VALUES OF SURFACE TEMPERATURE,
C     DEEP TEMPERATURE, SKIN SOIL WATER,
C     SURFACE SOIL MOISTURE (EXPRESSED AS A WATER
C     CONTENT), DEEP SOIL MOISTURE (IN COMPARABLE NUMBERS) AND SNOW
C     DEPTH (IN WATER EQUIVALENT). THIS IS DONE VIA A FORWARD TIME STEP
C     DAMPED WITH SOME IMPLICIT LINEAR CONSIDERATIONS: AS IF ALL FLUXES
C     THAT EXPLICITELY DEPEND ON THE VARIABLE HAD ONLY A LINEAR
C     VARIATION AROUND THE T-1 VALUE. FOR CONSISTENCY WITH THE
C     ATMOSPHERIC TREATMENT A TIME FILTER IS APPLIED ON ALL FIVE SURFACE
C     PRONOSTIC VARIABLES. HOWEVER ONLY THE FORWARD HALF OF THE TIME
C     FILTER IS USED TO AVOID KEEPING INDEFINITELY TRACES OF SNOW FOR
C     EXAMPLE. CLIMATIC TEMPERATURE AND MOISTURE, IN A THIRD, DEEPER
C     LAYER, ARE USED AS LOWER BOUNDARY CONDITIONS.
C
C**   INTERFACE.
C     ----------
C
C          *SURF* IS CALLED FROM *PHYSC*.
C          THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE:
C     TS,TD,WL,WS,WD,SN AT T-1,
C     SURFACE FLUXES COMPUTED IN OTHER PARTS OF
C     THE PHYSIC (THE LONG-WAVE RADIATIVE FLUX IS RECOMPUTED HERE FROM
C     THE EMISSIVITY), CLIMATIC T AND W AND LAND-SEA MASK. IT RETURNS
C     ITS OUTPUT TO THE SAME SPACE: SAME VARIABLES AT T+1 AND FILTERED
C     VALUES OF THE SAME VARIABLES AT T.
C
C     METHOD.
C     -------
C
C          STRAIGHTFORWARD ONCE THE DEFINITION OF THE CONSTANTS IS
C     UNDERSTOOD. FOR THIS REFER TO DOCUMENTATION. FOR THE TIME FILTER
C     SEE CORRESPONDING PART OF THE DOCUMENTATION OF THE ADIABATIC CODE.
C
C     EXTERNALS.
C     ----------
C
C          NONE.
C
C     REFERENCE.
C     ----------
C
C          SEE SOIL PROCESSES' PART OF THE MODEL'S DOCUMENTATION FOR
C     DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.
C
C
C
      SUBROUTINE SURF
     &   (KIDIA  , KFDIA  , KLON  , KLEVP1, KDEEP,
     &    NSTEP  , NSTART , IEXC  , WSMX  , RGCGN ,
     &    TLAMBDA, DLAMBDA, PORVOL, FCAP  , LWDIF ,
     &    TWODT  , LSURF  , LOGLAC, LOLAND,
C         -- SURF TEMP AND MOIST FIELDS  --
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
     &    TSL, TSLM1M,
     &    WS , WSM1M ,
     &    SN , SNM1M ,
     &    WL , WLM1M ,
     &    TD , TDM1M ,
CTS 250100
C         -- SURFACE FLUXES  --
     &    SRFL, THFL, QHFL, XHFL,
     &    RSFC, SSFC, RSFL, SSFL,
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
     &    ALSOL, ALBEDO,
CTS 250100
C         -- FLUX DERIVATIVES --
     &    DHFT, DHFQW, DHFQS,
C         -- EVAPORATION ---
     &    EVAP, EVAPM  ,
     &    VGRAT, CVS, WLMX, EMTER,
     &    DSNAC, RUNOFF, SNMEL  ,
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
     &    TDCL , TDCLM1M,
     &    TD3  , TD3M1M ,
     &    TD4  , TD4M1M ,
     &    TD5  , TD5M1M ,
     &    WI3  , WI3M1M ,
     &    WI4  , WI4M1M ,
     &    WI5  , WI5M1M ,
     &    WI   , WIM1M  ,
     &    WICL , WICLM1M,
     &    TSLIN,
     &    TSN  , TSNM1M ,
CTS 250100
     &    TSURF, VAROR  , DRAIN,
     &    BETA , WMINLOK, WMAXLOK,
CSH  *** - SOIL MOISTURE LAYERS, EVAP. FLUXES OVER LAND
     &    I5LAYER, WSI   , WSIM1M, DZR    ,
     &    FKSAT  , FMPOT , BCLAPP, VPOR   ,
     &    AETRANS, AEBSOIL, AESNOW, AESKIN, AERES,
     &    ZDEL   , DZRSI  , DZSI  , ZWSAT , ZWSFC)
C
      IMPLICIT NONE
C     
      INCLUDE "COMCON"
      INCLUDE "COMPH2"
      INCLUDE "COMVEG"
C not used      INCLUDE "faktinf.h"
      INCLUDE "higkon.h"
C-----------------------------------------------------------------------
C     Dummy Arguments
C
      INTEGER, INTENT(IN) :: KIDIA  , KFDIA  , KLON  , KLEVP1, KDEEP,
     &                       NSTEP  , NSTART , IEXC
      REAL, INTENT(INOUT) :: WSMX  (KLON), RGCGN  (KLON),
     &     TLAMBDA(KLON), DLAMBDA(KLON),
     &     PORVOL(KLON), FCAP   (KLON),
C          -- SURF TEMP AND MOIST FIELDS  --
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
     &     TSL(KLON), TSLM1M(KLON),
     &     WS(KLON) , WSM1M(KLON) ,
     &     SN(KLON) , SNM1M(KLON) ,
     &     WL(KLON) , WLM1M(KLON) ,
     &     TD(KLON) , TDM1M(KLON) ,
CTS 250100
C          -- SURFACE FLUXES  --
     &     SRFL(KLON), THFL(KLON), QHFL(KLON), XHFL(KLON),
     &     RSFC(KLON), SSFC(KLON), RSFL(KLON), SSFL(KLON),
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
     &     ALSOL(KLON), ALBEDO(KLON),
CTS 250100
C          -- FLUX DERIVATIVES --
     &     DHFT(KLON), DHFQW(KLON), DHFQS(KLON),
C          -- EVAPORATION AND VEGETATION ---
     &     EVAP(KLON), EVAPM(KLON), VGRAT(KLON),
     &     CVS(KLON) , WLMX(KLON) ,
C          --  EMISSIVITIES --
     &     EMTER (KLON,KLEVP1),
C          -- REMAINING ELEMENTS IN *SURF* --
     &     DSNAC(KLON), RUNOFF(KLON), SNMEL(KLON)  ,
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
     &     TDCL(KLON) , TDCLM1M(KLON),
     &     TD3(KLON)  , TD3M1M(KLON) ,
     &     TD4(KLON)  , TD4M1M(KLON) ,
     &     TD5(KLON)  , TD5M1M(KLON) ,
     &     WI3(KLON)  , WI3M1M(KLON) ,
     &     WI4(KLON)  , WI4M1M(KLON) ,
     &     WI5(KLON)  , WI5M1M(KLON) ,
     &     WI(KLON)   , WIM1M(KLON)  ,
     &     WICL(KLON) , WICLM1M(KLON),
     &     TSLIN(KLON),
     &     TSN(KLON)  , TSNM1M(KLON) ,
CTS 250100
     &     TSURF(KLON), VAROR(KLON) , DRAIN(KLON),
     &     BETA(KLON) , WMINLOK(KLON), WMAXLOK(KLON)
CSH  *** - SOIL MOISTURE LAYERS, EVAP. FLUXES OVER LAND
      REAL, INTENT(INOUT) :: WSI(KLON, KDEEP), WSIM1M(KLON, KDEEP),
     &     DZR(KLON)   ,
     &     FKSAT(KLON) , FMPOT(KLON) , BCLAPP(KLON), VPOR(KLON),
     &     AETRANS(KLON), AEBSOIL(KLON), AESNOW(KLON) , AESKIN(KLON),
     &     AERES(KLON),
     &     DZRSI(KLON,KDEEP), DZSI (KLON,KDEEP),
     &     ZWSAT(KLON,KDEEP), ZWSFC(KLON,KDEEP)
C          -- LAND AND ICE MASK      --
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
      LOGICAL, INTENT(INOUT) :: LOGLAC (KLON),  LSURF
CTS 250100
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
      LOGICAL, INTENT(INOUT) :: LOLAND(KLON)
CTS 250100
C-----------------------------------------------------------------------
C     Local Declarations
C
      LOGICAL :: LO
C     -------------------------------------------------------
C     DECLARATION OF REAL ONE-DIMENSIONAL WORKING ARRAYS
C
      REAL TEM1(KLON)  , TEM2(KLON)  , ZAIR(KLON)  , FSMELT(KLON),
     &     ZDAIR(KLON) , ZDRAIN(KLON),
     &     ZLAC(KLON)  ,
     &     ZPRFL(KLON) , ZPSFL(KLON) ,
     &     ZROS (KLON) , ZSNMEL(KLON), ZSNDP(KLON) ,
     &     ZTCOE(KLON) , ZTRFL(KLON) ,
     &     ZWL(KLON)   ,
CTS CHANGE FOR TEMPERATURE DEPENDENT HEAT CONDUCTIVITY, DENSITY
CTS OF SNOW, HEAT CONTENT OF THE UPPERMOST 10 CM OF SNOW PER
CTS KELVIN AND HEAT CONTENT OF THE UPPERMOST 10 CM OF SNOW PER
CTS KELVIN, NORMALIZED WITH THE TIME STEP
     &     ZALPHAS(KLON), ZRHOS(KLON) ,ZCPCONS(KLON), ZCPSDT(KLON),
C
C     -------------------------------------------------------
C     DECLARATION OF REAL TWO-DIMENSIONAL WORKING ARRAYS
C
     &     ZCFH (KLON,KDEEP), ZEB(KLON,KDEEP),
     &     ZSDIF(KLON,KDEEP), ZDEL(KDEEP)     , ZCGN (KLON,KDEEP),
     &     ZWI  (KLON,KDEEP), ZWWS(KDEEP)     ,
     &     ZWIS (KDEEP)     , ZTD (KLON,KDEEP), ZCONB1(KLON,KDEEP),
     &     ZDIFI(KLON,KDEEP), ZWQ (KLON)      , ZTDM1M(KLON,KDEEP)
C
CSH   *** IMPROVED ARNOSCHEME VARIABLES
C     *** NOTE: ZINFIL, ZVOL ELIMINATED
C     *** ZBETAG = GESAMT-BETAPARAMETER DES ARNOSCHEMAS
C     ***          OLD = OROGRAPHY-BETA  NEW: BETA + OROGRAPHY-BETA
      REAL ZBETAG(KLON)
C
CSH   *** SWITCH ZUM AUSSCHALTEN DER SOIL MOISTURE ABHAENGIGKEIT VON
C     *** SOIL DIFFUSIVITY AND HEAT CAPACITY: LWDIF
      LOGICAL LWDIF
C
CSH   *** VARIABLES NECESSARY FOR EVAPORATION FLUXES AND 5 SOIL LAYERS
      REAL ZFDUM(KLON)
      REAL ZINFIL(KLON)     ! INFILTRATION FOR TRANSFER TO SOILCHANGE
CSK
      REAL ZINFT            ! DUMMY ARGUMENT
CSK
      REAL REDEVAP(KLON)    ! REDUCTION OF EVAVOTRANSPIRATION FROM SOILCHANGE
      INTEGER :: I5LAYER       ! SWITCH FOR 5 LAYER CALCULATION OFF/ON
      INTEGER :: ISCH          ! SWITCH TO CONTROL SOILCHANGE
      INTEGER, PARAMETER :: ILOG=0      ! SWITCH FOR LOGOUTPUT IN SOILHYD
      INTEGER, PARAMETER :: JLLOG=1 ! GRIDBOX NUMBER FOR LOGOUTPUT IN SOILHYD
      REAL HLAMBDA, TWODT, ZBWS, ZCFHN, ZCONB2, ZCONS5, ZCONW2, ZCONW3, 
     &     ZCPCOEF, ZCPICE, ZCPSNOW, ZCPSOIL, ZCPWATE, ZDEFF, ZDHFQW, 
     &     ZDIAGS, ZDIAGT, ZDIAGW, ZDIFIZ, ZDISC
      REAL ZDQSDT, ZDQSNOW, ZDREXP, ZDRMAX, ZDRMIN, ZDSFLX, ZDSNOWFL, 
     &     ZDT, ZDTHFL, ZDTRFL, ZEMISS, ZFAC, ZIPRCP, ZMPRCP, ZORVARI, 
     &     ZORVARS, ZPSFR, ZQHFLW, ZRICI, ZROEFF
      REAL ZSA, ZSB, ZSFLX, ZSN, ZSNCRI, ZSNMELU, ZSNMLT, ZSNOWFL, 
     &     ZSODIF, ZSOFL, ZTCRITH, ZTCRITL, ZTHFLAD, ZTMST, ZTNEU, 
     &     ZTPFAC1, ZTPFAC2, ZTPFAC3, ZTPRCP, ZTRFLAD
      REAL ZVINTER, ZWDTR, ZWG, ZWIMEL, ZWMAX, ZWSLIM, ZWSUP, ZZDP, 
     &     ZZQSNFL
      INTEGER IWDIF, JK, JKK, JL, NSL
C-----------------------------------------------------------------------
C     *** TYPICAL TEST GRIDBOXES: JLLOG =  40*81+41 OR JLLOG = 81+67
               ! NO LOGOUTPUT FOR ILOG=0
C
C*    PHYSICAL CONSTANTS.
C     -------- ----------
C
C          *ZRGCG* IS THE PRODUCT OF DENSITY BY HEAT CAPACITY FOT THE
C     SOIL, *ZDIF* IS THE THERMAL DIFFUSIVITY OF THE SOIL, *ZD1* AND
C     *ZD2* ARE THE THICKNESSES OF THE TWO GROUND LAYERS (S,D) AND *ZTQ*
C     IS THE RATIO OF DIFFUSIVITIES FOR HEAT AND MOISTURE. *ZQSNCR* IS
C     IS THE INVERSE OF A CRITICAL VALUE FOR SNOW DEPTH (SEE *VDIFF*).
C
      ZSNCRI =0.01
C     M WATER EQUIVALENT  CRITICAL SNOW HEIGHT
C
C
      ZEMISS=0.996
      ZRICI=2.09E+06
      ZDIFIZ=12.E-07
      ZVINTER=CVINTER
      ZORVARI=10.
      ZORVARS=820.


CHG   CONSTANTS FOR MELT- AND FREEZINGPROCESSES
      ZTCRITH=TCRITH
      ZTCRITL=TCRITL
CKS   INITIALIZATION OF FSMELT (SECURITY)
      CALL SETRA(FSMELT,KLON,0.0)

CTS CHANGE FOR MODIFIED SNOW PARAMETRISATION
CTS SPECIFIC HEAT CAPACITY OF SNOW IN PIELKE: MESOSCALE
CTS METEOROLOGICAL MODELLING, ACADEMIC PRESS, 1984, PAGE 384:
CTS 2093 J/(KG*K)
      ZCPSNOW=2093.
CTS   SPECIFIC HEAT CAPACITY OF ICE IN PIELKE: MESOSCALE
CTS METEOROLOGICAL MODELLING, ACADEMIC PRESS, 1984, PAGE 384:
CTS 2093 J/(KG*K)
      ZCPICE = 2093.
C
C*    SECURITY PARAMETERS
C     -------- ----------
C
C
C
C*    COMPUTATIONAL CONSTANTS.
C     ------------- ----------
C
C
C     NUMBER OF VERTICAL LEVELS
C
      NSL = KDEEP
      ZTMST=TWODT
      IF(NSTEP.EQ.NSTART) ZTMST=0.5*TWODT
      ZDIAGT=0.5*TWODT
      ZDIAGW=ZDIAGT/RHOH2O
      ZDIAGS=ZDIAGT/ZTMST
C
      ZCONS5=ZTMST/RHOH2O
C
      ZPSFR=1.
      ZDRMIN=0.001/(3600.*1000.)
      ZDRMAX=0.1/(3600.*1000.)
      ZDREXP=1.5
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
C
C*                  SET CONSTANTS FOR TS SCHEME
C                   --- --------- --- -- ------
C
CSH   *** SOIL LAYER DEFINITION MOVED TO PHYORG
      ZCONB2 =1./(0.5*(ZDEL(1)+ZDEL(2)))
      ZTPFAC1=CVDIFTS
      ZTPFAC2=1./ZTPFAC1
      ZTPFAC3=1.-ZTPFAC2
CTS 250100
C
CSH   *** CALL SOILDEF MOVED TO PHYORG
C
C     ------------------------------------------------------------------
C
C*         1.     LOCATE SPACE.
C                 ------ ------
C
C***
      IF(LSURF) THEN
C***
C
C     ------------------------------------------------------------------
C
C*         2.     IF NO LAND POINT BY-PASS COMPUTATIONS.
C                 -- -- ---- ----- ------- -------------
C
C     ------------------------------------------------------------------
C
C         3.     SURFACE TEMPERATURE
C                ------- -----------
C
C         3.1.1. INPUT FOR ENERGY BALANCE, SNOW TEMPERATURE, SNOW MELT
C                ----- --- ------ -------- ---- ------------ ---- ----
C
CTS 311-LOOP REMOVED (CHANGE FOR MODIFIED SNOW PARAMETRISATION)
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
               IF(LOGLAC(JL)) THEN
                  ZZQSNFL=(EVAP(JL)-EVAPM(JL))/ZDIAGW
                  ZLAC(JL)=RSFL(JL)+RSFC(JL)+SSFL(JL)+SSFC(JL)+ZZQSNFL
               ELSE
                  ZLAC(JL)=0.
               ENDIF
            ENDIF
         ENDDO
C
C
C*        3.1.2.   CALCULATE SNOW TEMPERATURE OVER LAND
C                  --------- ---- ----------- ---- ----
C
!DIR$ IVDEP
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
               ZSNMEL(JL)=0.
C
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION AND FOR
CTS            MODIFIED SNOW PARAMETRISATION
CTS
               IF((SNM1M(JL).GE.ZSNCRI).AND.(.NOT.LOGLAC(JL))) THEN
CTS
CTS CALCULATION OF HEAT CONDUCTIVITY AND DENSITY OF SNOW:
CTS
CTS HEAT CONDUCTIVITY OF SNOW (ZALPHAS) BETWEEN 0.008 W/(M*K) AND
CTS 2.0 W/(M*K)
CTS (LINKE: METEOROLOGISCHES TASCHENBUCH IV, LEIPZIG 1939:
CTS BETWEEN 0.00002 CAL/(CM*S*K) AND 0.005 CAL/(CM*S*K))
CTS
CTS DENSITY OF SNOW (ZRHOS) BETWEEN 50 KG/M3 AND 900 KG/M3
CTS (LINKE: METEOROLOGISCHES TASCHENBUCH IV, LEIPZIG 1939)
CTS
CTS THREE POINTS DEFINE TWO STRAIGHT LINES BETWEEN THE MAXIMUM AND THE
CTS MINIMUM HEAT CONDUCTIVITY AND DENSITY OF SNOW:
CTS TSN = 263.16 K -> ZALPHAS = 0.100 W/(M*K), ZRHOS = 199.9 KG/M3
CTS TSN = 270.16 K -> ZALPHAS = 0.177 W/(M*K), ZRHOS = 300.0 KG/M3
CTS TSN = 273.16 K -> ZALPHAS = 0.300 W/(M*K), ZRHOS = 450.0 KG/M3
CTS
                  IF (TSNM1M(JL).GE.273.16) THEN
                     ZALPHAS(JL) = 0.3
                     ZRHOS  (JL) = 450.
                  ENDIF
                  IF (TSNM1M(JL).LE.263.16) THEN
                     ZALPHAS(JL) = 0.1
                     ZRHOS  (JL) = 199.9
                  ENDIF
                  IF ((TSNM1M(JL).GT.263.16).AND.
     &                (TSNM1M(JL).LT.270.16)) THEN
                     ZALPHAS(JL) = 0.011*(TSNM1M(JL)-273.16) + 0.21
                     ZRHOS(JL) = 14.3*(TSNM1M(JL)-273.16) + 342.9
                  ENDIF
                  IF ((TSNM1M(JL).GE.270.16).AND.
     &                (TSNM1M(JL).LT.273.16)) THEN
                     ZALPHAS(JL) = 0.041*(TSNM1M(JL)-273.16) + 0.3
                     ZRHOS(JL) = 50.*(TSNM1M(JL)-273.16) + 450.
                  ENDIF
                  ZSNDP(JL)=(RHOH2O/ZRHOS(JL))*SNM1M(JL)
                  ZCPCONS(JL) = ZRHOS(JL)*ZCPSNOW*AMIN1(0.1,ZSNDP(JL))
                  ZCPSDT(JL)=ZCPCONS(JL)/ZTMST
C
                  ZTRFLAD=EMTER(JL,KLEVP1)*STBO*TSLM1M(JL)**4
     &                 +4.*ZEMISS*STBO*TSLM1M(JL)**4
                  ZSOFL=(1.-ALSOL(JL))*SRFL(JL)/(1.-ALBEDO(JL))
                  ZTHFLAD=THFL(JL)-DHFT(JL)*TSLM1M(JL)
                  ZSNOWFL=ZALPHAS(JL)*TD3M1M(JL)/ZSNDP(JL)
                  ZDQSNOW=ZCPSDT(JL)*TSLM1M(JL)
                  ZSFLX=ZTRFLAD+ZTHFLAD+ZSOFL+ZSNOWFL+ZDQSNOW
C
C*    DERIVATIONS
C
                  ZDTRFL=-4.*ZEMISS*STBO*TSLM1M(JL)**3
                  ZDTHFL=DHFT(JL)
                  ZDSNOWFL=ZALPHAS(JL)/ZSNDP(JL)
                  ZDQSDT=ZCPSDT(JL)
                  ZDSFLX=ZDQSDT+ZDSNOWFL-ZDTRFL-ZDTHFL
C
                  TSL(JL)=ZSFLX/ZDSFLX
                  TEM1(JL)= ZALPHAS(JL)*(TSLM1M(JL)-TD3M1M(JL))
     &                 /ZSNDP(JL)
                  TEM2(JL)= -ZALPHAS(JL)/ZSNDP(JL)
C
CHG+KS SNOW MELTIN FROM TOP
               ELSE IF((SNM1M(JL).GT.0.0).AND.(.NOT.LOGLAC(JL))) THEN

CTS CHANGE FOR MODIFIED SNOW PARAMETRISATION
                  ZTRFL(JL)=(EMTER(JL,KLEVP1)*STBO*TSLM1M(JL)**3)
     &                 *TSLM1M(JL)
CSK ZSOFL IST COMPUTED AND USED IN SURFACE ENERGY BALANCE
                  ZSOFL=((1.-ALSOL(JL))*SRFL(JL)/(1.-ALBEDO(JL)))
                  ZAIR(JL)=-(ZSOFL+ZTRFL(JL)+THFL(JL))
                  ZDAIR(JL)=4.*ZEMISS*STBO*TSLM1M(JL)**3-DHFT(JL)

CHG+KS  FOR POSITIVE ENERGY FLUX USE 20% FOR SNOW MELTING
                  IF (ZAIR(JL).LT.0.0) THEN
                     TEM1(JL)=-0.8*ZAIR(JL)
                     FSMELT(JL)=-0.2*ZAIR(JL)
                  ELSE
                     TEM1(JL)=-ZAIR(JL)
                     FSMELT(JL)=0.0
                  ENDIF

                  TEM2(JL)=-ZDAIR(JL)

               ELSE
C
C              NO SNOW ON THE GROUND
C              -- ---- -- --- ------
C
CTS CHANGE FOR MODIFIED SNOW PARAMETRISATION
                  ZTRFL(JL)=(EMTER(JL,KLEVP1)*STBO*TSLM1M(JL)**3)
     &                 *TSLM1M(JL)
CSK ZSOFL IST COMPUTED AND USED IN SURFACE ENERGY BALANCE
                  ZSOFL=((1.-ALSOL(JL))*SRFL(JL)/(1.-ALBEDO(JL)))
                  ZAIR(JL)=-(ZSOFL+ZTRFL(JL)+THFL(JL))
                  ZDAIR(JL)=4.*ZEMISS*STBO*TSLM1M(JL)**3-DHFT(JL)
                  TEM1(JL)=-ZAIR(JL)
                  TEM2(JL)=-ZDAIR(JL)
               ENDIF
            ENDIF
C
         ENDDO
C
C
C*        3.1.3     SOLUTION OF SOIL DIFFUSION EQ. AS IN VDIFF
C                   -------- -- ---- --------- --- -- -- -----
C
CSH   *** SOIL MOISTURE DEPENDENCE OF DIFFUSIVITY/CAPACITY?
         IF (LWDIF) THEN
            IWDIF = 1
         ELSE
            IWDIF = 0
         ENDIF
C
C
         DO JL=KIDIA, KFDIA
            ZWI(JL,1) = WI3M1M(JL)
            ZWI(JL,2) = WI4M1M(JL)
            ZWI(JL,3) = WI5M1M(JL)
            ZWI(JL,4) = WIM1M(JL)
            ZWI(JL,5) = WICLM1M(JL)
         ENDDO
         DO JK = 1,NSL
            DO JL=KIDIA, KFDIA
               IF (LOLAND(JL)) THEN
                  ZWQ(JL)= WSM1M(JL)/WSMX(JL)*FCAP(JL)
C
CSH
                  IF (IWDIF.EQ.1) THEN
                     ZCGN(JL,JK) = RGCGN(JL)+RHOH2O*ZWQ(JL)*
     &                    (CLW*(1.-ZWI(JL,JK))+ZCPICE*ZWI(JL,JK))
                     IF (LOGLAC(JL)) THEN
                        ZCONB1(JL,JK) = 1./ZRICI
                     ELSE
                        ZCONB1(JL,JK) = 1./ZCGN(JL,JK)
                     ENDIF
                     HLAMBDA = (1.+0.35*DLAMBDA(JL)) /
     &                    (1.+1.95*DLAMBDA(JL))
                     HLAMBDA = (4. * ZWQ(JL) / PORVOL(JL) - 1.) 
     &                    * HLAMBDA
                     HLAMBDA = MIN((4. * ZWQ(JL) / PORVOL(JL)),
     &                              1.+HLAMBDA)
                     HLAMBDA = HLAMBDA * DLAMBDA(JL) *
     &                    (0.25 + 0.3 * DLAMBDA(JL) /
     &                    (1. + 0.75 * DLAMBDA(JL)))
                     HLAMBDA = TLAMBDA(JL) + HLAMBDA
                     ZSODIF = HLAMBDA * ZCONB1(JL,JK)
                     IF (LOGLAC(JL)) THEN
                        ZDIFI(JL,JK) = ZDIFIZ
                     ELSE
                        ZDIFI(JL,JK) = ZSODIF
                     ENDIF
                  ELSE
C                    *** REMO 5.0 BEHANDLUNG
                     ZCGN(JL,JK) = RGCGN(JL)
                     IF (LOGLAC(JL)) THEN
                        ZCONB1(JL,JK)=1./ZRICI
                        ZDIFI(JL,JK)=ZDIFIZ
                     ELSE
                        ZCONB1(JL,JK)=1./RGCGN(JL)
                        ZDIFI(JL,JK)=TLAMBDA(JL)
                     ENDIF
                  ENDIF
C
               ENDIF
            ENDDO
         ENDDO
C
         DO JL=KIDIA, KFDIA
CTS CHANGE FOR MODIFIED SNOW PARAMETRISATION
            IF(LOLAND(JL).AND.(SNM1M(JL).LT.ZSNCRI.OR.LOGLAC(JL)))THEN
               ZZDP=1./ZDEL(1)
               ZSA=ZTMST*ZCONB1(JL,1)*ZZDP
               ZSB=ZTMST*ZCONB2*ZDIFI(JL,1)*ZZDP
               ZDEFF=-ZAIR(JL)-TEM1(JL)/(1.-0.5*
     &              (ZSA*AMIN1(ZSB/ZSA,TEM2(JL))-ZSB))
               TSLIN(JL)=TSLIN(JL)+ZDIAGT*ZDEFF
            ENDIF
         ENDDO
C
C*          SETTING OF RIGHT HAND SIDE
C           --------------------------
C
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
               ZSDIF(JL,1)=ZTPFAC2*TD3M1M(JL)
               ZSDIF(JL,2)=ZTPFAC2*TD4M1M(JL)
               ZSDIF(JL,3)=ZTPFAC2*TD5M1M(JL)
               ZSDIF(JL,4)=ZTPFAC2*TDM1M(JL)
               ZSDIF(JL,5)=ZTPFAC2*TDCLM1M(JL)
            ENDIF
         ENDDO
C
C
C*            SET A(K), B=1+A+C, C(K)=A(K-1)
C             ------------------------------
C
         DO JK=1,NSL-1
            DO JL=KIDIA, KFDIA
               IF (LOLAND(JL)) THEN
                  ZCFHN=ZTPFAC1*2.*ZTMST/(ZDEL(JK)+ZDEL(JK+1))
                  ZCFH(JL,JK)=ZCFHN*(ZDIFI(JL,JK+1)*ZDEL(JK+1) +
     &                               ZDIFI(JL,JK)  *ZDEL(JK)  )/
     &                              (ZDEL(JK+1)+ZDEL(JK))
               ENDIF
            ENDDO
         ENDDO
C
C
C*           TOP LAYER ELIMINATION
C            --- ----- -----------
C
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
               ZTCOE(JL)=ZCFH(JL,1)
               ZZDP=1./ZDEL(1)
               ZDISC=1./(1.+ZCFH(JL,1)*ZZDP)
               ZEB(JL,1)=ZDISC*(ZCFH(JL,1)*ZZDP)
               ZSA=ZTMST*ZCONB1(JL,1)*ZZDP
               ZSB=ZTMST*ZCONB2*ZDIFI(JL,1)*ZZDP
               ZSDIF(JL,1)=ZDISC*(ZSDIF(JL,1)+TEM1(JL)
     &              *ZSA/(1.-0.5*(ZSA*AMIN1(ZSB/ZSA,TEM2(JL))-ZSB)))
            ENDIF
         ENDDO
C
C
C*           ELIMINATION FOR MIDDLE LAYERS
C            ----------- --- ------ ------
C
         DO JK=2,NSL-1
            DO JL=KIDIA, KFDIA
               IF (LOLAND(JL)) THEN
                  ZZDP=1./ZDEL(JK)
                  ZFAC=ZTCOE(JL)*ZZDP
                  ZTCOE(JL)=ZCFH(JL,JK)
                  ZDISC=1./(1.+ZFAC*(1.-ZEB(JL,JK-1))+ZCFH(JL,JK)
     &                 *ZZDP)
                  ZEB(JL,JK)=ZDISC*(ZCFH(JL,JK)*ZZDP)
                  ZSDIF(JL,JK)=ZDISC*(ZSDIF(JL,JK)+ZFAC*ZSDIF(JL,JK-1))
               ENDIF
            ENDDO
         ENDDO
C
C
C*          BOTTOM LAYER ELIMINATION
C           ------ ----- -----------
C
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
               ZZDP=1./ZDEL(NSL)
               ZFAC=ZTCOE(JL)*ZZDP
               ZDISC=1./(1.+ZFAC*(1.-ZEB(JL,NSL-1)))
               ZSDIF(JL,NSL)=ZDISC*(ZSDIF(JL,NSL)+
     &              ZFAC*ZSDIF(JL,NSL-1))
            ENDIF
         ENDDO
C
C
C*           BACK-SUBSTITUTION
C            -----------------
C
         DO JKK=1,NSL-1
            JK=NSL-JKK
            DO JL=KIDIA, KFDIA
               IF (LOLAND(JL)) THEN
                  ZSDIF(JL,JK)=ZSDIF(JL,JK)+ZEB(JL,JK)*ZSDIF(JL,JK+1)
               ENDIF
            ENDDO
         ENDDO
C
C
C*          RETURN TO TS VALUES
C           -------------------
C
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
               TD3(JL)=ZSDIF(JL,1)+ZTPFAC3*TD3M1M(JL)
               TD4(JL)=ZSDIF(JL,2)+ZTPFAC3*TD4M1M(JL)
               TD5(JL)=ZSDIF(JL,3)+ZTPFAC3*TD5M1M(JL)
               TD(JL) =ZSDIF(JL,4)+ZTPFAC3*TDM1M(JL)
               TDCL(JL)=ZSDIF(JL,5)+ZTPFAC3*TDCLM1M(JL)
CTS            CHANGE FOR MODIFIED SNOW PARAMETRISATION
               IF((SNM1M(JL).GE.ZSNCRI).AND.(.NOT.LOGLAC(JL))) THEN
C
C  TSN INTERPOLIEREN
C
                  TSN(JL) = TSL(JL) + ((TD3(JL) - TSL(JL))
     &                * (0.5 * ZSNDP(JL)) / (ZSNDP(JL) + 0.5 * ZDEL(1)))
               ELSE
                  TSL(JL) = TD3(JL)
                  TSN(JL) = TD3(JL)
               ENDIF
CTS
CTS            CHANGE FOR MODIFIED SOIL PARAMETRISATION
               ZTD(JL,1) = TD3(JL)
               ZTD(JL,2) = TD4(JL)
               ZTD(JL,3) = TD5(JL)
               ZTD(JL,4) = TD(JL)
               ZTD(JL,5) = TDCL(JL)
CHG
CHG            CHANGE FOR SOIL FREEZING AND MELTING
               ZTDM1M(JL,1) = TD3M1M(JL)
               ZTDM1M(JL,2) = TD4M1M(JL)
               ZTDM1M(JL,3) = TD5M1M(JL)
               ZTDM1M(JL,4) = TDM1M(JL)
               ZTDM1M(JL,5) = TDCLM1M(JL)
C
            ENDIF
CTS
         ENDDO

         DO JK=1,NSL
            DO JL=KIDIA, KFDIA

               IF (LOLAND(JL)) THEN

CHG            MODIFIED FREEZING AND MELTING
C
C              DESCRIPTION FOLLOWING LATER

                  ! MELTING AND FREEZING OF SOIL
                  IF (WSM1M(JL).GE.1.E-12) THEN ! CHANGED TO WSM1M
                     ! WATER/ICE
                     ! OLD BUCKET SCHEME
                     IF (I5LAYER.EQ.0) THEN
                        ZWWS(JK) = ZWQ(JL) *
     &                       (1. - ZWI(JL,JK)) * ZDEL(JK) ! LIQUID WATER
                        ZWIS(JK) = ZWQ(JL) * ZWI(JL,JK) * ZDEL(JK) ! ICE
                        ZWG      = ZWWS(JK) + ZWIS(JK) ! TOTAL WATER OR WSM1M
                     ! NEW BUCKET SCHEME
                     ELSE
                        ZWWS(JK) = (1. - ZWI(JL,JK)) * WSIM1M(JL,JK)
                        ZWIS(JK) = ZWI(JL,JK)        * WSIM1M(JL,JK)
                        ZWG      = ZWWS(JK) + ZWIS(JK)
                     ENDIF

                     ZWIMEL = 0.

                     ! HEAT CAPACITY
                     ! HEAT CAPACITY OF SOIL
                     ZCPSOIL = ZCGN(JL,JK)*ZDEL(JK)
                     ! HEAT CAPACITY OF WATER
                     ZCPWATE = RHOH2O*ALF
                     ! COEFFICIENT OF HEAT CAPACITY SOIL/WATER
                     ZCPCOEF = ZCPSOIL/ZCPWATE

                     ! TEMPERATURE CHANGES
                     ZDT=ZTD(JL,JK)-ZTDM1M(JL,JK)

                     ! UPDATE OF TEMPERATURE FOR MELTING AND FREEZING
                     ZTNEU = (
     &                   ( -ZWIS(JK)+ZDT*ZCPCOEF+ZTDM1M(JL,JK)*ZCPCOEF )
     &                   * (ZTCRITH-ZTCRITL) + ZTCRITH*ZWG
     &                   ) / ( ZWG+ZCPCOEF*(ZTCRITH-ZTCRITL) )

                     ! FREEZING OR MELTING OF WATER IN SOIL
                     ZWIMEL = ZDT * ZCPCOEF - ( ZTNEU-ZTDM1M(JL,JK) ) 
     &                    * ZCPCOEF

                     ! UPDATE
                     ZTD(JL,JK)=ZTNEU
                     ZWIS(JK) = ZWIS(JK) - ZWIMEL

                     ! SECURITY I
                     IF (ZWIS(JK).LT.0.) THEN
                        ZTD(JL,JK)=ZTD(JL,JK)-ZWIS(JK)/ZCPCOEF
                        ZWIS(JK)=0.
                     ENDIF

                     ! SECURITY II
                     IF (ZWIS(JK).GT.ZWG) THEN
                        ZTD(JL,JK)=ZTD(JL,JK)-(ZWIS(JK)-ZWG)/ZCPCOEF
                        ZWIS(JK)=ZWG
                     ENDIF

                     ! CALCULATION FRACTION OF FROZEN SOIL
                     IF (ZWG.GT.1E-4) THEN
                        ZWI(JL,JK) = ZWIS(JK) / ZWG
                        ZWI(JL,JK) = AMIN1(ZWI(JL,JK),1.)
                        ZWI(JL,JK) = AMAX1(ZWI(JL,JK),0.)
                     ELSE
                        ZWI(JL,JK) = 0.0
                     ENDIF

                  ENDIF         ! IF (WSM1M(JL).GE.1.E-12) THEN
               ENDIF            ! IF (LOLAND(JL)) THEN

            ENDDO
         ENDDO
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
               TD3(JL) = ZTD(JL,1)
               TD4(JL) = ZTD(JL,2)
               TD5(JL) = ZTD(JL,3)
               TD(JL)  = ZTD(JL,4)
               TDCL(JL)= ZTD(JL,5)
               WI3(JL) = ZWI(JL,1)
               WI4(JL) = ZWI(JL,2)
               WI5(JL) = ZWI(JL,3)
               WI(JL)  = ZWI(JL,4)
               WICL(JL)= ZWI(JL,5)
C
            ELSE
               WI3(JL) = 0.
               WI4(JL) = 0.
               WI5(JL) = 0.
               WI(JL)  = 0.
               WICL(JL)= 0.
            ENDIF
CTS
         ENDDO
C
C
C*         3.2     MOISTURE CHANGES.
C
C           EVAPORATION OF THE SKIN RESERVOIR
CSH
         CALL SETRA(ZFDUM,KLON,0.0)
C
!CDIR NOASSOC
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
               LO=DHFQW(JL).LT.0.
CTS 250100
               ZDHFQW=ZCONS5*DHFQW(JL)
               ZQHFLW=ZDHFQW*WLM1M(JL)
               IF (LO) THEN
                  ZWL(JL)=WLM1M(JL)+ZQHFLW
               ELSE
                  ZWL(JL)=WLM1M(JL)
               ENDIF
               ZWL(JL)=AMAX1(0.,ZWL(JL))
C              *** ZQHFLW IS ZERO IF NO SKIN EVAP., ELSE IT IS THE SKIN EVAP
C              *** FLUX OR THE SKIN RES. CONTENT IF IT IS EMPTIED.
               ZQHFLW=(ZWL(JL)-WLM1M(JL))/ZCONS5
               QHFL(JL)=QHFL(JL)-ZQHFLW
CSH
C              *** WL MUSS NOCH RAUFMULTIPLIZIERT WERDEN
C              ** EVAP SEPARATION WAS IN VDIFF - MUST BE ADJUSTED, HERE
C              AERES(JL) = AERES(JL) + ZQHFLW *
C              AERES(JL) = QHFL(JL) * ZDIAGW
               ZFDUM(JL) = DHFQW(JL) *WLM1M(JL) * ZDIAGW
C              *** ZFDUM = ACTUAL - DESIRED SKIN FLUX
C              *** ZFDUM > 0 IF DESIRED > SKIN RES. CONTENT
               ZFDUM(JL) = ZQHFLW  * ZDIAGW - ZFDUM(JL)
CSH
            ENDIF
         ENDDO
C
C
CSH      *** CHANGES IN AESKIN
         CALL SETRA(REDEVAP,KLON,0.0)
         DO  JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
               IF (DHFQW(JL).LT.0. .AND. ZFDUM(JL).GT.0) THEN
                  AESKIN(JL) = AESKIN(JL) + ZFDUM(JL)
                  REDEVAP(JL) = REDEVAP(JL) - ZFDUM(JL)
               ENDIF
            ENDIF
         ENDDO
CSH
C
C           COLLECTION OF DEW BY THE SKIN RESERVOIR
C
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
               ZQHFLW=ZCONS5*QHFL(JL)
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
               LO=ZQHFLW.GT.0.
CTS 250100
               IF (LO) THEN
CSH
                  IF (AETRANS(JL).GT.0) THEN
                     AESKIN(JL) = AESKIN(JL) + AETRANS(JL)
                     AETRANS(JL) = 0.
                  ENDIF
                  IF (AEBSOIL(JL).GT.0) THEN
                     AESKIN(JL) = AESKIN(JL) + AEBSOIL(JL)
                     AEBSOIL(JL) = 0.
                  ENDIF
                  ZFDUM(JL) = QHFL(JL)
CSH
                  ZWL(JL)=ZWL(JL)+ZQHFLW
                  QHFL(JL)=AMAX1(0.,ZWL(JL)-WLMX(JL))/ZCONS5
CSH
                  IF (QHFL(JL).GT.0) THEN
                     AEBSOIL(JL) = QHFL(JL) * ZDIAGW
                     AESKIN(JL) = AESKIN(JL) - (QHFL(JL) * ZDIAGW)
                  ENDIF
CSH
               ENDIF
               ZWL(JL)=AMIN1(ZWL(JL),WLMX(JL))
CSH
C      AERES(JL) = QHFL(JL) * ZDIAGW
CSH
            ENDIF
         ENDDO
C
C
C          INTERCEPTION OF PRECIPITATION BY THE VEGETATION
C
C          LARGE SCALE PRECIPITATION
C
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
               ZTPRCP=RSFL(JL)
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
               ZIPRCP=ZTPRCP*VGRAT(JL)*ZVINTER
CTS 250100
               ZMPRCP=AMIN1(ZIPRCP,(WLMX(JL)-ZWL(JL))/ZCONS5)
               ZWL(JL)=ZWL(JL)+ZMPRCP*ZCONS5
               RSFL(JL)=RSFL(JL)-ZMPRCP
            ENDIF
         ENDDO
C
C          CONVECTIVE PRECIPITATION
C
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
               ZTPRCP=RSFC(JL)
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
               ZIPRCP=ZTPRCP*VGRAT(JL)*ZVINTER
CTS 250100
               ZMPRCP=AMIN1(ZIPRCP,(WLMX(JL)-ZWL(JL))/ZCONS5)*ZPSFR
               ZWL(JL)=ZWL(JL)+ZMPRCP*ZCONS5
               RSFC(JL)=RSFC(JL)-ZMPRCP
            ENDIF
         ENDDO
C
C
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
CSH   *** OBSOLETE CONFUSING CODE WAS ELIMINATED
               ZPSFL(JL)=SSFL(JL)+SSFC(JL)+CVS(JL)*DHFQS(JL)
               ZROS(JL)=0.
CSH
               ZBETAG(JL)=0.
               ZDRAIN(JL)=0.
            ENDIF
         ENDDO
C
C
C*         3.3     SNOW CHANGES.
C
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
               IF(LOGLAC(JL)) THEN
                  SN(JL)=SNM1M(JL)
               ELSE
                  SN(JL)=SNM1M(JL)+ZCONS5*ZPSFL(JL)
                  SN(JL)=AMAX1(SN(JL),0.)
               ENDIF
            ELSE
               SN(JL)=0.
            ENDIF
         ENDDO
C
C     ------------------------------------------------------------------
C
C*         4.     SNOW MELT AND CORRECTIONS (INCLUDING RUN OFF).
C                 ---- ---- --- ----------- ---------- --- -----
C
C*         4.1     SNOW MELT AND SNOW CORRECTION.
C
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
CTS CHANGE FOR MODIFIED SNOW PARAMETRISATION
               IF((SN(JL).GT.0.).AND.(.NOT.LOGLAC(JL))) THEN
                  ZSNMEL(JL)=0.
                  IF (TD3(JL).GT.TMELT) THEN
                     ZSNMEL(JL)=(TD3(JL)-TMELT)*ZCGN(JL,1)*ZDEL(1)/
     &                    (RHOH2O*ALF)
                     TD3(JL) = TMELT
                     SN(JL) = SN(JL) - ZSNMEL(JL)
                  ENDIF
                  IF (SNM1M(JL).GE.ZSNCRI) THEN
                     IF (TSL(JL).GT.TMELT) THEN
                        ZSNMELU=(TSL(JL)-TMELT)*ZCPCONS(JL)/(RHOH2O*ALF)
                        ZSNMEL(JL)=ZSNMEL(JL)+ZSNMELU
                        TSL(JL) = TMELT
                        SN(JL) = SN(JL) - ZSNMELU
                     ENDIF
CHG+KS SNOW MELTING FROM TOP
                  ELSE
                     ZSNMELU=FSMELT(JL)*ZTMST/(RHOH2O*ALF)
                     ZSNMEL(JL)=ZSNMEL(JL)+ZSNMELU
                     SN(JL) = SN(JL) - ZSNMELU
                  ENDIF
                  IF (SN(JL).LT.0.) THEN
                     ZSNMEL(JL) = ZSNMEL(JL) + SN(JL)
                     TD3(JL) = TD3(JL)-SN(JL)*RHOH2O*ALF/
     &                    (ZCGN(JL,1)*ZDEL(1))
                     TSL(JL) = TD3(JL)
                     SN(JL) = 0.
                  ENDIF
               ENDIF
C
C
C          ALLOW FOR SOME WATER TO REMAIN ON THE GROUND OR ON THE LEAVES
C
               ZSNMLT=AMAX1(0.,ZSNMEL(JL))
               ZWL(JL)=ZWL(JL)+ZVINTER*ZSNMLT
               ZSNMLT=(1.-ZVINTER)*ZSNMLT+AMAX1(0.,ZWL(JL)-WLMX(JL))
               ZWL(JL)=AMIN1(WLMX(JL),ZWL(JL))
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
               IF (.NOT.LOGLAC(JL)) THEN
CTS 250100
                  ZPRFL(JL)=ZSNMLT*RHOH2O/ZTMST
               ELSE
                  ZPRFL(JL)=0.
               ENDIF
C
            ENDIF
         ENDDO
C
C
C*          4.2 COMPUTE SOIL HYDROLOGY
C
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
C
C      ACCOUNT FOR SURFACE RUNOFF DUE TO SLOPING TERRAIN
C
               ZROEFF=AMAX1(0.,SQRT(VAROR(JL))-ZORVARI)
     &                       /(SQRT(VAROR(JL))+ZORVARS)
               ZBWS=AMAX1(AMIN1(ZROEFF,0.5),0.01)
C
CSH  {
               IF (IEXC.EQ.1 .OR. IEXC.EQ.7) THEN
C              *** STANDARD ARNOSCHEME: IEXC = 1
                  ZBETAG(JL) = ZBWS
               ELSE IF (IEXC.EQ.5) THEN
C              *** IMPROVED ARNOSCHEME: IEXC = 5
                  ZBETAG(JL) = BETA(JL)
                  IF (ZBWS.GT.0.01) ZBETAG(JL)=ZBETAG(JL)+ZBWS
               ENDIF
CSH  }
C
               LO=QHFL(JL).GE.0.
               IF (LO) THEN
                  ZPRFL(JL)=ZPRFL(JL)+RSFL(JL)+RSFC(JL)
     &                 +QHFL(JL)
               ELSE
                  ZPRFL(JL)=ZPRFL(JL)+RSFL(JL)+RSFC(JL)
               ENDIF
CSH  {
C     *** UNIT OF ZPRFL IS KG/M^2S: ARNOSCHEME EXPECTS M/S
C     ***   --> DIVIDE BY DENSITY
               IF (IEXC.NE.7) ZPRFL(JL) = ZPRFL(JL) / RHOH2O
            ENDIF
         ENDDO
C
C*    COMPUTE INFILTRATION FROM RAINFALL AND RUNOFF
C     _____________________________________________
C
C     *** AUFRUF STANDARD ARNOSCHEME: IEXC = 1
C     *** AUFRUF IMPROVED ARNOSCHEME: IEXC = 5
C     *** AUFRUF OLD REMO SIMULACRUM ARNOSCHEME: IEXC = 7
         IF (IEXC.EQ.7) THEN
            CALL ASURF
     &          (KLON , KIDIA, KFDIA , WSM1M , WS    , ZROS, ZPRFL,
     &           ZTMST, WSMX , ZBETAG, LOLAND, LOGLAC)

         ELSE
            CALL ASCHEME
     &          (KLON  , KIDIA , KFDIA , WSM1M, WS     , ZROS   ,
     &           ZPRFL , IEXC  , ZTMST , WSMX , WMINLOK, WMAXLOK,
     &           ZBETAG, LOLAND, LOGLAC)
         ENDIF
C
         CALL SETRA(ZINFIL,KLON,0.0)
C
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
C
C              *** CALCULATE INFILTRATION FOR TRANSFER TO ROUTINE SOILCHANGE
C              *** UNIT: [M/S]
               ZINFIL(JL) = ZPRFL(JL) - ZROS(JL)
CSK
               ZINFT      = (1.-WI3(JL)) * ZINFIL(JL)
               ZROS(JL)   = ZROS(JL) + (ZINFIL(JL) - ZINFT)
               ZINFIL(JL) = ZINFT
CSK
C
C              *** UNIT OF ZROS IS M/S WHEN EXITING FROM ARNOSCHEME
C              ***   REMO EXPECTS M --> MULTIPLY BY TIMESTEP ZTMST
               IF (IEXC.NE.7) ZROS(JL) = ZROS(JL) * ZTMST
CSH  }
C
C*             SUBTRACT EVAPORATION
C              ____________________
C
               IF (QHFL(JL).LT.0.) THEN
                  WS(JL)=WS(JL)+QHFL(JL)*ZTMST/RHOH2O
               END IF
C
               WS(JL)=WS(JL)+XHFL(JL)*ZCONS5
CSH {
            ENDIF
         ENDDO
C
C     ******* DRAINAGE IF THE OLD BUCKET IS USED
         IF (I5LAYER.EQ.0) THEN
            DO JL=KIDIA, KFDIA
               IF (LOLAND(JL)) THEN
CSH }
C
C       *** TERMUMBELEGUNG AUS ANFANG DES ORIGINAL SOIL HYDROLOGY LOOPS
                  ZWMAX=WSMX(JL)
                  ZWDTR=0.90*ZWMAX
                  ZCONW2=ZWMAX-ZWDTR
                  ZCONW3=ZDRMAX-ZDRMIN
C
                  IF (WS(JL).LT.0.) THEN
                     WS(JL)=0.
                  ELSE
C
C*                SOIL WATER DEPLETION
C                 ____________________
C
                     IF (.NOT.LOGLAC(JL)) THEN !HGKS.AND.TD3(JL).GT.TMELT) THEN
C
                        ZDRAIN(JL)=ZDRMIN*WS(JL)/ZWMAX
                        IF (WS(JL).GT.ZWDTR) THEN
                           ZDRAIN(JL)=ZDRAIN(JL)+ZCONW3*
     &                          ((WS(JL)-ZWDTR)/ZCONW2)**ZDREXP
                        END IF
                        ZDRAIN(JL)=ZDRAIN(JL)*ZTMST
C
                        ZWSLIM=0.05*WSMX(JL)
                        ZDRAIN(JL)=AMIN1(ZDRAIN(JL),
     &                                   AMAX1(WS(JL)-ZWSLIM,0.))
                        WS(JL)=WS(JL)-ZDRAIN(JL)
C
                     ELSE
                        ZDRAIN(JL)=0.
                     ENDIF
C
                  ENDIF
C
                  ZWSUP=AMAX1(WS(JL)-WSMX(JL),0.)
                  WS(JL)=WS(JL)-ZWSUP
                  ZROS(JL)=ZROS(JL)+ZDRAIN(JL)+ZWSUP
C
               ENDIF
            ENDDO
C
C     *** NOT USED IN I5LAYER = 0
            CALL SETRA(AERES,KLON,0.0)
C
CSH { ******* CHANGE OF SOIL MOISTURE IN 5 LAYERS
C     ******* CHANGE BELONGS TO ONE TIME STEP = 0.5 * TWODT
         ELSE IF (I5LAYER.GE.1) THEN
C
C       *** INITILIZATION WSI
            DO JK=1,KDEEP
               DO JL=1,KLON
                  WSI(JL,JK)=WSIM1M(JL,JK)
               ENDDO
            ENDDO
C
            ISCH = 1            ! CHANGES BY EVAPORATION FLUXES
            CALL SOILCHANGE(KIDIA , KFDIA  , KLON   , KDEEP,
     &           ISCH   , ZDIAGT, ZDIAGS , LOLAND ,
     &           WSI    , WSIM1M,
     &           DZR    , DZRSI , DZSI   , ZWSAT  , ZWSFC,
     &           ZINFIL , ZDRAIN, AETRANS, AEBSOIL,
     &           REDEVAP)

            ISCH = 2            ! CHANGES BY INFILTRATION AND DRAINAGE
            CALL SOILCHANGE(KIDIA , KFDIA  , KLON   , KDEEP,
     &           ISCH   , ZDIAGT, ZDIAGS , LOLAND ,
     &           WSI    , WSIM1M,
     &           DZR    , DZRSI , DZSI   , ZWSAT  , ZWSFC,
     &           ZINFIL , ZDRAIN, AETRANS, AEBSOIL,
     &           REDEVAP)
C
C           *** DRAINAGE CALCULATION WITH 5 LAYER SOIL SCHEME
C
C
            CALL SOILHYD(KIDIA , KFDIA , KLON  , KDEEP ,
     &           ZDIAGT, ILOG  , JLLOG ,
     &           WSI   , DZSI  , ZWSFC ,
     &           FKSAT , VPOR  , BCLAPP,
     &           FMPOT , ZDRAIN)
C
            DO JL=1,KLON
               AERES(JL)=REDEVAP(JL)
            ENDDO
C
C           ***  WS MUST BE RE-DERIVED FROM WSI AFTER SOILHYD
            CALL SETRA(WS,KLON,0.0)
            DO JK=1, KDEEP
               DO JL=KIDIA, KFDIA
                  IF (LOLAND(JL)) THEN
                     IF (DZRSI(JL,JK).GE.ZDEL(JK)) THEN
                        WS(JL) = WS(JL) + WSI(JL,JK)
                     ELSE IF (DZRSI(JL,JK).GT.0.) THEN
                        WS(JL) = WS(JL) + 
     &                       WSI(JL,JK) * DZRSI(JL,JK)/DZSI(JL,JK)
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ELSE
C        NO CHANGES
            DO JK=1,KDEEP
               DO JL=1,KLON
                  WSI(JL,JK)=WSIM1M(JL,JK)
               ENDDO
            ENDDO
         ENDIF
CSH }
C
C*        4.3 ACCUMULATE SURFACE PARAMETERS FOR DIAGNOSTICS
C
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
               TSURF(JL)=TSURF(JL)+ZDIAGT*TSL(JL)
CSH {
               IF (I5LAYER.GE.1) THEN
C
C              *** TOTAL RUNOFF = RUNOFF + DRAINAGE
C              *** ZROS WURDE FUER ZTMST (TWODT) AUSGERECHNET, WAEHREND ZDRAIN
C              *** FUER ZDIAGT (=DT) AUSGERECHNET WURDE.
C
                  RUNOFF(JL)=RUNOFF(JL)+ZDIAGS*ZROS(JL)+ZDRAIN(JL)
                  DRAIN(JL)=DRAIN(JL) + ZDRAIN(JL)
               ELSE
C
C              *** ZROS UND ZDRAIN WURDEN FUER ZTMST (TWODT) AUSGERECHNET,
C              *** ZDRAIN IST SCHON IN ZROS ENTHALTEN.
                  RUNOFF(JL)=RUNOFF(JL)+ZDIAGS*ZROS(JL)
                  DRAIN(JL)=DRAIN(JL)+ZDIAGS*ZDRAIN(JL)
               ENDIF
CSH }
               SNMEL(JL)=SNMEL(JL)+ZDIAGS*ZSNMEL(JL)
               ZSN=(SN(JL)-SNM1M(JL))
               DSNAC(JL)=DSNAC(JL)+ZDIAGS*ZSN+ZDIAGW*ZLAC(JL)
C
C
            ELSE
               WS(JL)=WSM1M(JL)
               ZWL(JL)=0.
               TSN(JL)=TSNM1M(JL)
               TSL(JL)=TD3M1M(JL)
               TD3(JL)=TD3M1M(JL)
               TD4(JL)=TD4M1M(JL)
               TD5(JL)=TD5M1M(JL)
               TD(JL)=TDM1M(JL)
               TDCL(JL)=TDCLM1M(JL)
            ENDIF
         ENDDO
C
C
C     ------------------------------------------------------------------
C
C*         6.     NECESSARY COMPUTATIONS IF SUBROUTINE IS BY-PASSED.
C                 --------- ------------ -- ---------- -- ----------
C
      ELSE
C***
         DO JL=KIDIA, KFDIA
            TD(JL)=TDM1M(JL)
            ZWL(JL)=WLM1M(JL)
            WS(JL)=WSM1M(JL)
            SN(JL)=SNM1M(JL)
            TSN(JL)=TSNM1M(JL)
            TD3(JL)=TD3M1M(JL)
            TD4(JL)=TD4M1M(JL)
            TD5(JL)=TD5M1M(JL)
            TDCL(JL)=TDCLM1M(JL)
C
            REDEVAP(JL)=0.
            AETRANS(JL)=0.
            AEBSOIL(JL)=0.
            AESNOW(JL)=0.
            AESKIN(JL)=0.
            AERES(JL)=0.
C
         ENDDO
CSH
         ISCH = 0               ! NO CHANGES
         CALL SOILCHANGE(KIDIA , KFDIA  , KLON   , KDEEP,
     &        ISCH   , ZDIAGT, ZDIAGS , LOLAND ,
     &        WSI    , WSIM1M,
     &        DZR    , DZRSI , DZSI   , ZWSAT  , ZWSFC,
     &        ZINFIL , ZDRAIN, AETRANS, AEBSOIL,
     &        REDEVAP)
C***
      ENDIF
C***
C
C     ------------------------------------------------------------------
C
C*         7.     RETURN WORKSPACE.
C                 ------ ----------
C
      DO JL=KIDIA,KFDIA
         WL(JL)=ZWL(JL)
      ENDDO
C
      RETURN
      END SUBROUTINE SURF
