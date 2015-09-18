C
C
C**** *RADLSW* - RADIATION INTERFACE
C
C     PURPOSE.
C     --------
C           CONTROLS RADIATION COMPUTATIONS
C
C**   INTERFACE.
C     ----------
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C PMU0   : (KLON)              ; SOLAR ANGLE
C PALBS  : (KLON)              ; ALBEDO IN THE TWO INTERVALS
C                              ;       .25-.68 AND .68-4.
C******************************************************** ONLY 1]
C PCCO2  :                      ; CONCENTRATION IN CO2 (PA/PA)
C POZON  : (KLON ,KLEV)       ; CONCENTRATION IN OZONE (PA/PA)
C PTH    : (KLON ,KLEV+1)     ; HALF LEVEL TEMPERATURE
C PAPH   : (KLON  ,KLEV+1)     ; HALF LEVEL PRESSURE
C PSURP  : (KLON)               ; SURFACE PRESSURE
C PQ     : (KLON ,KLEV)       ; SPECIFIC HUMIDITY PA/PA
C PQW    : (KLON ,KLEV)       ; LIQUID WATER KG/KG
C PQS    : (KLON ,KLEV)       ; SATURATION SPECIFIC HUMIDITY (KG/KG)
C PCLFR  : (KLON ,KLEV)       ; CLOUD FRACTIONAL COVER
C PAER   : (KLON,KLEV,5)    ; AEROSOL AMOUNT IN PART/CM2 (NAER=11,12)
C                                AEROSOL OPTICAL THICKNESS  (NAER=5)
C     ==== OUTPUTS ===
C PFLT   : (KLON ,KLEV+1)     ; LONG WAVE FLUXES
C PFLS   : (KLON ,KLEV+1)     ; SHORT WAVE FLUXES
C PFLTC  : (KLON,2)            ; CLEAR SKY LW FLUXES (TOP, BOTTOM)
C PFLSC  : (KLON,2)            ; CLEAR SKY SW FLUXES (TOP, BOTTOM)
C                                 ONLY IF LSOLC=.TRUE.
C
C        IMPLICIT ARGUMENTS :   NONE
C        --------------------
C
C     METHOD.
C     -------
C        SEE DOCUMENTATION
C
C     EXTERNALS.
C     ----------
C
C
C     AUTHORS.
C     --------
C        J.-J. MORCRETTE         *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 88-02-04
C      REWRITTEN   U. SCHLESE  DKRZ - HAMBURG   JUL 93
C-----------------------------------------------------------------------
      SUBROUTINE RADLSW
     &    (KLON , KLEV  , KSTART, KSTOP, KMODE , KAER , KCFC,
     &     ZSCT , LALAND, LAGLAC, PAER , PALBS , PCLFR, PMU0,
     &     POZON, PSURP , PDP   , PQW  , PQ    , PQS  , PAPH,
     &     PT   , PTH   , PFLS  , PFLT , PFLSC , PFLTC, PQI)
C
      IMPLICIT NONE
C
      INCLUDE "COMCON"
      INCLUDE "YOMRDU"
      INCLUDE "YOMRDI"
      INCLUDE "YOMCFC"
C-----------------------------------------------------------------------
C
C     -----------------------------------------------------------------
C
C*       0.1   ARGUMENTS.
C              ----------
      INTEGER, INTENT(IN)   :: KLON , KLEV  , KSTART, KSTOP, 
     &                         KMODE , KAER , KCFC
C
      REAL,    INTENT(IN)   :: ZSCT
C
      REAL,    INTENT(IN)   ::
     &     PALBS(KLON)   , PMU0(KLON)      , POZON(KLON,KLEV) ,
     &     PT(KLON,KLEV) , PTH(KLON,KLEV+1), PAPH(KLON,KLEV+1),
     &     PSURP(KLON)   , PDP(KLON,KLEV)  , PQ(KLON,KLEV)    ,
     &     PQW(KLON,KLEV), PQS(KLON ,KLEV) , 
     &     PAER(KLON,KLEV,5)
      REAL,    INTENT(INOUT) :: 
     &     PCLFR(KLON,KLEV)
C      === COMPUTED IN RADLSW ===
      REAL,    INTENT(INOUT) :: 
     &     PFLS(KLON,KLEV+1), PFLT(KLON,KLEV+1),
     &     PFLSC(KLON,2)    , PFLTC(KLON,2),
     &     PQI(KLON,KLEV)
C
      LOGICAL, INTENT(IN) :: LALAND(KLON), LAGLAC(KLON)
C     -----------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS.
C              -------------
C     -----------------------------------------------------------------
C
C-- ARRAYS FOR LOCAL VARIABLES -----------------------------------------
C
      REAL ZALBSU(KLON,2)      , ZCG(KLON,2,KLEV)     ,
     &     ZCLDLW(KLON,KLEV)   , ZCLDSW(KLON,KLEV)    ,
     &     ZFLUX(KLON,2,KLEV+1), ZFLUXC(KLON,2,KLEV+1),
     &     ZFLWP(KLON)         , ZOMEGA(KLON,2,KLEV+1),
     &     ZPMB(KLON,KLEV+1)   , ZPSOL(KLON)          ,
     &     ZTAU (KLON,2,KLEV)  , ZTAVE (KLON,KLEV)    ,
     &     ZTL(KLON,KLEV+1)    , ZFSDWN(KLON,KLEV+1)  ,
     &     ZFSUP(KLON,KLEV+1)  ,
     &     ZRADIP(KLON)        , ZRADLP(KLON)         ,
     &     ZN1(KLON)           , ZN2(KLON)            ,
     &     ZKAP(KLON)          , ZCFCABS(4)           ,
     &     ZIWP(KLON)
C
      LOGICAL :: LO2
      INTEGER :: IKL, JCFC, JK, JKL, JKLP1, JL, JM, NEXP
      REAL    :: TLAB, ZASIC, ZASYMX1, ZASYMX2, ZCCO2, ZCDNC, ZEFFIR, 
     &           ZFCFC, ZFPI, ZG1I, ZG1L, ZG2I, ZG2L, ZHEY, ZHPBL, ZICE, 
     &           ZICEWP, ZIP1, ZIP2, ZEFFLR
      REAL    :: ZIP3, ZIP4, ZIWC, ZIWGKG, ZLIQWP, ZLOGI, ZLOGP2, 
     &           ZLP1, ZLP2, ZLP3, ZLWC, ZLWGKG, ZMACI, ZMACL, ZMEANR, 
     &           ZO1L, ZO2I, ZO2L, ZO1I, ZLOGP3
      REAL    :: ZOMGMX1, ZOMGMX2, ZPMBM, ZPRAT, ZREF, ZREX, ZRH2O, 
     &           ZTAUMX1, ZTAUMX2, ZTO1I, ZTO1L, ZTO2I, ZTO2L, ZWATER,
     &           ZRI0
C
!DIR$ NOBOUNDS ZCLDSW
C     ------------------------------------------------------------------
C
C*       0.3   FUNCTIONS
C              ---------
C
C
C
C
C      0.0  SET UP LOCAL CONSTANTS
C      ---------------------------
C
      ZFPI=1.2
      ZASIC=0.91
      ZHEY=1.
      ZRH2O=1000.
      ZHPBL=800.
      NEXP=2
      ZREX=1./3.
      ZREF=1.E6*(3.E-9/(4.*API*ZRH2O))**ZREX
C     MEASUREMENT TEMPERATURE CFC'S (296 K)
      TLAB=296.
C
C     ------------------------------------------------------------------
C
C*         1.     SET-UP INPUT QUANTITIES FOR RADIATION
C                 -------------------------------------
C
      ZCCO2 =ZCARDI
C
C
      DO JL = KSTART , KSTOP
         ZFLUX(JL,1,KLEV+1) = 0.
         ZFLUX(JL,2,KLEV+1) = 0.
         ZFLUXC(JL,1,KLEV+1) = 0.
         ZFLUXC(JL,2,KLEV+1) = 0.
         ZPSOL(JL) = PSURP(JL)
         ZPMB(JL,1) = PSURP(JL) / 100.
      ENDDO
C
C*         1.1    INITIALIZE VARIOUS FIELDS
C                 -------------------------
C
      DO JL = KSTART , KSTOP
         ZALBSU(JL,1)=PALBS(JL)
         ZALBSU(JL,2)=PALBS(JL)
         IF(LALAND(JL) .AND. (.NOT. LAGLAC(JL))) THEN
            ZN1(JL)=50.
            ZN2(JL)=220.
            ZKAP(JL)=1.143
         ELSE
            ZN1(JL)=50.
            ZN2(JL)=100.
            ZKAP(JL)=1.077
         ENDIF
      ENDDO
C
      DO JK = 1 , KLEV
         JKL = KLEV+ 1 - JK
         JKLP1 = JKL + 1
         DO JL = KSTART , KSTOP
            ZPMB(JL,JK+1)=PAPH(JL,JKL)/100.
            ZFLUX(JL,1,JK) = 0.
            ZFLUX(JL,2,JK) = 0.
            ZFLUXC(JL,1,JK) = 0.
            ZFLUXC(JL,2,JK) = 0.
         ENDDO
      ENDDO
C
C
      DO JK=1,KLEV
         JKL=KLEV+1-JK
         JKLP1=JKL+1
         DO JL=1,KLON
            ZTL(JL,JK)=PTH(JL,JKLP1)
            ZTAVE(JL,JK)=PT(JL,JKL)
         ENDDO
      ENDDO
      DO JL=1,KLON
         ZTL(JL,KLEV+1)= PTH(JL,1)
      ENDDO
C***
C
C     ------------------------------------------------------------------
C
C*         2.     CLOUD, AEROSOL AND CFC PARAMETERS
C                 ---------------------------------
C
      DO JK = 1 , KLEV

         IKL = KLEV + 1 - JK
C
C --- INITIALIZE OPTICAL PROPERTIES TO CLEAR SKY VALUES
         DO JL = KSTART , KSTOP
            ZCLDSW(JL,JK)  = 0.
            ZTAU(JL,1,JK)  = 0.
            ZTAU(JL,2,JK)  = 0.
            ZOMEGA(JL,1,JK)= 1.
            ZOMEGA(JL,2,JK)= 1.
            ZCG(JL,1,JK)   = 0.
            ZCG(JL,2,JK)   = 0.
            ZCLDLW(JL,JK)  = 0.
         ENDDO
C
         DO JL = KSTART , KSTOP
            ZPMBM=0.5*(ZPMB(JL,JK)+ZPMB(JL,JK+1))
            PCLFR(JL,IKL)=AMAX1(ZEPSC,AMIN1(PCLFR(JL,IKL),1.-ZEPSC))
            LO2=PCLFR(JL,IKL).GT.(2.*ZEPCLC)
C
C --- LIQUID WATER CONTENT (G.M-3) AND LIQUID WATER PATH (G.M-2)
            ZLWGKG=PQW(JL,IKL)*1000.
            ZIWGKG=PQI(JL,IKL)*1000.
            IF (LO2) THEN
               ZLWGKG=ZLWGKG/PCLFR(JL,IKL)
               ZIWGKG=ZIWGKG/PCLFR(JL,IKL)
            ELSE
               ZLWGKG=0.
               ZIWGKG=0.
            ENDIF
            ZFLWP(JL)= ZLWGKG*PDP(JL,IKL)/G
            ZLWC=ZLWGKG*100.*ZPMBM/(RD*ZTAVE(JL,JK))
            ZIWC=ZIWGKG*100.*ZPMBM/(RD*ZTAVE(JL,JK))
            ZIWP(JL)=ZIWGKG*PDP(JL,IKL)/G

C --- PARTITIONING BETWEEN LIQUID AND ICE WATER DROPLETS
CSP    VOR DER BERECHNUNG DER EFFEKTIVRADII ANSTATT DANACH.

C      ZTC=ZTAVE(JL,JK)-ZTMELT
C      IF(ZTAVE(JL,JK).LT.ZTMELT) THEN
C       ZFICE(JL)=1.-(ZCAA+(1.-ZCAA)*EXP(-ZCAB*ZTC**2))
C      ELSE
C       ZFICE(JL)=0.
C      END IF

CSP > ZICE, ZWATER: ENTSPRICHT DEM EIS- UND DEM FLUSSIFGEB ANTEIL DES
CSP   GESAMTWASSERGEHALTES ZLWC. GEHT EIN IN DIE BERECHNUNG DER EFFEKTIV-
CSP   RADII FUER EISKRISTALLE BZW. WASSERTROPFEN.

C      ZICE=ZLWC*ZFICE(JL)
C      ZWATER=ZLWC*(1.-ZFICE(JL))
            ZICE=ZIWC
            ZWATER=ZLWC
CSP

C
C --- EFFECTIVE RADIUS FOR WATER AND ICE PARTICLES (MICROMETER)
C     ACCORDING TO HEYMSFIELD (1977) AND MC FARLANE ET AL (1992)
C
            ZLOGI=ALOG10(AMAX1(1.E-4,ZHEY*ZICE))
            ZLOGP2=ZLOGI*ZLOGI
            ZLOGP3=ZLOGP2*ZLOGI
            IF(ZPMBM.LT.ZHPBL) THEN
               ZPRAT=(AMIN1(2.,ZHPBL/ZPMBM))**NEXP
               ZCDNC=ZN1(JL)+(ZN2(JL)-ZN1(JL))*(EXP(1.-ZPRAT))
            ELSE
               ZCDNC=ZN2(JL)
            ENDIF
            ZMEANR=0.001*(0.698+0.366*ZLOGI+0.122*ZLOGP2+0.0136*ZLOGP3)
            ZEFFIR=ZFPI*5640.*ZMEANR**0.786
            ZEFFLR=ZREF*ZKAP(JL)*(ZWATER/ZCDNC)**ZREX
            ZRADIP(JL)=AMIN1(80.,ZEFFIR)
            ZRADLP(JL)=AMAX1(4.,AMIN1(24.,ZEFFLR))
C
         ENDDO
C
C
C             SW OPTICAL PARAMETERS
C             -- ------- ----------
C         LIQUID WATER AND ICE (ROCKEL, 1993)
C         ------ ----- --- --- --------------
C
         DO JL = 1, KLON
            IF (ZFLWP(JL).NE.0.) THEN
               ZICEWP = ZIWP(JL)
               ZLIQWP = ZFLWP(JL)
               ZLP1=ALOG10(ZRADLP(JL))
               ZLP2=ZLP1*ZLP1
               ZLP3=ZLP2*ZLP1
               ZTO1L=ZLIQWP*1.8706058417*ZRADLP(JL)**(-1.0756364457)
               ZTO2L=ZLIQWP*1.9655460426*ZRADLP(JL)**(-1.0778999732)
               ZG1L=0.78756640717+0.10660598895*ZLP1-0.031012468401*ZLP2
               ZG2L=0.79208639802-0.044930076174*ZLP1+0.18980672305*ZLP2
     &              -0.082590933352*ZLP3
               ZO1L=0.99999999
               ZO2L=0.9854369057+0.013584242533*ZLP1-0.024856960461*ZLP2
     &              +0.0055918314369*ZLP3
C
               ZIP1=ALOG10(ZRADIP(JL))
               ZIP2=ZIP1*ZIP1
               ZIP3=ZIP2*ZIP1
               ZIP4=ZIP2*ZIP2
               ZTO1I=ZICEWP*1.9056067246*ZRADIP(JL)**(-1.0318784654)
               ZTO2I=ZICEWP*2.1666771102*ZRADIP(JL)**(-1.0634702711)
               ZG1I=0.7700034985+0.19598466851*ZIP1-0.11836420885*ZIP2
     &              +0.025209205131*ZIP3
               ZG2I=0.83631171237-0.19965998649*ZIP1+0.46130320487*ZIP2
     &              -0.29719270332*ZIP3+0.062554483594*ZIP4
               ZG1I=ZG1I*ZASIC
               ZG2I=ZG2I*ZASIC
               ZO1I=0.99999999
               ZO2I=0.98475089485+0.0053152066002*ZIP1
     &              -0.0061150583857*ZIP2-0.0032775655896*ZIP3
C
C  - MIX OF WATER AND ICE CLOUDS
               ZTAUMX1= ZTO1L + ZTO1I
               ZTAUMX2= ZTO2L + ZTO2I
               ZOMGMX1= ZTO1L*ZO1L + ZTO1I*ZO1I
               ZOMGMX2= ZTO2L*ZO2L + ZTO2I*ZO2I
               ZASYMX1= ZTO1L*ZO1L*ZG1L + ZTO1I*ZO1I*ZG1I
               ZASYMX2= ZTO2L*ZO2L*ZG2L + ZTO2I*ZO2I*ZG2I
C
               ZASYMX1= ZASYMX1/ZOMGMX1
               ZASYMX2= ZASYMX2/ZOMGMX2
               ZOMGMX1= ZOMGMX1/ZTAUMX1
               ZOMGMX2= ZOMGMX2/ZTAUMX2
C
C --- SW FINAL CLOUD OPTICAL PARAMETERS
C
               ZCLDSW(JL,JK)  = PCLFR(JL,IKL)
               ZTAU(JL,1,JK)  = ZTAUMX1
               ZTAU(JL,2,JK)  = ZTAUMX2
               ZOMEGA(JL,1,JK)= ZOMGMX1
               ZOMEGA(JL,2,JK)= ZOMGMX2
               ZCG(JL,1,JK)   = ZASYMX1
               ZCG(JL,2,JK)   = ZASYMX2
C
C
C             LW OPTICAL PARAMETERS
C             -- ------- ----------
C         LIQUID WATER AND ICE (ROCKEL, 1993)
C         ------ ----- --- --- --------------
C
               ZMACL=0.025520637+0.2854650784
     &              *EXP  (-0.088968393014*ZRADLP(JL))
               ZMACI=0.020219423+0.2058619832
     &              *EXP  (-0.067631070625*ZRADIP(JL))
               ZCLDLW(JL,JK)=PCLFR(JL,IKL)*(1.-EXP  (-ZMACL*ZLIQWP)*
     &                                         EXP  (-ZMACI*ZICEWP))
            ENDIF
         ENDDO
C
      ENDDO
C
      DO JL = KSTART , KSTOP
         ZPMB (JL,KLEV+1) = 0.0
      ENDDO
C
      DO JM=1,4
         ZCFCABS(JM)=0.
      ENDDO
C
C        CFC'S
C
      IF(KCFC.NE.0)THEN
         ZFCFC=RD*TLAB*1.E+03/101325.
         DO JM=1,4
            DO JCFC=1,NCFC
               ZCFCABS(JM)=ZCFCABS(JM)+ZCFC(JCFC)*CFCWMO(JM,JCFC)
            ENDDO
            ZCFCABS(JM)=ZCFCABS(JM)*ZFCFC
         ENDDO
      ENDIF
C
      ZRI0=ZSCT
C
C
C     ------------------------------------------------------------------
C
C*         3.     CALL LONGWAVE RADIATION CODE
C                 ----------------------------
C
      CALL LW(KLON, KLEV , KAER  , KCFC   , ZCCO2, ZCLDLW,
     &        PDP , ZPMB , POZON , ZTL    , PAER , ZTAVE ,
     &        PQ  , ZFLUX, ZFLUXC, ZCFCABS)
C
C
C     ------------------------------------------------------------------
C
C*         4.     CALL SHORTWAVE RADIATION CODE
C                 -----------------------------
C
      CALL SW(KLON , KLEV, ZRI0, ZCCO2, ZPSOL , ZALBSU,
     &        PQ   , PQS , PMU0, ZCG  , ZCLDSW, ZOMEGA,
     &        POZON, ZPMB, ZTAU, ZTAVE, PAER  , ZFSDWN,
     &        ZFSUP)
C
C     ------------------------------------------------------------------
C
C*         5.     FILL UP THE MODEL NET LW AND SW RADIATIVE FLUXES
C                 ------------------------------------------------
C
      DO JKL = 1 , KLEV+1
         JK = KLEV+1 + 1 - JKL
         DO JL = KSTART , KSTOP
            PFLS(JL,JKL) = ZFSDWN(JL,JK) - ZFSUP(JL,JK)
            PFLT(JL,JKL) = - ZFLUX(JL,1,JK) - ZFLUX(JL,2,JK)
         ENDDO
      ENDDO
C
C
C     LONG-WAVE CLEAR SKY FLUXES
C
      DO JL = KSTART , KSTOP
         PFLTC(JL,1)=-ZFLUXC(JL,1,KLEV+1)-ZFLUXC(JL,2,KLEV+1)
         PFLTC(JL,2)=-ZFLUXC(JL,1,1)-ZFLUXC(JL,2,1)
      ENDDO
C
C ---------------------------------------------------------------------
C
C*         6.  SHORT WAVE CLEAR SKY FLUXES
C              ----- ---- ----- --- ------
C
      IF (KMODE.EQ.1) THEN
         DO JK = 1 , KLEV
            DO JL = KSTART , KSTOP
               ZCLDSW(JL,JK)=0.
            ENDDO
         ENDDO
C
         CALL SW(KLON , KLEV, ZRI0, ZCCO2, ZPSOL , ZALBSU,
     &           PQ   , PQS , PMU0, ZCG  , ZCLDSW, ZOMEGA,
     &           POZON, ZPMB, ZTAU, ZTAVE, PAER  , ZFSDWN,
     &           ZFSUP)
C
         DO JL = KSTART , KSTOP
            PFLSC(JL,1)=ZFSDWN(JL,KLEV+1)-ZFSUP(JL,KLEV+1)
            PFLSC(JL,2)=ZFSDWN(JL,1)-ZFSUP(JL,1)
         ENDDO
      ELSE
         DO JL = KSTART , KSTOP
            PFLSC(JL,1)=0.
            PFLSC(JL,2)=0.
         ENDDO
      ENDIF
C
C  -----------------------------------------------------------
C
      RETURN
      END SUBROUTINE RADLSW
