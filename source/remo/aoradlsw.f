      SUBROUTINE AORADLSW
     &   (KLON , KLEV  , KSTART, KSTOP, KMODE , KAER , KCFC,
     &    KEWAER, KAERH,
     &    ZSCT , LALAND, LAGLAC, PAER , PALBS , PCLFR, PMU0,
     &    POZON, PSURP , PDP   , PQW  , PQ    , PQS  , PAPH,
     &    PT   , PTH   , PFLS  , PFLT , PFLSC , PFLTC,
     &    PSO4ALL, PSO4NAT)
!-----------------------------------------------------------------------
      IMPLICIT NONE
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
C PAER   : (KDLON,KFLEV,5+KEWAER); AEROSOL OPTICAL THICKNESS (1,..,5)
C                                  TANRE ET AL., 1984
C                                  AEROSOL MASS MIXING RATIO (KG/KG)
C                                  (6,...,5+KEWAER) COMPUTED IN ECHAM
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
C       MODIFIED: ROB VAN DORLAND, KNMI, 95-05-10
C-----------------------------------------------------------------------
C      IMPLICIT LOGICAL (L)
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
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN) :: KAER, KCFC, KEWAER, KLEV, KLON, KMODE
      INTEGER, INTENT(IN) :: KSTART, KSTOP, KAERH(KLON,KLEV)
      LOGICAL, INTENT(IN) :: LALAND(KLON), LAGLAC(KLON)
      REAL,    INTENT(IN) :: 
     &     PALBS(KLON), PMU0(KLON), POZON(KLON,KLEV),
     &     PT(KLON,KLEV) , PTH(KLON,KLEV+1), PAPH(KLON,KLEV+1),
     &     PSURP(KLON)   , PDP(KLON,KLEV)  , PQ(KLON,KLEV)    ,
     &     PQW(KLON,KLEV), PQS(KLON ,KLEV) , 
     &     PAER(KLON,KLEV,5+KEWAER),
     &     PSO4ALL(KLON,KLEV)              , PSO4NAT(KLON,KLEV),
     &     ZSCT
C
      REAL, INTENT(INOUT) :: PCLFR(KLON,KLEV)
C
C      === COMPUTED IN RADLSW ===
      REAL, INTENT(INOUT) :: PFLS(KLON,KLEV+1), PFLT(KLON,KLEV+1),
     &                       PFLSC(KLON,2)    , PFLTC(KLON,2)
C     -----------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS.
C              -------------
C     -----------------------------------------------------------------
C     Local Variables
C    
      REAL :: TLAB, ZASIC, ZASYMX1, ZASYMX2, ZCAA, ZCAB, ZCCO2, ZCDNA,
     &        ZCDNC, ZCDNLMN, ZCDNLMX, ZCDNN, ZCDNOMN, ZCDNOMX, ZEFFIR,
     &        ZEFFLR, ZFCFC, ZFPI, ZG1I, ZG1L
      REAL :: ZG2I, ZG2L, ZHEY, ZHPBL, ZICE, ZICEWP, ZIP1, ZIP2, ZIP3,
     &        ZIP4, ZLIQWP, ZLOGI, ZLOGP2, ZLOGP3, ZLP1, ZLP2, ZLP3,
     &        ZLWC, ZLWGKG, ZMACI
      REAL :: ZMACL, ZMEANR, ZO1I, ZO1L, ZO2I, ZO2L, ZOMGMX1, ZOMGMX2,
     &        ZPMBM, ZPRAT, ZRADMAX, ZRADMIN, ZREF, ZREX, ZRH2O, ZRI0,
     &        ZSO4ALLGM3, ZSO4NATGM3, ZTAUMX1
      REAL :: ZTAUMX2, ZTC, ZTMELT, ZTO1I, ZTO1L, ZTO2I, ZTO2L, ZWATER,
     &        ZZG1L, ZZG2L, ZZLP1, ZZLP2, ZZLP3, ZZO2L, ZZTO1L, ZZTO2L,
     &        ZZZG1L, ZZZG2L, ZZZLP1, ZZZLP2
      REAL :: ZZZLP3, ZZZO2L, ZZZTO1L, ZZZTO2L
C
      INTEGER :: IKL, JCFC, JK, JKL, JKLP1, JL, JM, NEXP
      LOGICAL :: LO2
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
     &     ZFSUP(KLON,KLEV+1)  , ZFICE(KLON)          ,
     &     ZRADIP(KLON)        , ZRADLP(KLON)         ,
     &     ZN1(KLON)           , ZN2(KLON)            ,
     &     ZKAP(KLON)          , ZCFCABS(4)           ,
     &     ZRADALL(KLON)       , ZRADNAT(KLON)
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
      ZCAA=0.0059
      ZCAB=0.003102
      ZTMELT=273.16
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
      ZRADMAX=24.
      ZRADMIN=4.
      ZCDNOMN=10.
      ZCDNOMX=500.
      ZCDNLMN=50.
      ZCDNLMX=2000.
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
         END IF
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
CRP   CVMGT ERSETZT
C     ZLWGKG= CVMGT( ZLWGKG/PCLFR(JL,IKL), 0. , LO2 )
            IF (LO2) THEN
               ZLWGKG=ZLWGKG/PCLFR(JL,IKL)
            ELSE
               ZLWGKG=0.
            ENDIF
            ZFLWP(JL)= ZLWGKG*PDP(JL,IKL)/G
            ZLWC=ZLWGKG*100.*ZPMBM/(RD*ZTAVE(JL,JK))

C --- PARTITIONING BETWEEN LIQUID AND ICE WATER DROPLETS
CSP    VOR DER BERECHNUNG DER EFFEKTIVRADII ANSTATT DANACH.

            ZTC=ZTAVE(JL,JK)-ZTMELT
            IF(ZTAVE(JL,JK).LT.ZTMELT) THEN
               ZFICE(JL)=1.-(ZCAA+(1.-ZCAA)*EXP(-ZCAB*ZTC**2))
            ELSE
               ZFICE(JL)=0.
            END IF

CSP > ZICE, ZWATER: ENTSPRICHT DEM EIS- UND DEM FLUSSIFGEB ANTEIL DES
CSP   GESAMTWASSERGEHALTES ZLWC. GEHT EIN IN DIE BERECHNUNG DER EFFEKTIV-
CSP   RADII FUER EISKRISTALLE BZW. WASSERTROPFEN.

            ZICE=ZLWC*ZFICE(JL)
            ZWATER=ZLWC*(1.-ZFICE(JL))
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
            END IF
            ZMEANR=0.001*(0.698+0.366*ZLOGI+0.122*ZLOGP2+0.0136*ZLOGP3)
            ZEFFIR=ZFPI*5640.*ZMEANR**0.786
            ZEFFLR=ZREF*ZKAP(JL)*(ZWATER/ZCDNC)**ZREX
            ZRADIP(JL)=AMIN1(80.,ZEFFIR)
            ZRADLP(JL)=AMAX1(4.,AMIN1(24.,ZEFFLR))
C
C
C      EFFECTIVE DROPLET RADIUS AS A FUNCTION OF SULFATE CONTENT
C      NAT: BASED ON NATURAL SULFUR SOURCES ONLY
C      ALL: BASED ON ALL (NATURAL + ANTHROPOGENIC) S-SOURCES
C
C
C     DEVIDE BY 3 BECAUSE "SULPHUR" IS ALREADY SULPHATE
C
            ZSO4NATGM3=MAX(PSO4NAT(JL,IKL),1.E-18)*3.55E9/3.
            ZSO4ALLGM3=MAX(PSO4ALL(JL,IKL),1.E-18)*3.55E9/3.
C
            IF (LALAND(JL)) THEN
               ZCDNN=174.*ZSO4NATGM3**0.257
               ZCDNA=174.*ZSO4ALLGM3**0.257
               ZCDNN=AMIN1(MAX(ZCDNN,ZCDNLMN),ZCDNLMX)
               ZCDNA=AMIN1(MAX(ZCDNA,ZCDNLMN),ZCDNLMX)
C
            ELSE
               ZCDNN=115.*ZSO4NATGM3**0.48
               ZCDNA=115.*ZSO4ALLGM3**0.48
               ZCDNN=AMIN1(MAX(ZCDNN,ZCDNOMN),ZCDNOMX)
               ZCDNA=AMIN1(MAX(ZCDNA,ZCDNOMN),ZCDNOMX)
            ENDIF
C
            ZRADNAT(JL)=ZREF*ZKAP(JL)*(ZLWC/ZCDNN)**ZREX
            ZRADALL(JL)=ZREF*ZKAP(JL)*(ZLWC/ZCDNA)**ZREX
            ZRADNAT(JL)=MAX(ZRADMIN,MIN(ZRADNAT(JL),ZRADMAX))
            ZRADALL(JL)=MAX(ZRADMIN,MIN(ZRADALL(JL),ZRADMAX))
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
               ZICEWP = ZFLWP(JL) * ZFICE(JL)
               ZLIQWP = ZFLWP(JL) * (1.-ZFICE(JL))
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
C     SHORTWAVE OPTICAL PROPERTIES FOR CLOUD DROPLETS
C     CALCULATED FROM NATURAL SULFUR SOURCES ONLY
C
               ZZLP1=ALOG10(ZRADNAT(JL))
               ZZLP2=ZZLP1*ZZLP1
               ZZLP3=ZZLP2*ZZLP1
               ZZTO1L=ZLIQWP*1.8706058417*ZRADNAT(JL)**(-1.0756364457)
               ZZTO2L=ZLIQWP*1.9655460426*ZRADNAT(JL)**(-1.0778999732)
               ZZG1L=0.78756640717+0.10660598895*ZZLP1-0.031012468401
     &              *ZZLP2
               ZZG2L=0.79208639802-0.044930076174*ZZLP1+0.18980672305
     &              *ZZLP2-0.082590933352*ZZLP3
               ZZO2L=0.9854369057+0.013584242533*ZZLP1-0.024856960461
     &              *ZZLP2+0.0055918314369*ZZLP3

C
C     SHORTWAVE OPTICAL PROPERTIES FOR CLOUD DROPLETS
C     CALCULATED FROM ALL S-SOURCES (NATURAL + ANTHROPOGENIC)
C
               ZZZLP1=ALOG10(ZRADALL(JL))
               ZZZLP2=ZZZLP1*ZZZLP1
               ZZZLP3=ZZZLP2*ZZZLP1
               ZZZTO1L=ZLIQWP*1.8706058417*ZRADALL(JL)**(-1.0756364457)
               ZZZTO2L=ZLIQWP*1.9655460426*ZRADALL(JL)**(-1.0778999732)
               ZZZG1L=0.78756640717+0.10660598895*ZZZLP1-0.031012468401
     &              *ZZZLP2
               ZZZG2L=0.79208639802-0.044930076174*ZZZLP1+0.18980672305
     &              *ZZZLP2-0.082590933352*ZZZLP3
               ZZZO2L=0.9854369057+0.013584242533*ZZZLP1-0.024856960461
     &              *ZZZLP2+0.0055918314369*ZZZLP3

C
C     ANTHROPOGENIC CORRECTION OF 'STANDARD'
C     SHORTWAVE OPTICAL PROPERTIES FOR CLOUD DROPLETS
C
C     (NO CORRECTION FOR ICE CLOUDS AND LONGWAVE PROPERTIES)
C
               ZTO1L=ZTO1L+ZZZTO1L-ZZTO1L
               ZTO2L=ZTO2L+ZZZTO2L-ZZTO2L
               ZG1L=ZG1L+ZZZG1L-ZZG1L
               ZG2L=ZG2L+ZZZG2L-ZZG2L
               ZO2L=ZO2L+ZZZO2L-ZZO2L


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
               ZO2I=0.98475089485+0.0053152066002*ZIP1-0.0061150583857
     &              *ZIP2-0.0032775655896*ZIP3
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
     &              EXP  (-ZMACI*ZICEWP))
            END IF
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
      CALL AOLW(KLON, KLEV , KEWAER, KAERH,  KAER  , KCFC  ,
     &        ZCCO2, ZCLDLW,
     &        PDP , ZPMB , POZON , ZTL    , PAER , ZTAVE ,
     &        PQ  , ZFLUX, ZFLUXC, ZCFCABS)
C
C
C     ------------------------------------------------------------------
C
C*         4.     CALL SHORTWAVE RADIATION CODE
C                 -----------------------------
C
         CALL AOSW(KLON  , KLEV, KEWAER, KAERH,
     &            ZRI0, ZCCO2, ZPSOL , ZALBSU,
     &            PQ   , PQS , PMU0, ZCG  , ZCLDSW, ZOMEGA,
     &            POZON, ZPMB, ZTAU, ZTAVE, PAER  , ZFSDWN,
     &            ZFSUP)
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
CDJ-041294   ZCLDSW(JL,1)=0.    ADD 604 LOOP
C
         CALL AOSW(KLON  , KLEV, KEWAER, KAERH,
     &           ZRI0, ZCCO2, ZPSOL , ZALBSU,
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
      END SUBROUTINE AORADLSW
