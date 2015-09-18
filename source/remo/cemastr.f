C
C      SUBROUTINE CEMASTR
C
C**** *CEMASTR*  MASTER ROUTINE FOR CUMULUS MASSFLUX-SCHEME
C
C     M.TIEDTKE      E.C.M.W.F.     1986/1987/1989
C
C
C     PURPOSE
C     -------
C
C          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE
C     PROGNOSTIC VARIABLES T,Q,U AND V DUE TO CONVECTIVE PROCESSES.
C     PROCESSES CONSIDERED ARE: CONVECTIVE FLUXES, FORMATION OF
C     PRECIPITATION, EVAPORATION OF FALLING RAIN BELOW CLOUD BASE,
C     SATURATED CUMULUS DOWNDRAFTS.
C
C**   INTERFACE.
C     ----------
C
C          *CEMASTR* IS CALLED FROM *CECALL*
C     THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE
C     T,Q,U,V,PHI AND P AND MOISTURE TENDENCIES.
C     IT RETURNS ITS OUTPUT TO THE SAME SPACE
C      1.MODIFIED TENDENCIES OF MODEL VARIABLES
C      2.RATES OF CONVECTIVE PRECIPITATION
C        (USED IN SUBROUTINE SURF)
C      3.CLOUD BASE, CLOUD TOP AND PRECIP FOR RADIATION
C        (USED IN SUBROUTINE CLOUD)
C
C     METHOD
C     -------
C
C     PARAMETERIZATION IS DONE USING A MASSFLUX-SCHEME.
C        (1) DEFINE CONSTANTS AND PARAMETERS
C        (2) SPECIFY VALUES (T,Q,QS...) AT HALF LEVELS AND
C            INITIALIZE UPDRAFT- AND DOWNDRAFT-VALUES IN 'CEINI'
C        (3) CALCULATE CLOUD BASE IN 'CEBASE'
C            AND SPECIFY CLOUD BASE MASSFLUX FROM PBL MOISTURE BUDGET
C        (4) DO CLOUD ASCENT IN 'CEASC' IN ABSENCE OF DOWNDRAFTS
C        (5) DO DOWNDRAFT CALCULATIONS:
C              (A) DETERMINE VALUES AT LFS IN 'CEDLFS'
C              (B) DETERMINE MOIST DESCENT IN 'CEDDRAF'
C              (C) RECALCULATE CLOUD BASE MASSFLUX CONSIDERING THE
C                  EFFECT OF CU-DOWNDRAFTS
C        (6) DO FINAL CLOUD ASCENT IN 'CEASC'
C        (7) DO FINAL ADJUSMENTS TO CONVECTIVE FLUXES IN 'CEFLX',
C            DO EVAPORATION IN SUBCLOUD LAYER
C        (8) CALCULATE INCREMENTS OF T AND Q IN 'CEDTDQ'
C        (9) CALCULATE INCREMENTS OF U AND V IN 'CEDUDV'
C
C     EXTERNALS.
C     ----------
C
C       CEINI:  INITIALIZES VALUES AT VERTICAL GRID USED IN CU-PARAMETR.
C       CEBASE: CLOUD BASE CALCULATION FOR PENETR.AND SHALLOW CONVECTION
C       CEASC:  CLOUD ASCENT FOR ENTRAINING PLUME
C       CEDLFS: DETERMINES VALUES AT LFS FOR DOWNDRAFTS
C       CEDDRAF:DOES MOIST DESCENT FOR CUMULUS DOWNDRAFTS
C       CEFLX:  FINAL ADJUSTMENTS TO CONVECTIVE FLUXES (ALSO IN PBL)
C       CUDQDT: UPDATES TENDENCIES FOR T AND Q
C       CEDUDV: UPDATES TENDENCIES FOR U AND V
C
C     SWITCHES.
C     --------
C
C          LMFPEN=.TRUE.   PENETRATIVE CONVECTION IS SWITCHED ON
C          LMFSCV=.TRUE.   SHALLOW CONVECTION IS SWITCHED ON
C          LMFMID=.TRUE.   MIDLEVEL CONVECTION IS SWITCHED ON
C          LMFDD=.TRUE.    CUMULUS DOWNDRAFTS SWITCHED ON
C          LMFDUDV=.TRUE.  CUMULUS FRICTION SWITCHED ON
C
C
C     MODEL PARAMETERS (DEFINED IN SUBROUTINE CUPARAM)
C     ------------------------------------------------
C     ENTRPEN    ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
C     ENTRSCV    ENTRAINMENT RATE FOR SHALLOW CONVECTION
C     ENTRMID    ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
C     ENTRDD     ENTRAINMENT RATE FOR CUMULUS DOWNDRAFTS
C     CMFCTOP    RELATIVE CLOUD MASSFLUX AT LEVEL ABOVE NONBUOYANCY LEVE
C     CMFCMAX    MAXIMUM MASSFLUX VALUE ALLOWED FOR
C     CMFCMIN    MINIMUM MASSFLUX VALUE (FOR SAFETY)
C     CMFDEPS    FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
C     CPRCON     COEFFICIENT FOR CONVERSION FROM CLOUD WATER TO RAIN
C
C     REFERENCE.
C     ----------
C
C          PAPER ON MASSFLUX SCHEME (TIEDTKE,1989)
C
C
      SUBROUTINE CEMASTR
     &    (KIDIA, KFDIA,
     &     KLON,     KLEV,
     &     KLEVP1,   KLEVM1,   ILAB,
     &     NSTEP ,   NSTART,   TWODT,
     &     PTEN,     PQEN,     PXEN,    PUEN,     PVEN,
     &     PVERV,    PQSEN,    PQHFL,
     &     PAP,      PAPH,     PGEO,     LDLAND,
     &     PTTE,     PQTE,     PVOM,     PVOL,
     &     PRSFC,    PSSFC,    PAPRC,    PAPRS,  PXTEC,
     &     LDCUM,    KTYPE,    KCBOT,    KCTOP,
     &     PTU,      PQU,      PLU,      PLUDE,
     &     PMFU,     PMFD,     PRAIN,    PCEVAPCU,
     &     PSRAIN,   PSEVAP,   PSHEAT,   PSDISS, PSMELT,
     &     PCAPE ,   LCOLDCONV)
C
      IMPLICIT NONE
C
      INCLUDE "COMCON"
      INCLUDE "COMCUMF"
      INCLUDE "COMPH2"
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN)   :: KIDIA, KFDIA, KLON, KLEV,
     &                         KLEVP1,   KLEVM1,
     &                         NSTEP ,   NSTART
C
      INTEGER, INTENT(INOUT)   :: ILAB(KLON,KLEV)
C
      REAL,    INTENT(IN)    :: TWODT
C
      REAL,    INTENT(INOUT) :: PSRAIN,   PSEVAP,   
     &                          PSHEAT,   PSDISS,   PSMELT
C
      REAL,    INTENT(IN)    ::   
     &         PTEN(KLON,KLEV),        PQEN(KLON,KLEV),
     &         PXEN(KLON,KLEV),
     &         PUEN(KLON,KLEV),        PVEN(KLON,KLEV),
     &         PQSEN(KLON,KLEV),       PGEO(KLON,KLEV),
     &         PAP(KLON,KLEV),         PAPH(KLON,KLEVP1),
     &         PVERV(KLON,KLEV),
     &         PQHFL(KLON)
      REAL,    INTENT(INOUT) ::
     &         PTTE(KLON,KLEV),        PQTE(KLON,KLEV),
     &         PVOM(KLON,KLEV),        PVOL(KLON,KLEV)
      REAL,    INTENT(INOUT)  ::    
     &         PTU(KLON,KLEV),         PQU(KLON,KLEV),
     &         PLU(KLON,KLEV),         PLUDE(KLON,KLEV),
     &         PMFU(KLON,KLEV),
     &         PAPRC(KLON),            PAPRS(KLON),
     &         PRSFC(KLON),            PSSFC(KLON),
     &         PRAIN(KLON)
C
      REAL,    INTENT(IN)    ::  PCEVAPCU(KLEV)
      REAL,    INTENT(INOUT) ::  PXTEC(KLON,KLEV)
      REAL,    INTENT(INOUT) ::  PMFD(KLON,KLEV),  PCAPE(KLON)
      INTEGER, INTENT(INOUT) ::  KCBOT(KLON),      KCTOP(KLON)
      INTEGER, INTENT(INOUT) ::  KTYPE(KLON)
      LOGICAL, INTENT(IN)    ::  LDLAND(KLON)
      LOGICAL, INTENT(INOUT) ::  LDCUM(KLON)      
      LOGICAL, INTENT(IN)    ::  LCOLDCONV(KLON)
C
C     Local Variables and Fields
C
      REAL    ::  
     &           ZTENH(KLON,KLEV),       ZQENH(KLON,KLEV),
     &           ZXENH(KLON,KLEV),
     &           ZGEOH(KLON,KLEV),       ZQSENH(KLON,KLEV),
     &           ZTD(KLON,KLEV),         ZQD(KLON,KLEV),
     &           ZMFUS(KLON,KLEV),       ZMFDS(KLON,KLEV),
     &           ZMFUQ(KLON,KLEV),       ZMFDQ(KLON,KLEV),
     &           ZDMFUP(KLON,KLEV),      ZDMFDP(KLON,KLEV),
     &           ZMFUL(KLON,KLEV),       ZRFL(KLON),
     &           ZUU(KLON,KLEV),         ZVU(KLON,KLEV),
     &           ZUD(KLON,KLEV),         ZVD(KLON,KLEV) 
      REAL    :: ZENTR(KLON),            ZHCBASE(KLON),
     &           ZMFUB(KLON),            ZMFUB1(KLON),
     &           ZDQPBL(KLON),           ZDQCV(KLON)
      REAL    :: ZSFL(KLON),             ZDPMEL(KLON,KLEV)
C   --  TEN94 -----
      REAL    :: ZCAPE(KLON),            ZHEAT(KLON)
      REAL    :: ZHMIN(KLON)
      REAL    :: ZHHATT(KLON,KLEV)
      INTEGER :: IHMIN(KLON)
C   ---------------
      INTEGER :: IDTOP(KLON), ICTOP0(KLON), ILWMIN(KLON)
      LOGICAL :: LODDRAF(KLON)
      LOGICAL :: LLO1
C
      INTEGER :: IKB, ITOPM2, JK, JL
      REAL    :: ZALVDCP, ZB, ZBI, ZCONS2, ZDEPTH, ZDHDZ,     
     &           ZDQMIN, ZDZ, ZEPS, ZFAC, ZGAM, ZHHAT, ZHSAT, ZMFMAX,
     &           ZPBMPT, ZQALV, ZQUMQE, ZRH, ZRO, ZTAU, ZTMST, ZZZ
C
C---------------------------------------------------------------------
C
C     1.           SPECIFY CONSTANTS AND PARAMETERS
C                  --------------------------------
C
      ZTMST=TWODT
      IF(NSTEP.EQ.NSTART) ZTMST=0.5*TWODT
      ZCONS2=1./(G*ZTMST)
C
C----------------------------------------------------------------------
C
C*    2.           INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CEINI'
C                  ---------------------------------------------------
C
      CALL CEINI
     &    (KIDIA, KFDIA,
     &     KLON,
     &     KLEV,     KLEVP1,   KLEVM1,
     &     PTEN,     PQEN,     PQSEN,    PXEN, PUEN, PVEN,
     &     PVERV,    PGEO,     PAPH,     ZGEOH,
     &     ZTENH,    ZQENH,    ZQSENH,   ZXENH,  ILWMIN,
     &     PTU,      PQU,      ZTD,      ZQD,
     &     ZUU,      ZVU,      ZUD,      ZVD,
     &     PMFU,     PMFD,     ZMFUS,    ZMFDS,
     &     ZMFUQ,    ZMFDQ,    ZDMFUP,   ZDMFDP,
     &     ZDPMEL,   PLU,      PLUDE,    ILAB)
C
C
C---------------------------------------------------------------------
C
C*    3.0          CLOUD BASE CALCULATIONS
C                  -----------------------
C
C
C*             (A) DETERMINE CLOUD BASE VALUES IN 'CEBASE'
C                  ---------------------------------------
C
      CALL CEBASE
     &    (KIDIA, KFDIA,
     &     KLON,     KLEV,     KLEVP1,   KLEVM1,
     &     ZTENH,    ZQENH,    ZGEOH,    PAPH,
     &     PTU,      PQU,      PLU,
     &     PUEN,     PVEN,     ZUU,      ZVU,
     &     LDCUM,    KCBOT,    ILAB)
C
C*             (B) DETERMINE TOTAL MOISTURE CONVERGENCE AND
C*                 THEN DECIDE ON TYPE OF CUMULUS CONVECTION
C                  -----------------------------------------
C
      JK=1
      DO JL=KIDIA,KFDIA
         ZDQCV(JL) =PQTE(JL,JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))
         ZDQPBL(JL)=0.0
         IDTOP(JL)=0
      ENDDO
      DO JK=2,KLEV
         DO JL=KIDIA,KFDIA
            ZDQCV(JL)=ZDQCV(JL)+PQTE(JL,JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))
            IF(JK.GE.KCBOT(JL)) ZDQPBL(JL)=ZDQPBL(JL)+PQTE(JL,JK)
     &           *(PAPH(JL,JK+1)-PAPH(JL,JK))
         ENDDO
      ENDDO
      DO JL=KIDIA,KFDIA
         IF (ZDQCV(JL).GT.MAX(0.,-1.1*PQHFL(JL)*G)) THEN
            KTYPE(JL)=1
         ELSE
            KTYPE(JL)=2
CSP
            IF (LCOLDCONV(JL)) THEN
               KTYPE(JL)=4
            ENDIF
CSP
         ENDIF
      ENDDO
C
C*             (C) DETERMINE MOISTURE SUPPLY FOR BOUNDARY LAYER
C*                 AND DETERMINE CLOUD BASE MASSFLUX IGNORING
C*                 THE EFFECTS OF DOWNDRAFTS AT THIS STAGE
C                  ------------------------------------------
C
!DIR$ IVDEP
      DO JL=KIDIA,KFDIA
         IKB=KCBOT(JL)
         ZQUMQE=PQU(JL,IKB)+PLU(JL,IKB)-ZQENH(JL,IKB)
         ZDQMIN=MAX(0.01*ZQENH(JL,IKB),1.E-10)
         LLO1=ZDQPBL(JL).GT.0..AND.ZQUMQE.GT.ZDQMIN.AND.LDCUM(JL)
         IF (LLO1) THEN
            ZMFUB(JL)=ZDQPBL(JL)/(G*MAX(ZQUMQE,ZDQMIN))
         ELSE
            ZMFUB(JL)=0.01
         ENDIF
         ZMFMAX=(PAPH(JL,IKB)-PAPH(JL,IKB-1))*ZCONS2
         ZMFUB(JL)=MIN(ZMFUB(JL),ZMFMAX)
         IF(.NOT.LLO1) LDCUM(JL)=.FALSE.
         ZCAPE(JL)=0.
         ZHEAT(JL)=0.
CSP      ADDED '.OR. KTYPE(JL) .EQ. 4'
         IF (KTYPE(JL).EQ.1 .OR. KTYPE(JL) .EQ. 4) THEN
            ZENTR(JL)=ENTRPEN
         ELSE
            ZENTR(JL)=ENTRSCV
         ENDIF
CSP
      ENDDO
C
C
C-----------------------------------------------------------------------
C
C*    4.0          DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
C                  -------------------------------------------
C
C*             (A) ESTIMATE CLOUD HEIGHT FOR ENTRAINMENT/DETRAINMENT
C*                 CALCULATIONS IN CEASC (MAX.POSSIBLE CLOUD HEIGHT
C*                 FOR NON-ENTRAINING PLUME, FOLLOWING A.-S.,1974)
C                  -------------------------------------------------
C
!DIR$ IVDEP
      DO JL=KIDIA,KFDIA
         IKB=KCBOT(JL)
         ZHCBASE(JL)=CPD*PTU(JL,IKB)+ZGEOH(JL,IKB)+ALV*PQU(JL,IKB)
         ICTOP0(JL)=KCBOT(JL)-1
      ENDDO
      ZALVDCP=ALV/CPD
      ZQALV=1./ALV
      DO JK=KLEVM1,3,-1
!DIR$ IVDEP
         DO JL=KIDIA,KFDIA
            ZHSAT=CPD*ZTENH(JL,JK)+ZGEOH(JL,JK)+ALV*ZQSENH(JL,JK)
            ZGAM=C5LES*ZALVDCP*ZQSENH(JL,JK)/
     &           ((1.-VTMPC1*ZQSENH(JL,JK))*(ZTENH(JL,JK)-C4LES)**2)
            ZZZ=CPD*ZTENH(JL,JK)*0.608
            ZHHAT=ZHSAT-(ZZZ+ZGAM*ZZZ)/(1.+ZGAM*ZZZ*ZQALV)*
     &           MAX(ZQSENH(JL,JK)-ZQENH(JL,JK),0.)
            ZHHATT(JL,JK)=ZHHAT
            IF(JK.LT.ICTOP0(JL).AND.ZHCBASE(JL).GT.ZHHAT) ICTOP0(JL)=JK
         ENDDO
      ENDDO
C
      DO JL=1,KLON
         JK=KCBOT(JL)
         ZHSAT=CPD*ZTENH(JL,JK)+ZGEOH(JL,JK)+ALV*ZQSENH(JL,JK)
         ZGAM=C5LES*ZALVDCP*ZQSENH(JL,JK)/
     &        ((1.-VTMPC1*ZQSENH(JL,JK))*(ZTENH(JL,JK)-C4LES)**2)
         ZZZ=CPD*ZTENH(JL,JK)*0.608
         ZHHAT=ZHSAT-(ZZZ+ZGAM*ZZZ)/(1.+ZGAM*ZZZ*ZQALV)*
     &        MAX(ZQSENH(JL,JK)-ZQENH(JL,JK),0.)
         ZHHATT(JL,JK)=ZHHAT
      ENDDO
C
C
C                  FIND LOWEST POSSIBLE ORG. DETRAINMENT LEVEL
C                  -------------------------------------------
C
      DO JL=1,KLON
         LLO1=LDCUM(JL).AND.KTYPE(JL).EQ.1
         IF(LLO1) THEN
            IKB=KCBOT(JL)
            ZHMIN(JL)=0.
            IHMIN(JL)=IKB
         ENDIF
      ENDDO
C
      ZB=25.
      ZBI=1./(ZB*G)
      DO JK=KLEV,1,-1
         DO JL=1,KLON
            LLO1=LDCUM(JL).AND.KTYPE(JL).EQ.1.AND.IHMIN(JL).EQ.KCBOT(JL)
            IF(LLO1.AND.JK.LT.KCBOT(JL).AND.JK.GE.ICTOP0(JL)) THEN
               IKB=KCBOT(JL)
               ZRO=PAPH(JL,JK)/(RD*ZTENH(JL,JK))
               ZDZ=(PAPH(JL,JK)-PAPH(JL,JK-1))/(G*ZRO)
               ZDHDZ=( CPD*(PTEN(JL,JK-1)-PTEN(JL,JK))+
     &              ALV*(PQEN(JL,JK-1)-PQEN(JL,JK))+
     &              (PGEO(JL,JK-1)-PGEO(JL,JK)) )*G/
     &              (PGEO(JL,JK-1)-PGEO(JL,JK))
               ZDEPTH=ZGEOH(JL,JK)-ZGEOH(JL,IKB)
               ZFAC=SQRT(1.+ZDEPTH*ZBI)
               ZHMIN(JL)=ZHMIN(JL) + ZDHDZ*ZFAC*ZDZ
               ZRH=-ALV*(ZQSENH(JL,JK)-ZQENH(JL,JK))*ZFAC
               IF(ZHMIN(JL).GT.ZRH) IHMIN(JL)=JK
            ENDIF
         ENDDO
      ENDDO
C
      DO JL=1,KLON
         IF(LDCUM(JL).AND.KTYPE(JL).EQ.1) THEN
            IF(IHMIN(JL).LT.ICTOP0(JL)) IHMIN(JL)=ICTOP0(JL)
         ENDIF
      ENDDO

C
C*             (B) DO ASCENT IN 'CEASC'IN ABSENCE OF DOWNDRAFTS
C                  --------------------------------------------
C
      CALL CEASC
     &    (KIDIA, KFDIA,
     &     KLON,
     &     KLEV,     KLEVP1,   KLEVM1,
     &     NSTEP,    NSTART,   TWODT,
     &     ZTENH,    ZQENH,    ZXENH,    PUEN,     PVEN,
     &     PTEN,     PQEN,     PQSEN,
     &     PGEO,     ZGEOH,    PAP,      PAPH,
     &     PVERV,    ILWMIN,
     &     LDLAND,   LDCUM,    KTYPE,    ILAB,
     &     PTU,      PQU,      PLU,      ZUU,      ZVU,
     &     PMFU,     ZMFUB,    ZENTR,
     &     ZMFUS,    ZMFUQ,
     &     ZMFUL,    PLUDE,    ZDMFUP,
     &     IHMIN,    ZHHATT,   ZHCBASE,   ZQSENH,
     &     KCBOT,    KCTOP,    ICTOP0)
C
C
C*         (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
C              CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
C              -----------------------------------------------------
C
C
      DO JL=KIDIA,KFDIA
         ZPBMPT=PAPH(JL,KCBOT(JL))-PAPH(JL,KCTOP(JL))
         IF(LDCUM(JL).AND.KTYPE(JL).EQ.1.AND.ZPBMPT.LT.2.E4) THEN
            KTYPE(JL)=2
CSP
            IF (LCOLDCONV(JL)) THEN
               KTYPE(JL)=4
            ENDIF
CSP
         ENDIF
CTEN94   IF(LDCUM(JL)) ICTOP0(JL)=KCTOP(JL)
         IF(KTYPE(JL).EQ.2) ZENTR(JL)=ENTRSCV
         ZRFL(JL)=ZDMFUP(JL,1)
      ENDDO
      DO JK=2,KLEV
         DO JL=KIDIA,KFDIA
            ZRFL(JL)=ZRFL(JL)+ZDMFUP(JL,JK)
         ENDDO
      ENDDO
C
C
C-----------------------------------------------------------------------
C
C*    5.0          CUMULUS DOWNDRAFT CALCULATIONS
C                  ------------------------------
C
C
      IF(LMFDD) THEN
C
C*             (A) DETERMINE LFS IN 'CEDLFS'
C                  -------------------------
C
         CALL CEDLFS
     &    (KIDIA, KFDIA,
     &     KLON,     KLEV,     KLEVP1,
     &     ZTENH,    ZQENH,    PUEN,     PVEN,
     &     ZGEOH,    PAPH,
     &     PTU,      PQU,      ZUU,      ZVU,
     &     LDCUM,    KCBOT,    KCTOP,    ZMFUB,    ZRFL,
     &     ZTD,      ZQD,      ZUD,      ZVD,
     &     PMFD,     ZMFDS,    ZMFDQ,    ZDMFDP,
     &     IDTOP,    LODDRAF)
C
C*            (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CEDDRAF'
C                  -----------------------------------------------
C
         CALL CEDDRAF
     &    (KIDIA, KFDIA,
     &     KLON,     KLEV,     KLEVP1,
     &     ZTENH,    ZQENH,    PUEN,     PVEN,
     &     ZGEOH,    PAPH,     ZRFL,
     &     ZTD,      ZQD,      ZUD,      ZVD,
     &     PMFD,     ZMFDS,    ZMFDQ,    ZDMFDP,
     &     LODDRAF)
C
C ENDIF FOR LMFDD
      END IF
C
C
C*    5.5          RECALCULATE CLOUD BASE MASSFLUX FROM A
C*                 CAPE CLOSURE FOR DEEP AND COLD CONVECTION (KTYPE=1, KTYPE=4)
C*                 AND BY PBL EQUILIBRUM TAKING DOWNDRAFTS INTO
C*                 ACCOUNT FOR SHALLOW CONVECTION (KTYPE=2)
C                  -------------------------------------------
C
C
      DO JL=1,KLON
         ZHEAT(JL)=0.
         ZCAPE(JL)=0.
         ZMFUB1(JL)=ZMFUB(JL)
      ENDDO
C
      DO JK=1,KLEV
         DO JL=1,KLON
CPL     LLO1=LDCUM(JL).AND.KTYPE(JL).EQ.1
            LLO1=LDCUM(JL)
            IF(LLO1.AND.JK.LE.KCBOT(JL).AND.JK.GE.KCTOP(JL)) THEN
               IKB=KCBOT(JL)
               ZRO=PAPH(JL,JK)/(RD*ZTENH(JL,JK))
               ZDZ=(PAPH(JL,JK)-PAPH(JL,JK-1))/(G*ZRO)
               ZHEAT(JL)=ZHEAT(JL) +
     &              (  (PTEN(JL,JK-1)-PTEN(JL,JK) 
     &              + G*ZDZ/CPD)/ZTENH(JL,JK)
     &              +  0.608*(PQEN(JL,JK-1)-PQEN(JL,JK))  ) *
     &              (G*(PMFU(JL,JK)+PMFD(JL,JK)))/ZRO
               ZCAPE(JL)=ZCAPE(JL) +
     &              (G*(PTU(JL,JK)-ZTENH(JL,JK))/ZTENH(JL,JK)
     &              +G*0.608*(PQU(JL,JK)-ZQENH(JL,JK))
     &              -G*PLU(JL,JK) ) * ZDZ
            ENDIF
         ENDDO
      ENDDO
C
      ZTAU=CTAU
      DO JL=1,KLON
         PCAPE(JL)=PCAPE(JL)+ZCAPE(JL)*TWODT*0.5
CSP   ADDED '.OR.KTYPE(JL).EQ.4'
         IF(LDCUM(JL).AND.
     &        (KTYPE(JL).EQ.1.OR.KTYPE(JL).EQ.4)) THEN
CSP
            IKB=KCBOT(JL)
            IF (ZHEAT(JL) .NE. 0.) THEN
               ZMFUB1(JL)=(ZCAPE(JL)*ZMFUB(JL))/(ZHEAT(JL)*ZTAU)
            ELSE
               ZMFUB1(JL)=0.001
            ENDIF
            ZMFUB1(JL)=MAX(ZMFUB1(JL),0.001)
            ZMFMAX=(PAPH(JL,IKB)-PAPH(JL,IKB-1))*ZCONS2
            ZMFUB1(JL)=MIN(ZMFUB1(JL),ZMFMAX)
         ENDIF
      ENDDO
C
C
C
C
C
C*                 RECALCULATE CONVECTIVE FLUXES DUE TO EFFECT OF
C*                 DOWNDRAFTS ON BOUNDARY LAYER MOISTURE BUDGET
C*                 FOR SHALLOW CONVECTION (KTYPE=2)
C                  --------------------------------------------
C
      DO JL=KIDIA,KFDIA
CSP      'ADDED .AND. KTYPE(JL).NE.4: MOISTURE BUDGET CLOSURE NOT APPLIED TO KTYPE 4
         IF(KTYPE(JL).NE.1.AND.KTYPE(JL).NE.4) THEN
CSP
            IKB=KCBOT(JL)
            LLO1=PMFD(JL,IKB).LT.0..AND.LODDRAF(JL)
            IF (LLO1) THEN
               ZEPS=CMFDEPS
            ELSE
               ZEPS=0.
            ENDIF
            ZQUMQE=PQU(JL,IKB)+PLU(JL,IKB)-
     &           ZEPS*ZQD(JL,IKB)-(1.-ZEPS)*ZQENH(JL,IKB)
            ZDQMIN=MAX(0.01*ZQENH(JL,IKB),1.E-10)
            ZMFMAX=(PAPH(JL,IKB)-PAPH(JL,IKB-1))*ZCONS2
            LLO1=ZDQPBL(JL).GT.0..AND.ZQUMQE.GT.ZDQMIN.AND.LDCUM(JL)
     &           .AND.ZMFUB(JL).LT.ZMFMAX
            IF (LLO1) THEN
               ZMFUB1(JL)=ZDQPBL(JL)/(G*MAX(ZQUMQE,ZDQMIN))
            ELSE
               ZMFUB1(JL)=ZMFUB(JL)
            ENDIF
CSP 'ADDED .OR. KTYPE(JL).EQ.4
            IF ((KTYPE(JL).EQ.1.OR.KTYPE(JL).EQ.2.OR.KTYPE(JL).EQ.4)
     &           .AND.(ABS(ZMFUB1(JL)-ZMFUB(JL)).LT.0.2*ZMFUB(JL))) THEN
CSP
               ZMFUB1(JL)=ZMFUB1(JL)
            ELSE
               ZMFUB1(JL)=ZMFUB(JL)
            ENDIF
         END IF
      ENDDO
      DO JK=1,KLEV
         DO JL=KIDIA,KFDIA
            IF(LDCUM(JL)) THEN
               ZFAC=ZMFUB1(JL)/MAX(ZMFUB(JL),1.E-10)
               PMFD(JL,JK)=PMFD(JL,JK)*ZFAC
               ZMFDS(JL,JK)=ZMFDS(JL,JK)*ZFAC
               ZMFDQ(JL,JK)=ZMFDQ(JL,JK)*ZFAC
               ZDMFDP(JL,JK)=ZDMFDP(JL,JK)*ZFAC
            ENDIF
         ENDDO
      ENDDO
C
C
C*                 NEW VALUES OF CLOUD BASE MASS FLUX
C                  ----------------------------------
C
      DO JL=KIDIA,KFDIA
         IF(LDCUM(JL)) ZMFUB(JL)=ZMFUB1(JL)
      ENDDO
C
C
C
C-----------------------------------------------------------------------
C
C*    6.0          DETERMINE FINAL CLOUD ASCENT FOR ENTRAINING PLUME
C*                 FOR PENETRATIVE CONVECTION (TYPE=1),
C*                 FOR SHALLOW TO MEDIUM CONVECTION (TYPE=2)
C*                 AND FOR MID-LEVEL CONVECTION (TYPE=3).
C                  -------------------------------------------------
C
      CALL CEASC
     &    (KIDIA, KFDIA,
     &     KLON,
     &     KLEV,     KLEVP1,   KLEVM1,
     &     NSTEP,    NSTART,   TWODT,
     &     ZTENH,    ZQENH,    ZXENH,    PUEN,     PVEN,
     &     PTEN,     PQEN,     PQSEN,
     &     PGEO,     ZGEOH,    PAP,      PAPH,
     &     PVERV,    ILWMIN,
     &     LDLAND,   LDCUM,    KTYPE,    ILAB,
     &     PTU,      PQU,      PLU,      ZUU,      ZVU,
     &     PMFU,     ZMFUB,    ZENTR,
     &     ZMFUS,    ZMFUQ,
     &     ZMFUL,    PLUDE,    ZDMFUP,
     &     IHMIN,    ZHHATT,   ZHCBASE,  ZQSENH,
     &     KCBOT,    KCTOP,    ICTOP0)
C
C
C-----------------------------------------------------------------------
C
C*    7.0          DETERMINE FINAL CONVECTIVE FLUXES IN 'CEFLX'
C                  ------------------------------------------
C
      CALL CEFLX
     &    (KIDIA, KFDIA,
     &     KLON,     KLEV,     KLEVP1,
     &     NSTEP,    NSTART,   TWODT,
     &     PQEN,     PQSEN,    ZTENH,    ZQENH,
     &     PAPH,     ZGEOH,    PCEVAPCU,
     &     KCBOT,    KCTOP,    IDTOP,
     &     KTYPE,    LODDRAF,  LDCUM,
     &     PMFU,     PMFD,     ZMFUS,    ZMFDS,
     &     ZMFUQ,    ZMFDQ,    ZMFUL,    PLUDE,
     &     ZDMFUP,   ZDMFDP,   ZRFL,     PRAIN,
     &     PTEN,     ZSFL,     ZDPMEL,   ITOPM2)
C
C----------------------------------------------------------------------
C
C*    8.0          UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CEDTDQ
C                  --------------------------------------------------
C
      CALL CEDTDQ
     &    (KIDIA, KFDIA,
     &     KLON,     KLEV,     KLEVP1,
     &     TWODT,
     &     ITOPM2,   PAPH,
     &     LDCUM,    PTEN,     PTTE,     PQTE,
     &     PXTEC,
     &     ZMFUS,    ZMFDS,    ZMFUQ,    ZMFDQ,
     &     ZMFUL,    ZDMFUP,   ZDMFDP,   PLUDE,
     &     ZDPMEL,   PRAIN,    ZRFL,     ZSFL,
     &     PSRAIN,   PSEVAP,   PSHEAT,   PSMELT,
     &     PRSFC,    PSSFC,    PAPRC,    PAPRS)
C
C
C----------------------------------------------------------------------
C
C*    9.0          UPDATE TENDENCIES FOR U AND U IN SUBROUTINE CEDUDV
C                  --------------------------------------------------
C
      IF(LMFDUDV) THEN
         CALL CEDUDV
     &    (KIDIA, KFDIA,
     &     KLON,     KLEV,     KLEVP1,
     &     ITOPM2,   KTYPE,    KCBOT,    PAPH,     LDCUM,
     &     PUEN,     PVEN,     PVOM,     PVOL,
     &     ZUU,      ZUD,      ZVU,      ZVD,
     &     PMFU,     PMFD,     PSDISS)
C
      ENDIF
C
      RETURN
      END SUBROUTINE CEMASTR
