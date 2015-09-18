      SUBROUTINE CEASC
     &    (KIDIA, KFDIA, KLON,
     &     KLEV,     KLEVP1,   KLEVM1,
     &     NSTEP,    NSTART,   TWODT,
     &     PTENH,    PQENH,    PXENH,   PUEN,     PVEN,
     &     PTEN,     PQEN,     PQSEN,
     &     PGEO,     PGEOH,    PAP,      PAPH,
     &     PVERV,    KLWMIN,
     &     LDLAND,   LDCUM,    KTYPE,    KLAB,
     &     PTU,      PQU,      PLU,      PUU,      PVU,
     &     PMFU,     PMFUB,    PENTR,
     &     PMFUS,    PMFUQ,
     &     PMFUL,    PLUDE,    PDMFUP,
     &     KHMIN,    PHHATT,   PHCBASE,  PQSENH,
     &     KCBOT,    KCTOP,    KCTOP0 )
C
      IMPLICIT NONE
C
C          THIS ROUTINE DOES THE CALCULATIONS FOR CLOUD ASCENTS
C          FOR CUMULUS PARAMETERIZATION
C
C          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
C
C          PURPOSE.
C          --------
C          TO PRODUCE CLOUD ASCENTS FOR CU-PARAMETRIZATION
C          (VERTICAL PROFILES OF T,Q,L,U AND V AND CORRESPONDING
C           FLUXES AS WELL AS PRECIPITATION RATES)
C
C          INTERFACE
C          ---------
C
C          THIS ROUTINE IS CALLED FROM *CEMASTR*.
C
C          METHOD.
C          --------
C          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
C          AND THEN CALCULATE MOIST ASCENT FOR
C          ENTRAINING/DETRAINING PLUME.
C          ENTRAINMENT AND DETRAINMENT RATES DIFFER FOR
C          SHALLOW AND DEEP CUMULUS CONVECTION.
C          IN CASE THERE IS NO PENETRATIVE OR SHALLOW CONVECTION
C          CHECK FOR POSSIBILITY OF MID LEVEL CONVECTION
C          (CLOUD BASE VALUES CALCULATED IN *CEBASMC*)
C
C          EXTERNALS
C          ---------
C          *CEADJTQ* ADJUST T AND Q DUE TO CONDENSATION IN ASCENT
C          *CEENTR*  CALCULATE ENTRAINMENT/DETRAINMENT RATES
C          *CEBASMC* CALCULATE CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
C
C          REFERENCE
C          ---------
C          (TIEDTKE,1989)
C
      INCLUDE "COMCON"
      INCLUDE "COMCUMF"
      INCLUDE "comconv.h"
C
      INTEGER, INTENT(IN)   :: KIDIA, KFDIA, KLON,
     &                         KLEV,  KLEVP1, KLEVM1,
     &                         NSTEP, NSTART
C
      REAL,    INTENT(IN)   :: TWODT
C
      REAL,    INTENT(IN)   :: 
     &         PTENH(KLON,KLEV),       PQENH(KLON,KLEV),
     &         PXENH(KLON,KLEV),
     &         PUEN(KLON,KLEV),        PVEN(KLON,KLEV),
     &         PTEN(KLON,KLEV),        PQEN(KLON,KLEV),
     &         PGEO(KLON,KLEV),        PGEOH(KLON,KLEV),
     &         PAP(KLON,KLEV),         PAPH(KLON,KLEVP1),
     &         PQSEN(KLON,KLEV),
     &         PVERV(KLON,KLEV)
C
      REAL,    INTENT(INOUT) :: 
     &         PTU(KLON,KLEV),         PQU(KLON,KLEV),
     &         PUU(KLON,KLEV),         PVU(KLON,KLEV),
     &         PMFU(KLON,KLEV),
     &         PMFUB(KLON),            PENTR(KLON),
     &         PMFUS(KLON,KLEV),       PMFUQ(KLON,KLEV),
     &         PLU(KLON,KLEV),         PLUDE(KLON,KLEV),
     &         PMFUL(KLON,KLEV),       PDMFUP(KLON,KLEV)
C
      INTEGER, INTENT(IN)    ::        KLWMIN(KLON)
C
      INTEGER, INTENT(INOUT) ::        KTYPE(KLON),
     &         KLAB(KLON,KLEV),        KCBOT(KLON),
     &         KCTOP(KLON),            KCTOP0(KLON)
CDJ MODIFICATION DUE TO T.E.NORDENG SCHEME
      INTEGER, INTENT(IN)    ::        KHMIN(KLON)
      REAL,    INTENT(IN)    ::        PHHATT(KLON,KLEV)
      REAL,    INTENT(IN)    ::        PHCBASE(KLON)
      REAL,    INTENT(IN)    ::        PQSENH(KLON,KLEV)
      LOGICAL, INTENT(IN)    ::        LDLAND(KLON)
      LOGICAL, INTENT(INOUT) ::        LDCUM(KLON)
C
C     Local Variables
C
      REAL     :: ZDMFEN(KLON), ZDMFDE(KLON),
     &            ZMFUU(KLON),  ZMFUV(KLON),
     &            ZPBASE(KLON), ZQOLD(KLON)
      REAL     :: ZDLAND(KLON)
      REAL     :: ZPH(KLON)
CDJ  T.E.NORDENG
      REAL     :: ZODETR(KLON,KLEV)
      REAL     :: ZOENTR(KLON,KLEV)
      REAL     :: ZBUOY(KLON)
      LOGICAL  :: LOFLAG(KLON)
C
      REAL    :: ZDPHI, ZDNOPRC, ZDMFEU, ZDMFDU, ZCONS2, ZBUOYZ, ZBUO
      REAL    :: ZMFULK, ZMFTEST, ZMFMAX, ZLNEW, ZGA, ZFAC, ZDZ, ZDT
      REAL    :: ZDRODZ, ZDPRHO, ZSCDE, ZQUDE, ZQEEN, ZQCOD, ZPRCON
      REAL    :: ZODMAX, ZNEVN, ZMSE, ZMFUSK, ZMFUQK, ZZDMF, ZZ, ZTMST
      REAL    :: ZTGLACE, ZSEEN, ZSCOD, FACCO
      INTEGER :: IS,JK,JL,IKT,IKB,IK,ICALL
C
C----------------------------------------------------------------------
C
C*    1.           SPECIFY PARAMETERS
C                  ------------------
C
      ZTMST=TWODT
      IF(NSTEP.EQ.NSTART) ZTMST=0.5*TWODT
      ZCONS2=1./(G*ZTMST)
      ZTGLACE=TMELT-13.
C
C
C----------------------------------------------------------------------
C
C     2.           SET DEFAULT VALUES
C                  ------------------
C
      DO JL=KIDIA,KFDIA
         ZMFUU(JL)=0.
         ZMFUV(JL)=0.
         IF(.NOT.LDCUM(JL)) KTYPE(JL)=0
      ENDDO
      DO JK=1,KLEV
         DO JL=KIDIA,KFDIA
            PLU(JL,JK)=0.
            PMFU(JL,JK)=0.
            PMFUS(JL,JK)=0.
            PMFUQ(JL,JK)=0.
            PMFUL(JL,JK)=0.
            PLUDE(JL,JK)=0.
            PDMFUP(JL,JK)=0.
            IF(.NOT.LDCUM(JL).OR.KTYPE(JL).EQ.3) KLAB(JL,JK)=0
            IF(.NOT.LDCUM(JL).AND.PAPH(JL,JK).LT.4.E4) KCTOP0(JL)=JK
         ENDDO
      ENDDO
!DIR$ IVDEP
      DO JL=KIDIA,KFDIA
         IF(LDLAND(JL)) THEN
            ZDLAND(JL)=DLAND
            ZDPHI=PGEOH(JL,KCTOP0(JL))-PGEOH(JL,KCBOT(JL))
            IF(PTU(JL,KCTOP0(JL)).GE.ZTGLACE) ZDLAND(JL)=ZDPHI
            ZDLAND(JL)=MAX(3.0E4,ZDLAND(JL))
            ZDLAND(JL)=MIN(5.0E4,ZDLAND(JL))
            ZDLAND(JL)=DLAND
         ENDIF
      ENDDO
      DO JK=1,KLEV
         DO JL=1,KLON
            ZOENTR(JL,JK)=0.
            ZODETR(JL,JK)=0.
         ENDDO
      ENDDO
C
C
C----------------------------------------------------------------------
C
C     3.0          INITIALIZE VALUES AT LIFTING LEVEL
C                  ----------------------------------
C
      DO JL=KIDIA,KFDIA
         KCTOP(JL)=KLEVM1
         IF(.NOT.LDCUM(JL)) THEN
            KCBOT(JL)=KLEVM1
            PMFUB(JL)=0.
            PQU(JL,KLEV)=0.
         ENDIF
         PMFU(JL,KLEV)=PMFUB(JL)
         PMFUS(JL,KLEV)=PMFUB(JL)*(CPD*PTU(JL,KLEV)+PGEOH(JL,KLEV))
         PMFUQ(JL,KLEV)=PMFUB(JL)*PQU(JL,KLEV)
         IF(LMFDUDV) THEN
            ZMFUU(JL)=PMFUB(JL)*PUU(JL,KLEV)
            ZMFUV(JL)=PMFUB(JL)*PVU(JL,KLEV)
         ENDIF
      ENDDO
C
      DO JL=KIDIA,KFDIA
         LDCUM(JL)=.FALSE.
      ENDDO
C
C
C
C
C T.E.NORDENG SCHEME
C----------------------------------------------------------------------
C
C     3.5          FIND ORGANIZED ENTRAINMENT AT CLOUD BASE
C                  ----------------------------------------
C
      DO JL=1,KLON
         IF(KTYPE(JL).EQ.1) THEN
            IKB=KCBOT(JL)
            ZBUOY(JL)=G*(PTU(JL,IKB)-PTENH(JL,IKB))/PTENH(JL,IKB) +
     &          G*0.608*(PQU(JL,IKB)-PQENH(JL,IKB))
            IF(ZBUOY(JL).GT.0.) THEN
               ZDZ=(PGEO(JL,IKB-1)-PGEO(JL,IKB))/G
               ZDRODZ=-ALOG(PTEN(JL,IKB-1)/PTEN(JL,IKB))/ZDZ
     &              -G/(RD*PTENH(JL,IKB))
         ! NB ZOENTR IS HERE A FRACTIONAL VALUE
               ZOENTR(JL,IKB-1)=ZBUOY(JL)*0.5/(1.+ZBUOY(JL)*ZDZ)
     &              + ZDRODZ
               IF(ZOENTR(JL,IKB-1).GT.1.E-3) ZOENTR(JL,IKB-1)=1.E-3
               IF(ZOENTR(JL,IKB-1).LT.0.0) ZOENTR(JL,IKB-1)=0.0
            ENDIF
         ENDIF
      ENDDO
C
C----------------------------------------------------------------------
C
C     4.           DO ASCENT: SUBCLOUD LAYER (KLAB=1) ,CLOUDS (KLAB=2)
C                  BY DOING FIRST DRY-ADIABATIC ASCENT AND THEN
C                  BY ADJUSTING T,Q AND L ACCORDINGLY IN *CEADJTQ*,
C                  THEN CHECK FOR BUOYANCY AND SET FLAGS ACCORDINGLY
C                  -------------------------------------------------
C
      DO JK=KLEVM1,2,-1
C
C                  SPECIFY CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
C                  IN *CEBASMC* IN CASE THERE IS NOT ALREADY CONVECTION
C                  ----------------------------------------------------
C
         IK=JK
         CALL CEBASMC
     &    (KIDIA, KFDIA,
     &     KLON,     KLEV,     KLEVM1,   IK,
     &     PTEN,     PQEN,     PQSEN,    PUEN,     PVEN,
     &     PVERV,    PGEO,     PGEOH,    LDCUM,    KTYPE,    KLAB,
     &     PMFU,     PMFUB,    PENTR,    KCBOT,
     &     PTU,      PQU,      PLU,      PUU,      PVU,
     &     PMFUS,    PMFUQ,    PMFUL,    PDMFUP,   ZMFUU,    ZMFUV)
C
         IS=0
         DO JL=KIDIA,KFDIA
            IS=IS+KLAB(JL,JK+1)
            IF(KLAB(JL,JK+1).EQ.0) KLAB(JL,JK)=0
C     LOFLAG(JL)=LCVMGT(.TRUE.,.FALSE.,KLAB(JL,JK+1).GT.0)
CMB SUBSTITUTED BY
            IF (KLAB(JL,JK+1).GT.0) THEN
               LOFLAG(JL)=.TRUE.
            ELSE
               LOFLAG(JL)=.FALSE.
            ENDIF
C
            ZPH(JL)=PAPH(JL,JK)
            IF(KTYPE(JL).EQ.3.AND.JK.EQ.KCBOT(JL)) THEN
               ZMFMAX=(PAPH(JL,JK)-PAPH(JL,JK-1))*ZCONS2
               IF(PMFUB(JL).GT.ZMFMAX) THEN
                  ZFAC=PMFUB(JL)/ZMFMAX
                  PMFU(JL,JK+1)=PMFU(JL,JK+1)*ZFAC
                  PMFUS(JL,JK+1)=PMFUS(JL,JK+1)*ZFAC
                  PMFUQ(JL,JK+1)=PMFUQ(JL,JK+1)*ZFAC
                  ZMFUU(JL)=ZMFUU(JL)*ZFAC
                  ZMFUV(JL)=ZMFUV(JL)*ZFAC
                  PMFUB(JL)=ZMFMAX
               ENDIF
            ENDIF
         ENDDO
         IF(IS.EQ.0) CYCLE
C
C
C*                 SPECIFY TURBULENT ENTRAINMENT AND DETRAINMENTS
C                  RATES PLUS ORGANIZED DETRAINMENT RATES IN *CEENTR*
C                   -------------------------------------
C
         IK=JK
         CALL CEENTR
     &    (KIDIA, KFDIA,
     &     KLON,     KLEV,     KLEVP1,   IK,
     &     PTENH,    PAPH,     PAP,
     &     KLWMIN,   LDCUM,    KTYPE,    KCBOT,    KCTOP0,
     &     ZPBASE,   PMFU,     PENTR,    ZODETR,
     &     KHMIN,    PGEOH,
     &     ZDMFEN,   ZDMFDE)
C
C
C
C                  DO ADIABATIC ASCENT FOR ENTRAINING/DETRAINING PLUME
C                  THE CLOUD ENSEMBLE ENTRAINS ENVIRONMENTAL VALUES
C                  IN TURBULENT DETRAINMENT CLOUD ENSEMBLE VALUES
C                  ARE DETRAINED
C                  IN ORGANIZED DETRAINMENT THE DRY STATIC ENERGY AND
C                  MOISTURE THAT ARE NEUTRAL COMPARED TO THE
C                  ENVIRONMENTAL AIR ARE DETRAINED
C                  ---------------------------------------------------
C
         DO JL=KIDIA,KFDIA
            IF(LOFLAG(JL)) THEN
               IF(JK.LT.KCBOT(JL)) THEN
                  ZMFTEST=PMFU(JL,JK+1)+ZDMFEN(JL)-ZDMFDE(JL)
                  ZMFMAX=MIN(ZMFTEST,(PAPH(JL,JK)-PAPH(JL,JK-1))*ZCONS2)
                  ZDMFEN(JL)=MAX(ZDMFEN(JL)-MAX(ZMFTEST-ZMFMAX,0.),0.)
               END IF
               ZDMFDE(JL)=MIN(ZDMFDE(JL),0.75*PMFU(JL,JK+1))
               PMFU(JL,JK)=PMFU(JL,JK+1)+ZDMFEN(JL)-ZDMFDE(JL)
               IF(JK.LT.KCBOT(JL)) THEN
                  ZDPRHO=(PGEOH(JL,JK)-PGEOH(JL,JK+1))/G
                  ZOENTR(JL,JK)=ZOENTR(JL,JK)*ZDPRHO*PMFU(JL,JK+1)
                  ZMFTEST=PMFU(JL,JK)+ZOENTR(JL,JK)-ZODETR(JL,JK)
                  ZMFMAX=MIN(ZMFTEST,(PAPH(JL,JK)-PAPH(JL,JK-1))*ZCONS2)
                  ZOENTR(JL,JK)=MAX(ZOENTR(JL,JK)-
     &                 MAX(ZMFTEST-ZMFMAX,0.),0.)
               ENDIF
               IF(KTYPE(JL).EQ.1.AND.JK.LT.KCBOT(JL)
     &                          .AND.JK.LE.KHMIN(JL)) THEN
         ! LIMIT ORGANIZED DETRAINMENT TO NOT ALLOWING FOR TOO
         ! DEEP CLOUDS
                  ZMSE =CPD*PTU(JL,JK+1)+ALV*PQU(JL,JK+1)+PGEOH(JL,JK+1)
                  IKT=KCTOP0(JL)
                  ZNEVN=(PGEOH(JL,IKT)-PGEOH(JL,JK+1))*
     &                 (ZMSE-PHHATT(JL,JK+1))/G
                  IF(ZNEVN.LE.0.) ZNEVN=1.
                  ZDPRHO=(PGEOH(JL,JK)-PGEOH(JL,JK+1))/G
                  ZODMAX=((PHCBASE(JL)-ZMSE)/ZNEVN)*ZDPRHO*PMFU(JL,JK+1)
                  ZODMAX=MAX(ZODMAX,0.)
                  ZODETR(JL,JK)=MIN(ZODETR(JL,JK),ZODMAX)
               ENDIF
               ZODETR(JL,JK)=MIN(ZODETR(JL,JK),0.75*PMFU(JL,JK))
               PMFU(JL,JK)=PMFU(JL,JK)+ZOENTR(JL,JK)-ZODETR(JL,JK)
               ZQEEN=PQENH(JL,JK+1)*ZDMFEN(JL)
               ZQEEN=ZQEEN+
     &              PQENH(JL,JK+1)*ZOENTR(JL,JK)
               ZSEEN=(CPD*PTENH(JL,JK+1)+PGEOH(JL,JK+1))*ZDMFEN(JL)
               ZSEEN=ZSEEN+
     &              (CPD*PTENH(JL,JK+1)+PGEOH(JL,JK+1))*ZOENTR(JL,JK)
               ZSCDE=(CPD*PTU(JL,JK+1)+PGEOH(JL,JK+1))*ZDMFDE(JL)
         ! FIND MOIST STATIC ENERGY THAT GIVE NONBUOYANT AIR
               ZGA=ALV*PQSENH(JL,JK+1)/(RV*(PTENH(JL,JK+1)**2))
               ZDT=(PLU(JL,JK+1)-0.608*(PQSENH(JL,JK+1)-PQENH(JL,JK+1)))
     &              /(1./PTENH(JL,JK+1) + 0.608*ZGA)
               ZSCOD=CPD*PTENH(JL,JK+1)+PGEOH(JL,JK+1)+CPD*ZDT
               ZSCDE=ZSCDE+ZODETR(JL,JK)*ZSCOD
               ZQUDE=PQU(JL,JK+1)*ZDMFDE(JL)
               ZQCOD=PQSENH(JL,JK+1)+ZGA*ZDT
               ZQUDE=ZQUDE+ZODETR(JL,JK)*ZQCOD
               PLUDE(JL,JK)=PLU(JL,JK+1)*ZDMFDE(JL)
               PLUDE(JL,JK)=PLUDE(JL,JK)+PLU(JL,JK+1)*ZODETR(JL,JK)
               ZMFUSK=PMFUS(JL,JK+1)+ZSEEN-ZSCDE
               ZMFUQK=PMFUQ(JL,JK+1)+ZQEEN-ZQUDE
               ZMFULK=PMFUL(JL,JK+1)    -PLUDE(JL,JK)
               PLU(JL,JK)=ZMFULK*(1./MAX(CMFCMIN,PMFU(JL,JK)))
               PQU(JL,JK)=ZMFUQK*(1./MAX(CMFCMIN,PMFU(JL,JK)))
               PTU(JL,JK)=(ZMFUSK*(1./MAX(CMFCMIN,PMFU(JL,JK)))-
     &              PGEOH(JL,JK))*RCPD
               PTU(JL,JK)=MAX(100.,PTU(JL,JK))
               PTU(JL,JK)=MIN(400.,PTU(JL,JK))
               ZQOLD(JL)=PQU(JL,JK)
            ENDIF
         ENDDO
C
C
C                  DO CORRECTIONS FOR MOIST ASCENT
C                  BY ADJUSTING T,Q AND L IN *CEADJTQ*
C                  -----------------------------------
C
         IK=JK
         ICALL=1
         CALL CEADJTQ
     &        (KIDIA, KFDIA,
     &        KLON,     KLEV,     IK,
     &        ZPH,      PTU,      PQU,      LOFLAG,   ICALL)
C
!DIR$ IVDEP
         DO JL=KIDIA,KFDIA
            IF(LOFLAG(JL)) THEN
               IF(PQU(JL,JK).NE.ZQOLD(JL)) THEN
                  KLAB(JL,JK)=2
                  PLU(JL,JK)=PLU(JL,JK)+ZQOLD(JL)-PQU(JL,JK)
                  ZBUO=PTU(JL,JK)*(1.+VTMPC1*PQU(JL,JK)-PLU(JL,JK))-
     &               PTENH(JL,JK)*(1.+VTMPC1*PQENH(JL,JK)-PXENH(JL,JK))
                  IF(KLAB(JL,JK+1).EQ.1) ZBUO=ZBUO+0.5
                  IF(ZBUO.GT.0..AND.PMFU(JL,JK).GE.0.01*PMFUB(JL).AND.
     &                 JK.GE.KCTOP0(JL)) THEN
                     KCTOP(JL)=MAX(JK,4)
                     LDCUM(JL)=.TRUE.
CRP   CVMGT ERSETZT
C           ZDNOPRC=CVMGT(ZDLAND(JL),1.5E4,LDLAND(JL))
C           ZPRCON=CVMGT(0.,CPRCON,ZPBASE(JL)-PAPH(JL,JK).LT.ZDNOPRC)
                     IF (LDLAND(JL)) THEN
                        ZDNOPRC=ZDLAND(JL)
                     ELSE
                        ZDNOPRC=DNOPRC
                     ENDIF
CSP                  HIER WIRD FUER KALTE KONVEKTION (KTYPE=4)
C                    NIEDERSCHLAGSBILDUNG ERLAUBTT
                     IF (LCOLDC) THEN
                        IF ((ZPBASE(JL)-PAPH(JL,JK).LT.ZDNOPRC) .AND. 
     &                       (KTYPE(JL) .NE. 4)) THEN
                           ZPRCON=0.
                        ELSE
                           ZPRCON=CPRCON
CSP TEMPERATURABHAENGIGE KONVERSIONSRATE FUER KALTE KONVEKTION:
                           IF (KTYPE(JL) .EQ. 4) THEN
                              FACCO=1. + 0.5 * 
     &                             SQRT(MAX(268.26-PTU(JL,JK),8.26))
                              IF (PTU(JL,JK) .GT. 268.26) THEN
                                 FACCO=1.
                              ENDIF
                              ZPRCON=CPRCON*FACCO
                           ENDIF
                        ENDIF
CSP
                     ELSE
                        IF (ZPBASE(JL)-PAPH(JL,JK).LT.ZDNOPRC) THEN
                           ZPRCON=0.
                        ELSE
                           ZPRCON=CPRCON
                        ENDIF
                     ENDIF
                     ZLNEW=PLU(JL,JK)/
     &                    (1.+ZPRCON*(PGEOH(JL,JK)-PGEOH(JL,JK+1)))
                     PDMFUP(JL,JK)=
     &                    MAX(0.,(PLU(JL,JK)-ZLNEW)*PMFU(JL,JK))
                     PLU(JL,JK)=ZLNEW
                  ELSE
                     KLAB(JL,JK)=0
                     PMFU(JL,JK)=0.
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         DO JL=KIDIA,KFDIA
            IF(LOFLAG(JL)) THEN
               PMFUL(JL,JK)=PLU(JL,JK)*PMFU(JL,JK)
               PMFUS(JL,JK)=(CPD*PTU(JL,JK)+PGEOH(JL,JK))*PMFU(JL,JK)
               PMFUQ(JL,JK)=PQU(JL,JK)*PMFU(JL,JK)
            ENDIF
         ENDDO
C
         IF(LMFDUDV) THEN
            DO JL=1,KLON
               ZDMFEN(JL)=ZDMFEN(JL)+ZOENTR(JL,JK)
               ZDMFDE(JL)=ZDMFDE(JL)+ZODETR(JL,JK)
            ENDDO
            DO JL=KIDIA,KFDIA
               IF(LOFLAG(JL)) THEN
                  IF (KTYPE(JL).EQ.1.OR.KTYPE(JL).EQ.3) THEN
CRP   CVMGT ERSETZT
C                ZZ=CVMGT(3.,2.,ZDMFEN(JL).EQ.0.)
                     IF (ZDMFEN(JL).EQ.0.) THEN
                        ZZ=3.
                     ELSE
                        ZZ=2.
                     ENDIF
                  ELSE
C                ZZ=CVMGT(1.,0.,ZDMFEN(JL).EQ.0.)
                     IF (ZDMFEN(JL).EQ.0.) THEN
                        ZZ=1.
                     ELSE
                        ZZ=0.
                     ENDIF
                  END IF
                  ZDMFEU=ZDMFEN(JL)+ZZ*ZDMFDE(JL)
                  ZDMFDU=ZDMFDE(JL)+ZZ*ZDMFDE(JL)
                  ZDMFDU=MIN(ZDMFDU,0.75*PMFU(JL,JK+1))
                  ZMFUU(JL)=ZMFUU(JL)+
     &                 ZDMFEU*PUEN(JL,JK)-ZDMFDU*PUU(JL,JK+1)
                  ZMFUV(JL)=ZMFUV(JL)+
     &                 ZDMFEU*PVEN(JL,JK)-ZDMFDU*PVU(JL,JK+1)
                  IF(PMFU(JL,JK).GT.0.) THEN
                     PUU(JL,JK)=ZMFUU(JL)*(1./PMFU(JL,JK))
                     PVU(JL,JK)=ZMFUV(JL)*(1./PMFU(JL,JK))
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
C
C
C
C                  COMPUTE ORGANIZED ENTRAINMENT
C                  FOR USE AT NEXT LEVEL
C                  ------------------------------
C
         DO JL=1,KLON
            IF(LOFLAG(JL).AND.KTYPE(JL).EQ.1) THEN
               ZBUOYZ=G*(PTU(JL,JK)-PTENH(JL,JK))/PTENH(JL,JK) +
     &              G*0.608*(PQU(JL,JK)-PQENH(JL,JK))-G*PLU(JL,JK)
               ZBUOYZ=MAX(ZBUOYZ,0.0)
               ZDZ=(PGEO(JL,JK-1)-PGEO(JL,JK))/G
               ZDRODZ=-ALOG(PTEN(JL,JK-1)/PTEN(JL,JK))/ZDZ
     &              -G/(RD*PTENH(JL,JK))
               ZBUOY(JL)=ZBUOY(JL)+ZBUOYZ*ZDZ
               ZOENTR(JL,JK-1)=ZBUOYZ*0.5/(1.+ZBUOY(JL))
     &              + ZDRODZ
               IF(ZOENTR(JL,JK-1).GT.1.E-3) ZOENTR(JL,JK-1)=1.E-3
               IF(ZOENTR(JL,JK-1).LT.0.0)   ZOENTR(JL,JK-1)=0.
C
            ENDIF
         ENDDO
C
C
      ENDDO
C
C
C----------------------------------------------------------------------
C
C     5.           DETERMINE CONVECTIVE FLUXES ABOVE NON-BUOYANCY LEVEL
C                  ----------------------------------------------------
C                  (NOTE: CLOUD VARIABLES LIKE T,Q AND L ARE NOT
C                         AFFECTED BY DETRAINMENT AND ARE ALREADY KNOWN
C                         FROM PREVIOUS CALCULATIONS ABOVE)
C
      DO JL=KIDIA,KFDIA
         IF(KCTOP(JL).EQ.KLEVM1) LDCUM(JL)=.FALSE.
         KCBOT(JL)=MAX(KCBOT(JL),KCTOP(JL))
      ENDDO
      IS=0
      DO JL=KIDIA,KFDIA
CRP   CVMGT ERSETZT
C     IS=IS+CVMGT(1,0,LDCUM(JL))
         IF (LDCUM(JL)) IS=IS+1
      ENDDO
      IF(IS.NE.0) THEN
!DIR$ IVDEP
         DO JL=KIDIA,KFDIA
            IF(LDCUM(JL)) THEN
               JK=KCTOP(JL)-1
               ZZDMF=CMFCTOP
               ZDMFDE(JL)=(1.-ZZDMF)*PMFU(JL,JK+1)
               PLUDE(JL,JK)=ZDMFDE(JL)*PLU(JL,JK+1)
               PMFU(JL,JK)=PMFU(JL,JK+1)-ZDMFDE(JL)
               ZLNEW=PLU(JL,JK)
               PDMFUP(JL,JK)=MAX(0.,(PLU(JL,JK)-ZLNEW)*PMFU(JL,JK))
               PLU(JL,JK)=ZLNEW
               PMFUS(JL,JK)=(CPD*PTU(JL,JK)+PGEOH(JL,JK))*PMFU(JL,JK)
               PMFUQ(JL,JK)=PQU(JL,JK)*PMFU(JL,JK)
               PMFUL(JL,JK)=PLU(JL,JK)*PMFU(JL,JK)
               PLUDE(JL,JK-1)=PMFUL(JL,JK)
            ENDIF
         ENDDO
C
         IF(LMFDUDV) THEN
!DIR$ IVDEP
            DO JL=KIDIA,KFDIA
               IF(LDCUM(JL)) THEN
                  JK=KCTOP(JL)-1
                  PUU(JL,JK)=PUU(JL,JK+1)
                  PVU(JL,JK)=PVU(JL,JK+1)
               ENDIF
            ENDDO
         ENDIF
C
      ENDIF
C
      RETURN
      END SUBROUTINE CEASC
