C
C           SUBROUTINE CECALL
C
C          *CECALL* - MASTER ROUTINE - PROVIDES INTERFACE FOR:
C                     *CEMASTR* (CUMULUS PARAMETERIZATION)
C                     *CUCCDIA* (CUMULUS CLOUDS FOR RADIATION)
C                     *STRATCU* (PBL_STRATOCUMULUS)
C
C           M.TIEDTKE      E.C.M.W.F.     12/1989
C
C**   PURPOSE.
C     --------
C
C          *CECALL* - INTERFACE FOR *CEMASTR*,*CUCCDIA* AND *STRATCU*:
C                     PROVIDES INPUT FOR CEMASTR, CUCCDIA AND STRATCU.
C                     RECEIVES UPDATED TENDENCIES, PRECIPITATION
C                     AND CONVECTIVE CLOUD PARAMETERS FOR RADIATION.
C
C**   INTERFACE.
C     ----------
C
C          *CECALL* IS CALLED FROM *PHECHAM*
C
C     EXTERNALS.
C     ----------
C
C          CEMASTR
C          CUCCDIA
C          STRATCU
C
      SUBROUTINE CECALL
     &   (KIDIA , KFDIA  , KLON  , KLEV    , KLEVP1,
     &    KLEVM1, ILAB   , NSTEP , NSTART  , TWODT ,
     &    PTM1  , PQM1   , PUM1  , PVM1    , PXM1  ,
     &    PTTE  , PQTE   , PVOM  , PVOL    , PXTE  ,
     &    PVERV , PQHFL  , PXTEC , PCEVAPCU,
     &    PAP   , PAPH   , PGEO  , LDLAND  ,
     &    PRSFC , PSSFC  , PAPRC , PAPRS   ,
     &    KTYPE , PTOPMAX, PSRAIN, PSEVAP  ,
     &    PSHEAT, PSDISS , PSMELT, PCAPE   , PXIM1 ,
     &    PXITE , PTSM1M)
C
      IMPLICIT NONE

      INCLUDE "COMCON"
      INCLUDE "YOTLUC"
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: KIDIA, KFDIA , KLON, KLEV, KLEVP1,
     &                          KLEVM1, NSTEP, NSTART
C
      REAL,    INTENT(IN)    :: TWODT
C
      REAL,    INTENT(IN)    :: 
     &     PTM1(KLON,KLEV) , PQM1(KLON,KLEV), PUM1(KLON,KLEV)  ,
     &     PVM1(KLON,KLEV) , PTTE(KLON,KLEV), PQTE(KLON,KLEV)  ,
     &     PVOM(KLON,KLEV) , PVOL(KLON,KLEV), PVERV(KLON,KLEV) ,
     &     PGEO(KLON,KLEV) , PAP(KLON,KLEV) , PAPH(KLON,KLEVP1),
     &     PQHFL(KLON)     , PAPRC(KLON)    , PAPRS(KLON)      ,
     &     PRSFC(KLON)     , PSSFC(KLON)    , 
     &     PXTEC(KLON,KLEV), PCEVAPCU(KLEV) , PXM1(KLON,KLEV)  ,
     &     PXTE(KLON,KLEV) , PCAPE(KLON)
C
      INTEGER, INTENT(IN)    :: KTYPE(KLON)  ,   ILAB(KLON,KLEV)
C
      LOGICAL, INTENT(IN)    :: LDLAND(KLON)

      
      REAL,    INTENT(IN)    :: PXIM1(KLON,KLEV), PXITE(KLON,KLEV),
     &                          PTSM1M(KLON)
C
      REAL,    INTENT(INOUT) :: PTOPMAX(KLON)
C
C     Local Variables
C
      REAL ZTP1(KLON,KLEV), ZQP1(KLON,KLEV), ZXP1(KLON,KLEV) ,
     &     ZUP1(KLON,KLEV), ZVP1(KLON,KLEV), ZTU(KLON,KLEV)  ,
     &     ZQU(KLON,KLEV) , ZLU(KLON,KLEV) , ZLUDE(KLON,KLEV),
     &     ZMFU(KLON,KLEV), ZMFD(KLON,KLEV), ZQSAT(KLON,KLEV),
     &     ZRAIN(KLON),     ZXIP1(KLON,KLEV)
C 
      INTEGER :: ITOPEC2(KLON), ICBOT(KLON), ICTOP(KLON)
      LOGICAL :: LOCUM(KLON)
CSP
      LOGICAL :: LCOLDCONV(KLON)
C
      INTEGER :: JL,JK,IT
      REAL    :: ZTOPMAX, ZTMST, ZCOLDC, PSRAIN
      REAL    :: PSEVAP, PSHEAT, PSDISS, PSMELT
CSP
CSP   CALCULATE TEMPERATURE DIFFERENCE LOWEST LEVEL AND SURFACE
CSP   FOR USE IN COLD CONVECTION PARAMETERIZATION (KTYPE = 4)
C
      DO JL=1,KLON
         LCOLDCONV(JL)=.FALSE.
      ENDDO
C
      IF (LCOLDC) THEN
         ZCOLDC=10.
         DO JL=1,KLON
            IF ((PTM1(JL,KLEV) .LT. TMELT-10.) .AND. ((PTSM1M(JL) -
     &           PTM1(JL,KLEV)) .GT. ZCOLDC)) THEN
               LCOLDCONV(JL)=.TRUE.
            ENDIF
         ENDDO
      ENDIF
CSP
C
C
C-----------------------------------------------------------------------
C
C*    1.           CALCULATE T,Q AND QS AT MAIN LEVELS
C*                 -----------------------------------
C
      ZTMST=TWODT
      IF(NSTEP.EQ.NSTART) ZTMST=0.5*TWODT
      DO JK=1,KLEV
         DO JL=KIDIA,KFDIA
            ZTP1(JL,JK)=PTM1(JL,JK)+PTTE(JL,JK)*ZTMST
            ZQP1(JL,JK)=PQM1(JL,JK)+PQTE(JL,JK)*ZTMST
            ZXP1(JL,JK)=PXM1(JL,JK)+PXTE(JL,JK)*ZTMST
CSP
            ZXIP1(JL,JK)=PXIM1(JL,JK)+PXITE(JL,JK)*ZTMST
            ZXP1(JL,JK)=ZXP1(JL,JK)+ZXIP1(JL,JK)
CSP
            ZUP1(JL,JK)=PUM1(JL,JK)+PVOM(JL,JK)*ZTMST
            ZVP1(JL,JK)=PVM1(JL,JK)+PVOL(JL,JK)*ZTMST
C     ZTP1(JL,JK)=PTM1(JL,JK)
C     ZQP1(JL,JK)=PQM1(JL,JK)
C     ZXP1(JL,JK)=PXM1(JL,JK)
C     ZUP1(JL,JK)=PUM1(JL,JK)
C     ZVP1(JL,JK)=PVM1(JL,JK)
            IT=INT(ZTP1(JL,JK)*1000.)
            ZQSAT(JL,JK)=TLUCUA(IT)/PAP(JL,JK)
            ZQSAT(JL,JK)=MIN(0.5,ZQSAT(JL,JK))
            ZQSAT(JL,JK)=ZQSAT(JL,JK)/(1.-VTMPC1*ZQSAT(JL,JK))
         ENDDO
      ENDDO
      DO JL=KIDIA,KFDIA
         ZRAIN(JL)=0.
         LOCUM(JL)=.FALSE.
      ENDDO
C
C
C-----------------------------------------------------------------------
C
C*    2.     CALL 'CEMASTR'(MASTER-ROUTINE FOR CUMULUS PARAMETERIZATION)
C*           -----------------------------------------------------------
C
C
CSP
      CALL CEMASTR
     &    (KIDIA , KFDIA , KLON  , KLEV    , KLEVP1,
     &     KLEVM1, ILAB  , NSTEP , NSTART  , TWODT ,
     &     ZTP1  , ZQP1  , ZXP1  , ZUP1    , ZVP1  ,
     &     PVERV , ZQSAT , PQHFL ,
     &     PAP   , PAPH  , PGEO  , LDLAND  ,
     &     PTTE  , PQTE  , PVOM  , PVOL    ,
     &     PRSFC , PSSFC , PAPRC , PAPRS   , PXTEC ,
     &     LOCUM , KTYPE , ICBOT , ICTOP   ,
     &     ZTU   , ZQU   , ZLU   , ZLUDE   ,
     &     ZMFU  , ZMFD  , ZRAIN , PCEVAPCU,
     &     PSRAIN, PSEVAP, PSHEAT, PSDISS  , PSMELT,
     &     PCAPE , LCOLDCONV)
C
C
C ---------------------------------------------------------------------
C
C*    2.1 GEOPOTENTIAL OF CONVECTIVE CLOUD TOPS (KUO0)
C
C
      DO JL=KIDIA,KFDIA
         ITOPEC2(JL)=KLEVP1
      ENDDO
C
      DO JK=1,KLEV-4
         DO JL=KIDIA,KFDIA
            IF(ILAB(JL,JK).EQ.2 .AND. ITOPEC2(JL).EQ.KLEVP1) THEN
               ITOPEC2(JL)=JK
            ENDIF
         ENDDO
      ENDDO
C
      DO JL=KIDIA,KFDIA
         IF(ITOPEC2(JL).EQ.1) THEN
            ZTOPMAX=PAP(JL,1)
         ELSE IF(ITOPEC2(JL).NE.KLEVP1) THEN
            ZTOPMAX=PAPH(JL,ITOPEC2(JL))
         ELSE
            ZTOPMAX=99999.
         END IF
         PTOPMAX(JL)=AMIN1(PTOPMAX(JL),ZTOPMAX)
      ENDDO
C
C
C---------------------------------------------------------------------
C
      RETURN
      END SUBROUTINE CECALL
