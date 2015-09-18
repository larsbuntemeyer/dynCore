C
C          SUBROUTINE CEENTR
C
C          M.TIEDTKE         E.C.M.W.F.     12/89
C
C          PURPOSE.
C          --------
C          THIS ROUTINE CALCULATES ENTRAINMENT/DETRAINMENT RATES
C          FOR UPDRAFTS IN CUMULUS PARAMETERIZATION
C
C          INTERFACE
C          ---------
C
C          THIS ROUTINE IS CALLED FROM *CEASC*.
C          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
C          AND UPDRAFT VALUES T,Q ETC
C          IT RETURNS ENTRAINMENT/DETRAINMENT RATES
C
C          METHOD.
C          --------
C          S. TIEDTKE (1989)
C
C          EXTERNALS
C          ---------
C          NONE
C
      SUBROUTINE CEENTR
     &    (KIDIA, KFDIA,
     &     KLON,     KLEV,     KLEVP1,   KK,
     &     PTENH,    PAPH,     PAP,
     &     KLWMIN,   LDCUM,    KTYPE,    KCBOT,    KCTOP0,
     &     PPBASE,   PMFU,     PENTR,    PODETR,
     &     KHMIN,    PGEOH,
     &     PDMFEN,   PDMFDE)
C
      IMPLICIT NONE
C
      INCLUDE "COMCON"
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: KIDIA, KFDIA,
     &                          KLON, KLEV, KLEVP1, KK
C          
      INTEGER, INTENT(IN)    :: KHMIN (KLON),
     &                          KLWMIN(KLON),      KTYPE(KLON),
     &                          KCBOT(KLON),       KCTOP0(KLON)
C
      REAL,    INTENT(IN)    :: PTENH(KLON,KLEV),
     &                          PAP(KLON,KLEV),    PAPH(KLON,KLEVP1),
     &                          PMFU(KLON,KLEV),
     &                          PENTR(KLON),       PGEOH (KLON,KLEV)    
C 
      LOGICAL, INTENT(IN)    :: LDCUM(KLON)
C      
      REAL,    INTENT(INOUT) :: PODETR(KLON,KLEV), PPBASE(KLON),
     &                          PDMFEN(KLON),      PDMFDE(KLON)
C
C     Local Variables
C
      LOGICAL :: LLO1,LLO2
      REAL    :: ZZMZK,ZTMZK,ZRRHO,ZRG,ZPMID,ZORGDE,ZENTR,ZDPRHO,ARG
      INTEGER :: JL,IKT,IKLWMIN,IKH
C
C
C----------------------------------------------------------------------
C
C*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
C                  -------------------------------------------
C
C
C*    1.1          SPECIFY ENTRAINMENT RATES FOR SHALLOW CLOUDS
C                  --------------------------------------------
C
C
C*    1.2          SPECIFY ENTRAINMENT RATES FOR DEEP CLOUDS
C                  -----------------------------------------
C
      ZRG=1./G
      DO JL=KIDIA,KFDIA
         PPBASE(JL)=PAPH(JL,KCBOT(JL))
         ZRRHO=(RD*PTENH(JL,KK+1))/PAPH(JL,KK+1)
         ZDPRHO=(PAPH(JL,KK+1)-PAPH(JL,KK))*ZRG
         ZPMID=0.5*(PPBASE(JL)+PAPH(JL,KCTOP0(JL)))
         ZENTR=PENTR(JL)*PMFU(JL,KK+1)*ZDPRHO*ZRRHO
         LLO1=KK.LT.KCBOT(JL).AND.LDCUM(JL)
         IF (LLO1) THEN
            PDMFDE(JL)=ZENTR
         ELSE
            PDMFDE(JL)=0.
         ENDIF
         LLO2=LLO1.AND.(KTYPE(JL).EQ.2).AND.
     &        (PPBASE(JL)-PAPH(JL,KK).LT.0.2E5.OR.
     &        PAPH(JL,KK).GT.ZPMID)
         IF (LLO2) THEN
            PDMFEN(JL)=ZENTR
         ELSE
            PDMFEN(JL)=0.
         ENDIF
         IKLWMIN=MAX(KLWMIN(JL),KCTOP0(JL)+2)
         LLO2=LLO1.AND.KTYPE(JL).EQ.3.AND.
     &        (KK.GE.IKLWMIN.OR.PAP(JL,KK).GT.ZPMID)
         IF(LLO2) PDMFEN(JL)=ZENTR
CSP      ADDED '.OR.KTYPE(JL).EQ.4'
         LLO2=LLO1.AND.(KTYPE(JL).EQ.1.OR.KTYPE(JL).EQ.4)
         IF(LLO2) PDMFEN(JL)=ZENTR !TURBULENT ENTRAINMENT
C        NO ORGANIZED DETRAINMENT FOR COLD CONVECTION (KTYPE 4)
         LLO2=LLO1.AND.KTYPE(JL).EQ.1.
C        ORGANIZED DETRAINMENT, DETRAINMENT STARTS AT KHMIN
CSP
         PODETR(JL,KK)=0.
         IF(LLO2.AND.KK.LE.KHMIN(JL).AND.KK.GE.KCTOP0(JL)) THEN
            IKT=KCTOP0(JL)
            IKH=KHMIN(JL)
            IF(IKH.GT.IKT) THEN
               ZZMZK  =-(PGEOH(JL,IKH)-PGEOH(JL,KK))*ZRG
               ZTMZK  =-(PGEOH(JL,IKH)-PGEOH(JL,IKT))*ZRG
               ARG=3.1415*(ZZMZK/ZTMZK)*0.5
               ZORGDE=TAN(ARG)*3.1415*0.5/ZTMZK
               ZDPRHO=(PAPH(JL,KK+1)-PAPH(JL,KK))*(ZRG*ZRRHO)
               PODETR(JL,KK)=MIN(ZORGDE,1.E-3)*PMFU(JL,KK+1)*ZDPRHO
            ENDIF
         ENDIF
      ENDDO
C
      RETURN
      END SUBROUTINE CEENTR
