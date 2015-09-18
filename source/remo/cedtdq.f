C      SUBROUTINE CEDTDQ
C
C**** *CEDTDQ* - UPDATES T AND Q TENDENCIES, PRECIPITATION RATES
C                DOES GLOBAL DIAGNOSTICS
C
C          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
C
C**   INTERFACE.
C     ----------
C
C          *CEDTDQ* IS CALLED FROM *CEMASTR*
C
      SUBROUTINE CEDTDQ
     &    (KIDIA, KFDIA,
     &     KLON,     KLEV,     KLEVP1,
     &     TWODT,
     &     KTOPM2,   PAPH,
     &     LDCUM,    PTEN,     PTTE,     PQTE,
     &     PXTEC,
     &     PMFUS,    PMFDS,    PMFUQ,    PMFDQ,
     &     PMFUL,    PDMFUP,   PDMFDP,   PLUDE,
     &     PDPMEL,   PRAIN,    PRFL,     PSFL,
     &     PSRAIN,   PSEVAP,   PSHEAT,   PSMELT,
     &     PRSFC,    PSSFC,    PAPRC,    PAPRS)
C
      IMPLICIT NONE
C
      INCLUDE "COMCON"
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: KIDIA, KFDIA,
     &                          KLON,  KLEV,  KLEVP1, KTOPM2
C
      REAL,    INTENT(IN)    :: TWODT
C
      REAL,    INTENT(IN)    :: PTEN(KLON,KLEV),   PAPH(KLON,KLEVP1),
     &                          PMFUS(KLON,KLEV),  PMFDS(KLON,KLEV),
     &                          PMFUQ(KLON,KLEV),  PMFDQ(KLON,KLEV),
     &                          PMFUL(KLON,KLEV),  PLUDE(KLON,KLEV),
     &                          PDMFUP(KLON,KLEV), PDMFDP(KLON,KLEV),
     &                          PRFL(KLON),        PRAIN(KLON),
     &                          PDPMEL(KLON,KLEV), PSFL(KLON)
C
      LOGICAL, INTENT(IN)    :: LDCUM(KLON)
C
      REAL,    INTENT(INOUT) :: PTTE(KLON,KLEV),  PQTE(KLON,KLEV),
     &                          PAPRC(KLON),      PAPRS(KLON),
     &                          PRSFC(KLON),      PSSFC(KLON),
     &                          PXTEC(KLON,KLEV)
C
      REAL,    INTENT(INOUT) :: PSRAIN, PSEVAP, PSHEAT, PSMELT
C
C     Local Variables
C
      REAL    ::  ZMELT(KLON), ZSHEAT(KLON)
      INTEGER ::  JL, JK
      REAL    ::  ZDTDT, ZDQDT, ZDIAGT, ZALV, ZDIAGW
      LOGICAL ::  LLO1
C
C
C----------------------------------------------------------------------
C
C*    1.0          SPECIFY PARAMETERS
C                  ------------------
C
      ZDIAGT=0.5*TWODT
      ZDIAGW=ZDIAGT/RHOH2O
C
C
C----------------------------------------------------------------------
C
C*    2.0          INCREMENTATION OF T AND Q TENDENCIES
C                  ------------------------------------
C
      DO JL=KIDIA,KFDIA
         ZMELT(JL)=0.
         ZSHEAT(JL)=0.
      ENDDO
C
      DO JK=KTOPM2,KLEV
C
         IF(JK.LT.KLEV) THEN
            DO JL=KIDIA,KFDIA
               IF(LDCUM(JL)) THEN
                  LLO1=(PTEN(JL,JK)-TMELT).GT.0.
                  IF (LLO1) THEN
                     ZALV=ALV
                  ELSE
                     ZALV=ALS
                  ENDIF
                  ZDTDT=(G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*RCPD*
     &                 (PMFUS(JL,JK+1)-PMFUS(JL,JK)+
     &                 PMFDS(JL,JK+1)-PMFDS(JL,JK)
     &                 -ALF*PDPMEL(JL,JK)
     &                 -ZALV*(PMFUL(JL,JK+1)-PMFUL(JL,JK)-
     &                 PLUDE(JL,JK)-
     &                 (PDMFUP(JL,JK)+PDMFDP(JL,JK))))
                  PTTE(JL,JK)=PTTE(JL,JK)+ZDTDT
                  ZDQDT=(G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*
     &                 (PMFUQ(JL,JK+1)-PMFUQ(JL,JK)+
     &                 PMFDQ(JL,JK+1)-PMFDQ(JL,JK)+
     &                 PMFUL(JL,JK+1)-PMFUL(JL,JK)-
     &                 PLUDE(JL,JK)-
     &                 (PDMFUP(JL,JK)+PDMFDP(JL,JK)))
                  PQTE(JL,JK)=PQTE(JL,JK)+ZDQDT
                  PXTEC(JL,JK)=(G/(PAPH(JL,JK+1)-PAPH(JL,JK)))
     &                 *PLUDE(JL,JK)
                  ZSHEAT(JL)=ZSHEAT(JL)+ZALV*(PDMFUP(JL,JK)
     &                 +PDMFDP(JL,JK))
                  ZMELT(JL)=ZMELT(JL)+PDPMEL(JL,JK)
               ENDIF
            ENDDO
C
         ELSE
            DO JL=KIDIA,KFDIA
               IF(LDCUM(JL)) THEN
                  LLO1=(PTEN(JL,JK)-TMELT).GT.0.
                  IF (LLO1) THEN
                     ZALV=ALV
                  ELSE
                     ZALV=ALS
                  ENDIF
                  ZDTDT=-(G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*RCPD*
     &                 (PMFUS(JL,JK)+PMFDS(JL,JK)+ALF*PDPMEL(JL,JK)-ZALV
     &                 *(PMFUL(JL,JK)+PDMFUP(JL,JK)
     &                 +PDMFDP(JL,JK)+PLUDE(JL,JK)))
                  PTTE(JL,JK)=PTTE(JL,JK)+ZDTDT
                  ZDQDT=-(G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*
     &                 (PMFUQ(JL,JK)+PMFDQ(JL,JK)+
     &                 PLUDE(JL,JK)+
     &                 (PMFUL(JL,JK)+PDMFUP(JL,JK)+PDMFDP(JL,JK)))
                  PQTE(JL,JK)=PQTE(JL,JK)+ZDQDT
                  PXTEC(JL,JK)=(G/(PAPH(JL,JK+1)-PAPH(JL,JK)))
     &                 *PLUDE(JL,JK)
                  ZSHEAT(JL)=ZSHEAT(JL)+ZALV*(PDMFUP(JL,JK)
     &                 +PDMFDP(JL,JK))
                  ZMELT(JL)=ZMELT(JL)+PDPMEL(JL,JK)
               ENDIF
            ENDDO
C
         ENDIF
C
      ENDDO
C
C
C---------------------------------------------------------------------
C
C      3.          UPDATE SURFACE FIELDS AND DO GLOBAL BUDGETS
C                  -------------------------------------------
C
      DO JL=KIDIA,KFDIA
         PRSFC(JL)=PRFL(JL)
         PSSFC(JL)=PSFL(JL)
         PAPRC(JL)=PAPRC(JL)+ZDIAGW*(PRFL(JL)+PSFL(JL))
         PAPRS(JL)=PAPRS(JL)+ZDIAGW*PSFL(JL)
         PSHEAT=PSHEAT+ZSHEAT(JL)
         PSRAIN=PSRAIN+PRAIN(JL)
         PSEVAP=PSEVAP-(PRFL(JL)+PSFL(JL))
         PSMELT=PSMELT+ZMELT(JL)
      ENDDO
C
      PSEVAP=PSEVAP+PSRAIN
C
C
      RETURN
      END SUBROUTINE CEDTDQ
