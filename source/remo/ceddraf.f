C
C          SUBROUTINE CEDDRAF
C
C          THIS ROUTINE CALCULATES CUMULUS DOWNDRAFT DESCENT
C
C          M.TIEDTKE         E.C.M.W.F.    12/86 MODIF. 12/89
C
C          PURPOSE.
C          --------
C          TO PRODUCE THE VERTICAL PROFILES FOR CUMULUS DOWNDRAFTS
C          (I.E. T,Q,U AND V AND FLUXES)
C
C          INTERFACE
C          ---------
C
C          THIS ROUTINE IS CALLED FROM *CEMASTR*.
C          INPUT IS T,Q,P,PHI,U,V AT HALF LEVELS.
C          IT RETURNS FLUXES OF S,Q AND EVAPORATION RATE
C          AND U,V AT LEVELS WHERE DOWNDRAFT OCCURS
C
C          METHOD.
C          --------
C          CALCULATE MOIST DESCENT FOR ENTRAINING/DETRAINING PLUME BY
C          A) MOVING AIR DRY-ADIABATICALLY TO NEXT LEVEL BELOW AND
C          B) CORRECTING FOR EVAPORATION TO OBTAIN SATURATED STATE.
C
C          EXTERNALS
C          ---------
C          *CEADJTQ* FOR ADJUSTING T AND Q DUE TO EVAPORATION IN
C          SATURATED DESCENT
C
C          REFERENCE
C          ---------
C          (TIEDTKE,1989)
C
      SUBROUTINE CEDDRAF
     &    (KIDIA, KFDIA,
     &     KLON,     KLEV,     KLEVP1,
     &     PTENH,    PQENH,    PUEN,     PVEN,
     &     PGEOH,    PAPH,     PRFL,
     &     PTD,      PQD,      PUD,      PVD,
     &     PMFD,     PMFDS,    PMFDQ,    PDMFDP,
     &     LDDRAF)
C
      IMPLICIT NONE
C
      INCLUDE "COMCON"
      INCLUDE "COMCUMF"
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: KIDIA, KFDIA,
     &                          KLON, KLEV, KLEVP1
C
      REAL,    INTENT(IN)    :: PTENH(KLON,KLEV), PQENH(KLON,KLEV),
     &                          PUEN(KLON,KLEV),  PVEN(KLON,KLEV),
     &                          PGEOH(KLON,KLEV), PAPH(KLON,KLEVP1)
C
      REAL,    INTENT(INOUT) :: PTD(KLON,KLEV),   PQD(KLON,KLEV),
     &                          PUD(KLON,KLEV),   PVD(KLON,KLEV),
     &                          PMFD(KLON,KLEV),  PMFDS(KLON,KLEV),
     &                          PMFDQ(KLON,KLEV), PDMFDP(KLON,KLEV),
     &                          PRFL(KLON)
      LOGICAL, INTENT(IN)    :: LDDRAF(KLON)
C
C     Local Variables
C
      REAL     :: ZDMFEN(KLON), ZDMFDE(KLON), ZCOND(KLON)
      REAL     :: ZPH(KLON)
      LOGICAL  :: LLO2(KLON)
      LOGICAL  :: LLO1
C
      INTEGER  :: ICALL,IK,IS,ITOPDE,JK,JL
      REAL     :: ZSEEN,ZSDDE,ZQEEN,ZQDDE,ZMFDVK,ZMFDUK,ZMFDSK,ZMFDQK
      REAL     :: ZENTR,ZDMFDP,ZBUO
C
C
C----------------------------------------------------------------------
C
C     1.           CALCULATE MOIST DESCENT FOR CUMULUS DOWNDRAFT BY
C                     (A) CALCULATING ENTRAINMENT RATES, ASSUMING
C                         LINEAR DECREASE OF MASSFLUX IN PBL
C                     (B) DOING MOIST DESCENT - EVAPORATIVE COOLING
C                         AND MOISTENING IS CALCULATED IN *CEADJTQ*
C                     (C) CHECKING FOR NEGATIVE BUOYANCY AND
C                         SPECIFYING FINAL T,Q,U,V AND DOWNWARD FLUXES
C                    -------------------------------------------------
C
      DO JK=3,KLEV
         IS=0
         DO JL=KIDIA,KFDIA
            ZPH(JL)=PAPH(JL,JK)
            LLO2(JL)=LDDRAF(JL).AND.PMFD(JL,JK-1).LT.0.
            IF (LLO2(JL)) IS=IS+1
         ENDDO
         IF(IS.EQ.0) CYCLE
         DO JL=KIDIA,KFDIA
            IF(LLO2(JL)) THEN
               ZENTR=ENTRDD*PMFD(JL,JK-1)*RD*PTENH(JL,JK-1)/
     &              (G*PAPH(JL,JK-1))*(PAPH(JL,JK)-PAPH(JL,JK-1))
               ZDMFEN(JL)=ZENTR
               ZDMFDE(JL)=ZENTR
            ENDIF
         ENDDO
         ITOPDE=KLEV-2
         IF(JK.GT.ITOPDE) THEN
            DO JL=KIDIA,KFDIA
               IF(LLO2(JL)) THEN
                  ZDMFEN(JL)=0.
                  ZDMFDE(JL)=PMFD(JL,ITOPDE)*
     &                 (PAPH(JL,JK)-PAPH(JL,JK-1))/
     &                 (PAPH(JL,KLEVP1)-PAPH(JL,ITOPDE))
               ENDIF
            ENDDO
         ENDIF
C
         DO JL=KIDIA,KFDIA
            IF(LLO2(JL)) THEN
               PMFD(JL,JK)=PMFD(JL,JK-1)+ZDMFEN(JL)-ZDMFDE(JL)
               ZSEEN=(CPD*PTENH(JL,JK-1)+PGEOH(JL,JK-1))*ZDMFEN(JL)
               ZQEEN=PQENH(JL,JK-1)*ZDMFEN(JL)
               ZSDDE=(CPD*PTD(JL,JK-1)+PGEOH(JL,JK-1))*ZDMFDE(JL)
               ZQDDE=PQD(JL,JK-1)*ZDMFDE(JL)
               ZMFDSK=PMFDS(JL,JK-1)+ZSEEN-ZSDDE
               ZMFDQK=PMFDQ(JL,JK-1)+ZQEEN-ZQDDE
               PQD(JL,JK)=ZMFDQK*(1./MIN(-CMFCMIN,PMFD(JL,JK)))
               PTD(JL,JK)=(ZMFDSK*(1./MIN(-CMFCMIN,PMFD(JL,JK)))-
     &              PGEOH(JL,JK))*RCPD
               PTD(JL,JK)=MIN(400.,PTD(JL,JK))
               PTD(JL,JK)=MAX(100.,PTD(JL,JK))
               ZCOND(JL)=PQD(JL,JK)
            ENDIF
         ENDDO
C
         IK=JK
         ICALL=2
         CALL CEADJTQ
     &        (KIDIA, KFDIA,
     &        KLON,     KLEV,     IK,
     &        ZPH,      PTD,      PQD,      LLO2,     ICALL)
C
         DO JL=KIDIA,KFDIA
            IF(LLO2(JL)) THEN
               ZCOND(JL)=ZCOND(JL)-PQD(JL,JK)
               ZBUO=PTD(JL,JK)*(1.+VTMPC1*PQD(JL,JK))-
     &              PTENH(JL,JK)*(1.+VTMPC1*PQENH(JL,JK))
               LLO1=ZBUO.LT.0..AND.
     &              (PRFL(JL)-PMFD(JL,JK)*ZCOND(JL).GT.0.)
               IF (.NOT.LLO1) PMFD(JL,JK)=0.
               PMFDS(JL,JK)=(CPD*PTD(JL,JK)+PGEOH(JL,JK))*PMFD(JL,JK)
               PMFDQ(JL,JK)=PQD(JL,JK)*PMFD(JL,JK)
               ZDMFDP=-PMFD(JL,JK)*ZCOND(JL)
               PDMFDP(JL,JK-1)=ZDMFDP
               PRFL(JL)=PRFL(JL)+ZDMFDP
            ENDIF
         ENDDO
C
         IF(LMFDUDV) THEN
            DO JL=KIDIA,KFDIA
               IF(LLO2(JL).AND.PMFD(JL,JK).LT.0.) THEN
                  ZMFDUK=PMFD(JL,JK-1)*PUD(JL,JK-1)+
     &                 ZDMFEN(JL)*PUEN(JL,JK-1)-ZDMFDE(JL)*PUD(JL,JK-1)
                  ZMFDVK=PMFD(JL,JK-1)*PVD(JL,JK-1)+
     &                 ZDMFEN(JL)*PVEN(JL,JK-1)-ZDMFDE(JL)*PVD(JL,JK-1)
                  PUD(JL,JK)=ZMFDUK*(1./MIN(-CMFCMIN,PMFD(JL,JK)))
                  PVD(JL,JK)=ZMFDVK*(1./MIN(-CMFCMIN,PMFD(JL,JK)))
               ENDIF
            ENDDO
         ENDIF
C
      ENDDO
C
      RETURN
      END SUBROUTINE CEDDRAF
