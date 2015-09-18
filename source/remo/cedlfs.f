C
C          SUBROUTINE CEDLFS
C
C          THIS ROUTINE CALCULATES LEVEL OF FREE SINKING FOR
C          CUMULUS DOWNDRAFTS AND SPECIFIES T,Q,U AND V VALUES
C
C          M.TIEDTKE         E.C.M.W.F.    12/86 MODIF. 12/89
C
C          PURPOSE.
C          --------
C          TO PRODUCE LFS-VALUES FOR CUMULUS DOWNDRAFTS
C          FOR MASSFLUX CUMULUS PARAMETERIZATION
C
C          INTERFACE
C          ---------
C          THIS ROUTINE IS CALLED FROM *CEMASTR*.
C          INPUT ARE ENVIRONMENTAL VALUES OF T,Q,U,V,P,PHI
C          AND UPDRAFT VALUES T,Q,U AND V AND ALSO
C          CLOUD BASE MASSFLUX AND CU-PRECIPITATION RATE.
C          IT RETURNS T,Q,U AND V VALUES AND MASSFLUX AT LFS.
C
C          METHOD.
C          --------
C          CHECK FOR NEGATIVE BUOYANCY OF AIR OF EQUAL PARTS OF
C          MOIST ENVIRONMENTAL AIR AND CLOUD AIR.
C
C          EXTERNALS
C          ---------
C          *CEADJTQ* FOR CALCULATING WET BULB T AND Q AT LFS
C
C
      SUBROUTINE CEDLFS
     &    (KIDIA,    KFDIA,
     &     KLON,     KLEV,     KLEVP1,
     &     PTENH,    PQENH,    PUEN,     PVEN,
     &     PGEOH,    PAPH,
     &     PTU,      PQU,      PUU,      PVU,
     &     LDCUM,    KCBOT,    KCTOP,    PMFUB,    PRFL,
     &     PTD,      PQD,      PUD,      PVD,
     &     PMFD,     PMFDS,    PMFDQ,    PDMFDP,
     &     KDTOP,    LDDRAF)
C
      IMPLICIT NONE
C
      INCLUDE "COMCON"
      INCLUDE "COMCUMF"
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN)     :: KIDIA, KFDIA,
     &                           KLON,  KLEV,  KLEVP1
C
      REAL,    INTENT(IN)     :: PTENH(KLON,KLEV), PQENH(KLON,KLEV),
     &                           PUEN(KLON,KLEV),  PVEN(KLON,KLEV),
     &                           PGEOH(KLON,KLEV), PAPH(KLON,KLEVP1),
     &                           PTU(KLON,KLEV),   PQU(KLON,KLEV),
     &                           PUU(KLON,KLEV),   PVU(KLON,KLEV),
     &                           PMFUB(KLON)          
C
      REAL,    INTENT(INOUT)  :: PTD(KLON,KLEV),   PQD(KLON,KLEV),
     &                           PUD(KLON,KLEV),   PVD(KLON,KLEV),
     &                           PMFD(KLON,KLEV),  PMFDS(KLON,KLEV),
     &                           PMFDQ(KLON,KLEV), PDMFDP(KLON,KLEV), 
     &                           PRFL(KLON)
C
      INTEGER, INTENT(IN)     :: KCBOT(KLON),      KCTOP(KLON)
      INTEGER, INTENT(INOUT)  :: KDTOP(KLON)
      LOGICAL, INTENT(IN)     :: LDCUM(KLON)
      LOGICAL, INTENT(INOUT)  :: LDDRAF(KLON)
C
C     Local Variables
C
      REAL    ::   ZTENWB(KLON,KLEV), ZQENWB(KLON,KLEV)
      REAL    ::   ZCOND(KLON),       ZPH(KLON)
      LOGICAL ::   LLO2(KLON)
      REAL    ::   ZTTEST,ZQTEST,ZBUO,ZMFTOP
      INTEGER ::   KE,JL,JK,IS,IK,ICALL
C
C
C----------------------------------------------------------------------
C
C     1.           SET DEFAULT VALUES FOR DOWNDRAFTS
C                  ---------------------------------
C
      DO JL=KIDIA,KFDIA
         LDDRAF(JL)=.FALSE.
         KDTOP(JL)=KLEVP1
      ENDDO
C
      IF(LMFDD) THEN
C
C
C----------------------------------------------------------------------
C
C     2.           DETERMINE LEVEL OF FREE SINKING BY
C                  DOING A SCAN FROM TOP TO BASE OF CUMULUS CLOUDS
C
C                  FOR EVERY POINT AND PROCEED AS FOLLOWS:
C
C                    (1) DETEMINE WET BULB ENVIRONMENTAL T AND Q
C                    (2) DO MIXING WITH CUMULUS CLOUD AIR
C                    (3) CHECK FOR NEGATIVE BUOYANCY
C
C                  THE ASSUMPTION IS THAT AIR OF DOWNDRAFTS IS MIXTURE
C                  OF 50% CLOUD AIR + 50% ENVIRONMENTAL AIR AT WET BULB
C                  TEMPERATURE (I.E. WHICH BECAME SATURATED DUE TO
C                  EVAPORATION OF RAIN AND CLOUD WATER)
C                  ----------------------------------------------------
C
         KE=KLEV-3
         DO JK=3,KE
C
C
C     2.1          CALCULATE WET-BULB TEMPERATURE AND MOISTURE
C                  FOR ENVIRONMENTAL AIR IN *CEADJTQ*
C                  -------------------------------------------
C
            IS=0
            DO JL=KIDIA,KFDIA
               ZTENWB(JL,JK)=PTENH(JL,JK)
               ZQENWB(JL,JK)=PQENH(JL,JK)
               ZPH(JL)=PAPH(JL,JK)
               LLO2(JL)=LDCUM(JL).AND.PRFL(JL).GT.0.
     &              .AND..NOT.LDDRAF(JL).AND.
     &              (JK.LT.KCBOT(JL).AND.JK.GT.KCTOP(JL))
               IF (LLO2(JL)) IS=IS+1
            ENDDO
            IF(IS.EQ.0) CYCLE
C
            IK=JK
            ICALL=2
            CALL CEADJTQ
     &           (KIDIA, KFDIA,
     &            KLON,     KLEV,     IK,
     &            ZPH,      ZTENWB,   ZQENWB,   LLO2,     ICALL)
C
C
C     2.2          DO MIXING OF CUMULUS AND ENVIRONMENTAL AIR
C                  AND CHECK FOR NEGATIVE BUOYANCY.
C                  THEN SET VALUES FOR DOWNDRAFT AT LFS.
C                  ----------------------------------------
C
!DIR$ IVDEP
            DO JL=KIDIA,KFDIA
               IF(LLO2(JL)) THEN
                  ZTTEST=0.5*(PTU(JL,JK)+ZTENWB(JL,JK))
                  ZQTEST=0.5*(PQU(JL,JK)+ZQENWB(JL,JK))
                  ZBUO=ZTTEST*(1.+VTMPC1*ZQTEST)-
     &                 PTENH(JL,JK)*(1.+VTMPC1*PQENH(JL,JK))
                  ZCOND(JL)=PQENH(JL,JK)-ZQENWB(JL,JK)
                  ZMFTOP=-CMFDEPS*PMFUB(JL)
                  IF(ZBUO.LT.0.
     &               .AND.PRFL(JL).GT.10.*ZMFTOP*ZCOND(JL)) THEN
                     KDTOP(JL)=JK
                     LDDRAF(JL)=.TRUE.
                     PTD(JL,JK)=ZTTEST
                     PQD(JL,JK)=ZQTEST
                     PMFD(JL,JK)=ZMFTOP
                     PMFDS(JL,JK)=PMFD(JL,JK)*(CPD*PTD(JL,JK)
     &                    +PGEOH(JL,JK))
                     PMFDQ(JL,JK)=PMFD(JL,JK)*PQD(JL,JK)
                     PDMFDP(JL,JK-1)=-0.5*PMFD(JL,JK)*ZCOND(JL)
                     PRFL(JL)=PRFL(JL)+PDMFDP(JL,JK-1)
                  ENDIF
               ENDIF
            ENDDO
C
            IF(LMFDUDV) THEN
               DO JL=KIDIA,KFDIA
                  IF(PMFD(JL,JK).LT.0.) THEN
                     PUD(JL,JK)=0.5*(PUU(JL,JK)+PUEN(JL,JK-1))
                     PVD(JL,JK)=0.5*(PVU(JL,JK)+PVEN(JL,JK-1))
                  ENDIF
               ENDDO
            ENDIF
C
         ENDDO
C
      ENDIF !IF(LMFDD) THEN
      RETURN
      END SUBROUTINE CEDLFS
