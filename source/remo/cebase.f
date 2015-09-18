      SUBROUTINE CEBASE
     &    (KIDIA, KFDIA,
     &     KLON,     KLEV,     KLEVP1,   KLEVM1,
     &     PTENH,    PQENH,    PGEOH,    PAPH,
     &     PTU,      PQU,      PLU,
     &     PUEN,     PVEN,     PUU,      PVU,
     &     LDCUM,    KCBOT,    KLAB)
C
      IMPLICIT NONE
C
C          THIS ROUTINE CALCULATES CLOUD BASE VALUES (T AND Q)
C          FOR CUMULUS PARAMETERIZATION
C
C          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
C
C          PURPOSE.
C          --------
C          TO PRODUCE CLOUD BASE VALUES FOR CU-PARAMETRIZATION
C
C          INTERFACE
C          ---------
C          THIS ROUTINE IS CALLED FROM *CEMASTR*.
C          INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
C          IT RETURNS CLOUD BASE VALUES AND FLAGS AS FOLLOWS;
C                 KLAB=1 FOR SUBCLOUD LEVELS
C                 KLAB=2 FOR CONDENSATION LEVEL
C
C          METHOD.
C          --------
C          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
C          (NON ENTRAINING PLUME,I.E.CONSTANT MASSFLUX)
C
C          EXTERNALS
C          ---------
C          *CEADJTQ* FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT
C
      INCLUDE "COMCON"
      INCLUDE "COMCUMF"
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: KIDIA, KFDIA, KLON, KLEV, KLEVP1, KLEVM1
C
      REAL,    INTENT(IN)    :: PTENH(KLON,KLEV), PQENH(KLON,KLEV),
     &                          PGEOH(KLON,KLEV), PAPH(KLON,KLEVP1)
C
      REAL,    INTENT(IN)    :: PUEN(KLON,KLEV),  PVEN(KLON,KLEV)
      REAL,    INTENT(INOUT) :: PTU(KLON,KLEV),   PQU(KLON,KLEV),
     &                          PLU(KLON,KLEV)
      REAL,    INTENT(INOUT) :: PUU(KLON,KLEV),   PVU(KLON,KLEV)
      INTEGER, INTENT(INOUT) :: KLAB(KLON,KLEV),  KCBOT(KLON)
      LOGICAL, INTENT(INOUT) :: LDCUM(KLON)
C
C     Local Variables
C
      REAL,    DIMENSION(KLON) :: ZQOLD, ZPH
      LOGICAL, DIMENSION(KLON) :: LOFLAG
      INTEGER                  :: ICALL,IK,IKB,IS,JK,JL
      REAL                     :: ZZ,ZBUO
C
C
C----------------------------------------------------------------------
C
C     1.           INITIALIZE VALUES AT LIFTING LEVEL
C                  ----------------------------------
C
      DO JL=KIDIA,KFDIA
         KLAB(JL,KLEV)=1
         KCBOT(JL)=KLEVM1
         LDCUM(JL)=.FALSE.
         PUU(JL,KLEV)=PUEN(JL,KLEV)*(PAPH(JL,KLEVP1)-PAPH(JL,KLEV))
         PVU(JL,KLEV)=PVEN(JL,KLEV)*(PAPH(JL,KLEVP1)-PAPH(JL,KLEV))
      ENDDO
C
C
C----------------------------------------------------------------------
C
C     2.0          DO ASCENT IN SUBCLOUD LAYER,
C                  CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
C                  ADJUST T,Q AND L ACCORDINGLY IN *CEADJTQ*,
C                  CHECK FOR BUOYANCY AND SET FLAGS
C                  -------------------------------------
C
      DO JK=KLEVM1,2,-1
         IS=0
         DO JL=KIDIA,KFDIA
            IF (KLAB(JL,JK+1).EQ.1) THEN
               IS=IS+1
               LOFLAG(JL)=.TRUE.
            ELSE
               LOFLAG(JL)=.FALSE.
            ENDIF
            ZPH(JL)=PAPH(JL,JK)
         ENDDO
         IF(IS.EQ.0) CYCLE
         DO JL=KIDIA,KFDIA
            IF(LOFLAG(JL)) THEN
               PQU(JL,JK)=PQU(JL,JK+1)
               PTU(JL,JK)=(CPD*PTU(JL,JK+1)+PGEOH(JL,JK+1)
     &              -PGEOH(JL,JK))*RCPD
               ZBUO=PTU(JL,JK)*(1.+VTMPC1*PQU(JL,JK))-
     &              PTENH(JL,JK)*(1.+VTMPC1*PQENH(JL,JK))+0.5
               IF(ZBUO.GT.0.) KLAB(JL,JK)=1
               ZQOLD(JL)=PQU(JL,JK)
            ENDIF
         ENDDO
C
         IK=JK
         ICALL=1
         CALL CEADJTQ
     &    (KIDIA, KFDIA,
     &     KLON,     KLEV,     IK,
     &     ZPH,      PTU,      PQU,      LOFLAG,   ICALL)
C
!DIR$ IVDEP
         DO JL=KIDIA,KFDIA
            IF(LOFLAG(JL).AND.PQU(JL,JK).NE.ZQOLD(JL)) THEN
               KLAB(JL,JK)=2
               PLU(JL,JK)=PLU(JL,JK)+ZQOLD(JL)-PQU(JL,JK)
               ZBUO=PTU(JL,JK)*(1.+VTMPC1*PQU(JL,JK))-
     &              PTENH(JL,JK)*(1.+VTMPC1*PQENH(JL,JK))+0.5
               IF(ZBUO.GT.0.) THEN
                  KCBOT(JL)=JK
                  LDCUM(JL)=.TRUE.
               ENDIF
            ENDIF
         ENDDO
C
C             CALCULATE AVERAGES OF U AND V FOR SUBCLOUD ARA,.
C             THE VALUES WILL BE USED TO DEFINE CLOUD BASE VALUES.
C
         IF(LMFDUDV) THEN
            DO JL=KIDIA,KFDIA
               IF(JK.GE.KCBOT(JL)) THEN
                  PUU(JL,KLEV)=PUU(JL,KLEV)+
     &                 PUEN(JL,JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))
                  PVU(JL,KLEV)=PVU(JL,KLEV)+
     &                 PVEN(JL,JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))
               ENDIF
            ENDDO
         ENDIF
C
      ENDDO
C
C
      IF(LMFDUDV) THEN
         DO JL=KIDIA,KFDIA
            IF(LDCUM(JL)) THEN
               IKB=KCBOT(JL)
               ZZ=1./(PAPH(JL,KLEVP1)-PAPH(JL,IKB))
               PUU(JL,KLEV)=PUU(JL,KLEV)*ZZ
               PVU(JL,KLEV)=PVU(JL,KLEV)*ZZ
            ELSE
               PUU(JL,KLEV)=PUEN(JL,KLEVM1)
               PVU(JL,KLEV)=PVEN(JL,KLEVM1)
            ENDIF
         ENDDO
      ENDIF
C
      RETURN
      END SUBROUTINE CEBASE
