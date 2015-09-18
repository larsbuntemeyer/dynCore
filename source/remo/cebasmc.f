      SUBROUTINE CEBASMC
     &    (KIDIA, KFDIA,
     &     KLON,     KLEV,     KLEVM1,   KK,
     &     PTEN,     PQEN,     PQSEN,    PUEN,     PVEN,
     &     PVERV,    PGEO,     PGEOH,    LDCUM,    KTYPE,    KLAB,
     &     PMFU,     PMFUB,    PENTR,    KCBOT,
     &     PTU,      PQU,      PLU,      PUU,      PVU,
     &     PMFUS,    PMFUQ,    PMFUL,    PDMFUP,   PMFUU,    PMFUV)
C
      IMPLICIT NONE
C
C          M.TIEDTKE         E.C.M.W.F.     12/89
C
C          PURPOSE.
C          --------
C          THIS ROUTINE CALCULATES CLOUD BASE VALUES
C          FOR MIDLEVEL CONVECTION
C
C          INTERFACE
C          ---------
C
C          THIS ROUTINE IS CALLED FROM *CEASC*.
C          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
C          IT RETURNS CLOUDBASE VALUES FOR MIDLEVEL CONVECTION
C
C          METHOD.
C          --------
C          S. TIEDTKE (1989)
C
C          EXTERNALS
C          ---------
C          NONE
C
      INCLUDE "COMCON"
      INCLUDE "COMCUMF"
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: KIDIA, KFDIA, KLON, KLEV, KLEVM1, KK
C
      REAL,    INTENT(IN)    ::  PTEN(KLON,KLEV),   PQEN(KLON,KLEV),
     &                           PUEN(KLON,KLEV),   PVEN(KLON,KLEV),
     &                           PQSEN(KLON,KLEV),  PVERV(KLON,KLEV),
     &                           PGEO(KLON,KLEV),   PGEOH(KLON,KLEV)
C
      LOGICAL, INTENT(IN)    ::  LDCUM(KLON)
C
      REAL,    INTENT(INOUT) ::  PTU(KLON,KLEV),    PQU(KLON,KLEV),
     &                           PUU(KLON,KLEV),    PVU(KLON,KLEV),
     &                           PLU(KLON,KLEV),    PMFU(KLON,KLEV),
     &                           PMFUB(KLON),       PENTR(KLON),
     &                           PMFUS(KLON,KLEV),  PMFUQ(KLON,KLEV),
     &                           PMFUL(KLON,KLEV),  PDMFUP(KLON,KLEV),
     &                           PMFUU(KLON),       PMFUV(KLON)
C
      INTEGER, INTENT(INOUT) ::  KTYPE(KLON),       KCBOT(KLON),
     &                           KLAB(KLON,KLEV)
C
C     Local Variables
C
      INTEGER  ::  JL
      REAL     ::  ZZZMB
C
C----------------------------------------------------------------------
C
C*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
C                  -------------------------------------------
C
      IF(LMFMID.AND.KK.LT.KLEVM1.AND.KK.GT.KLEV/2) THEN
!DIR$ IVDEP
         DO JL=KIDIA,KFDIA
            IF(.NOT.LDCUM(JL).AND.KLAB(JL,KK+1).EQ.0.
     &         .AND.PQEN(JL,KK).GT.0.80*PQSEN(JL,KK)) THEN
               PTU(JL,KK+1)=(CPD*PTEN(JL,KK)+PGEO(JL,KK)-PGEOH(JL,KK+1))
     &              *RCPD
               PQU(JL,KK+1)=PQEN(JL,KK)
               PLU(JL,KK+1)=0.
               ZZZMB=MAX(CMFCMIN,-PVERV(JL,KK)/G)
               ZZZMB=MIN(ZZZMB,CMFCMAX)
               PMFUB(JL)=ZZZMB
               PMFU(JL,KK+1)=PMFUB(JL)
               PMFUS(JL,KK+1)=PMFUB(JL)
     &              *(CPD*PTU(JL,KK+1)+PGEOH(JL,KK+1))
               PMFUQ(JL,KK+1)=PMFUB(JL)*PQU(JL,KK+1)
               PMFUL(JL,KK+1)=0.
               PDMFUP(JL,KK+1)=0.
               KCBOT(JL)=KK
               KLAB(JL,KK+1)=1
               KTYPE(JL)=3
               PENTR(JL)=ENTRMID
               IF(LMFDUDV) THEN
                  PUU(JL,KK+1)=PUEN(JL,KK)
                  PVU(JL,KK+1)=PVEN(JL,KK)
                  PMFUU(JL)=PMFUB(JL)*PUU(JL,KK+1)
                  PMFUV(JL)=PMFUB(JL)*PVU(JL,KK+1)
               ENDIF
            ENDIF
         ENDDO
C
      ENDIF
C
      RETURN
      END SUBROUTINE CEBASMC
