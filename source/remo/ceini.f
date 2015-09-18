C
C          SUBROUTINE CEINI
C
C          M.TIEDTKE         E.C.M.W.F.     12/89
C
C          PURPOSE
C          -------
C
C          THIS ROUTINE INTERPOLATES LARGE-SCALE FIELDS OF T,Q ETC.
C          TO HALF LEVELS (I.E. GRID FOR MASSFLUX SCHEME),
C          DETERMINES LEVEL OF MAXIMUM VERTICAL VELOCITY
C          AND INITIALIZES VALUES FOR UPDRAFTS AND DOWNDRAFTS
C
C          INTERFACE
C          ---------
C          THIS ROUTINE IS CALLED FROM *CEMASTR*.
C
C          METHOD.
C          --------
C          FOR EXTRAPOLATION TO HALF LEVELS SEE TIEDTKE(1989)
C
C          EXTERNALS
C          ---------
C          *CEADJTQ* TO SPECIFY QS AT HALF LEVELS
C
      SUBROUTINE CEINI
     &    (KIDIA, KFDIA,
     &     KLON,
     &     KLEV,     KLEVP1,   KLEVM1,
     &     PTEN,     PQEN,     PQSEN,    PXEN, PUEN, PVEN,
     &     PVERV,    PGEO,     PAPH,     PGEOH,
     &     PTENH,    PQENH,    PQSENH,   PXENH,  KLWMIN,
     &     PTU,      PQU,      PTD,      PQD,
     &     PUU,      PVU,      PUD,      PVD,
     &     PMFU,     PMFD,     PMFUS,    PMFDS,
     &     PMFUQ,    PMFDQ,    PDMFUP,   PDMFDP,
     &     PDPMEL,   PLU,      PLUDE,    KLAB)
C
      IMPLICIT NONE
C
      INCLUDE "COMCON"
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN) :: KIDIA, KFDIA,  KLON,
     &                        KLEV,  KLEVP1, KLEVM1
C
      REAL, INTENT(IN)    ::  
     &         PTEN(KLON,KLEV),        PQEN(KLON,KLEV),
     &         PUEN(KLON,KLEV),        PVEN(KLON,KLEV),
     &         PQSEN(KLON,KLEV),       PVERV(KLON,KLEV),
     &         PGEO(KLON,KLEV),        
     &         PAPH(KLON,KLEVP1),      
     &         PXEN(KLON,KLEV)
     
      REAL, INTENT(INOUT) ::
     &         PTENH(KLON,KLEV),       PXENH(KLON,KLEV),
     &         PQENH(KLON,KLEV),       PQSENH(KLON,KLEV),
     &         PGEOH(KLON,KLEV)
C
      REAL, INTENT(INOUT) ::
     &         PTU(KLON,KLEV),         PQU(KLON,KLEV),
     &         PTD(KLON,KLEV),         PQD(KLON,KLEV),
     &         PUU(KLON,KLEV),         PUD(KLON,KLEV),
     &         PVU(KLON,KLEV),         PVD(KLON,KLEV),
     &         PMFU(KLON,KLEV),        PMFD(KLON,KLEV),
     &         PMFUS(KLON,KLEV),       PMFDS(KLON,KLEV),
     &         PMFUQ(KLON,KLEV),       PMFDQ(KLON,KLEV),
     &         PDMFUP(KLON,KLEV),      PDMFDP(KLON,KLEV),
     &         PLU(KLON,KLEV),         PLUDE(KLON,KLEV)
C
      REAL, INTENT(INOUT)    :: PDPMEL(KLON,KLEV)
      INTEGER, INTENT(INOUT) :: KLAB(KLON,KLEV), KLWMIN(KLON)
C
C     Local Variables
C
      REAL     ZWMAX(KLON)
      REAL     ZPH(KLON)
      LOGICAL  LOFLAG(KLON)
      INTEGER :: ICALL,IK,JK,JL
      REAL    :: ZDP,ZZS
C
C
C----------------------------------------------------------------------
C
C*    1.           SPECIFY LARGE SCALE PARAMETERS AT HALF LEVELS
C*                 ADJUST TEMPERATURE FIELDS IF STATICLY UNSTABLE
C*                 FIND LEVEL OF MAXIMUM VERTICAL VELOCITY
C                  ----------------------------------------------
C
      ZDP=0.5
      DO JK=2,KLEV
         DO JL=KIDIA,KFDIA
            PGEOH(JL,JK)=PGEO(JL,JK)+(PGEO(JL,JK-1)-PGEO(JL,JK))*ZDP
            PTENH(JL,JK)=(MAX(CPD*PTEN(JL,JK-1)+PGEO(JL,JK-1),
     &           CPD*PTEN(JL,JK)+PGEO(JL,JK))-PGEOH(JL,JK))*RCPD
            PQSENH(JL,JK)=PQSEN(JL,JK-1)
            ZPH(JL)=PAPH(JL,JK)
            LOFLAG(JL)=.TRUE.
         ENDDO
C
         IK=JK
         ICALL=0
         CALL CEADJTQ
     &        (KIDIA, KFDIA,
     &         KLON,     KLEV,     IK,
     &         ZPH,      PTENH,    PQSENH,   LOFLAG,   ICALL)
C
         DO JL=KIDIA,KFDIA
            PXENH(JL,JK)=(PXEN(JL,JK)+PXEN(JL,JK-1))*ZDP
            PQENH(JL,JK)=MIN(PQEN(JL,JK-1),PQSEN(JL,JK-1))
     &           +(PQSENH(JL,JK)-PQSEN(JL,JK-1))
            PQENH(JL,JK)=MAX(PQENH(JL,JK),0.)
         ENDDO
      ENDDO
C
      DO JL=KIDIA,KFDIA
         PTENH(JL,KLEV)=(CPD*PTEN(JL,KLEV)+PGEO(JL,KLEV)-
     &        PGEOH(JL,KLEV))*RCPD
         PXENH(JL,KLEV)=PXEN(JL,KLEV)
         PQENH(JL,KLEV)=PQEN(JL,KLEV)
         PTENH(JL,1)=PTEN(JL,1)
         PXENH(JL,1)=PXEN(JL,1)
         PQENH(JL,1)=PQEN(JL,1)
         PGEOH(JL,1)=PGEO(JL,1)
         KLWMIN(JL)=KLEV
         ZWMAX(JL)=0.
      ENDDO
C
      DO JK=KLEVM1,2,-1
         DO JL=KIDIA,KFDIA
            ZZS=MAX(CPD*PTENH(JL,JK)+PGEOH(JL,JK),
     &           CPD*PTENH(JL,JK+1)+PGEOH(JL,JK+1))
            PTENH(JL,JK)=(ZZS-PGEOH(JL,JK))*RCPD
         ENDDO
      ENDDO
C
      DO JK=KLEV,3,-1
!DIR$ IVDEP
         DO JL=KIDIA,KFDIA
            IF(PVERV(JL,JK).LT.ZWMAX(JL)) THEN
               ZWMAX(JL)=PVERV(JL,JK)
               KLWMIN(JL)=JK
            ENDIF
         ENDDO
      ENDDO
C
C
C-----------------------------------------------------------------------
C
C*    2.0          INITIALIZE VALUES FOR UPDRAFTS AND DOWNDRAFTS
C*                 ---------------------------------------------
C
      DO JK=1,KLEV
         IK=JK-1
         IF(JK.EQ.1) IK=1
         DO JL=KIDIA,KFDIA
            PTU(JL,JK)=PTENH(JL,JK)
            PTD(JL,JK)=PTENH(JL,JK)
            PQU(JL,JK)=PQENH(JL,JK)
            PQD(JL,JK)=PQENH(JL,JK)
            PLU(JL,JK)=0.
            PUU(JL,JK)=PUEN(JL,IK)
            PUD(JL,JK)=PUEN(JL,IK)
            PVU(JL,JK)=PVEN(JL,IK)
            PVD(JL,JK)=PVEN(JL,IK)
            PMFU(JL,JK)=0.
            PMFD(JL,JK)=0.
            PMFUS(JL,JK)=0.
            PMFDS(JL,JK)=0.
            PMFUQ(JL,JK)=0.
            PMFDQ(JL,JK)=0.
            PDMFUP(JL,JK)=0.
            PDMFDP(JL,JK)=0.
            PDPMEL(JL,JK)=0.
            PLUDE(JL,JK)=0.
            KLAB(JL,JK)=0
         ENDDO
C
      ENDDO
C
      RETURN
      END SUBROUTINE CEINI
