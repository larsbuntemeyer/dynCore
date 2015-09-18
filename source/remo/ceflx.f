C
C          SUBROUTINE CEFLX
C
C          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
C
C          PURPOSE
C          -------
C
C          THIS ROUTINE DOES THE FINAL CALCULATION OF CONVECTIVE
C          FLUXES IN THE CLOUD LAYER AND IN THE SUBCLOUD LAYER
C
C          INTERFACE
C          ---------
C          THIS ROUTINE IS CALLED FROM *CEMASTR*.
C
C          EXTERNALS
C          ---------
C          NONE
C
      SUBROUTINE CEFLX
     &    (KIDIA, KFDIA,
     &     KLON,     KLEV,     KLEVP1,
     &     NSTEP,    NSTART,   TWODT,
     &     PQEN,     PQSEN,    PTENH,    PQENH,
     &     PAPH,     PGEOH,    PCEVAPCU,
     &     KCBOT,    KCTOP,    KDTOP,
     &     KTYPE,    LDDRAF,   LDCUM,
     &     PMFU,     PMFD,     PMFUS,    PMFDS,
     &     PMFUQ,    PMFDQ,    PMFUL,    PLUDE,
     &     PDMFUP,   PDMFDP,   PRFL,     PRAIN,
     &     PTEN,     PSFL,     PDPMEL,   KTOPM2)
C
      IMPLICIT NONE
C
      INCLUDE "COMCON"
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: KIDIA, KFDIA, KLON, 
     &                          KLEV, KLEVP1, NSTEP, NSTART
      REAL,    INTENT(IN)    :: TWODT
      INTEGER, INTENT(INOUT) :: KTOPM2
C
      REAL,    INTENT(IN)    ::   
     &         PQEN(KLON,KLEV),        PQSEN(KLON,KLEV),
     &         PTENH(KLON,KLEV),       PQENH(KLON,KLEV),
     &         PAPH(KLON,KLEVP1),      PGEOH(KLON,KLEV),
     &         PCEVAPCU(KLEV),         PTEN(KLON,KLEV)
C
      REAL,    INTENT(INOUT) ::    
     &         PMFU(KLON,KLEV),        PMFD(KLON,KLEV),
     &         PMFUS(KLON,KLEV),       PMFDS(KLON,KLEV),
     &         PMFUQ(KLON,KLEV),       PMFDQ(KLON,KLEV),
     &         PDMFUP(KLON,KLEV),      PDMFDP(KLON,KLEV),
     &         PMFUL(KLON,KLEV),       PLUDE(KLON,KLEV),
     &         PRFL(KLON),             PRAIN(KLON),
     &         PDPMEL(KLON,KLEV),      PSFL(KLON)
C 
      INTEGER, INTENT(IN)    ::  
     &         KCBOT(KLON),            KCTOP(KLON),
     &         KDTOP(KLON)
C
      INTEGER, INTENT(INOUT) :: KTYPE(KLON)
C
      LOGICAL, INTENT(IN)    :: LDCUM(KLON)
      LOGICAL, INTENT(INOUT) :: LDDRAF(KLON)
C
C     Local Variables
C
      INTEGER :: IKB,ITOP,JK,JL
      REAL    :: ZPSUBCL(KLON)
      REAL    :: ZCONS1,ZCONS2,ZCUCOV,ZZP,ZTMST,ZSNMLT,ZRSUM,
     &           ZRNEW,ZRMIN,ZRFLN,ZRFL,ZFAC,ZTMELP2,ZDRFL,ZDPEVAP
C
C*             SPECIFY CONSTANTS
C
      ZTMST=TWODT
      IF(NSTEP.EQ.NSTART) ZTMST=0.5*TWODT
      ZCONS1=CPD/(ALF*G*ZTMST)
      ZCONS2=1./(G*ZTMST)
      ZCUCOV=0.05
      ZTMELP2=TMELT+2.
C
C
C*    1.0          DETERMINE FINAL CONVECTIVE FLUXES
C                  ---------------------------------
C
      ITOP=KLEV
      DO JL=KIDIA,KFDIA
         ITOP=MIN(ITOP,KCTOP(JL))
         IF(.NOT.LDCUM(JL).OR.KDTOP(JL).LT.KCTOP(JL)) LDDRAF(JL)=.FALSE.
         IF(.NOT.LDCUM(JL)) KTYPE(JL)=0
      ENDDO
      KTOPM2=ITOP-2
      DO JK=KTOPM2,KLEV
!DIR$ IVDEP
         DO JL=KIDIA,KFDIA
            IF(LDCUM(JL).AND.JK.GE.KCTOP(JL)-1) THEN
               PMFUS(JL,JK)=PMFUS(JL,JK)-PMFU(JL,JK)*
     &              (CPD*PTENH(JL,JK)+PGEOH(JL,JK))
               PMFUQ(JL,JK)=PMFUQ(JL,JK)-PMFU(JL,JK)*PQENH(JL,JK)
               IF(LDDRAF(JL).AND.JK.GE.KDTOP(JL)) THEN
                  PMFDS(JL,JK)=PMFDS(JL,JK)-PMFD(JL,JK)*
     &                 (CPD*PTENH(JL,JK)+PGEOH(JL,JK))
                  PMFDQ(JL,JK)=PMFDQ(JL,JK)-PMFD(JL,JK)*PQENH(JL,JK)
               ELSE
                  PMFD(JL,JK)=0.
                  PMFDS(JL,JK)=0.
                  PMFDQ(JL,JK)=0.
                  PDMFDP(JL,JK-1)=0.
               ENDIF
            ELSE
               PMFU(JL,JK)=0.
               PMFD(JL,JK)=0.
               PMFUS(JL,JK)=0.
               PMFDS(JL,JK)=0.
               PMFUQ(JL,JK)=0.
               PMFDQ(JL,JK)=0.
               PMFUL(JL,JK)=0.
               PDMFUP(JL,JK-1)=0.
               PDMFDP(JL,JK-1)=0.
               PLUDE(JL,JK-1)=0.
            ENDIF
         ENDDO
C
      ENDDO
      DO JK=KTOPM2,KLEV
!DIR$ IVDEP
         DO JL=KIDIA,KFDIA
            IF(LDCUM(JL).AND.JK.GT.KCBOT(JL)) THEN
               IKB=KCBOT(JL)
               ZZP=((PAPH(JL,KLEVP1)-PAPH(JL,JK))/
     &              (PAPH(JL,KLEVP1)-PAPH(JL,IKB)))
               IF (KTYPE(JL).EQ.3) ZZP=ZZP**2
               PMFU(JL,JK)=PMFU(JL,IKB)*ZZP
               PMFUS(JL,JK)=PMFUS(JL,IKB)*ZZP
               PMFUQ(JL,JK)=PMFUQ(JL,IKB)*ZZP
               PMFUL(JL,JK)=PMFUL(JL,IKB)*ZZP
            ENDIF
         ENDDO
C
      ENDDO
C
C
C*    2.            CALCULATE RAIN/SNOW FALL RATES
C*                  CALCULATE MELTING OF SNOW
C*                  CALCULATE EVAPORATION OF PRECIP
C                   -------------------------------
C
      DO JL=KIDIA,KFDIA
        PRFL(JL)=0.
        PSFL(JL)=0.
        PRAIN(JL)=0.
      ENDDO
      DO JK=KTOPM2,KLEV
        DO JL=KIDIA,KFDIA
          IF(LDCUM(JL)) THEN
            PRAIN(JL)=PRAIN(JL)+PDMFUP(JL,JK)
            IF(PTEN(JL,JK).GT.TMELT) THEN
              PRFL(JL)=PRFL(JL)+PDMFUP(JL,JK)+PDMFDP(JL,JK)
              IF(PSFL(JL).GT.0..AND.PTEN(JL,JK).GT.ZTMELP2) THEN
                ZFAC=ZCONS1*(PAPH(JL,JK+1)-PAPH(JL,JK))
                ZSNMLT=MIN(PSFL(JL),ZFAC*(PTEN(JL,JK)-ZTMELP2))
                PDPMEL(JL,JK)=ZSNMLT
                PSFL(JL)=PSFL(JL)-ZSNMLT
                PRFL(JL)=PRFL(JL)+ZSNMLT
              ENDIF
            ELSE
              PSFL(JL)=PSFL(JL)+PDMFUP(JL,JK)+PDMFDP(JL,JK)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      DO JL=KIDIA,KFDIA
        PRFL(JL)=MAX(PRFL(JL),0.)
        PSFL(JL)=MAX(PSFL(JL),0.)
        ZPSUBCL(JL)=PRFL(JL)+PSFL(JL)
      ENDDO
      DO JK=KTOPM2,KLEV
        DO JL=KIDIA,KFDIA
          IF(LDCUM(JL).AND.JK.GE.KCBOT(JL).AND.
     &         ZPSUBCL(JL).GT.1.E-20) THEN
            ZRFL=ZPSUBCL(JL)
            ZRNEW=(MAX(0.,SQRT(ZRFL/ZCUCOV)-
     &           PCEVAPCU(JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))*
     &           MAX(0.,PQSEN(JL,JK)-PQEN(JL,JK))))**2*ZCUCOV
            ZRMIN=ZRFL-ZCUCOV*MAX(0.,0.8*PQSEN(JL,JK)-PQEN(JL,JK))
     &           *ZCONS2*(PAPH(JL,JK+1)-PAPH(JL,JK))
            ZRNEW=MAX(ZRNEW,ZRMIN)
            ZRFLN=MAX(ZRNEW,0.)
            ZDRFL=MIN(0.,ZRFLN-ZRFL)
            PDMFUP(JL,JK)=PDMFUP(JL,JK)+ZDRFL
            ZPSUBCL(JL)=ZRFLN
          ENDIF
        ENDDO
      ENDDO
      DO JL=KIDIA,KFDIA
        ZRSUM=PRFL(JL)+PSFL(JL)
        ZDPEVAP=ZPSUBCL(JL)-ZRSUM
        PRFL(JL)=PRFL(JL)+ZDPEVAP*PRFL(JL)*
     &       (1./MAX(1.E-20,ZRSUM))
        PSFL(JL)=PSFL(JL)+ZDPEVAP*PSFL(JL)*
     &       (1./MAX(1.E-20,ZRSUM))
      ENDDO
C
      RETURN
      END SUBROUTINE CEFLX
