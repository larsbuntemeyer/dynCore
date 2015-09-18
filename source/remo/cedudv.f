C
C      SUBROUTINE CEDUDV
C
C**** *CEDUDV* - UPDATES U AND V TENDENCIES,
C                DOES GLOBAL DIAGNOSTIC OF DISSIPATION
C
C          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
C
C     THE CRAY INTRINSIC 'SSUM' HAS BEEN REPLACED BY THE FORTRAN90
C     INTRINSIC 'SUM', INCLUDING A LOGICAL MASK.
C
C          T.DIEHL           DKRZ           09/98
C
C**   INTERFACE.
C     ----------
C
C          *CEDUDV* IS CALLED FROM *CEMASTR*
C
      SUBROUTINE CEDUDV
     &    (KIDIA,    KFDIA,
     &     KLON,     KLEV,     KLEVP1,
     &     KTOPM2,   KTYPE,    KCBOT,    PAPH,     LDCUM,
     &     PUEN,     PVEN,     PVOM,     PVOL,
     &     PUU,      PUD,      PVU,      PVD,
     &     PMFU,     PMFD,     PSDISS)
C
      IMPLICIT NONE
C
      INCLUDE "COMCON"
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN) :: KIDIA,    KFDIA,
     &                       KLON,     KLEV,     KLEVP1,
     &                       KTOPM2
C
      INTEGER, INTENT(IN) :: KTYPE(KLON), KCBOT(KLON)
C
      REAL, INTENT(IN)    :: PUEN(KLON,KLEV), PVEN(KLON,KLEV),
     &                       PAPH(KLON,KLEVP1),
     &                       PUU(KLON,KLEV),         PUD(KLON,KLEV),
     &                       PVU(KLON,KLEV),         PVD(KLON,KLEV),
     &                       PMFU(KLON,KLEV),        PMFD(KLON,KLEV)
C
      LOGICAL, INTENT(IN) :: LDCUM(KLON)
C
      REAL, INTENT(INOUT) :: PVOL(KLON,KLEV), PVOM(KLON,KLEV)
C
      REAL, INTENT(INOUT) :: PSDISS
C
C     Local Variables
C
      REAL    ::  ZMFUU(KLON,KLEV),       ZMFDU(KLON,KLEV),
     &            ZMFUV(KLON,KLEV),       ZMFDV(KLON,KLEV),
     &            ZDISS(KLON)
      REAL    ::  ZZP,ZSUM,ZDVDT,ZDUDT
      INTEGER ::  IK,IKB,JK,JL
C
C
C----------------------------------------------------------------------
C
C*    1.0          CALCULATE FLUXES AND UPDATE U AND V TENDENCIES
C                  ----------------------------------------------
C
      IF (KTOPM2.EQ.1) THEN
         DO JK=2,KLEV
            IK=JK-1
            DO JL=KIDIA,KFDIA
               IF(LDCUM(JL)) THEN
                  ZMFUU(JL,JK)=PMFU(JL,JK)*(PUU(JL,JK)-PUEN(JL,IK))
                  ZMFUV(JL,JK)=PMFU(JL,JK)*(PVU(JL,JK)-PVEN(JL,IK))
                  ZMFDU(JL,JK)=PMFD(JL,JK)*(PUD(JL,JK)-PUEN(JL,IK))
                  ZMFDV(JL,JK)=PMFD(JL,JK)*(PVD(JL,JK)-PVEN(JL,IK))
               ENDIF
            ENDDO
         ENDDO
         DO JL=KIDIA,KFDIA
            IF(LDCUM(JL)) THEN
               ZMFUU(JL,1)=ZMFUU(JL,2)
               ZMFUV(JL,1)=ZMFUV(JL,2)
               ZMFDU(JL,1)=ZMFDU(JL,2)
               ZMFDV(JL,1)=ZMFDV(JL,2)
            ENDIF
         ENDDO
      ELSE
         DO JK=KTOPM2,KLEV
            IK=JK-1
            DO JL=KIDIA,KFDIA
               IF(LDCUM(JL)) THEN
                  ZMFUU(JL,JK)=PMFU(JL,JK)*(PUU(JL,JK)-PUEN(JL,IK))
                  ZMFUV(JL,JK)=PMFU(JL,JK)*(PVU(JL,JK)-PVEN(JL,IK))
                  ZMFDU(JL,JK)=PMFD(JL,JK)*(PUD(JL,JK)-PUEN(JL,IK))
                  ZMFDV(JL,JK)=PMFD(JL,JK)*(PVD(JL,JK)-PVEN(JL,IK))
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      DO JK=KTOPM2,KLEV
!DIR$ IVDEP
         DO JL=KIDIA,KFDIA
            IF(LDCUM(JL).AND.JK.GT.KCBOT(JL)) THEN
               IKB=KCBOT(JL)
               ZZP=((PAPH(JL,KLEVP1)-PAPH(JL,JK))/
     &              (PAPH(JL,KLEVP1)-PAPH(JL,IKB)))
               IF (KTYPE(JL).EQ.3) ZZP=ZZP**2
               ZMFUU(JL,JK)=ZMFUU(JL,IKB)*ZZP
               ZMFUV(JL,JK)=ZMFUV(JL,IKB)*ZZP
               ZMFDU(JL,JK)=ZMFDU(JL,IKB)*ZZP
               ZMFDV(JL,JK)=ZMFDV(JL,IKB)*ZZP
            ENDIF
         ENDDO
      ENDDO
C
      DO JL=KIDIA,KFDIA
         ZDISS(JL)=0.
      ENDDO
C
      DO JK=KTOPM2,KLEV
C
         IF(JK.LT.KLEV) THEN
            DO JL=KIDIA,KFDIA
               IF(LDCUM(JL)) THEN
                  ZDUDT=(G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*
     &                    (ZMFUU(JL,JK+1)-ZMFUU(JL,JK)+
     &                     ZMFDU(JL,JK+1)-ZMFDU(JL,JK))
                  ZDVDT=(G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*
     &                    (ZMFUV(JL,JK+1)-ZMFUV(JL,JK)+
     &                     ZMFDV(JL,JK+1)-ZMFDV(JL,JK))
                  ZDISS(JL)=ZDISS(JL)+
     &                 PUEN(JL,JK)*(ZMFUU(JL,JK+1)-ZMFUU(JL,JK)+
     &                              ZMFDU(JL,JK+1)-ZMFDU(JL,JK))+
     &                 PVEN(JL,JK)*(ZMFUV(JL,JK+1)-ZMFUV(JL,JK)+
     &                              ZMFDV(JL,JK+1)-ZMFDV(JL,JK))
                  PVOM(JL,JK)=PVOM(JL,JK)+ZDUDT
                  PVOL(JL,JK)=PVOL(JL,JK)+ZDVDT
               ENDIF
            ENDDO
C
         ELSE
            DO JL=KIDIA,KFDIA
               IF(LDCUM(JL)) THEN
                  ZDUDT=-(G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*
     &                 (ZMFUU(JL,JK)+ZMFDU(JL,JK))
                  ZDVDT=-(G/(PAPH(JL,JK+1)-PAPH(JL,JK)))*
     &                 (ZMFUV(JL,JK)+ZMFDV(JL,JK))
                  ZDISS(JL)=ZDISS(JL)-
     &                 (PUEN(JL,JK)*(ZMFUU(JL,JK)+ZMFDU(JL,JK))+
     &                  PVEN(JL,JK)*(ZMFUV(JL,JK)+ZMFDV(JL,JK)))
                  PVOM(JL,JK)=PVOM(JL,JK)+ZDUDT
                  PVOL(JL,JK)=PVOL(JL,JK)+ZDVDT
               ENDIF
            ENDDO
         ENDIF
C
      ENDDO
C
C      ZSUM=SSUM(KLON,ZDISS(1),1)
      ZSUM=SUM(ZDISS)
      PSDISS=PSDISS+ZSUM
C
      RETURN
      END SUBROUTINE CEDUDV
