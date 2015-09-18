      SUBROUTINE ASURF
     &        (KLON , KIDIA, KFDIA, WSM1M , WS    , ZROS, ZPRFL,
     &         ZTMST, WSMX , BETA , LOLAND, LOGLAC)
C
      IMPLICIT NONE
C
C****************************************************************************
C
C     ******* ROUTINE, DIE ARNOSCHEMA AUS REMO-OROGINAL NACHBILDET
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: KLON, KIDIA, KFDIA
      REAL,    INTENT(IN)    :: ZTMST, WSM1M(KLON), ZPRFL(KLON)
      REAL,    INTENT(IN)    :: WSMX(KLON), BETA(KLON)
      LOGICAL, INTENT(IN)    :: LOLAND(KLON), LOGLAC(KLON)
      REAL,    INTENT(INOUT) :: WS(KLON), ZROS(KLON)
C
C     Local Varibles
C
      REAL    :: ZINFIL(KLON), ZVOL(KLON)
      REAL    :: ZB1, ZBM, ZCONW1, ZWMAX, ZBWS, RHOH2O, ZLYEPS, ZLYSIC
      INTEGER :: JL
C
C     *********** ORIGINAL REMO-PART
      RHOH2O = 1000.0
C
      DO JL=KIDIA, KFDIA
         IF (LOLAND(JL)) THEN
            ZINFIL(JL) = 0.
            ZVOL(JL) = 0.
C
C      ACCOUNT FOR SURFACE RUNOFF DUE TO SLOPING TERRAIN
C
            ZWMAX=WSMX(JL)
            ZBWS = BETA(JL)
            ZB1=1.+ZBWS
            ZBM=1./ZB1
            ZCONW1=ZWMAX*ZB1
C
C*    COMPUTE INFILTRATION FROM RAINFALL AND RUNOFF
C     _____________________________________________
C
            ZLYEPS=0.
CHGKS ALLOW INFILTRATION FOR ANY TEMPERATURE
C      IF(TD3(JL).LT.TMELT) THEN
C        ZROS(JL)=ZPRFL(JL)*ZCONS5
C        WS(JL)=WSM1M(JL)
C      ELSE
            IF(ZPRFL(JL).GT.0.) THEN
               IF(WSM1M(JL).GT.ZWMAX) THEN
                  ZLYEPS=WSM1M(JL)-ZWMAX
               ELSE
                  ZLYEPS=0.
               ENDIF
               ZLYSIC=(WSM1M(JL)-ZLYEPS)/ZWMAX
               ZLYSIC=AMIN1(ZLYSIC,1.)
               ZVOL(JL)=(1.-ZLYSIC)**ZBM
     &              -ZTMST*ZPRFL(JL)/(RHOH2O*ZCONW1)
               ZROS(JL)=ZPRFL(JL)*ZTMST/RHOH2O-(ZWMAX-WSM1M(JL))
               IF (ZVOL(JL).GT.0.) THEN
                  ZROS(JL)=ZROS(JL)+ZWMAX*ZVOL(JL)**ZB1
               END IF
               ZROS(JL)=AMAX1(ZROS(JL),0.)
               ZINFIL(JL)=ZPRFL(JL)*ZTMST/RHOH2O-ZROS(JL)
            END IF
            WS(JL)=WSM1M(JL)+ZINFIL(JL)
CHGKS      ENDIF
C
            IF(LOGLAC(JL)) THEN
               ZROS(JL)=0.
               WS(JL)=WSM1M(JL)
            ENDIF
         ENDIF
      ENDDO
C
C     ***********THE END
C
      RETURN
      END SUBROUTINE ASURF
