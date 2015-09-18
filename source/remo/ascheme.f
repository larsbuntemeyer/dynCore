      SUBROUTINE ASCHEME
     &        (KLON , KIDIA , KFDIA , WSM1M, WS     , ARUN   ,
     &         FTHRO, IEXC  , DT    , WSMX , WMINLOK, WMAXLOK,
     &         BETA , LOLAND, LOGLAC)
C
      IMPLICIT NONE
C
C****************************************************************************
C
C     ******** VERSION 2.0 - FEBRUAR 2005
C              PROGRAMMIERUNG UND ENTWICKLUNG: STEFAN HAGEMANN
C
C              ARNO SCHEME LAEUFT NUN ALS SUBROUTINE VON EXCESS ALS VORBEREITUNG
C              ZUR IMPLEMENTIERUNG IN REMO. AUCH TAGES-SUBZEITSCHRITTE SIND
C              NUN ZUM TESTEN MOEGLICH.
C              WATER STORAGE AND CAPACITIES: [M] STATT [M^3]
C              WATER FLUXES IN [M/S] STATT [M^3/S] PER DAY
C              ZEITSCHRITT IN [S] STATT [DAY]
C              ROUTINE WURDE ETWAS OPTIMIERT HINSICHTLICH GESCHWINDIGKEIT
C
C              VS. 2.1 - SEPT. 2007 - STEFAM HAGEMANN
C              ABFRAGE VON WSMX AUF 1E-10. WIRD GEAENDERT WEGEN 
C              RECHENUNGENAUIGKEIT DER NEC AUF WSMX-1.E-10 < ZEPS (1E-15)
C
C     ***      DT = ZEITINTERVALL DER ZEITREIHE [S]
C     ***    KLON = FELDGROESSE (:= NL * NB)
C     ***   KIDIA = INDEX AB DEM IN DEN BODENFELDERN DIE VERARBEITUNG BEGINNT
C     ***   KFDIA = INDEX BEI DEM IN DEN BODENFELDERN DIE VERARBEITUNG ENDET
C     ***
C     ***   WSM1M = SPEICHER UEBER DEM EXCESS RUNOFF STATTFINDET [M]
C     ***           = ACTUAL SOIL MOISTURE CONTENT OF LAST TIMESTEP
C     ***      WS = NEW SOIL MOISTURE CONTENT [M]
C     ***           = WSM1M + INFILTRATION
C     ***    ARUN = AKTUELLER EXCESS RUNOFF
C     ***   FTHRO = THROUGHFALL = RAIN + SNOWMELT
C     ***    CINF = SOIL INFILTRATION CAPACITY FOR IEXC = 1
C     ***  CINMAX = MAXIMUM SOIL INFILTRATION CAPACITY FOR IEXC = 1
C     ***    WSMX = TOTAL SOIL MAXIMUM WATER CAPACITY [M] = WCAP
C     ***
C     *** WMINLOK = MINIMUM SUBGRID WCAP ON MAIN GRID [M]
C     *** WMAXLOK = MAXIMUM SUBGRID WCAP ON MAIN GRID [M]
C     ***    BETA = B-PARAMETER
C     ***
C     ***  LOGLAC = GLACIER JA/NEIN MASK = (1/0)
C     ***  LOLAND = LAND > 0:  JA/NEIN MASK
C     ***    TD3 = SOIL TEMPERATURE OF UPPER LAYER
C     ***  TMELT = METLING TEMPERATURE OF SOIL
C     ***
C     ***    IEXC = ART DER EXCESS RUNOFF BERECHNUNG
C     ***           1 = DUEMENIL/TODINI BZW. ECHAM4
C     ***           5 = IMPROVED ARNO SCHEME NACH HAGEMANN & DUMENIL GATES, 2003
C
CRP   INTEGER IEXC, KLON, KIDIA, KFDIA
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: KLON, KIDIA, KFDIA, IEXC
      REAL,    INTENT(IN)    :: DT
      REAL,    INTENT(IN)    :: WSM1M(KLON), FTHRO(KLON)
      REAL,    INTENT(IN)    :: WSMX(KLON), WMINLOK(KLON), 
     &                          WMAXLOK(KLON), BETA(KLON)
      LOGICAL, INTENT(IN)    :: LOLAND(KLON), LOGLAC(KLON)
      REAL,    INTENT(INOUT) :: ARUN(KLON), WS(KLON)     
C
C     Local Variables
C
      REAL, PARAMETER  :: ZEPS=1.E-15
      REAL             :: CINF(KLON), CINMAX(KLON)
      REAL             :: FINI, FDUM, XDUM
      REAL             :: ZD1, ZWA, ZWAKT, ZWS
      INTEGER          :: JL
CKS      REAL TD3(KLON), TMELT
C
C     ******* ORIGINAL ARNO-SCHEMA
      IF (IEXC.EQ.1) THEN
C
C     ******* BERECHNUNG DER INFILTRATIONSKAPAZITAET
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
               CINMAX(JL) = (1. + BETA(JL)) * WSMX(JL)

               IF (WSM1M(JL).GE.WSMX(JL)) THEN
                  CINF(JL) = CINMAX(JL)
               ELSE IF (WSMX(JL)-1.E-10.GT.ZEPS) THEN
C
C              *** EXPONENT: DUEMENIL 1/(1+BB)
                  CINF(JL) = CINMAX(JL) * (1 - (1-WSM1M(JL)
     &                 /WSMX(JL))**(1/(1+BETA(JL))) )
               ELSE
                  CINF(JL)=0.
               ENDIF
            ENDIF
         ENDDO
C
C     ******* BERECHNUNG DES EXCESS RUNOFFS *********************************
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
C
CHGKS ALLOW INFILTRATION FOR ANY TEMPERATURE
               IF (FTHRO(JL).LE.0. !HGKS.OR. TD3(JL).LT.TMELT
     &              .OR. WSMX(JL)-1.E-10.LE.ZEPS) THEN
C
C           *** SPERREN BEI CRITISCHEM PUNKT ODER BEI INFILT. CAPACITY =0.
                  ARUN(JL) = FTHRO(JL)
               ELSE
                  IF (FTHRO(JL)*DT .GE.
     &                 CINMAX(JL)-CINF(JL)) THEN
                     ARUN(JL) = FTHRO(JL) + (WSM1M(JL)
     &                    - WSMX(JL)) / DT
                  ELSE
                     ARUN(JL) = FTHRO(JL) + (WSM1M(JL)
     &                    - WSMX(JL) + WSMX(JL) * ( 1. -
     &                    (CINF(JL)+FTHRO(JL)*DT)
     &                    /CINMAX(JL) )**(1+BETA(JL)) )
     &                    / DT
                  ENDIF
               ENDIF
C
               IF (ARUN(JL).LT.0.) ARUN(JL) = 0.
            ENDIF
         ENDDO
      ENDIF
C
C     ******* IMPROVED ARNO SCHEME  *********************************
C
      IF (IEXC.EQ.5) THEN
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
C
CHGKS ALLOW INFILTRATION FOR ANY TEMPERATURE
               IF (FTHRO(JL).LE.0. !HGKS.OR. TD3(JL).LT.TMELT
     &              .OR. WSMX(JL)-1.E-10.LE.ZEPS) THEN
C
C           *** SPERREN BEI CRITISCHEM PUNKT ODER BEI INFILT. CAPACITY =0.
                  ARUN(JL) = FTHRO(JL)
               ELSE
C
                  ZWS = WSM1M(JL) + FTHRO(JL)*DT
                  ZD1 = 0.
                  IF (WSM1M(JL).LE.WMINLOK(JL)) THEN
                     ZWAKT = WSM1M(JL)
                  ELSE IF (WSM1M(JL).GE.WSMX(JL)) THEN
                     ZWAKT = WMAXLOK(JL)
                  ELSE
                     FINI = ( 1. - (WSM1M(JL)-WMINLOK(JL)) /
     &                    (WSMX(JL)-WMINLOK(JL)) )
                     IF (FINI.LE.0) THEN
                        ZWAKT = WMAXLOK(JL)
                     ELSE
                        ZD1 = FINI**(1./(1.+BETA(JL)))
                        ZWAKT = WMAXLOK(JL) - 
     &                       ((WMAXLOK(JL)-WMINLOK(JL)) * ZD1)
                     ENDIF
                  ENDIF
                  ZWA = ZWAKT + FTHRO(JL)*DT
C
C
                  IF (WSM1M(JL) .LT. WMINLOK(JL)) THEN
                     IF (ZWS .LE. WMINLOK(JL)) THEN
C
                        ARUN(JL) = 0.
                     ELSE
                        IF (ZWA.GT.WMINLOK(JL) .AND.
     &                      ZWA.LT.WMAXLOK(JL)) THEN
                           FINI = (WMAXLOK(JL) - ZWA) /
     &                          (WMAXLOK(JL)-WMINLOK(JL))
                           XDUM = AMAX1(FINI, 1.E-20)
                           FDUM = (1.+BETA(JL)) * LOG10(XDUM)
                           IF (FDUM.GT.-30.) THEN
                              XDUM = FINI ** (1.+BETA(JL))
                           ELSE
                              XDUM = 0.
                           ENDIF
                           ARUN(JL) = FTHRO(JL)
     &                          - (WMINLOK(JL) - WSM1M(JL)) / DT
     &                          - ( (WMAXLOK(JL)-WMINLOK(JL))
     &                          / (1.+BETA(JL)) * (1. - XDUM ) ) / DT
C
                        ELSE
                           ARUN(JL) = FTHRO(JL) +
     &                          (WSM1M(JL)-WSMX(JL))/DT
CCC                 IF (ZWS.LT.WSMX(JL))
CCC     &              WRITE(*,*) " EXCESS:", JL," WARNING ZWS = ",
CCC     &              ZWS, " < ", WSMX(JL), " = WSMX"
                        ENDIF
                     ENDIF
C
                  ELSE
                     IF (ZWA.LT.WMAXLOK(JL)) THEN
                        FINI = (WMAXLOK(JL) - ZWAKT) /
     &                       (WMAXLOK(JL)-WMINLOK(JL))
                        XDUM = AMAX1(FINI, 1.E-20)
                        FDUM = (1.+BETA(JL)) * LOG10(XDUM)
                        IF (FDUM.GT.-30.) THEN
                           XDUM = FINI ** (1.+BETA(JL))
                        ELSE
                           XDUM = 0.
                        ENDIF
C
                        FINI = (WMAXLOK(JL) - ZWA) /
     &                       (WMAXLOK(JL)-WMINLOK(JL))
                        FDUM = AMAX1(FINI, 1.E-20)
                        FDUM = (1.+BETA(JL)) * LOG10(FDUM)
                        IF (FDUM.GT.-30.) THEN
                           FDUM = FINI ** (1.+BETA(JL))
                        ELSE
                           FDUM = 0.
                        ENDIF
C
                        ARUN(JL) = FTHRO(JL) -
     &                       (WMAXLOK(JL)-WMINLOK(JL)) / (1.+BETA(JL))
     &                       / DT * ( XDUM - FDUM )
C
                     ELSE
                        ARUN(JL) = FTHRO(JL) +
     &                       (WSM1M(JL)-WSMX(JL))/DT
CCC               IF (ZWS.LT.WSMX(JL))
CCC     &           WRITE(*,*) "2. EXCESS:", JL," WARNING ZWS = ",
CCC     &           ZWS, " < ", WSMX(JL), " = WSMX"
                     ENDIF
                  ENDIF
               ENDIF
C
               IF (ARUN(JL).LT.0.) ARUN(JL) = 0.
            ENDIF
         ENDDO
      ENDIF
C
C     *** NEW SOIL MOISTURE CONTENT = OLD + INFILTRATION
      DO JL=KIDIA, KFDIA
         IF (LOLAND(JL)) THEN
            WS(JL)= WSM1M(JL) + ( FTHRO(JL)-ARUN(JL) )*DT
            IF (LOGLAC(JL)) THEN
               WS(JL)= WSM1M(JL)
               ARUN(JL) = 0.
            ENDIF
         ENDIF
      ENDDO
C
C        *** THE END
      RETURN
      END SUBROUTINE ASCHEME
