C
C     SUBROUTINE SOILCHANGE
C
C     INTERMEDIATE ROUTINE FOR CHANGING 5 SOIL MOISTURE LAYERS DUE TO
C     EVAPORATION FLUXES, INFILTRATION AND DRAINAGE THAT ARE CALCULATED
C     WITH THE OLD BUCKET SCHEME. THE SOIL MOISTURE CHANGES ARE CALCULATED
C     FOR ONE TIME STEP, NOT TWO AS PREVIOUSLY DONE FOR WS USING THE
C     3 TIME STEP SCHEME
C
C     *** VS. 1.0 -- JULI 2006 -- STEFAN HAGEMANN
C
C     *** VS 1.1 AUGUST 2006 -- STEFAN HAGEMANN
C     *** SICHERHEITSCHECK IN ISCH =1 SCHLEIFE VERSCHOBEN
C
C     *** VS 1.2 -- JANUARY 2007 -- STEFAN HAGEMANN
C     *** IWORK NACH ISCH UMBENANNT.
C
C     *** VS 1.3 -- SEPTEMBER 2007 -- STEFAN HAGEMANN
C     *** SMALL ERROR IN TRANSPIRATION LOOP ELIMINATED THAT OCCURED IF TRANSPIRATION
C     ***   WAS TAKEN FROM A LAYER WHERE WSI < WILTING POINT. THE ERROR LEAD TO THE
C     ***   THE RESULT THAT THE WSI OF THE CORRESPONDING LAYER WAS SET TO THE WILTING POINT.
C     ***   THE ERROR ONLY CAUSES PROBLEMS IN DRY REGION GRIDBOXES WITH SMALL WSMX THAT ARE
C     ***   CLOSE TO THE BARE SOIL EVAPORATION WATER DEPTH OF 10 CM, AND PROBABLY
C     ***   ONLY IF SMALL TIME STEPS ARE USED, E.G. AT 0.25 DEGREE.
C     ***   RESOLUTION.
C     *** INTRODUCTION OF ZEPS AS SMALL INCREMENT TO BE USED INSTEAD OF ZERO IN SOME
C     ***   COMPARISON-IF STATEMENTS. SMALL FLUX: 0.01 MM/DAY = 1.16 E-10 M/S
C     ***   ZEPS = 1E-15 << AGAINST A SMALL FLUX.
C     *** IF WSI < WWILT FOR ONE LAYER, TRANSPIRATION AMOUNT WILL BE REDUCED
C     ***   INSTEAD OF RE-DISTRIBUTION TO OTHER LAYERS.
C     ***   THUS, AETRANS IS CHANGED IN SOILCHANGE, BUT ALSO EVAP AND EVAPL
C     ***   MUST BE CHANGED.
C     *** BARE SOIL EVAPORATION WILL ONLY BE TAKEN FROM THE MOST UPPER LAYER.
C     ***   OVERSHOOTING AMOUNTS WILL BE DISREGARDED AND AEBSOIL WILL BE REDUCED.
C     ***   THUS, AEBSOIL IS CHANGED IN SOILCHANGE, BUT ALSO EVAP AND EVAPL
C     ***   MUST BE CHANGED.
C     ***  --> NEW ARRAY SUBMITTED TO ROUTINE: REDEVAP
C     ***
C     *** VS 1.4 -- MARCH 2008 -- STEFAN HAGEMANN
C     *** REDUCTION OF EVAPORATION IS TAKEN OUT AGAIN AS THERE WAS NO FEEDBACK
C     *** (TECHNICAL) TO THE ATMOSPHERE. AN IMPLEMENTATION OF THIS
C     *** FEEDBACK WOULD
C     *** REQUIRE MORE COMPLEX CHANGES TO VDIFF AND COMPREHENSIVE TESTING. IN
C     *** MY OPINION A REDUCTION OF TRANSPIRATION IN VDIFF WOULD BE DESIRABLE.
C     *** FELD REDEVAP BLEIBT ERHALTEN, UM ALS DIAGNOSTIC ZU FUNKTIONIEREN.
C
C     ******** INPUTFELDER
C     ***AEBSOIL = BARE SOIL EVAPORATION IN [M], ACCUMULATED OVER TWODT/2
C     ***AETRANS = TRANSPIRATION IN [M], ACCUMULATED OVER TWODT/2
C     ***  DRAIN = DRAINAGE IN [M], ACCUMULATED OVER TWODT
C     ***          --> CONVERSION AUF TWODT/2 NOTWENDIG
C     *** FINFIL = INFILTRATION IN [M/S] --> MULTIPLICATION WITH DT REQUIRED.
C
C     ***    DZR = ROOTING DEPTH IN [M]
C     ***  DZRSI = ROOTED DEPTH PER SOIL LAYER (UNTIL ROOTING DEPTH DZR)
C     ***   DZSI = DEPTH THAT CAN BE SATURATED WITH WATER PER SOIL LAYER
C     ***          (UNTIL BEDROCK DZS)
C     ***  ZWSAT = MAXIMALE SPEICHERGEHALT DER SOIL LAYER I (POROSITAET) [M]
C     ***          WICHTIG FUER AUFNAHME DER INFILTRATION IN DEN BODENSCHICHTEN.
C     ***          IST DEFINIERT BIS ZUR TIEFE DES BODENS DZS,
C     ***          WO DER BEDROCK BEGINNT.
C     ***  ZWSFC = FELD KAPAZITAET DER SOIL LAYER I (ANALOG ZU ZWSAT) [M]
C
C     ***    DT = TIME STEP [USUALLY IN SEC]
C     *** ISCH  = SWITCH WHICH VALUES WILL CHANGE THE SOIL MOISTURE
C     ***         0 = NO CHANGE: SET WSI = WSIM1M
C     ***         1 = CHANGE BY EVAPORATION FLUXES
C     ***         2 = CHANGE BY INFILTRATION
C     ***         3 = CHANGE BY INFILTRATION AND DRAINAGE
C
C
C     ---------------------
      SUBROUTINE SOILCHANGE
     &   (KIDIA , KFDIA  , KLON  , KDEEP,  ISCH ,
     &    DT, ZDIAGS, LOLAND,
C
C         *** - SOIL MOISTURE LAYERS,
     &    WSI   , WSIM1M,
     &    DZR, DZRSI  , DZSI, ZWSAT, ZWSFC,
C
C         ***   VERTICAL WATER FLUXES  -- EVAP. FLUXES OVER LAND
     &    FINFIL, DRAIN, AETRANS, AEBSOIL,
C         ***   REDUCTION OF EVAPORATION IF NECESSARY
     &    REDEVAP)
C
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     Dummy Arguments
C
      INTEGER, INTENT(IN) :: KIDIA , KFDIA  , KLON  , KDEEP,  ISCH
      REAL,    INTENT(IN) :: DT, ZDIAGS

      REAL,    INTENT(IN) :: DRAIN(KLON)
      REAL,    INTENT(IN) :: DZR(KLON), WSIM1M(KLON, KDEEP),
     &                       AETRANS(KLON), AEBSOIL(KLON),
     &                       FINFIL(KLON)
C
      REAL, INTENT(INOUT) :: WSI(KLON, KDEEP)
      REAL, INTENT(IN)    :: DZRSI(KLON,KDEEP), DZSI(KLON,KDEEP)
      REAL, INTENT(IN)    :: ZWSAT(KLON,KDEEP), ZWSFC(KLON,KDEEP)
      REAL, INTENT(INOUT) :: REDEVAP(KLON)
C
C          -- LAND AND ICE MASK      --
      LOGICAL, INTENT(IN) :: LOLAND(KLON)
C-----------------------------------------------------------------------
C     Local Declarations
C
      REAL            :: ZREBSOIL(KLON), ZRETRANS(KLON)
      REAL, PARAMETER :: ZEPS=1.E-15
      REAL            :: ZDUM, ZDUM2, ZREL
      REAL            :: ZOVER(KLON)
      INTEGER         :: JK, JL
      REAL            :: ZWDMIN, ZWILT
C-----------------------------------------------------------------------
C     *** KONSTANTEN
      ZWDMIN=0.05
      ZWILT=0.35                ! SPAETER WAVA EINLESEN
C
C     ******* NO CHANGE = SET WSI TO WSIM1M
      IF (ISCH.EQ.0) THEN
         WSI(1:KLON,1:KDEEP) = WSIM1M(1:KLON,1:KDEEP)
         RETURN
      ENDIF
C
C     ******* CHANGE BY EVAPORATION FLUXES
      IF (ISCH.EQ.1) THEN
         ZREBSOIL(1:KLON) = 0.
         ZRETRANS(1:KLON) = 0.
         ZOVER(1:KLON) = AEBSOIL(1:KLON)
C
C     *** BARE SOIL EVAPORATION LOOP
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
C        *** TAU WIRD IN SURF IN NIEDERSCHLAG GEPACKT
C        ***     --> IST IN INFILTRATION ENTHALTEN!
               IF (ZOVER(JL).LT.-ZEPS) THEN
                  WSI(JL,1) = WSI(JL,1) + ZOVER(JL)
                  IF (WSI(JL,1).LT.0.) THEN
                     ZOVER(JL) = WSI(JL,1)
                     WSI(JL,1) = 0.
                  ELSE
                     ZOVER(JL)=0.
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
C
C       *** REDUCTION OF AEBSOIL NECESSARY? --> DIAGNOSTIC
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
               IF (ZOVER(JL).LT.-ZEPS) THEN
                  ZREBSOIL(JL) = ZOVER(JL)
CCC            AEBSOIL(JL) = AMIN1(AEBSOIL(JL) - ZOVER(JL), 0.)
               ENDIF
            ENDIF
         ENDDO
C
C       *** TRANSPIRATION LOOP
C
         ZOVER(1:KLON) = 0.
         DO JK=1, KDEEP
            DO JL=KIDIA, KFDIA
               IF (LOLAND(JL)) THEN
                  IF (AETRANS(JL).LT.-ZEPS) THEN
C
C         *** TRANSPIRATION BER WURZELTIEFE VERTEILEN
C         *** HIERBEI ERRECHNUNG DURCHWURZELTE TIEFE NOTWENDIG#####
                     IF (DZR(JL).LE.ZEPS) THEN
C
C           *** SONDERFALL VERSIEGELTE GITTERBOX
                        IF (JK.EQ.1) THEN
                           ZRETRANS(JL) = ZRETRANS(JL) + AETRANS(JL)
                        ENDIF
C
                     ELSE IF (DZRSI(JL,JK).GT.0) THEN
C           *** NOTE THAT TRANSPIRATION SHOULD ONLY ACCESS WATER UNTIL THE
C           *** ROOTING DEPTH --> PAY REGARD WHERE THE RELATIVE TERM OF THE
C           *** ROOTED DEPTH WITHIN THE LAYER HAS TO BE APPLIED
                        ZREL = AMIN1(DZRSI(JL,JK)/DZSI(JL,JK), 1.)

                        ZDUM2 = ZWSFC(JL,JK)*ZWILT * ZREL
                        IF (WSI(JL,JK)*ZREL.GE.ZDUM2) THEN
                           ZDUM = WSI(JL,JK) * ZREL +
     &                          AETRANS(JL) * DZRSI(JL,JK)/DZR(JL)
                           IF (ZDUM.LT.ZDUM2) THEN
                              ZOVER(JL) = ZDUM - ZDUM2
                              WSI(JL,JK) = ZDUM2 + WSI(JL,JK) 
     &                             * AMAX1(1.-ZREL, 0.)
                           ELSE
                              WSI(JL,JK) = ZDUM + WSI(JL,JK) 
     &                             * AMAX1(1.-ZREL, 0.)
                              ZOVER(JL)=0.
                           ENDIF
                        ELSE
                           ZOVER(JL) = AETRANS(JL) 
     &                          * DZRSI(JL,JK)/DZR(JL)
                        ENDIF
                        ZRETRANS(JL) = ZRETRANS(JL) + ZOVER(JL)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
C
C       *** ADDING ZRETRANS AND ZREBSOIL TO REDEVAP
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL)) THEN
CCC          IF (ZRETRANS(JL).LT.-ZEPS) THEN
CCC            AETRANS(JL) = AMIN1(AETRANS(JL) - REDEVAP(JL), 0.)
CCC          ENDIF
               REDEVAP(JL) = REDEVAP(JL) + ZREBSOIL(JL) + ZRETRANS(JL)
            ENDIF
         ENDDO
C
C       *** IF REDEVAP NOT ZERO SOIL MOISTURE MUST BE CHANGED
C       *** TO CLOSE THE LAND SURFACE WATER BALANCE.
C       *** REDUCTION WILL BE RELATIVELY DISTRIBUTED TO THE WETNESS
C       *** OF EACH LAYER.
         DO JK=1, KDEEP
            DO JL=KIDIA, KFDIA
               IF (LOLAND(JL)) THEN
                  IF (REDEVAP(JL).LT.-ZEPS) THEN
                     IF (ZWSFC(JL,JK).GT.0) THEN
                        ZDUM = WSI(JL,JK) +
     &                       REDEVAP(JL) * WSI(JL,JK)/ZWSFC(JL,JK)
                        WSI(JL,JK) = AMAX1(ZDUM, 0.)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDIF
C
C     ******* CHANGE BY INFILTRATION
      IF (ISCH.EQ.2 .OR. ISCH.EQ.3) THEN
         ZOVER(1:KLON) = FINFIL(1:KLON)*DT
C
         DO JK = 1, KDEEP
            DO JL=KIDIA, KFDIA
               IF (LOLAND(JL)) THEN
C
C         *** INFILTRATION-SCHLEIFE AUS 5 LAYER ROUTINE
                  WSI(JL,JK) = WSI(JL,JK) + ZOVER(JL)
                  IF (WSI(JL,JK).GT.ZWSAT(JL,JK)) THEN
                     ZOVER(JL) = WSI(JL,JK) - ZWSAT(JL,JK)
                     WSI(JL,JK) = ZWSAT(JL,JK)
                  ELSE
                     ZOVER(JL)=0.
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDIF
C
C
C     ******* CHANGE BY DRAINAGE (UNIT M IN SURF), USED ONLY FOR TESTING
      IF (ISCH.EQ.3) THEN
C
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL).AND.DRAIN(JL).GT.0.) THEN
               ZOVER(JL) = -DRAIN(JL) * ZDIAGS
               DO JK=KDEEP,1,-1
C           *** START DURCHWSSERTE TIEFE
                  ZDUM = ZWSFC(JL,JK)*ZWDMIN
                  IF (ZWSAT(JL,JK).GT.0 .AND.
     &                 WSI(JL,JK)-ZDUM .GT.0) THEN
                     WSI(JL,JK) =  WSI(JL,JK) + ZOVER(JL)
                     IF (WSI(JL,JK)-ZDUM .LT.0) THEN
                        ZOVER(JL) = WSI(JL,JK)-ZDUM
                        WSI(JL,JK) = ZDUM
                     ELSE
                        ZOVER(JL) =0.
                        EXIT    ! AUSSTEIGEN AUS SOIL LAYER LOOP
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF
C
C     ****** AVOID NEGATIVE VALUES - SHOULDN'T HAPPEN
      WSI(1:KLON,1:KDEEP) = AMAX1( WSI(1:KLON,1:KDEEP), 0.)
C
C     *** THE END
      RETURN
      END SUBROUTINE SOILCHANGE
