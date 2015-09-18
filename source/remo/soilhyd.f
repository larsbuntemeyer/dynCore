C
C     SUBROUTINE SOILHYD
C
C     **************************************************************************
C
C     ******* SOIL HYDROLOGY - ROUTINE FUER SPAETEREN EINBAU IN REMO
C     ******* PROGRAMM, WELCHE DIE SOIL HYDROLOGY MIT 5 SCHICHTEN BERECHNET.
C     ******* SOIL HYDROLOGY = PERCOLATION UND VERTICAL DIFFUSION OF WATER
C     ******* DIFFUSIONSSART:  RICHTMEYR/MORTON
C
C     *** PROGRAMMIERUNG UND ENTWICKLUNG: STEFAN HAGEMANN
C     ****
C     *** VS. 1.0 - VERSION ZUM TESTEN DER PROZESSE IM BODEN OHNE RUNOFF UND SO.
C              PERKOLATIONSFORMULIERUNG IN ABHAENGIGKEIT VON ART IPERC
C              IN REMO WIRD IN [M] GERECHNET, NUR AUSGABE IST
C              IN [MM]/AUSGABEZEITSCHRITT.
C     ***
C     *** VS. 1.1 - VERSION FUER REMO
C     ***           UNITS MUSSTEN ANGEPASST WERDEN.
C
C     *** VS. 1.2 - BESCHRAENKUNG AUF IPERC=2, IDIFF=2
C     ***           --> ALTERNATIVEN RAUSGENOMMEN
C
C     *** VS. 1.3 - LAST LOOP SEPARIERT UND GEDREHT, DAMIT LAENGSTE SCHLEIFE
C                   INNEN
C
C     *** VS. 1.4 - SEPTEMBER 2007 - HAG
C     *** LOOP AUS NUMERICAL RECIPES UMPROGRAMMIERT, SO DASS
C     ***     LAENGSTE SCHLEIFE INNEN LIEGT. DAZU WURDEN ZTRI(KDEEP) UND ZDC
C     ***     IN FELDER ZTRI(KLON,KDEEP) UND ZDC(KLON,KDEEP) UMGEWANDELT.
C     *** INCONSISTENCY FOR PERCOLATION OVERFLOW OF LAYES ELIMINATED
C     ***   NOW OVERFLOW EXISTS IF WSI > ZWSFC (INSTEAD OF ZWSAT AS USED BEFORE)
C     *** ABFRAGEN AUF 1E-10. WERDEN GEAENDERT WEGEN
C     ***   RECHENUNGENAUIGKEIT DER NEC AUF WSMX-1.E-10 < ZEPS (1E-15)
C     ***
C     ***
C     *** ES WERDEN IMMER FLUESSE BERECHNET, DIE FUER 2 ZEITSCHRITTE IM
C     *** VORAUS GELTEN.
C     *** DAHER MUESSEN FLUESSE BEI AUFAKKUMULIERUNG DURCH 2 GETEILT WERDEN.
C     *** IN DER BERECHNUNG SELBST WERDEN VOLUMINA BERECHNET, NICHT FLUESSE.
C     *** DER FLUSS ENTSTEHT DADURCH, DASS SPAETER DIE EINHEIT DES FLUSSES
C     *** DAS AKKUMULIERTE VOLUMEN PRO ZEITINTERVALL IST.
C     ***
C     *** VS. 1.5 - DECEMBER 2007 - HAG
C     ***           USE OF ZWSFC INSTEAD OF ZWSAT FOR PERCOLATION
C
C     SOILHYD
C     ***  FKSAT = SATURATED HYDRAULIC CONDUCTIVITY: UNIT [M/S] -> 
C     ***          MUSS IN ECMGBDT.F SICHERGESTELLT WERDEN.
C     ***  FMPOT = MATRIX POTENTIAL [M]
C     *** BCLAPP = EXPONENT B IN CLAPP AND HORNBERGER
C     ***   VPOR = VOLUMETRIC SOIL POROSITY [M/M]   FOLLOWING COSBY ET AL.
C
C     *** ZDIAGT = ZEITSCHRITT: UNIT [S]   ##### ORG. WAR [DAY]
C
C     *** ZDRAIN = DRAINAGE: VOLUMEN IM ZEITSCHRITT: UNIT [M]
C     ***  ZWSAT = MAXIMALE SPEICHERGEHALT DER SOIL LAYER I (POROSITAET) [M]
C     ***          WICHTIG FUER AUFNAHME DER INFILTRATION IN DEN BODENSCHICHTEN.
C     ***          IST DEFINIERT BIS ZUR TIEFE DES BODENS DZS,
C     ***          WO DER BEDROCK BEGINNT.
C     ***  ZWSFC = FELD KAPAZITAET DER SOIL LAYER I (ANALOG ZU ZWSAT) [M]
C     ***
C     ***  ZPERC = DOWNWARD FLUX AUS SOIL LAYER I: ERST IN M/DAY, DANN IN [M].
C     ***          ZUERST: DUE TO GRAVITATIONAL DRAINAGE
C     ***  ZDIFF = DIFFUSION OF WATER BETWEEN SOIL LAYER I AND
C     ***          LAYER I+1 [M^2/DAY]
C     ***
C     ***  NSOIL - NUMBER OF SOIL LAYERS
C     *** WSI(I) = SOIL MOISTURE CONTENT IN I. SOIL LAYER [M]
C     ***          WIRD UPGEDATED DURCH DIE IN SOILHYD BERECHNETEN FLUESSE.
C     ***  ILOG =  1 --> LOGOUTPUT NACH STDOUT MIT KEYWORD SOILFLUX
C     ***                FOR GRIDBOX JLLOG
C     ***                5 * DIFF. FLUXES, 5 PERC. FLUXES, ZDRAIN
C     ***          2 --> LOGOUTPUT NACH STDOUT MIT KEYWORD SOILHYD1 FOR
C     ***                GRIDBOX JLLOG
C
C     ***  IPERC = PERCOLATIONSART (NACH TESTEN IN REMO RAUSGENOMMEN)
C     ***          0  KEINE
C     ***          1  ECHAM4 DRAINAGE
C     ***          2  CLAPP & HORNBERGER  --> GETESTET & ERWAEHLT
C     ***  ANMERKUNG: IN PROGRAMM TESTSHYD.F GIBT ES AUCH NOCH
C     ***   3  BROOKS & COREY,  4  VAN GENUCHTEN,
C     ***   5  VAN GENUCHTEN METHODE 2 NACH DISSE
C     ***   DIESE BENOETIGEN ABER DEN PORE SIZE DISTRIBUTION INDEX PSIBC,
C     ***   DER IN REMO NICHT VORHANDEN IST.
C     ***
C     ***   DIE FORMULIERUNGEN FUER 2,3,4,5 SIND NUR DANN VALID, WENN
C     ***   DER ZEITSCHRITT 1 STUNDE ODER KUERZER IST. 1 IST AUCH MIT
C     ***   TAGESZEITSCHRITTEN LAUFFAEHIG. JEDOCH BESCHREIBT 1 EINEN LATERALEN
C     ***   PROZESS, UND KEINEN VERTIKALEN WIE 2-5.
C     ***
C
      SUBROUTINE SOILHYD(KIDIA, KFDIA, KLON, KDEEP,
     &   ZDIAGT, ILOG, JLLOG,
     &   WSI, DZSI, ZWSFC,
     &   FKSAT, VPOR, BCLAPP, FMPOT,
     &   ZDRAIN
     &   )
C
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: KIDIA, KFDIA, KLON, KDEEP,
     &                          ILOG, JLLOG
      REAL,    INTENT(IN)    :: ZDIAGT
      REAL,    INTENT(INOUT) :: WSI(KLON, KDEEP)
      REAL,    INTENT(IN)    :: DZSI(KLON, KDEEP)
      REAL,    INTENT(IN)    :: FKSAT(KLON)
      REAL,    INTENT(IN)    :: VPOR(KLON), BCLAPP(KLON), FMPOT(KLON)
      REAL,    INTENT(INOUT) :: ZDRAIN(KLON)
      REAL,    INTENT(IN)    :: ZWSFC(KLON,KDEEP)
C-----------------------------------------------------------------------
C     Local Declarations
C
      REAL, PARAMETER :: ZEPS=1.E-15
      REAL :: ZDUM, ZDUM2, ZDUM3
      REAL :: ZPERC(KLON,KDEEP), ZDUM4
      REAL :: ZWDMIN
      REAL :: ZDIFLOG(10), WSLOG(KDEEP)
      INTEGER :: I,JL
C
C     *** DIFFUSION
      REAL :: ZDA(KLON,KDEEP), ZDB(KLON,KDEEP), ZDIFF(KLON,KDEEP-1)
      REAL :: ZTRI(KLON,KDEEP), ZDC(KLON,KDEEP)
C-----------------------------------------------------------------------
      ZWDMIN=0.05
C
C     ***** NULL-INIS UND LOKALE FELDER
      ZDRAIN(1:KLON)=0.
      ZPERC(1:KLON,1:KDEEP)=0.
C
C     *** PERCOLATION = GRAVITATIONAL DRAINAGE - VON UNTEN NACH OBEN BERECHNET.
C     *** MAXIMUM PERCOPLATION WEGEN NUMERISCHER INSTABILITY BEI GROSSEN
C     *** ZEITSCHRITTEN UND KLEINEN RAEUMLICHEN ENTFERNUNGEN.
C
      DO I = 1, KDEEP
         DO JL=KIDIA, KFDIA
CCC        IF (DZSI(JL,I).GT.0) THEN
            IF (ZWSFC(JL,I).GT.0) THEN
               ZDUM = WSI(JL,I) / ZWSFC(JL,I)
               IF (ZDUM-1.E-10.GT.ZEPS) THEN
                  ZPERC(JL,I) = FKSAT(JL) *
     &                 ZDUM ** (2*BCLAPP(JL)+3)
               ELSE
                  ZPERC(JL,I) = 0.
               ENDIF
            ENDIF
         ENDDO
      ENDDO
C
C     *** DIFFUSION OF WATER *******************************************
C
C     *** NUMERICAL RECIPES/RICHTMYER & -MORTON DIFFUSION  *************
C
      DO I = 1, KDEEP-1
         DO JL=KIDIA, KFDIA
            IF (DZSI(JL,I).GT.0) THEN
C
C              ***  SOIL MOISTURE DIFFUSIVITY [M^2/DAY]
C              ***  METHOD C: DIFFUSIVITY OF BOTH LEVELS IS AVERAGED
C              ***            TO CALCULATE DIFFUSIVITY BETWEEN LAYERS
               IF (WSI(JL,I+1).GT.0) THEN
                  ZDUM = WSI(JL,I+1)/DZSI(JL,I+1)
C
C                 *** CALCULATING THE DIFFUSIVITY OF LAYER I+1
                  ZDUM4 = ZDUM / VPOR(JL)
                  IF (ZDUM4-1.E-10.GT.ZEPS) THEN
                     ZDUM2 = BCLAPP(JL) * FKSAT(JL) * FMPOT(JL) /
     &                    ZDUM * ZDUM4 ** (BCLAPP(JL)+3.)
                  ELSE
                     ZDUM2 = 0.
                  ENDIF
               ELSE
                  ZDUM2=0.
               ENDIF
               IF (WSI(JL,I).GT.0) THEN
                  ZDUM = WSI(JL,I)/DZSI(JL,I)
C
C                 *** CALCULATING THE DIFFUSIVITY OF LAYER I
                  ZDUM4 = ZDUM / VPOR(JL)
                  IF (ZDUM4-1.E-10.GT.ZEPS) THEN
                     ZDUM3 = BCLAPP(JL) * FKSAT(JL) * FMPOT(JL) /
     &                    ZDUM *  ZDUM4 ** (BCLAPP(JL)+3.)
                  ELSE
                     ZDUM3 = 0.
                  ENDIF
               ELSE
                  ZDUM3=0.
               ENDIF
C
C              *** CALCULATING THE DIFFUSIVITY AT THE BOTTOM OF LAYER I
               ZDIFF(JL,I) = (ZDUM2+ZDUM3)*0.5
            ENDIF
         ENDDO
      ENDDO
C
C     *** BERECHNEN DER DIFFUSIONSKOEFFIZIENTEN
C     *** DEEPEST LAYER KDEEP = 5
      DO JL=KIDIA, KFDIA
         ZDA(JL, KDEEP) = 0.
         IF (DZSI(JL,KDEEP).GT.0) THEN
C
C           *** [ZDIAGT] = S
            ZDB(JL, KDEEP) = ZDIFF(JL,KDEEP-1) *ZDIAGT /DZSI(JL,KDEEP)
     &           / (DZSI(JL,KDEEP)+DZSI(JL,KDEEP-1)) * 2.
         ELSE
            ZDB(JL, KDEEP) = 0.
         ENDIF
      ENDDO
C
C     *** LAYERS I = 4,3,2
      DO I=KDEEP-1, 2,-1
         DO JL=KIDIA, KFDIA
            IF (DZSI(JL,I).GT.0) THEN
               IF (DZSI(JL,I+1).GT.0) THEN
                  ZDA(JL,I) = ZDIFF(JL,I) * ZDIAGT / DZSI(JL,I) /
     &                 (DZSI(JL,I)+DZSI(JL,I+1)) * 2.
               ELSE
                  ZDA(JL,I) = 0.
               ENDIF
               ZDB(JL,I) = ZDIFF(JL,I-1) * ZDIAGT / DZSI(JL,I) /
     &              (DZSI(JL,I)+DZSI(JL,I-1)) * 2.
            ELSE
               ZDA(JL, I) = 0.
               ZDB(JL, I) = 0.
            ENDIF
         ENDDO
      ENDDO
C
C     *** LAYER 1
      DO JL=KIDIA, KFDIA
         ZDB(JL, 1) = 0.
         IF (DZSI(JL,1).GT.0) THEN
            IF (DZSI(JL,2).GT.0) THEN
               ZDA(JL,1) = ZDIFF(JL,1) * ZDIAGT / DZSI(JL,1) /
     &              (DZSI(JL,1)+DZSI(JL,2)) * 2.
            ELSE
               ZDA(JL,1) = 0.
            ENDIF
         ELSE
            ZDA(JL, 1) = 0.
         ENDIF
      ENDDO
C
C     *** ZWISCHENSPEICHER FUER LOGAUSGABE
      IF (ILOG.EQ.1 .OR. ILOG.EQ.3) THEN
         DO I=1, KDEEP
            WSLOG(I) = WSI(JLLOG, I)
         ENDDO
      ENDIF
C
C     *** ROUTINE TRIDAG AUS NUMERICAL RECIPES, P. 43
C     *** -ZDA = CI, -ZDB=AI, ZDA+ZDB+1=BI, WSI(T)=RI, WSI(T+1)=UI
      DO JL=KIDIA, KFDIA
         IF (DZSI(JL,1).GT.0) THEN
            ZDC(JL,1) = ZDA(JL,1)+ZDB(JL,1)+1.
            WSI(JL,1) = WSI(JL,1) / ZDC(JL,1)
         ENDIF
      ENDDO
C     *** DECOMPOSITION AND FORWARD SUBSTITUTION
      DO I=2, KDEEP
         DO JL=KIDIA, KFDIA
            IF (DZSI(JL,1).GT.0) THEN
C
               IF (DZSI(JL,I).GT.0) THEN
                  ZTRI(JL,I) = -ZDA(JL,I-1) / ZDC(JL,I-1)
                  ZDC(JL,I) = ZDA(JL,I)+ZDB(JL,I)+1. +
     &                 ZDB(JL,I)*ZTRI(JL,I)
                  WSI(JL,I) = (WSI(JL,I)/DZSI(JL,I) + ZDB(JL,I)
     &                 * WSI(JL,I-1)/DZSI(JL,I-1) )
     &                 / ZDC(JL,I) * DZSI(JL,I)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
C
      DO I=KDEEP-1,1,-1         ! BACKSUBSTITUTION
         DO JL=KIDIA, KFDIA
            IF (DZSI(JL,1).GT.0) THEN
C
               IF (DZSI(JL,I+1).GT.0) THEN
                  WSI(JL,I) = WSI(JL,I) - ZTRI(JL,I+1)*WSI(JL,I+1)*
     &                 (DZSI(JL,I)/DZSI(JL,I+1))
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!
!     *** LOGAUSGABE DER FLUXES ?
      IF (ILOG.EQ.1 .OR. ILOG.EQ.3) THEN
         DO I=1, KDEEP
            ZDIFLOG(I) = WSLOG(I) - WSI(JLLOG, I)
         ENDDO
      ENDIF
C
C     *** URSPRUENGLICHE PERKOLATIONS-FLUSSEINHEIT (=UNIT OF FKSAT)
C     *** WAR M/DAY
C     *** NUN IST ES ACCUMULATED M OVER ZDIAGT (UNIT VON ZDIAGT=[S])
C     *** --> ZEITSCHRITT BERUECKSICHTIGEN: * ZDIAGT
      ZPERC(1:KLON,1:KDEEP) = ZPERC(1:KLON,1:KDEEP) * ZDIAGT
C
      IF (ILOG.EQ.2 .OR. ILOG.EQ.3) THEN
         WRITE(*,'(A10, 11(E16.9E2) )')
     &        "SOILHYD1: ", FKSAT(JLLOG)*1000., WSI(JLLOG,1:KDEEP)*1000.
     &        ,ZWSFC(JLLOG,1:KDEEP)*1000.
      ENDIF
C
C     ******* NEW SOIL LAYER SOIL MOISTURE CONTENTS *********************
C
C     *** PERKOLATION PRO ZEITSCHRITT MUSS KLEINER ALS USEABLE SOIL LAYER
C     *** CONTENT SEIN. MAXIMUM PERCOPLATION WEGEN NUMERISCHER INSTABILITY
C     *** BEI GROSSEN ZEITSCHRITTEN UND KLEINEN RAEUMLICHEN ENTFERNUNGEN:
C
      DO I = 1, KDEEP
         DO JL=KIDIA, KFDIA
            IF (DZSI(JL,I).GT.0 .AND. ZPERC(JL,I).GT.0) THEN
C
               ZDUM = ZWDMIN*ZWSFC(JL,I)
               IF (WSI(JL,I)-ZPERC(JL,I) .LT. ZDUM) THEN
                  IF (WSI(JL,I).GT.ZDUM) THEN
                     ZPERC(JL,I) = WSI(JL,I) - ZDUM
                     WSI(JL,I) = ZDUM
                  ELSE
                     ZPERC(JL,I) = 0.
                  ENDIF
               ELSE
                  WSI(JL,I) = WSI(JL,I) - ZPERC(JL,I)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
C
C     *** DRAINAGE CALCULATION FROM PERCOLATION OVERFLOW
C
      DO JL=KIDIA, KFDIA
         IF (DZSI(JL,1).GT.0) THEN
            IF (WSI(JL,1).GT.ZWSFC(JL,1)) THEN
               ZDRAIN(JL) = ZDRAIN(JL) + WSI(JL,1) - ZWSFC(JL,1)
               WSI(JL,1) = ZWSFC(JL,1)
            ENDIF
         ENDIF
      ENDDO
C
      DO I = 1, KDEEP-1
         DO JL=KIDIA, KFDIA
            IF (DZSI(JL,1).GT.0) THEN
               IF (DZSI(JL,I+1).GT.0) THEN
                  WSI(JL,I+1) = WSI(JL,I+1) + ZPERC(JL,I)
                  IF (WSI(JL,I+1).GT.ZWSFC(JL,I+1)) THEN
                     ZDRAIN(JL) = ZDRAIN(JL) + WSI(JL,I+1) -
     &                    ZWSFC(JL,I+1)
                     ZPERC(JL,I) =
     &                MAX( ZPERC(JL,I) - WSI(JL,I+1) + ZWSFC(JL,I+1)
     &                    , 0.)
                     WSI(JL,I+1) = ZWSFC(JL,I+1)
                  ENDIF
               ELSE
                  ZDRAIN(JL) = ZDRAIN(JL) + ZPERC(JL,I)
                  ZPERC(JL,I) = 0.
               ENDIF
            ENDIF
         ENDDO
      ENDDO
C
      DO JL=KIDIA, KFDIA
         IF (DZSI(JL,1).GT.0) THEN
            ZDRAIN(JL) = ZDRAIN(JL) + ZPERC(JL,KDEEP)
         ENDIF
      ENDDO
C
C     *** LOGAUSGABE DER FLUXES ?
      IF (ILOG.EQ.1 .OR. ILOG.EQ.3) THEN
         DO I=1, KDEEP
            ZDIFLOG(KDEEP+I) = ZPERC(JLLOG,I)
         ENDDO
         WRITE(*, '(A10, 11(E15.6E3) )') "SOILFLUX",
     &        ZDIFLOG(1:10)*1000., ZDRAIN(JLLOG)*1000.
      ENDIF
C
C     *** LOGAUSGABE DER WSI!
      IF (ILOG.EQ.2 .OR. ILOG.EQ.3) THEN
         WRITE(*,'(A10, 10(E16.9E2) )' )
     &        "SOILHYD2: ", WSI(JLLOG,1:KDEEP)*1000.,
     &        DZSI(JLLOG,1:KDEEP)*1000.
      ENDIF
C
      RETURN
      END SUBROUTINE SOILHYD
C
