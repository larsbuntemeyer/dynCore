C
C     SUBROUTINE SOILDEF
C     ***********************************************************************
C
C     *** ROUTINE, IN DER DIVERSE FELDER FUERS 5-LAYER SCHEME DEFINIERT WERDEN.
C
C     ******** VERSION 1.0 - JULI 2006
C              PROGRAMMIERUNG UND ENTWICKLUNG: STEFAN HAGEMANN
C     ******** VERSION 1.1 - JUNI 2007
C              LAUFZEITOPTIMIERT : RALF PODZUN
C     ******** VERSION 1.2 - SEPTEMBER 2007
C              EIN LOOP OPTIMIERT : STEFAN HAGEMANN
C     ******** VERSION 1.3 - OKTOBER 2007
C              MINOR TECHNICAL CORRECTION WITHOUT EFFECT: STEFAN HAGEMANN
C     ******** VERSION 1.4 - JANUARY 2008
C              MINOR TECHNICAL CORRECTION WITHOUT EFFECT: STEFAN HAGEMANN
C     ******** VERSION 1.5 - FEBRUARY 2008 - STEFAN HAGEMANN
C              MINOR TECHNICAL CORRECTION THAT MAY EFFECT ONLY 1-2 GRIDBOXES
C              WHERE DZR=0. (NEAR LAKES)
C
C     ***   ZDEL = DICKE DER 5 BODENSCHICHTEN [M]
C     ***  DZRSI = ROOTED DEPTH PER SOIL LAYER (UNTIL ROOTING DEPTH DZR) [M]
C     ***   DZSI = DEPTH THAT CAN BE SATURATED WITH WATER PER SOIL LAYER
C     ***          (UNTIL BEDROCK DZS) [M]
C     ***    DZR = ROOTING DEPTH IN [M]
C     ***    DZS = SOIL DEPTH DERIVED FROM TEXTURES IN [M]
C     ***   WSMX = FELDKAPAZITAET DES GESAMTEN BODENS IN [M]
C     ***  KDEEP = NUMBER OF SOIL LAYERS
C     ***  ZWSAT = MAXIMALE SPEICHERGEHALT DER SOIL LAYER I (POROSITAET) [M]
C     ***          WICHTIG FUER AUFNAHME DER INFILTRATION IN DEN BODENSCHICHTEN.
C     ***          IST DEFINIERT BIS ZUR TIEFE DES BODENS DZS, WO DER BEDROCK BEGINNT.
C     ***  ZWSFC = FELD KAPAZITAET DER SOIL LAYER I (ANALOG ZU ZWSAT) [M]
C
C-----------------------------------------------------------------------
      SUBROUTINE SOILDEF
     &   (KIDIA, KFDIA, KLON, KDEEP,
     &    WSMX, DZR, DZS, VPOR, LOLAND, ZDEL,
     &    DZRSI, DZSI, ZWSAT, ZWSFC)
C
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: KIDIA, KFDIA, KLON, KDEEP
      REAL,    INTENT(INOUT) :: DZRSI(KLON, KDEEP), DZSI(KLON, KDEEP)
      REAL,    INTENT(INOUT) :: ZWSAT(KLON,KDEEP), ZWSFC(KLON,KDEEP)
      REAL,    INTENT(IN)    :: WSMX(KLON), DZS(KLON), 
     &                          DZR(KLON),  VPOR(KLON)
      REAL,    INTENT(IN)    :: ZDEL(KDEEP)
      LOGICAL, INTENT(IN)    :: LOLAND(KLON)
C-----------------------------------------------------------------------
C     Local Declarations
C
      REAL, PARAMETER :: ZEPS=1.E-15
      INTEGER :: I,JL
      LOGICAL :: LOFLAG1(KLON), LOFLAG2(KLON)
      REAL    :: ZVFC(KLON),  ZDELS(KDEEP), ZDELSM1(KDEEP)
C
C-----------------------------------------------------------------------
C     ************ FELDINI -
      DO JL=1,KLON
         LOFLAG1(JL)=.TRUE.
         LOFLAG2(JL)=.TRUE.
      ENDDO

      CALL SETRA(ZVFC,KLON,0.)
      CALL SETRA(DZSI,KLON*KDEEP,0.)
      CALL SETRA(DZRSI,KLON*KDEEP,0.)
C
!CDIR NOVECTOR
      DO I=2,KDEEP
         ZDELS(I)=SUM(ZDEL(1:I))
         ZDELSM1(I)=SUM(ZDEL(1:I-1))
      ENDDO
C
      ZDELS(1)=ZDEL(1)
      ZDELSM1(1)=0.
C
      DO JL=KIDIA, KFDIA
         IF (LOLAND(JL)) THEN
C
C       #### DEF.: ZVFC: WAS IST BEI WSMX > 0, DZS > 0, DZR = 0?
            IF (DZR(JL).GT.ZEPS) THEN
               ZVFC(JL) = WSMX(JL)/DZR(JL)
            ELSE IF(DZS(JL).GT.ZEPS) THEN
               ZVFC(JL) = WSMX(JL)/DZS(JL)
            ENDIF
         ENDIF
      ENDDO
C
      DO I=1,KDEEP
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL).AND.LOFLAG1(JL)) THEN
               IF (DZR(JL).GT.ZDELS(I)) THEN
                  DZRSI(JL,I) = ZDEL(I)
               ELSE
                  DZRSI(JL,I) = DZR(JL) - ZDELSM1(I)
                  LOFLAG1(JL) = .FALSE.
               ENDIF
            ENDIF
         ENDDO
      ENDDO
C
C       *** SATURATABLE DEPTH
      DO I=1, KDEEP
         DO JL=KIDIA, KFDIA
            IF (LOLAND(JL).AND.LOFLAG2(JL)) THEN
               IF (DZS(JL).GT.ZDELS(I)) THEN
                  DZSI(JL,I) = ZDEL(I)
               ELSE
                  DZSI(JL,I) = DZS(JL) - ZDELSM1(I)
                  LOFLAG2(JL) = .FALSE.
               ENDIF
            ENDIF
C
         ENDDO
      ENDDO
C
      DO I=1, KDEEP
         DO JL=1,KLON
            ZWSAT(JL,I) = VPOR(JL) * DZSI(JL,I)
            ZWSFC(JL,I) = ZVFC(JL) * DZSI(JL,I)
         ENDDO
      ENDDO
C
      RETURN
      END SUBROUTINE SOILDEF
