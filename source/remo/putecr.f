      SUBROUTINE PUTECR
     &  (YTYP  , UFELD  , IEJE   , KE, KO , KU , IVLOC, YVARN ,
     &   LVARDA, NZV    , NZFELD , UR     , VR , TR   , QDR   ,
     &   QWR   , TSWECHR, TSIECHR, SEAICER, PSR, QDBLR, SICEDR,
     &   NRD1  , NRD2)

C
C**** PUTECR   -   UP:BESETZEN DER FELDER IM LANGZEIT-SPEICHER; NUR
C****                 FELDER, DEREN DIMENSIONIERUNG VON IEJE ABHAENGT.
C**   AUFRUF   :   CALL PUTECR
C**  &  (YTYP  , UFELD  , IEJE   , KO     , KU , IVLOC, YVARN ,
C**  &   LVARDA, NZV    , NZFELD , UR     , VR , TR   , QDR   ,
C**  &   QWR   , TSWECHR, TSIECHR, SEAICER, PSR, QDBLR, SICEDR,
C**  &   NRD1  , NRD2)
C**
C**   ENTRIES  :   KEINE
C**   ZWECK    :   BESETZEN DER FELDER IM LANGZEIT-SPEICHER
C**   VERSIONS-
C**   DATUM    :   06.05.04
C**
C**   EXTERNALS:   KEINE
C**
C**   EINGABE-
C**   PARAMETER:   YTYP:   TYP DER DATEN:
C**                        R1: 1.RANDDATENSATZ    R2: 2.RANDDATENSATZ
C**                UFELD:  AUSGEPACKTES DATENFELD  UFELD(IEJE)
C**                IEJE:   DIMENSION VON UFELD:    UFELD(IEJE)
C**                KO:     OBERER  K-INDEX FUER MULTI-LEVEL-FELD
C**                KU:     UNTERER K-INDEX FUER MULTI-LEVEL-FELD
C**                IVLOC:  POSITION DES FELDES IN DER TABELLE DER EM-
C**                        VARIABLEN-NAMEN
C**                YVARN:  NAMEN DER ABZUSPEICHERNDEN A- ODER R-FELDER
C**                NZV:    DIMENSIONIERUNG VON YVARN UND LVARDA
C**   AUSGABE-
C**   PARAMETER:   LVARDA: CHECK-VEKTOR, ZUR PRUEFUNG, OB ALLE FELDER
C**                        GESCHRIEBEN SIND
C**                NZFELD: ANZAHL DER GESCHRIEBENEN FELDER, ZUR PRUE-
C**                        FUNG, OB ALLE FELDER GESCHRIEBEN SIND
C**   COMMON-
C**   BLOECKE  :   CORG, ORG, EMGBCH, EMGBRI, GRIB, HIGKON
C**
C**   METHODE  :
C**
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   R.PODZUN
C     ------------------------------------------------------------------
      IMPLICIT NONE
C
C not used      INCLUDE "corg.h"
      INCLUDE "emgbch.h"
      INCLUDE "emgbri.h"
      INCLUDE "grib.h"
C not used      INCLUDE "higkon.h"
C     ------------------------------------------------------------------
C     Dummy Arguments
C
      INTEGER,   INTENT(IN)    :: IEJE, KE, KO, KU, IVLOC,
     &                            NZV, NRD1, NRD2
      REAL,      INTENT(IN)    :: UFELD(IEJE)
      CHARACTER, INTENT(IN)    :: YVARN (NZV)*8, YTYP*(*)
      INTEGER,   INTENT(INOUT) :: NZFELD     
C
      REAL,      INTENT(INOUT) ::
     & UR     (IEJE,KE,2), VR     (IEJE,KE,2),
     & TR     (IEJE,KE,2), QDR    (IEJE,KE,2),
     & QWR    (IEJE,KE,2), TSWECHR(IEJE,   2),
     & TSIECHR(IEJE,   2), SEAICER(IEJE,   2),
     & PSR    (IEJE,   2), QDBLR  (IEJE,   2),
     & SICEDR (IEJE,   2)
C     ------------------------------------------------------------------
C     Local Declarations
C
      LOGICAL :: LVARDA(NZV)
      INTEGER :: IJ, K, NTAB, NTLEV
      REAL    :: ZBIAS, ZFAK
C     ------------------------------------------------------------------
C     RANDFELDER ABSPEICHERN
          IF(YTYP(2:2).EQ.'1') THEN
              NTLEV = NRD1
          ELSE
              NTLEV = NRD2
          ENDIF

C     ENTSPRECHENDES FELD IM EM-LANGZEITSPEICHER SUCHEN
          DO 40 NTAB = 1,NZV
          IF(YEMNAME(IVLOC).EQ.YVARN(NTAB)) THEN

C     BEI ZEITLICH INTEGRIERTEN GROESSEN (TIME RANGE INDIKATOR = 3)
C     DEN SKALIERUNGSFAKTOR SO AENDERN, DASS MIT VV(H) MULTIPLIZIERT
C     WIRD
              IF(NGBTRI(IVLOC).EQ.3) THEN
                  IF(IPDB(16).EQ.0) ZFAK = FLOAT(IPDB(18))/
     &                             (EMGBFK(IVLOC)*60.0)
                  IF(IPDB(16).EQ.1) ZFAK = FLOAT(IPDB(18))/
     &                              EMGBFK(IVLOC)
              ELSE
                  ZFAK = 1.0/EMGBFK(IVLOC)
              ENDIF

C     FELDER, VON DENEN NUR EIN NIVEAU VORLIEGT (SINGLE LEVEL FIELDS)
              IF(KO.EQ.KU) THEN
                  IF(IPDB(8).NE.1) THEN
                      RETURN
                  ENDIF
C     FELDER, VON DENEN VIELE NIVEAUS VORLIEGEN (MULTI LEVEL FIELDS)
              ELSE
                  IF(IPDB(8).EQ.110) THEN
                      K = IPDB( 9)
                  ELSE IF(IPDB(8).EQ.109) THEN
                      K = IPDB(10)
                  ELSE
                      RETURN
                  ENDIF
              ENDIF

C     FELD IN DER LISTE ABHAKEN, ANZAHL DER GELESENEN FELDER UM 1 ERHOE-
C     HEN.
              LVARDA(NTAB) = .TRUE.
              NZFELD       = NZFELD + 1
C
CRP   BEGINN DER ERSETZUNG VON LOCATE DES MEMORYMANAGERS
C
              ZBIAS=-EMGBBS(IVLOC)*EMGBFK(IVLOC)
C
              IF (YEMNAME(IVLOC).EQ.'U       ') THEN
                  DO 420 IJ=1,IEJE
                  UR(IJ,K,NTLEV)=(UFELD(IJ)+ZBIAS)*ZFAK
 420              CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'V       ') THEN
                  DO 425 IJ=1,IEJE
                  VR(IJ,K,NTLEV)=(UFELD(IJ)+ZBIAS)*ZFAK
 425              CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'T       ') THEN
                  DO 430 IJ=1,IEJE
                  TR(IJ,K,NTLEV)=(UFELD(IJ)+ZBIAS)*ZFAK
 430              CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'QD      ') THEN
                  DO 435 IJ=1,IEJE
                  QDR(IJ,K,NTLEV)=(UFELD(IJ)+ZBIAS)*ZFAK
 435              CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'QW      ') THEN
                  DO 440 IJ=1,IEJE
                  QWR(IJ,K,NTLEV)=(UFELD(IJ)+ZBIAS)*ZFAK
 440              CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'PS      ') THEN
                  DO 445 IJ=1,IEJE
                  PSR(IJ,NTLEV)=(UFELD(IJ)+ZBIAS)*ZFAK
 445              CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TSWECH  ') THEN
                  DO 450 IJ=1,IEJE
                  TSWECHR(IJ,NTLEV)=(UFELD(IJ)+ZBIAS)*ZFAK
 450              CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TSIECH  ') THEN
                  DO 455 IJ=1,IEJE
                  TSIECHR(IJ,NTLEV)=(UFELD(IJ)+ZBIAS)*ZFAK
 455              CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'SEAICE  ') THEN
                  DO 460 IJ=1,IEJE
                  SEAICER(IJ,NTLEV)=(UFELD(IJ)+ZBIAS)*ZFAK
 460              CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'QDBL    ') THEN
                  DO 465 IJ=1,IEJE
                  QDBLR(IJ,NTLEV)=(UFELD(IJ)+ZBIAS)*ZFAK
 465              CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'SICED   ') THEN
                  DO 470 IJ=1,IEJE
                  SICEDR(IJ,NTLEV)=(UFELD(IJ)+ZBIAS)*ZFAK
 470              CONTINUE
              ENDIF
C
              RETURN
          ENDIF
C
 40       CONTINUE
C
      RETURN
      END SUBROUTINE PUTECR
