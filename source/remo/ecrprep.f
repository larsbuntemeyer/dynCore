C
C     SUBROUTINE ECRPREP
C
C**** ECRPREP   -   UP:BEREITSTELLUNG DER EM-RANDDATEN IM EM-LANGZEIT-
C****                 SPEICHER
C**   AUFRUF   :   CALL ECRPREP(NZEIT, NAKTION) IN UP *EC4ORG*
C**   ENTRIES  :   KEINE
C**   ZWECK    :   BEREITSTELLUNG DER EM-RANDDATEN FUER DIE EM-PROGNOSE.
C**                HOLEN UND LESEN EINER RANDDATEI BZW. KOPIE DER AN-
C**                FANGSDATEN IN DIE 1.RANDDATEI FUER NZEIT=0.
C**   VERSIONS-
C**   DATUM    :   21.06.01
C**                2007
C**
C**   EXTERNALS:   MAKEPN , GETD   , GNZFLD , DATEIN , GRIDCHK, DATCHK ,
C**                GETLOC , PUTECR ,
C**                DCOMPLT, PRCV   , PRIV   , ECAMIXR
C**
C**   EINGABE-
C**   PARAMETER:   NZEIT:  ZU HOLENDER RANDZEITPUNKT (ZEITSCHRITTE BZGL.
C**                        YADAT)
C**                NAKTION:AUTRAGSNUMMER; 1: 1.ZEITPUNKT;
C**                                       2: FOLGEZEITPUNKTE
C**   AUSGABE-
C**   PARAMETER:   KEINE
C**
C**   COMMON-
C**   BLOECKE  :   ORG, CORG, COMDYN, COMPRI, EMGBCH, EMGBRI,
C**                UNITCH, UNITNR, COMPHY
C**
C**   METHODE  :   SEQUENTIELLER AUFRUF DER UNTERPROGRAMME ZUM HOLEN
C**                LESEN, PRUEFEN UND VERARBEITEN DER DATEN
C**   FEHLERBE-
C**   HANDLUNG :   STOP/ABORT IM FEHLERFALL
C**   VERFASSER:   R.PODZUN
C
      SUBROUTINE ECRPREP(NZEIT , NAKTION,
     &     UR   , VR     , TR     , QDR    , QWR   , PSR   ,
     &     QDBLR, TSWECHR, TSIECHR, SEAICER, U     , V     , T     ,
     &     QD   , QW     , PS     , QDBL   , TSWECH, TSIECH, SEAICE,
     &     SICED, SICEDR)
C
      IMPLICIT NONE
C
      INTEGER, PARAMETER :: KEMAX=50
C
      INCLUDE "parorg.h"
      INCLUDE "org.h"
      INCLUDE "corg.h"
      INCLUDE "emgbch.h"
      INCLUDE "emgbri.h"
      INCLUDE "unitnr.h"

      INCLUDE "comdyn.h"
      INCLUDE "comphy.h"
      INCLUDE "unitch.h"
      INCLUDE "grib.h"

      LOGICAL ::   LVARDA(NZMXVR)
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN) :: NZEIT, NAKTION
C
C     PUTECR
C
      REAL, INTENT(INOUT) :: UR     (IEJEKE,2), VR     (IEJEKE,2),
     &                       TR     (IEJEKE,2), QDR    (IEJEKE,2),
     &                       QWR    (IEJEKE,2),
     &                       PSR    (IEJE,2)  , QDBLR  (IEJE,2) ,
     &                       TSWECHR(IEJE,2)  , TSIECHR(IEJE,2) ,
     &                       SEAICER(IEJE,2)  , SICEDR (IEJE,2)
C
C     RDAUSAD, RDQDW, ADTORD
C
      REAL, INTENT(IN)   ::  U     (IEJEKE,3), V     (IEJEKE,3),
     &                       T     (IEJEKE,3), QD    (IEJEKE,3),
     &                       QW    (IEJEKE,3),
     &                       PS    (IEJE,3)  , QDBL  (IEJE,3)  ,
     &                       TSWECH(IEJE,3)  , TSIECH(IEJE,3)  ,
     &                       SEAICE(IEJE  )  , SICED (IEJE  )
C
C     LOKALE FELDER
C
      REAL      ::   UFELD(MOIEJE)
      REAL      ::   VFELD(IEJE+38)
      REAL      ::   ZAK(KEMAX), ZBK(KEMAX)
      CHARACTER ::   YTYP*3
C
C     Local Variables
C
      INTEGER   ::   ISTAT,NZRFELD,NZGFELD,NTLEV,NTAB,IVLOC,IJ,I
C
      ISTAT=0
      TAGCOUNT=0
      TAGTABLE(1)=1
      TAGTABLE(2)=2

C     AKTUELLE ZAHL DER RANDVARIABLEN; WENN QW IN DEN RANDFELDERN,
C     NZVR ERHOEHEN.
      IF (LQWR) THEN
         NZVR = 10
         YRVARN (NZVR) = 'QW     '
      ELSE
         NZVR = 9
      ENDIF
C
      IF (LSICED) THEN
         NZVR = NZVR + 1
         YRVARN(NZVR) = 'SICED   '
      ENDIF
C
      IF (NAKTION.EQ.1) THEN
         YTYP = 'R1'
      ELSE IF (NAKTION.EQ.2) THEN
         YTYP = 'R2'
      ENDIF
C
C     WENN NZT=0 (BEGINN DER VORHERSAGE) 1.RANDDATENSATZ ALS
C     KOPIE DER ANFANGSDATEN, WENN ANALYSEN ALS RANDWERTE BENUTZT WERDEN
      IF (NZT.EQ.0 .AND.NAKTION.EQ.1) THEN
C
         IF (YTYP(2:2).EQ.'1') THEN
            NTLEV = NRD1
         ELSE
            NTLEV = NRD2
         ENDIF
C
CRP      ERSETZEN VON CALL RDAUSAD SOWIE CALL ADTORD IND RDAUSAD
C
         DO IJ=1,IEJEKE
            UR (IJ,NTLEV)=U (IJ,1)
            VR (IJ,NTLEV)=V (IJ,1)
            TR (IJ,NTLEV)=T (IJ,1)
            QDR(IJ,NTLEV)=QD(IJ,1)
         ENDDO
C
         IF (LQWR) THEN
            DO IJ=1,IEJEKE
               QWR(IJ,NTLEV)=QW(IJ,1)
            ENDDO
         ELSE
            DO IJ=1,IEJEKE
               QWR(IJ,NTLEV)=0.
            ENDDO
         ENDIF
C
         DO IJ=1,IEJE
            PSR    (IJ,NTLEV)=PS    (IJ,1)
            TSWECHR(IJ,NTLEV)=TSWECH(IJ,1)
            TSIECHR(IJ,NTLEV)=TSIECH(IJ,1)
            SEAICER(IJ,NTLEV)=SEAICE(IJ  )
            QDBLR  (IJ,NTLEV)=QDBL  (IJ,1)
         ENDDO
C
         IF (LSICED) THEN
            DO IJ=1,IEJE
               SICEDR(IJ,NTLEV)=SICED(IJ)
            ENDDO
         ENDIF
C
         RETURN
C
      ENDIF
C
C     CHECK-VEKTOR FUER DIE FELDER MIT .FALSE. VORBESETZEN
      DO NTAB  = 1,NZVR
         LVARDA(NTAB) = .FALSE.
      ENDDO
C
C     GESAMTZAHL DER ZU LESENDEN FLAECHEN BERECHNEN
      CALL GNZFLD('R', NZGFELD)
      NZRFELD = 0
C
C
C     DATEINAMEN ERZEUGEN, DATEINAME WIRD IN 'YRDNAM'
C     (CB *CORG*) GESCHRIEBEN
      IF (MYID .EQ. 0) THEN
         CALL MAKEPN('R', NZEIT)
C
C     DATEI MIT DEN RANDDATEN HOLEN
         CALL GETD(NURDAT, YRDNAM, YRDCAT)
C
C     FELDER EINLESEN
         REWIND(NURDAT)
      ENDIF

C     LOOP OVER INPUT FIELDS
      DO

         IF (MYID .EQ. 0) THEN
            CALL DATEIN(NURDAT, UFELD, MOIEJE, ZAK, ZBK, KEMAX, ISTAT)
         ENDIF

         IF (ISTAT.EQ.0) THEN
            IF (MYID .EQ. 0) THEN

C
C           IE, JE, KE AUS GDB VERGLEICHEN MIT DEN WERTEN IN CB *PARAM*
C           POLPHI, POLLAM, DPHI, DLAM, PHILU, RLALU AUS GDB VERGLEICHEN
C           MIT DEN WERTEN IN CB *GRID*
               CALL GRIDCHK(MOIE,MOJE,KE,ISTAT)
               IF (ISTAT.NE.0) THEN
                  CALL REMARK( 'ECRPREP: FALSCHE RANDDATEN' )
                  STOP
               ENDIF
C
C              DATUM UND VERSIONSNUMMER AUS PDB VERGLEICHEN MIT DEN WERTEN IN
C              CB *CORG*, *ORG*
               CALL DATCHK('R', NZEIT, ISTAT)
               IF (ISTAT.NE.0) THEN
                  CALL REMARK( 'ECRPREP: FALSCHE RANDDATEN' )
                  STOP
               ENDIF
C
C              PLATZ DES FELDES IN DER TABELLE DER EM-VARIABLEN (CB *EMGRIB*)
C              SUCHEN, WENN NICHT BENOETIGT, NAECHSTES FELD EINLESEN
               CALL GETLOC(IVLOC)
               IF ( IVLOC.EQ. -1 ) CYCLE

               CALL SENDSUBGRID(IVLOC,IPDB,UFELD,VFELD)

            ELSE

               TAGCOUNT = 2
               CALL PTEST
               IF (TYPE .EQ. 1) THEN
                  CALL PRECVR(VFELD(1))
               ELSE IF (TYPE .EQ. 2) THEN
                  CALL PRECVI(ISTAT)
               ENDIF
               IF (ISTAT .NE. 0) CYCLE
            ENDIF

C           SYNCHRONISATION ZWISCHEN DEM EINLESEN DER EINZELNEN FLAECHEN
            CALL PSTOP

            IF (MYID .NE. 0) THEN
               IVLOC=INT(VFELD(IEJE+38))
               DO I=1,37
                  IPDB(I)=INT(VFELD(IEJE+I))
               ENDDO
            ENDIF
C
C           FELD AUF EM-EINHEITEN SKALIEREN UND IN ENTSPRECHENDES EM-FELD
C           KOPIEREN.
            CALL PUTECR
     &  (YTYP   , VFELD , IEJE, KE , KAKE(IVLOC,1), KAKE(IVLOC,2),
     &   IVLOC  , YRVARN, LVARDA, NZVR  , NZRFELD, UR           ,
     &   VR     , TR    , QDR   , QWR   , TSWECHR, TSIECHR      ,
     &   SEAICER, PSR   , QDBLR , SICEDR, NRD1   , NRD2)

C     NAECHSTES FELD EINLESEN
            CYCLE
C
C     DATEI IST FERTIG GELESEN
         ELSE
            EXIT
         ENDIF
      ENDDO

      IF (MYID .EQ. 0) THEN
         COUNT=1
         TYPE=2
         CALL PSENDALL2I(ISTAT)
      ENDIF
      CALL PSTOP
C
C     UEBERPRUEFEN, OB ALLE ZUR RECHNUNG NOETIGEN FELDER DA SIND
      CALL DCOMPLT('R', YRVARN, LVARDA, NZVR, NZGFELD, NZRFELD)
C
C     LOKALEN SPEICHER ZURUECKGEBEN
C
C     DATEI NURDAT ZURUECKGEBEN
      IF (MYID .EQ. 0) THEN
         CLOSE(NURDAT)
      ENDIF
C
      IF (YTYP(2:2).EQ.'1') THEN
         NTLEV = NRD1
      ELSE
         NTLEV = NRD2
      ENDIF
C
C     WENN QW NICHT IN DEN RANDFELDERN GEGEBEN IST, QWR=0. SETZEN
C
      IF (.NOT.LQWR) THEN
         DO IJ=1,IEJEKE
            QWR(IJ,NTLEV)=0.
         ENDDO
      ENDIF
C
      RETURN
      END SUBROUTINE ECRPREP
