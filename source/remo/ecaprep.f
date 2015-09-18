C
C     SUBROUTINE ECAPREP
C
C**** ECAPREP  -   UP:BEREITSTELLUNG DER ANFANGSDATEN IM LANGZEIT-
C****                 SPEICHER
C**   AUFRUF   :   CALL ECAPREP IN UP *EC4ORG*
C**   ENTRIES  :   KEINE
C**   ZWECK    :   BEREITSTELLUNG DER ANFANGSDATEN FUER DIE PROGNOSE;
C**                HOLEN UND LESEN DER ANFANGSDATEI UND BERECHNUNG
C**                DER KONSTANTEN FELDER.
C**   VERSIONS-
C**   DATUM    :   03.03.05
C**                2007
C**
C**   EXTERNALS:   MAKEPN , GETD   , GNZFLD , DATEIN , GRIDCHK, DATCHK ,
C**                GETAKBK, GETLOC , PUTECA ,
C**                DCOMPLT, ADJKAKE, GEOPOT ,
C**                GKONST , PRCV   , PRIV
C**
C**   EINGABE-
C**   PARAMETER:   KEINE
C**   AUSGABE-
C**   PARAMETER:   KEINE
C**
C**   COMMON-
C**   BLOECKE  :   ORG, CORG, COMDYN, EMGBCH, EMGBRI,
C**                COMPHY, UNITCH, UNITNR
C**
C**   METHODE  :   SEQUENTIELLER AUFRUF DER UNTERPROGRAMME ZUM HOLEN
C**                LESEN, PRUEFEN UND VERARBEITEN DER DATEN
C**   FEHLERBE-
C**   HANDLUNG :   STOP/ABORT IM FEHLERFALL
C**   VERFASSER:   R.PODZUN
C
      SUBROUTINE ECAPREP
     &    (
     &     U      , V      , T      , QD    , QW    , PS    , SN    ,
     &     TSLECH , TSWECH , TSIECH , USTRL , USTRW , USTRI , VSTRL ,
     &     VSTRW  , VSTRI  , EVAPL  , EVAPW , EVAPI ,
     &     AHFSL  , AHFSW  , AHFSI  , AZ0L  , AZ0W  , AZ0I  ,
     &     ALSOL  , ALSOW  , ALSOI  , AHFICE, QRES  ,
     &     TMCHL  , TMCHW  , TMCHI  , BFLHSL, BFLHSW, BFLHSI,
     &     BFLQDSL, BFLQDSW, BFLQDSI, QDBL  , QDBW  , QDBI  , WSECH ,
     &     TSECH  , TD     , TDCL   , WL    , TD3   , TD4   , TD5   ,
     &     TSN    , QDB    , TKE    , BLA   , AZ0   , FIB   , GLAC  ,
     &     VGRAT  , FOREST , ALBECH , WSMX  , VLT   , FAO   , BFLHS ,
     &     BFLQDS , RUNOFF , TMCM   , TMCH  , VERVEL, EMTEF , TRSOF ,
     &     APRL   , APRC   , APRS   , ACLC  , ACLCV , ALBEDO, EVAP  ,
     &     VAROR  , T2MAX  , T2MIN  , TSMAX , TSMIN , WIMAX , SEAICE,
     &     SICED  , TOPMAX , RGCGN  , EVAPM , DSNAC , RLA   ,
     &     PHI    , AHFS   , AHFL   , VDIS  , USTAR3, ACLCOV, U10   ,
     &     V10    , TEMP2  , DEW2   , TSURF , WIND10, SRADS , TRADS ,
     &     USTR   , SRAD0  , TRAD0  , VSTR  , SRAFS , TRAFS , USTRGW,
     &     SRAF0  , TRAF0  , VSTRGW , SCLFS , TCLFS , VDISGW, SCLF0 ,
     &     TCLF0  , SRAD0U , SRADSU , TRADSU, SNMEL , TSLIN , ACLCAC,
     &     QVI    , ALWCVI , FI     , TEFF  , EMTER , TRSOL , FC    ,
     &     CPHI   , ACPHIR , RMY    , AK    , BK    , AKH   , BKH   ,
     &     DAK    , DBK    , SISTM  , SIVMT , SIVMTI, SICQ  , SIGAM ,
     &     SITAU  , SINUE  , A1T    , A2T   , VVFH  , TRIGSI,
     &     TRIGSJ , RZ1I   , RZ2I   , RZ1J  , RZ2J  , IFAXI , IFAXJ ,
     &     TMKVMH , DRAIN  , SRFL   , QDBOXS, QWBOXS, EKBOXS, FHBOXS,
     &     FIBOXS , TLAMBDA, DLAMBDA, PORVOL, FCAP  , WI3   , WI4   ,
     &     WI5    , WI     , WICL   , BETA  ,WMINLOK,WMAXLOK, VBM10M,
     &     CAPE   , WS1    , WS2    , WS3   , WS4   , WS5   , DZR   ,
     &     DZS    , FKSAT  , FMPOT  , BCLAPP, VPOR  , ETRANS, EBSOIL,
     &     ESNOW  , ESKIN  , ERES   , QI    , QIVI  , QIBOXS,
     &     GCPHI  , GACPHIR,
     &     PINT   , DWDT   , W      , RPRAC)
C
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
      INCLUDE "org.h"
      INCLUDE "corg.h"
      INCLUDE "emgbch.h"
      INCLUDE "emgbri.h"
      INCLUDE "unitch.h"
      INCLUDE "unitnr.h"
      INCLUDE "comecphy.h"

      INCLUDE "comdyn.h"
      INCLUDE "comphy.h"
      INCLUDE "grib.h"
C
      INTEGER, PARAMETER :: KEMAX=50
C
C     PUTECA
C
      REAL ::
     & U     (IEJE,KE,3), V     (IEJE,KE,3), T     (IEJE,KE,3),
     & QD    (IEJE,KE,3), QW    (IEJE,KE,3), PS    (IEJE,   3),
     & SN    (IEJE,   3), WSECH (IEJE,   3), TSECH (IEJE,   3),
     & TSLECH(IEJE,   3), TSWECH(IEJE,   3), TSIECH(IEJE,   3),
     & TD    (IEJE,   3), TDCL  (IEJE,   3), WL    (IEJE,   3),
     & TD3   (IEJE,   3), TD4   (IEJE,   3), TD5   (IEJE,   3),
     & TSN   (IEJE,   3), QDB   (IEJE,   3), TKE   (IEJE,KE,3),
     & QDBL  (IEJE,   3), QDBW  (IEJE,   3), QDBI  (IEJE,   3),
     & USTRL (IEJE     ), USTRW (IEJE     ), USTRI (IEJE     ),
     & VSTRL (IEJE     ), VSTRW (IEJE     ), VSTRI (IEJE     ),
     & EVAPL (IEJE     ), EVAPW (IEJE     ), EVAPI (IEJE     ),
     & AHFSL (IEJE     ), AHFSW (IEJE     ), AHFSI (IEJE     ),
     & AZ0L  (IEJE     ), AZ0W  (IEJE     ), AZ0I  (IEJE     ),
     & ALSOL (IEJE     ), ALSOW (IEJE     ), ALSOI (IEJE     ),
     & AHFICE(IEJE     ), QRES  (IEJE     ), SRFL  (IEJE     ),
     & TMCHL (IEJE     ), TMCHW (IEJE     ), TMCHI (IEJE     ),
     & BFLHSL(IEJE     ), BFLHSW(IEJE     ), BFLHSI(IEJE     ),
     & BFLQDSL(IEJE    ), BFLQDSW(IEJE    ), BFLQDSI(IEJE    ),
     & BLA   (IEJE     ), AZ0   (IEJE     ), FIB   (IEJE     ),
     & GLAC  (IEJE     ), VGRAT (IEJE     ), FOREST(IEJE     ),
     & ALBECH(IEJE     ), WSMX  (IEJE     ), VLT   (IEJE     ),
     & FAO   (IEJE     ), BFLHS (IEJE     ), BFLQDS(IEJE     ),
     & RUNOFF(IEJE     ), TMCM  (IEJE     ), TMCH  (IEJE     ),
     & VERVEL(IEJE,KE  ), EMTEF (IEJE,2   ), TRSOF (IEJE,2   ),
     & APRL  (IEJE     ), APRC  (IEJE     ), APRS  (IEJE     ),
     & ACLC  (IEJE,KE  ), ACLCV (IEJE     ), ALBEDO(IEJE     ),
     & EVAP  (IEJE     ), VAROR (IEJE     ), T2MAX (IEJE     ),
     & T2MIN (IEJE     ), TSMAX (IEJE     ), TSMIN (IEJE     ),
     & WIMAX (IEJE     ), SEAICE(IEJE     ), SICED (IEJE     ),
     & TOPMAX(IEJE     ), RGCGN (IEJE     )
      REAL ::
     & EVAPM (IEJE     ), DSNAC (IEJE     ), QDBOXS(IEJE     ),
     & QWBOXS(IEJE     ), EKBOXS(IEJE     ), FHBOXS(IEJE     ),
     & FIBOXS(IEJE     ), TLAMBDA(IEJE    ), DLAMBDA(IEJE     ),
     & PORVOL(IEJE     ), FCAP  (IEJE     ), WI3   (IEJE,   3),
     & WI4   (IEJE,   3), WI5   (IEJE,   3), WI    (IEJE,   3),
     & WICL  (IEJE,   3), BETA  (IEJE     ), WMINLOK(IEJE     ),
     & WMAXLOK(IEJE    ), VBM10M(IEJE     ), CAPE  (IEJE     )
      REAL ::
     & RLA   (IEJE     ), PHI   (IEJE     ), AHFS  (IEJE     ),
     & AHFL  (IEJE     ), VDIS  (IEJE     ), USTAR3(IEJE     ),
     & ACLCOV(IEJE     ), U10   (IEJE     ), V10   (IEJE     ),
     & TEMP2 (IEJE     ), DEW2  (IEJE     ), TSURF (IEJE     ),
     & WIND10(IEJE     ), SRADS (IEJE     ), TRADS (IEJE     ),
     & USTR  (IEJE     ), SRAD0 (IEJE     ), TRAD0 (IEJE     ),
     & VSTR  (IEJE     ), SRAFS (IEJE     ), TRAFS (IEJE     ),
     & USTRGW(IEJE     ), SRAF0 (IEJE     ), TRAF0 (IEJE     ),
     & VSTRGW(IEJE     ), SCLFS (IEJE     ), TCLFS (IEJE     ),
     & VDISGW(IEJE     ), SCLF0 (IEJE     ), TCLF0 (IEJE     ),
     & SRAD0U(IEJE     ), SRADSU(IEJE     ), TRADSU(IEJE     ),
     & SNMEL (IEJE     ), TSLIN (IEJE     ), ACLCAC(IEJE,KE  ),
     & QVI   (IEJE     ), ALWCVI(IEJE     ), FI    (IEJE,KE,2),
     & TEFF  (IEJE     ), EMTER (IEJE,KE1 ),
     & TRSOL (IEJE,KE1 ), DRAIN (IEJE     )
      REAL ::
     & QI    (IEJE,KE,3), QIVI  (IEJE     ), QIBOXS(IEJE     ),
     & RPRAC (IEJE,KE)
C
      REAL ::
     & WS1   (IEJE     ), WS2   (IEJE     ), WS3   (IEJE     ),
     & WS4   (IEJE     ), WS5   (IEJE     ), DZR   (IEJE     ),
     & DZS   (IEJE     ), FKSAT (IEJE     ), FMPOT (IEJE     ),
     & BCLAPP(IEJE     ), VPOR  (IEJE     ), ETRANS(IEJE     ),
     & EBSOIL(IEJE     ), ESNOW (IEJE     ), ESKIN (IEJE     ),
     & ERES  (IEJE     )

      REAL ::   PINT(IEJE,KE1,3), DWDT(IEJE,KE ,3), 
     &          W   (IEJE,KE1,3)
C
C     GKONST
C
      REAL ::
     &          FC    (IE,JE), CPHI  (JE, 2),
     &          ACPHIR(JE, 2), RMY   (IE,JE,3)

      REAL ::   GCPHI(MOJE,2), GACPHIR(MOJE,2)

      REAL ::
     &          AK    (KE1  ), BK    (KE1  ),
     &          AKH   (KE   ), BKH   (KE   ),
     &          DAK   (KE   ), DBK   (KE   ),
     &          SISTM (KE,KE), SIVMT (KE,KE),
     &          SIVMTI(KE,KE), SICQ  (KE   ),
     &          SIGAM (KE,KE), SITAU (KE,KE),
     &          SINUE (KE   )

      REAL ::   A1T   (KE1  ), A2T   (KE1  )

      REAL ::   VVFH  (KE   )

      REAL ::
     &          TRIGSI(5*(MOIE-1)/2), TRIGSJ(5*(MOJE-1)/2),
     &          RZ1I(6)  , RZ2I(6),
     &          RZ1J(6)  , RZ2J(6)
      INTEGER ::IFAXI(11), IFAXJ(11)
C
      REAL ::   TMKVMH(IE*(KE-1),JE,2)
C
C     LOKALE FELDER
C
      REAL ::
     &          FTKVM (IEJE,KE)     , FTKVH (IEJE,KE)     ,
     &          ZFTKVM(IE*(KE-1),JE), ZFTKVH(IE*(KE-1),JE)
C
      REAL ::   ZAK(KEMAX), ZBK(KEMAX)
C
      LOGICAL ::LVARDA(NZMXVA), LAKBK
C
      REAL ::   UFELD(MOIEJE)
      REAL ::   VFELD(IEJE+38)
C
      INTEGER :: I,IJ,ISTAT,J,K,NTAB,NZAFELD,NZGFELD,IVLOC
C
      LAKBK=.FALSE.
      CALL SETRA(FTKVM ,IEJEKE     ,0.0)
      CALL SETRA(FTKVH ,IEJEKE     ,0.0)
      CALL SETRA(ZFTKVM,IEJE*(KE-1),0.0)
      CALL SETRA(ZFTKVH,IEJE*(KE-1),0.0)
      CALL SETRA(ZAK   ,KEMAX      ,0.0)
      CALL SETRA(ZBK   ,KEMAX      ,0.0)
C
      ISTAT=0
      TAGCOUNT=0
      TAGTABLE(1)=1
      TAGTABLE(2)=2
C
C     AKTUELLE ANZAHL DER ANFANGSVARIABLEN  (NZVA):
      IF (NZT.EQ.0) THEN
         IF (L5LAY) THEN
            NZVA=45
            IF (LSICED) THEN
               NZVA=46
               YAVARN(NZVA)='SICED   '
            ENDIF
         ELSE
            NZVA=34
            IF (LSICED) THEN
               NZVA=35
               YAVARN(NZVA)='SICED   '
            ENDIF
         ENDIF
      ELSE
         NZVA=NZMXVA
      ENDIF
C
C     TABELLE KAKE (CB *EMGRIB*) KORRIGIEREN, WENN KE.NE.20
      IF (KE.NE.20) CALL ADJKAKE(KE)
C
C     CHECK-VEKTOR FUER DIE FELDER MIT .FALSE. VORBESETZEN
      DO NTAB  = 1,NZVA
         LVARDA(NTAB) = .FALSE.
      ENDDO
C
C     GESAMTZAHL DER ZU LESENDEN FLAECHEN BERECHNEN
      CALL GNZFLD('A', NZGFELD)
      NZAFELD = 0
C     WENN NZT=0, BERUECKSICHTIGEN, DASS EVT. KEIN WOLKENWASSER (QW)
C     IN DEN ANFANGSDATEN; (DANN WIRD MIT QW=0 BEGONNEN)
      IF (NZT.EQ.0) THEN
         IF (.NOT.LQWR) THEN
            DO NTAB = 1, NZVA
               IF (YAVARN(NTAB).EQ.'QW') THEN
                  LVARDA(NTAB) = .TRUE.
               ENDIF
            ENDDO
            NZGFELD=NZGFELD-KE
         ENDIF
      ENDIF
C
C     DATEINAMEN ERZEUGEN, DATEINAME WIRD IN 'YADNAM' UND 'YUADAT'
C     (CB *CORG*) GESCHRIEBEN
      IF (MYID .EQ. 0) THEN
         IF (NZT.GT.0) THEN
C
            CALL MAKEPN('F', NZT+1)
C
C           DATEI MIT DEN FORTSETZUNGSDATEN HOLEN
            CALL GETD(NUFDAT, YFDNAM, YFDCAT)
C
C           FELDER EINLESEN
            REWIND(NUFDAT)
C
         ELSE
C
            CALL MAKEPN('A', NZT)
C
C           DATEI MIT DEN ANFANGSDATEN HOLEN
            CALL GETD(NUADAT, YADNAM, YADCAT)
C
C           FELDER EINLESEN
            REWIND(NUADAT)
C
         ENDIF
      ENDIF

C     LOOP OVER INPUT FIELDS
      DO

         IF (MYID .EQ. 0) THEN
            IF (NZT.GT.0) THEN
               CALL DATEINDBL(NUFDAT, UFELD, MOIEJE, ZAK, ZBK, KEMAX,
     &              ISTAT)
            ELSE
               CALL DATEIN(NUADAT, UFELD, MOIEJE, ZAK, ZBK, KEMAX,
     &              ISTAT)
            ENDIF
         ENDIF
         IF (ISTAT.EQ.0) THEN
            IF (MYID .EQ. 0) THEN
C
C              IE, JE, KE AUS GDB VERGLEICHEN MIT DEN WERTEN IN CB *PARAM*
C              POLPHI, POLLAM, DPHI, DLAM, PHILU, RLALU AUS GDB VERGLEICHEN
C              MIT DEN WERTEN IN CB *GRID*
               CALL GRIDCHK(MOIE,MOJE,KE,ISTAT)
               IF (ISTAT.NE.0) THEN
                  CALL REMARK( 'ECAPREP: FALSCHE ANFANGSDATEN' )
                  STOP
               ENDIF
C
C              DATUM UND VERSIONSNUMMER AUS PDB VERGLEICHEN MIT DEN WERTEN IN
C              CB *CORG*, *ORG*
               IF (NZT.GT.0) THEN
                  CALL DATCHK('F', NZT, ISTAT)
               ELSE
                  CALL DATCHK('A', NZT, ISTAT)
               ENDIF
               IF (ISTAT.NE.0) THEN
                  CALL REMARK( 'ECAPREP: FALSCHE ANFANGSDATEN' )
                  STOP
               ENDIF
C
C              VERTIKALKOORDINATEN-PARAMETER (AK, BK) AUS GDB HOLEN UND AUF
C              DATEI 'YUAUFTR' AUSGEBEN
               IF (.NOT.LAKBK) THEN
                  LAKBK = .TRUE.
                  CALL GETAKBK(AK,BK,KE1,ZAK,ZBK,KEMAX)
               ENDIF
C
C              PLATZ DES FELDES IN DER TABELLE DER VARIABLEN (CB *EMGRIB*)
C              SUCHEN, WENN NICHT BENOETIGT, NEUES FELD EINLESEN
               CALL GETLOC(IVLOC)
               IF ( IVLOC.EQ. -1 ) CYCLE
C
               CALL SENDSUBGRID(IVLOC,IPDB,UFELD,VFELD)

            ELSE                !IF (MYID .EQ. 0)

               TAGCOUNT = 2
               CALL PTEST
               IF (TYPE .EQ. 1) THEN
                  CALL PRECVR(VFELD(1))
               ELSE IF (TYPE .EQ. 2) THEN
                  CALL PRECVI(ISTAT)
                  IF (ISTAT .NE. 0) CYCLE
                  IF (ISTAT .EQ. 0) THEN
                     CALL REMARK('SYNCHRONISATIONSFEHLER IN ECAPREP !')
                     STOP
                  ENDIF
               ENDIF
            ENDIF               !IF (MYID .EQ. 0)

C           SYNCHRONISATION ZWISCHEN DEM EINLESEN DER EINZELNEN FLAECHEN
            CALL PSTOP

            IF (MYID .NE. 0) THEN
               IVLOC=INT(VFELD(IEJE+38))
               DO I=1,37
                  IPDB(I)=INT(VFELD(IEJE+I))
               ENDDO
            ENDIF

C
C           FELD AUF REMO-EINHEITEN SKALIEREN UND IN ENTSPRECHENDES REMO-
C           FELD KOPIEREN.
            CALL PUTECA
     &       ('A1', VFELD, IEJE, KAKE(IVLOC,1), KAKE(IVLOC,2), KE, KE1,
     &        IVLOC  , YAVARN , LVARDA, NZVA  , NZAFELD,NZT    ,
     &        U      , V      , T     , QD    , QW    , PS     ,
     &        SN     , WSECH  , TSECH , TD    , TDCL  , WL     ,
     &        TSLECH , TSWECH , TSIECH, USTRL , USTRW , USTRI  ,
     &        VSTRL  , VSTRW  , VSTRI , EVAPL , EVAPW , EVAPI  ,
     &        AHFSL  , AHFSW  , AHFSI , AZ0L  , AZ0W  , AZ0I   ,
     &        ALSOL  , ALSOW  , ALSOI , AHFICE, QRES  , TMCHL  ,
     &        TMCHW  , TMCHI  , BFLHSL, BFLHSW, BFLHSI, BFLQDSL,
     &        BFLQDSW, BFLQDSI, QDBL  , QDBW  , QDBI  , SRFL   ,
     &        TD3    , TD4    , TD5   , TSN   , QDB   , TKE    ,
     &        BLA    , AZ0    , FIB   , GLAC  , VGRAT , FOREST ,
     &        ALBECH , WSMX   , VLT   , FAO   , BFLHS , BFLQDS ,
     &        RUNOFF , TMCM   , TMCH  , VERVEL, EMTEF , TRSOF  ,
     &        APRL   , APRC   , APRS  , ACLC  , ACLCV , ALBEDO ,
     &        EVAP   , VAROR  , T2MAX , T2MIN , TSMAX , TSMIN  ,
     &        WIMAX  , SEAICE , SICED , TOPMAX, RGCGN ,
     &        EVAPM  , DSNAC  , RLA   , PHI   , AHFS  , AHFL   ,
     &        VDIS   , USTAR3 , ACLCOV, U10   , V10   , DRAIN  ,
     &        TEMP2  , DEW2   , TSURF , WIND10, SRADS , TRADS  ,
     &        USTR   , SRAD0  , TRAD0 , VSTR  , SRAFS , TRAFS  ,
     &        USTRGW , SRAF0  , TRAF0 , VSTRGW, SCLFS , TCLFS  ,
     &        VDISGW , SCLF0  , TCLF0 , SRAD0U, SRADSU, TRADSU ,
     &        SNMEL  , TSLIN  , ACLCAC, QVI   , ALWCVI, FI     ,
     &        TEFF   , FTKVM  , FTKVH , EMTER , TRSOL , QDBOXS ,
     &        QWBOXS , EKBOXS , FHBOXS, FIBOXS,TLAMBDA, DLAMBDA,
     &        PORVOL , FCAP   , WI3   , WI4   , WI5   , WI     ,
     &        WICL   , BETA   ,WMINLOK,WMAXLOK, VBM10M, CAPE   ,
     &        WS1    , WS2    , WS3   , WS4   , WS5   , DZR    ,
     &        DZS    , FKSAT  , FMPOT , BCLAPP, VPOR  , ETRANS ,
     &        EBSOIL , ESNOW  , ESKIN , ERES  , QI    , QIVI   ,
     &        QIBOXS ,
     &        PINT   , DWDT   , W     , RPRAC)
C
C        NAECHSTES FELD EINLESEN
            CYCLE
C
C     DATEI IST FERTIG GELESEN
         ELSE
            EXIT
         ENDIF
      ENDDO
C
C
      IF (MYID .EQ. 0) THEN
         COUNT=1
         TYPE=2
         CALL PSENDALL2I(ISTAT)
      ENDIF
      CALL PSTOP

      COUNT=KE1*1
      CALL PSENDALLR(AK)
      CALL PSTOP

      COUNT=KE1*1
      CALL PSENDALLR(BK)
      CALL PSTOP
C
C     UEBERPRUEFEN, OB ALLE ZUR RECHNUNG NOETIGEN FELDER DA SIND
      CALL DCOMPLT('A', YAVARN, LVARDA, NZVA, NZGFELD, NZAFELD)
C
C     GEOPOTENTIAL FI BERECHNEN FUER ZEITPUNKT NA (=1)
      IF (NZT.EQ.0) THEN
        DO K  = 1, KE1
          DO IJ = 1, IEJE
            PINT(IJ,K,1) = AK(K) + BK(K) * PS(IJ,1)
            PINT(IJ,K,2) = AK(K) + BK(K) * PS(IJ,2)
            PINT(IJ,K,3) = PINT(IJ,K,2)
          ENDDO
        ENDDO
        CALL GEOPOT(1, 1,
     &              1,IE, 1, JE ,
     &              T  , QD, QW  , FI    , FIB ,
     &              QI , PINT    , DWDT)
      ENDIF
C
      IF (MYID .EQ. 0) THEN
         IF (NZT.GT.0) THEN
C
C           DATEI NUFDAT ZURUECKGEBEN
            CLOSE(NUFDAT)
         ELSE
C
C           DATEI NUADAT ZURUECKGEBEN
            CLOSE(NUADAT)
         ENDIF
      ENDIF
C
C-----FORTSETZUNGSLAUF--------------------------------------------------
C
C     WENN FORTSETZUNGSLAUF (NZT>0), 2.DATENSATZ HOLEN.
      IF (NZT.GT.0) THEN
C
C        CHECK-VEKTOR FUER DIE FELDER MIT .FALSE. VORBESETZEN
         DO NTAB  = 1,NZVA
            LVARDA(NTAB) = .FALSE.
         ENDDO
         NZAFELD = 0
C
C        DATEINAMEN ERZEUGEN, DATEINAME WIRD IN 'YADNAM' UND 'YUADAT'
C        (CB *CORG*) GESCHRIEBEN
         YFDNAM(8:8) = 'g'
C
         IF (MYID .EQ. 0) THEN
C           DATEI MIT DEN ANFANGSDATEN HOLEN
            CALL GETD(NUFDAT, YFDNAM, YFDCAT)
C
C           FELDER EINLESEN
            REWIND(NUFDAT)
         ENDIF

         ISTAT=0
C        LOOP OVER INPUT FIELDS
         DO

            IF (MYID .EQ. 0) THEN
               CALL DATEINDBL(NUFDAT,UFELD, MOIEJE, ZAK, ZBK, KEMAX,
     &              ISTAT)
            ENDIF

            IF (ISTAT.EQ.0) THEN
               IF (MYID .EQ. 0) THEN
C
C                 IE, JE, KE AUS GDB VERGLEICHEN MIT DEN WERTEN IN CB *PARAM*
C                 POLPHI, POLLAM, DPHI, DLAM, PHILU, RLALU AUS GDB VERGLEICHEN
C                 MIT DEN WERTEN IN CB *GRID*
                  CALL GRIDCHK(MOIE,MOJE,KE,ISTAT)
                  IF (ISTAT.NE.0) THEN
                     CALL REMARK( 'ECAPREP: FALSCHE ANFANGSDATEN' )
                     STOP
                  ENDIF
C
C                 DATUM UND VERSIONSNUMMER AUS PDB VERGLEICHEN MIT DEN WERTEN IN
C                 CB *CORG*, *ORG*
                  CALL DATCHK('G', NZT+1, ISTAT)
                  IF (ISTAT.NE.0) THEN
                     CALL REMARK( 'ECAPREP: FALSCHE ANFANGSDATEN' )
                     STOP
                  ENDIF
C
C                 PLATZ DES FELDES IN DER TABELLE DER EM-VARIABLEN (CB *EMGRIB*)
C                 SUCHEN
                  CALL GETLOC(IVLOC)
                  IF ( IVLOC.EQ. -1 ) CYCLE

                  CALL SENDSUBGRID(IVLOC,IPDB,UFELD,VFELD)

               ELSE

                  TAGCOUNT = 2
                  CALL PTEST
                  IF (TYPE .EQ. 1) THEN
                     CALL PRECVR(VFELD)
                  ELSE IF (TYPE .EQ. 2) THEN
                     CALL PRECVI(ISTAT)
                     IF (ISTAT .NE. 0) CYCLE
                     IF (ISTAT .EQ. 0) THEN
                        CALL REMARK(
     &                       'SYNCHRONISATIONSFEHLER IN ECAPREP !' )
                        STOP
                     ENDIF
                  ENDIF
               ENDIF

C           SYNCHRONISATION
               CALL PSTOP

               IF (MYID .NE. 0) THEN
                  IVLOC=INT(VFELD(IEJE+38))
                  DO I=1,37
                     IPDB(I)=INT(VFELD(IEJE+I))
                  ENDDO
               ENDIF
C
C
C              FELD AUF REMO-EINHEITEN SKALIEREN UND IN ENTSPRECHENDES REMO
C              FELD KOPIEREN.
               CALL PUTECA
     &       ('A2', VFELD, IEJE, KAKE(IVLOC,1), KAKE(IVLOC,2), KE, KE1,
     &        IVLOC  , YAVARN , LVARDA, NZVA  , NZAFELD,NZT    ,
     &        U      , V      , T     , QD    , QW    , PS     ,
     &        SN     , WSECH  , TSECH , TD    , TDCL  , WL     ,
     &        TSLECH , TSWECH , TSIECH, USTRL , USTRW , USTRI  ,
     &        VSTRL  , VSTRW  , VSTRI , EVAPL , EVAPW , EVAPI  ,
     &        AHFSL  , AHFSW  , AHFSI , AZ0L  , AZ0W  , AZ0I   ,
     &        ALSOL  , ALSOW  , ALSOI , AHFICE, QRES  , TMCHL  ,
     &        TMCHW  , TMCHI  , BFLHSL, BFLHSW, BFLHSI, BFLQDSL,
     &        BFLQDSW, BFLQDSI, QDBL  , QDBW  , QDBI  , SRFL   ,
     &        TD3    , TD4    , TD5   , TSN   , QDB   , TKE    ,
     &        BLA    , AZ0    , FIB   , GLAC  , VGRAT , FOREST ,
     &        ALBECH , WSMX   , VLT   , FAO   , BFLHS , BFLQDS ,
     &        RUNOFF , TMCM   , TMCH  , VERVEL, EMTEF , TRSOF  ,
     &        APRL   , APRC   , APRS  , ACLC  , ACLCV , ALBEDO ,
     &        EVAP   , VAROR  , T2MAX , T2MIN , TSMAX , TSMIN  ,
     &        WIMAX  , SEAICE , SICED , TOPMAX, RGCGN ,
     &        EVAPM  , DSNAC  , RLA   , PHI   , AHFS  , AHFL   ,
     &        VDIS   , USTAR3 , ACLCOV, U10   , V10   , DRAIN  ,
     &        TEMP2  , DEW2   , TSURF , WIND10, SRADS , TRADS  ,
     &        USTR   , SRAD0  , TRAD0 , VSTR  , SRAFS , TRAFS  ,
     &        USTRGW , SRAF0  , TRAF0 , VSTRGW, SCLFS , TCLFS  ,
     &        VDISGW , SCLF0  , TCLF0 , SRAD0U, SRADSU, TRADSU ,
     &        SNMEL  , TSLIN  , ACLCAC, QVI   , ALWCVI, FI     ,
     &        TEFF   , FTKVM  , FTKVH , EMTER , TRSOL , QDBOXS ,
     &        QWBOXS , EKBOXS , FHBOXS, FIBOXS,TLAMBDA, DLAMBDA,
     &        PORVOL , FCAP   , WI3   , WI4   , WI5   , WI     ,
     &        WICL   , BETA   ,WMINLOK,WMAXLOK, VBM10M, CAPE   ,
     &        WS1    , WS2    , WS3   , WS4   , WS5   , DZR    ,
     &        DZS    , FKSAT  , FMPOT , BCLAPP, VPOR  , ETRANS ,
     &        EBSOIL , ESNOW  , ESKIN , ERES  , QI    , QIVI   ,
     &        QIBOXS ,
     &        PINT   , DWDT   , W     , RPRAC)
C
C           NAECHSTES FELD EINLESEN
               CYCLE
C
C        DATEI IST FERTIG GELESEN
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
C        UEBERPRUEFEN, OB ALLE ZUR RECHNUNG NOETIGEN FELDER DA SIND
         CALL DCOMPLT('A', YAVARN, LVARDA, NZVA, NZGFELD, NZAFELD)
C
C
C        ZEITSCHRITT UM 1 ERHOEHEN, WEIL NZT UND NZT+1 ALS FORTSETZUNGS-
C        DATEN VORLIEGEN
         NZT = NZT + 1
         DO K     = 2,KE
            DO J     = 1,JE
               DO I     = 1,IE
                  IJ           = I + (J-1)*IE
                  ZFTKVM(I+(K-2)*IE,J) = FTKVM(IJ,K)
                  ZFTKVH(I+(K-2)*IE,J) = FTKVH(IJ,K)
               ENDDO
            ENDDO
         ENDDO
         DO J = 1,JE
            CALL COPYRE(ZFTKVM(1,J),TMKVMH(1,J,1),IE*(KE-1))
            CALL COPYRE(ZFTKVH(1,J),TMKVMH(1,J,2),IE*(KE-1))
         ENDDO
C
         IF (MYID .EQ. 0) THEN
C           DATEI NUADAT ZURUECKGEBEN
            CLOSE(NUFDAT)
C
         ENDIF
      ENDIF

C
C-----ENDE FORTSETZUNGSLAUF---------------------------------------------
C
C     GEOPOTENTIAL BERECHNEN FUER ZEITPUNKT NJ (=2)
      IF (NZT.EQ.0)
     &     CALL GEOPOT(2, 2,
     &                 1, IE, 1, JE ,
     &                 T  , QD, QW  , FI    , FIB ,
     &                 QI , PINT,     DWDT)
C
C
C     KONSTANTE FELDER BERECHNEN:
C     FC:     CORIOLISPARAMETER
C     RLA:    GEOGRAPHISCHE LAENGE DER EM-GITTERPUNKTE
C     PHI:    GEOGRAPHISCHE BREITE DER EM-GITTERPUNKTE
C     RMY:    RANDRELAXATIONSFAKTOR NACH DAVIES/KALLBERG
C     SISTM:  STRUKTUR-MATRIX(TRANSPONIERT) FUER SI-VERFAHREN
C     SIVMT:  MATRIX DER EIGENVEKTOREN VON SISTM
C     SIVMTI: INVERSE MATRIX DER EIGENVEKTOREN VON SISTM
C     SICQ:   VEKTOR DES QUADRATES DER PHASENGESCHWINDIGKEITEN

      CALL GKONST
     &  (PHI   , RLA  ,
     &   FC    , CPHI , ACPHIR, RMY  , AK    , BK    , AKH   ,
     &   BKH   , DAK  , DBK   , SISTM, SIVMT , SIVMTI, SICQ  ,
     &   SIGAM , SITAU, SINUE , A1T  , A2T   , VVFH  , TRIGSI,
     &   TRIGSJ, RZ1I , RZ2I  , RZ1J , RZ2J  , IFAXI , IFAXJ ,
     &   GCPHI , GACPHIR)

      END SUBROUTINE ECAPREP
