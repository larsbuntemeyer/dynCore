      SUBROUTINE PUTEC4
     & (NUNIT , YTYP  , IVLOC , NX    ,
     &  AK    , BK    ,
     &  U     , V     , T     , QD    , QW    , FI    , EMTER  , TRSOL ,
     &  VERVEL, FTKVM , FTKVH , ACLC  , ACLCAC, TKE   , EMTEF  , TRSOF ,
     &  PS    , QDB   , TSECH , WSECH , SN    , TD    , TDCL   , WL    ,
     &  TSLECH, TSWECH, TSIECH, USTRL , USTRW , USTRI , VSTRL  , VSTRW ,
     &  VSTRI , EVAPL , EVAPW , EVAPI , AHFSL ,
     &  AHFSW , AHFSI , AZ0L  , AZ0W  , AZ0I  , ALSOL , ALSOW  , ALSOI ,
     &  AHFICE, QRES  , TMCHL , TMCHW , TMCHI ,
     &  QDBL  , QDBW  , QDBI  , BFLHSL, BFLHSW, BFLHSI, BFLQDSL, SRFL  ,
     &  TSN   , TD3   , TD4   , TD5   , SRADS , TRADS , SRAD0  , TRAD0 ,
     &  APRL  , APRC  , APRS  , VDIS  , AHFS  , FIB   , BFLQDSI, BLA   ,
     &  USTAR3, RUNOFF, ACLCV , ACLCOV, TMCM  , TMCH  , BFLQDSW, AHFL  ,
     &  PHI   , RLA   , BFLHS , BFLQDS, U10   , V10   , TEMP2  , DEW2  ,
     &  TSURF , WIND10, AZ0   , ALBECH, ALBEDO, USTR  , VSTR   , EVAP  ,
     &  EVAPM , SRAFS , TRAFS , SRAF0 , TRAF0 , SCLFS , TCLFS  , SCLF0 ,
     &  TCLF0 , USTRGW, VSTRGW, VDISGW, VGRAT , VAROR , VLT    , T2MAX ,
     &  SRAD0U, SRADSU, TRADSU, T2MIN , SEAICE, SICED , FOREST , TEFF  ,
     &  TSMAX , TSMIN , WIMAX , TOPMAX, SNMEL , TSLIN , DSNAC  , FAO   ,
     &  RGCGN , WSMX  , QVI   , ALWCVI, GLAC  , DRAIN , QDBOXS ,
     &  QWBOXS, EKBOXS, FHBOXS, FIBOXS,TLAMBDA,DLAMBDA, PORVOL , FCAP  ,
     &  WI3   , WI4   , WI5   , WI    , WICL  , PSRED , U10ER  , V10ER ,
     &  GHPBL , BETA  ,WMINLOK,WMAXLOK, VBM10M, CAPE  , WS1    , WS2   ,
     &  WS3   , WS4   , WS5   , DZR   , DZS   , FKSAT , FMPOT  , BCLAPP,
     &  VPOR  , ETRANS, EBSOIL, ESNOW , ESKIN , ERES  , QI     , QIVI  ,
     &  QIBOXS, PINT  , DWDT  , W     , RPRAC)
C
C
C**** PUTEC4   - UP:AUSGABE DER FELDER AUF ERGEBNISDATEI (ECHAM4-PHYSIK)
C****
C**   AUFRUF   :   CALL PUTEC4
C**
C**   ENTRIES  :   KEINE
C**   ZWECK    :   AUSGABE DER FELDER AUF ERGEBNISDATEI (ECHAM4-PHYSIK)
C**
C**   VERSIONS-
C**   DATUM    :   05.10.04
C**                2007
C**
C**   EXTERNALS:   SETRA, WRITEGB
C**
C**   EINGABE-
C**   PARAMETER:   NUNIT: UNIT-NUMMER DER AUSGABE-DATEI
C**                YTYP : TYP DER ERGEBNISDATEN:
C**                      'E': NORMALE ERGEBNISDATEN FUER GESAMTGEBIET
C**                      'D':         ERGEBNISDATEN FUER TEILGEBIET
C**                      'F': FORTSETZUNGSDATEN     FUER GESAMTGEBIET
C**                      'T': TRAJEKTORIEN-DATEN    FUER GESAMTGEBIET
C**                IVLOC: PLATZ DER VARIABLEN IN DER TABELLE *YEMNAME*
C**                NX   : GEWUENSCHTE ZEITEBENE FUER FELDER MIT DREI
C**                       ZEITEBENEN
C**   AUSGABE-
C**   PARAMETER:   KEINE
C**
C**   COMMON-
C**   BLOECKE  :   ORG, GRIB, EMGBCH, EMGBRI
C**
C**   METHODE  :   GEFUNDENES FELD AUSGEBEN
C**
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   R. PODZUN
C
C     ------------------------------------------------------------------
      IMPLICIT NONE
C     
      INCLUDE "parorg.h"
      INCLUDE "org.h"
      INCLUDE "grib.h"
      INCLUDE "emgbch.h"
      INCLUDE "emgbri.h"
      INCLUDE "comdyn.h"
C     ------------------------------------------------------------------
C
C     Dummy Arguments and Arrays
C
      INTEGER,   INTENT(IN) :: NUNIT, IVLOC, NX
      CHARACTER, INTENT(IN) :: YTYP*(*)
C
      REAL, INTENT(IN) ::
     & U     (IEJE,KE,3), V     (IEJE,KE,3), T     (IEJE,KE,3),
     & QD    (IEJE,KE,3), QW    (IEJE,KE,3), FI    (IEJE,KE,2),
     & VERVEL(IEJE,KE  ), FTKVM (IEJE,KE  ), FTKVH (IEJE,KE  ),
     & ACLC  (IEJE,KE  ), ACLCAC(IEJE,KE  ), TKE   (IEJE,KE,3),
     & EMTER (IEJE,KE1 ), TRSOL (IEJE,KE1 ), EMTEF (IEJE,2   ),
     & TRSOF (IEJE,2   )
C
      REAL, INTENT(IN) ::
     & PS    (IEJE,3), QDB   (IEJE,3), TSECH (IEJE,3), WSECH(IEJE,3),
     & TSLECH(IEJE,3), TSWECH(IEJE,3), TSIECH(IEJE,3), QDBL (IEJE,3),
     & QDBW  (IEJE,3), QDBI  (IEJE,3), SN    (IEJE,3), TD   (IEJE,3),
     & TDCL  (IEJE,3), WL    (IEJE,3), TSN   (IEJE,3), TD3  (IEJE,3),
     & TD4   (IEJE,3), TD5   (IEJE,3)
C
      REAL, INTENT(IN) ::
     & USTRL  (IEJE), USTRW (IEJE), USTRI (IEJE), VSTRL (IEJE),
     & VSTRW  (IEJE), VSTRI (IEJE), EVAPL (IEJE), EVAPW (IEJE),
     & EVAPI  (IEJE), AHFSL (IEJE), AHFSW (IEJE), AHFSI (IEJE),
     & AZ0L   (IEJE), AZ0W  (IEJE), AZ0I  (IEJE), ALSOL (IEJE),
     & ALSOW  (IEJE), ALSOI (IEJE), AHFICE(IEJE), QRES  (IEJE),
     & TMCHL  (IEJE), TMCHW (IEJE), WSMX  (IEJE), QVI   (IEJE),
     & TMCHI  (IEJE), BFLHSL(IEJE), BFLHSW(IEJE), BFLHSI(IEJE),
     & BFLQDSW(IEJE), TRADS (IEJE), SRAD0 (IEJE), TRAD0 (IEJE),
     & APRL   (IEJE), APRC  (IEJE), APRS  (IEJE), VDIS  (IEJE),
     & BFLQDSL(IEJE), FIB   (IEJE), BLA   (IEJE), AHFL  (IEJE),
     & BFLQDSI(IEJE), RUNOFF(IEJE), ACLCV (IEJE), ACLCOV(IEJE),
     & TMCM   (IEJE), TMCH  (IEJE), SRADS (IEJE), USTAR3(IEJE),
     & PHI    (IEJE), RLA   (IEJE), BFLHS (IEJE), BFLQDS(IEJE),
     & U10    (IEJE), V10   (IEJE), TEMP2 (IEJE), DEW2  (IEJE),
     & TSURF  (IEJE), WIND10(IEJE), AZ0   (IEJE), ALBECH(IEJE),
     & ALBEDO (IEJE), USTR  (IEJE), VSTR  (IEJE), EVAP  (IEJE),
     & EVAPM  (IEJE), SRAFS (IEJE), TRAFS (IEJE), SRAF0 (IEJE),
     & TRAF0  (IEJE), SCLFS (IEJE), TCLFS (IEJE), SCLF0 (IEJE),
     & TCLF0  (IEJE), USTRGW(IEJE), VSTRGW(IEJE), VDISGW(IEJE),
     & VGRAT  (IEJE), VAROR (IEJE), VLT   (IEJE), T2MAX (IEJE),
     & SRAD0U (IEJE), SRADSU(IEJE), TRADSU(IEJE), T2MIN (IEJE),
     & SEAICE (IEJE), SICED (IEJE), FOREST(IEJE), TEFF  (IEJE),
     & TSMAX  (IEJE), TSMIN (IEJE), WIMAX (IEJE), TOPMAX(IEJE),
     & SNMEL  (IEJE), TSLIN (IEJE), DSNAC (IEJE), FAO   (IEJE),
     & RGCGN  (IEJE), ALWCVI(IEJE), GLAC  (IEJE),
     & DRAIN  (IEJE), AHFS  (IEJE), SRFL  (IEJE), QDBOXS(IEJE),
     & QWBOXS (IEJE), EKBOXS(IEJE), FHBOXS(IEJE), FIBOXS(IEJE),
     & TLAMBDA(IEJE),DLAMBDA(IEJE), PORVOL(IEJE), FCAP  (IEJE),
     & WI3  (IEJE,3), WI4 (IEJE,3), WI5 (IEJE,3), WI  (IEJE,3),
     & WICL (IEJE,3), PSRED (IEJE), U10ER (IEJE), V10ER (IEJE),
     & GHPBL  (IEJE), BETA  (IEJE),WMINLOK(IEJE),WMAXLOK(IEJE),
     & VBM10M (IEJE), CAPE  (IEJE)
C
      REAL, INTENT(IN) ::
     & WS1    (IEJE), WS2   (IEJE), WS3   (IEJE), WS4   (IEJE),
     & WS5    (IEJE), DZR   (IEJE), DZS   (IEJE), FKSAT (IEJE),
     & FMPOT  (IEJE), BCLAPP(IEJE), VPOR  (IEJE), ETRANS(IEJE),
     & EBSOIL (IEJE), ESNOW (IEJE), ESKIN (IEJE), ERES  (IEJE)

      REAL ::   PINT(IEJE,KE1,3), DWDT(IEJE,KE ,3), 
     &          W   (IEJE,KE1,3)
C
      REAL, INTENT(IN) ::
     & AK(KE1), BK(KE1)
      REAL, INTENT(IN) ::
     & QI(IEJE,KE,3), QIVI (IEJE), QIBOXS (IEJE),
     & RPRAC(IEJE,KE)
C     ------------------------------------------------------------------
C
C     Local Declarations
C
      REAL    :: UFELD(IEJE)
      INTEGER :: IJ,K,NX2
C
C     ------------------------------------------------------------------
C     PDB FERTIG BESETZEN (ELEMENT-NUMMER UND LEVEL-TYP)
      IPDB( 7) = NEMGBNR(IVLOC)
      IF (NEMGBLT(IVLOC).EQ.1 ) THEN
          IPDB( 8) =   1
          IPDB( 9) =   0
          IPDB(10) =   0
C
      ELSE IF (NEMGBLT(IVLOC).EQ.102) THEN
          IPDB( 8) = 102
          IPDB( 9) =   0
          IPDB(10) =   0
C
      ELSE IF (NEMGBLT(IVLOC).EQ.105) THEN
          IPDB( 8) = 105
          IPDB( 9) =   0
          IF (YEMNAME(IVLOC).EQ. 'TEMP2  ' .OR.
     &        YEMNAME(IVLOC).EQ. 'DEW2   ' .OR.
     &        YEMNAME(IVLOC).EQ. 'T2MIN  ' .OR.
     &        YEMNAME(IVLOC).EQ. 'T2MAX  ' )   THEN
               IPDB(10) =   2
          ELSE IF (YEMNAME(IVLOC).EQ. 'U10    ' .OR.
     &             YEMNAME(IVLOC).EQ. 'V10    ' .OR.
     &             YEMNAME(IVLOC).EQ. 'U10ER  ' .OR.
     &             YEMNAME(IVLOC).EQ. 'V10ER  ' .OR.
     &             YEMNAME(IVLOC).EQ. 'WIND10 ' .OR.
     &             YEMNAME(IVLOC).EQ. 'VBM10M ' .OR.
     &             YEMNAME(IVLOC).EQ. 'WIMAX  ' )   THEN
               IPDB(10) =  10
          ENDIF
C
      ELSE IF (NEMGBLT(IVLOC).EQ.109) THEN
          IPDB( 8) = 109
          IPDB( 9) =   0
          IPDB(10) =   1
C
      ELSE IF (NEMGBLT(IVLOC).EQ.110) THEN
          IPDB( 8) = 110
      ENDIF
C
C     ZEITEBENE FESTSTELLEN FUER FELDER MIT ZWEI ZEITEBENEN
      IF (NX.EQ.NE) THEN
          NX2 = NA2
      ELSE
          NX2 = NJ2
      ENDIF
C
C     FELD AUS LANGZEIT-SPEICHER HOLEN; UNTERSCHEIDUNG BZGL. MULTI-
C     LEVEL, SINGLE-LEVEL, LANDPUNKT-FELDERN, ANZAHL DER ZEITEBENEN
C
C     1. SINGLE LEVEL FELDER
      IF (KAKE(IVLOC,1).EQ.KAKE(IVLOC,2)) THEN

C     3 ZEITEBENEN
C
          IF (NZEMTLV(IVLOC).EQ.3) THEN
C
          IF (YEMNAME(IVLOC).EQ.'PS       ') THEN
                  DO 110 IJ  = 1,IEJE
                  UFELD(IJ) = PS   (IJ,NX)
 110              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'QDB      ') THEN
                  DO 115 IJ  = 1,IEJE
                  UFELD(IJ) = QDB  (IJ,NX)
 115              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'QDBL     ') THEN
                  DO 116 IJ  = 1,IEJE
                  UFELD(IJ) = QDBL (IJ,NX)
 116              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'QDBW     ') THEN
                  DO 117 IJ  = 1,IEJE
                  UFELD(IJ) = QDBW (IJ,NX)
 117              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'QDBI     ') THEN
                  DO 118 IJ  = 1,IEJE
                  UFELD(IJ) = QDBI (IJ,NX)
 118              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TSECH    ') THEN
                  DO 120 IJ  = 1,IEJE
                  UFELD(IJ) = TSECH(IJ,NX)
 120              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TSLECH   ') THEN
                  DO 121 IJ  = 1,IEJE
                  UFELD(IJ) = TSLECH(IJ,NX)
 121              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TSWECH   ') THEN
                  DO 122 IJ  = 1,IEJE
                  UFELD(IJ) = TSWECH(IJ,NX)
 122              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TSIECH   ') THEN
                  DO 123 IJ  = 1,IEJE
                  UFELD(IJ) = TSIECH(IJ,NX)
 123              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'WSECH    ') THEN
                  DO 125 IJ  = 1,IEJE
                  UFELD(IJ) = WSECH(IJ,NX)
 125              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'WI3      ') THEN
                  DO 126 IJ  = 1,IEJE
                  UFELD(IJ) = WI3  (IJ,NX)
 126              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'WI4      ') THEN
                  DO 127 IJ  = 1,IEJE
                  UFELD(IJ) = WI4  (IJ,NX)
 127              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'WI5      ') THEN
                  DO 128 IJ  = 1,IEJE
                  UFELD(IJ) = WI5  (IJ,NX)
 128              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'WI       ') THEN
                  DO 129 IJ  = 1,IEJE
                  UFELD(IJ) = WI   (IJ,NX)
 129              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'WICL     ') THEN
                  DO 131 IJ  = 1,IEJE
                  UFELD(IJ) = WICL (IJ,NX)
 131              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'SN       ') THEN
                  DO 130 IJ  = 1,IEJE
                  UFELD(IJ) = SN   (IJ,NX)
 130              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TD       ') THEN
                  DO 135 IJ  = 1,IEJE
                  UFELD(IJ) = TD   (IJ,NX)
 135              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'WL       ') THEN
                  DO 140 IJ  = 1,IEJE
                  UFELD(IJ) = WL   (IJ,NX)
 140              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TDCL     ') THEN
                  DO 145 IJ  = 1,IEJE
                  UFELD(IJ) = TDCL (IJ,NX)
 145              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TSN      ') THEN
                  DO 150 IJ  = 1,IEJE
                  UFELD(IJ) = TSN  (IJ,NX)
 150              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TD3      ') THEN
                  DO 155 IJ  = 1,IEJE
                  UFELD(IJ) = TD3  (IJ,NX)
 155              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TD4      ') THEN
                  DO 160 IJ  = 1,IEJE
                  UFELD(IJ) = TD4  (IJ,NX)
 160              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TD5      ') THEN
                  DO 165 IJ  = 1,IEJE
                  UFELD(IJ) = TD5  (IJ,NX)
 165              CONTINUE
          ENDIF
C
C     2 ZEITEBENEN
          ELSE IF (NZEMTLV(IVLOC).EQ.2) THEN
CRP
C     *** MOMENTAN KEINES VERHANDEN ***
C
      IF (MYID .EQ. 0) THEN
      PRINT *,'ERROR PUTEC4: VERSUCH 1D-FELD MIT 2 ZEITEN ZU SCHREIBEN'
      ENDIF
      CALL EMABORT
C             CALL LOCATE (IZSLF2, YEMNAME(IVLOC), ITYPE)
C                 DO 80 IJ = 1,IEJE
C                 UFELD(IJ) = ZSLF2(IJ,NX2)
C  80             CONTINUE
C
C     1 ZEITEBENE
          ELSE IF (NZEMTLV(IVLOC).EQ.1) THEN
C
          IF (YEMNAME(IVLOC).EQ.'SRADS    ') THEN
                  DO 170 IJ  = 1,IEJE
                  UFELD(IJ) = SRADS (IJ)
 170              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TRADS    ') THEN
                  DO 175 IJ  = 1,IEJE
                  UFELD(IJ) = TRADS (IJ)
 175              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'SRAD0    ') THEN
                  DO 180 IJ  = 1,IEJE
                  UFELD(IJ) = SRAD0 (IJ)
 180              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TRAD0    ') THEN
                  DO 185 IJ  = 1,IEJE
                  UFELD(IJ) = TRAD0 (IJ)
 185              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'SRAFS    ') THEN
                  DO 190 IJ  = 1,IEJE
                  UFELD(IJ) = SRAFS (IJ)
 190              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TRAFS    ') THEN
                  DO 195 IJ  = 1,IEJE
                  UFELD(IJ) = TRAFS (IJ)
 195              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'SRAF0    ') THEN
                  DO 200 IJ  = 1,IEJE
                  UFELD(IJ) = SRAF0 (IJ)
 200              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TRAF0    ') THEN
                  DO 205 IJ  = 1,IEJE
                  UFELD(IJ) = TRAF0 (IJ)
 205              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'SCLFS    ') THEN
                  DO 210 IJ  = 1,IEJE
                  UFELD(IJ) = SCLFS (IJ)
 210              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TCLFS    ') THEN
                  DO 215 IJ  = 1,IEJE
                  UFELD(IJ) = TCLFS (IJ)
 215              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'SCLF0    ') THEN
                  DO 220 IJ  = 1,IEJE
                  UFELD(IJ) = SCLF0 (IJ)
 220              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TCLF0    ') THEN
                  DO 225 IJ  = 1,IEJE
                  UFELD(IJ) = TCLF0 (IJ)
 225              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'SEAICE   ') THEN
                  DO 230 IJ  = 1,IEJE
                  UFELD(IJ) = SEAICE(IJ)
 230              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'SICED    ') THEN
                  DO 235 IJ  = 1,IEJE
                  UFELD(IJ) = SICED (IJ)
 235              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'FOREST   ') THEN
                  DO 240 IJ  = 1,IEJE
                  UFELD(IJ) = FOREST(IJ)
 240              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TEFF     ') THEN
                  DO 245 IJ  = 1,IEJE
                  UFELD(IJ) = TEFF  (IJ)
 245              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TSMAX    ') THEN
                  DO 250 IJ  = 1,IEJE
                  UFELD(IJ) = TSMAX (IJ)
 250              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TSMIN    ') THEN
                  DO 255 IJ  = 1,IEJE
                  UFELD(IJ) = TSMIN (IJ)
 255              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'WIMAX    ') THEN
                  DO 260 IJ  = 1,IEJE
                  UFELD(IJ) = WIMAX (IJ)
 260              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TOPMAX   ') THEN
                  DO 265 IJ  = 1,IEJE
                  UFELD(IJ) = TOPMAX(IJ)
 265              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'SNMEL    ') THEN
                  DO 270 IJ  = 1,IEJE
                  UFELD(IJ) = SNMEL (IJ)
 270              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TSLIN    ') THEN
                  DO 275 IJ  = 1,IEJE
                  UFELD(IJ) = TSLIN (IJ)
 275              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'APRL     ') THEN
                  DO 280 IJ  = 1,IEJE
                  UFELD(IJ) = APRL  (IJ)
 280              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'APRC     ') THEN
                  DO 285 IJ  = 1,IEJE
                  UFELD(IJ) = APRC  (IJ)
 285              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'APRS     ') THEN
                  DO 290 IJ  = 1,IEJE
                  UFELD(IJ) = APRS  (IJ)
 290              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'VDIS     ') THEN
                  DO 295 IJ  = 1,IEJE
                  UFELD(IJ) = VDIS  (IJ)
 295              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'AHFS     ') THEN
                  DO 300 IJ  = 1,IEJE
                  UFELD(IJ) = AHFS  (IJ)
 300              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'FIB      ') THEN
                  DO 305 IJ  = 1,IEJE
                  UFELD(IJ) = FIB   (IJ)
 305              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'BLA      ') THEN
                  DO 310 IJ  = 1,IEJE
                  UFELD(IJ) = BLA   (IJ)
 310              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'AHFL     ') THEN
                  DO 315 IJ  = 1,IEJE
                  UFELD(IJ) = AHFL  (IJ)
 315              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'USTAR3   ') THEN
                  DO 320 IJ  = 1,IEJE
                  UFELD(IJ) = USTAR3(IJ)
 320              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'RUNOFF   ') THEN
                  DO 325 IJ  = 1,IEJE
                  UFELD(IJ) = RUNOFF(IJ)
 325              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'ACLCV    ') THEN
                  DO 330 IJ  = 1,IEJE
                  UFELD(IJ) = ACLCV (IJ)
 330              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'ACLCOV   ') THEN
                  DO 335 IJ  = 1,IEJE
                  UFELD(IJ) = ACLCOV(IJ)
 335              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TMCM     ') THEN
                  DO 340 IJ  = 1,IEJE
                  UFELD(IJ) = TMCM  (IJ)
 340              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TMCH     ') THEN
                  DO 345 IJ  = 1,IEJE
                  UFELD(IJ) = TMCH  (IJ)
 345              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'U10ER    ') THEN
                  DO 350 IJ  = 1,IEJE
                  UFELD(IJ) = U10ER (IJ)
 350              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'V10ER    ') THEN
                  DO 355 IJ  = 1,IEJE
                  UFELD(IJ) = V10ER (IJ)
 355              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'U10      ') THEN
                  DO 360 IJ  = 1,IEJE
                  UFELD(IJ) = U10   (IJ)
 360              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'V10      ') THEN
                  DO 365 IJ  = 1,IEJE
                  UFELD(IJ) = V10   (IJ)
 365              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'BFLHS    ') THEN
                  DO 370 IJ  = 1,IEJE
                  UFELD(IJ) = BFLHS (IJ)
 370              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'BFLQDS   ') THEN
                  DO 375 IJ  = 1,IEJE
                  UFELD(IJ) = BFLQDS(IJ)
 375              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TEMP2    ') THEN
                  DO 380 IJ  = 1,IEJE
                  UFELD(IJ) = TEMP2 (IJ)
 380              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'DEW2     ') THEN
                  DO 385 IJ  = 1,IEJE
                  UFELD(IJ) = DEW2  (IJ)
 385              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TSURF    ') THEN
                  DO 390 IJ  = 1,IEJE
                  UFELD(IJ) = TSURF (IJ)
 390              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'WIND10   ') THEN
                  DO 395 IJ  = 1,IEJE
                  UFELD(IJ) = WIND10(IJ)
 395              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'ALBECH   ') THEN
                  DO 400 IJ  = 1,IEJE
                  UFELD(IJ) = ALBECH(IJ)
 400              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'PHI      ') THEN
                  DO 405 IJ  = 1,IEJE
                  UFELD(IJ) = PHI   (IJ)
 405              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'RLA      ') THEN
                  DO 410 IJ  = 1,IEJE
                  UFELD(IJ) = RLA   (IJ)
 410              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'AZ0      ') THEN
                  DO 415 IJ  = 1,IEJE
                  UFELD(IJ) = AZ0   (IJ)
 415              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'ALBEDO   ') THEN
                  DO 420 IJ  = 1,IEJE
                  UFELD(IJ) = ALBEDO(IJ)
 420              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'USTR     ') THEN
                  DO 425 IJ  = 1,IEJE
                  UFELD(IJ) = USTR  (IJ)
 425              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'VSTR     ') THEN
                  DO 430 IJ  = 1,IEJE
                  UFELD(IJ) = VSTR  (IJ)
 430              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'EVAP     ') THEN
                  DO 435 IJ  = 1,IEJE
                  UFELD(IJ) = EVAP  (IJ)
 435              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'EVAPM    ') THEN
                  DO 440 IJ  = 1,IEJE
                  UFELD(IJ) = EVAPM (IJ)
 440              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'USTRGW   ') THEN
                  DO 445 IJ  = 1,IEJE
                  UFELD(IJ) = USTRGW(IJ)
 445              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'VSTRGW   ') THEN
                  DO 450 IJ  = 1,IEJE
                  UFELD(IJ) = VSTRGW(IJ)
 450              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'VDISGW   ') THEN
                  DO 455 IJ  = 1,IEJE
                  UFELD(IJ) = VDISGW(IJ)
 455              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'VGRAT    ') THEN
                  DO 460 IJ  = 1,IEJE
                  UFELD(IJ) = VGRAT (IJ)
 460              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'VAROR    ') THEN
                  DO 465 IJ  = 1,IEJE
                  UFELD(IJ) = VAROR (IJ)
 465              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'VLT      ') THEN
                  DO 470 IJ  = 1,IEJE
                  UFELD(IJ) = VLT   (IJ)
 470              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'T2MAX    ') THEN
                  DO 475 IJ  = 1,IEJE
                  UFELD(IJ) = T2MAX (IJ)
 475              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'T2MIN    ') THEN
                  DO 480 IJ  = 1,IEJE
                  UFELD(IJ) = T2MIN (IJ)
 480              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'SRAD0U   ') THEN
                  DO 485 IJ  = 1,IEJE
                  UFELD(IJ) = SRAD0U(IJ)
 485              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'SRADSU   ') THEN
                  DO 490 IJ  = 1,IEJE
                  UFELD(IJ) = SRADSU(IJ)
 490              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TRADSU   ') THEN
                  DO 495 IJ  = 1,IEJE
                  UFELD(IJ) = TRADSU(IJ)
 495              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'DSNAC    ') THEN
                  DO 500 IJ  = 1,IEJE
                  UFELD(IJ) = DSNAC (IJ)
 500              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'FAO      ') THEN
                  DO 505 IJ  = 1,IEJE
                  UFELD(IJ) = FAO   (IJ)
 505              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'RGCGN    ') THEN
                  DO 510 IJ  = 1,IEJE
                  UFELD(IJ) = RGCGN (IJ)
 510              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TLAMBDA  ') THEN
                  DO 511 IJ  = 1,IEJE
                  UFELD(IJ) =TLAMBDA(IJ)
 511              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'DLAMBDA  ') THEN
                  DO 512 IJ  = 1,IEJE
                  UFELD(IJ) =DLAMBDA(IJ)
 512              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'PORVOL   ') THEN
                  DO 513 IJ  = 1,IEJE
                  UFELD(IJ) = PORVOL(IJ)
 513              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'FCAP   ') THEN
                  DO 515 IJ  = 1,IEJE
                  UFELD(IJ) = FCAP  (IJ)
 515              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'WSMX     ') THEN
                  DO 520 IJ  = 1,IEJE
                  UFELD(IJ) = WSMX  (IJ)
 520              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'QVI      ') THEN
                  DO 525 IJ  = 1,IEJE
                  UFELD(IJ) = QVI   (IJ)
 525              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'ALWCVI   ') THEN
                  DO 530 IJ  = 1,IEJE
                  UFELD(IJ) = ALWCVI(IJ)
 530              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'GLAC     ') THEN
                  DO 535 IJ  = 1,IEJE
                  UFELD(IJ) = GLAC  (IJ)
 535              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'DRAIN    ') THEN
                  DO 537 IJ  = 1,IEJE
                  UFELD(IJ) = DRAIN (IJ)
 537              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'SRFL     ') THEN
                  DO 538 IJ  = 1,IEJE
                  UFELD(IJ) = SRFL  (IJ)
 538              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'PSRED    ') THEN
                  DO 539 IJ  = 1,IEJE
                  UFELD(IJ) = PSRED (IJ)
CHG                  UFELD(IJ) = 0
 539              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'USTRL    ') THEN
                  DO 700 IJ  = 1,IEJE
                  UFELD(IJ) = USTRL (IJ)
 700              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'USTRW    ') THEN
                  DO 705 IJ  = 1,IEJE
                  UFELD(IJ) = USTRW (IJ)
 705              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'USTRI    ') THEN
                  DO 710 IJ  = 1,IEJE
                  UFELD(IJ) = USTRI (IJ)
 710              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'VSTRL    ') THEN
                  DO 715 IJ  = 1,IEJE
                  UFELD(IJ) = VSTRL (IJ)
 715              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'VSTRW    ') THEN
                  DO 720 IJ  = 1,IEJE
                  UFELD(IJ) = VSTRW (IJ)
 720              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'VSTRI    ') THEN
                  DO 725 IJ  = 1,IEJE
                  UFELD(IJ) = VSTRI (IJ)
 725              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'EVAPL    ') THEN
                  DO 730 IJ  = 1,IEJE
                  UFELD(IJ) = EVAPL (IJ)
 730              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'EVAPW    ') THEN
                  DO 735 IJ  = 1,IEJE
                  UFELD(IJ) = EVAPW (IJ)
 735              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'EVAPI    ') THEN
                  DO 740 IJ  = 1,IEJE
                  UFELD(IJ) = EVAPI (IJ)
 740              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'AHFSL    ') THEN
                  DO 760 IJ  = 1,IEJE
                  UFELD(IJ) = AHFSL (IJ)
 760              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'AHFSW    ') THEN
                  DO 765 IJ  = 1,IEJE
                  UFELD(IJ) = AHFSW (IJ)
 765              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'AHFSI    ') THEN
                  DO 770 IJ  = 1,IEJE
                  UFELD(IJ) = AHFSI (IJ)
 770              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'AZ0L     ') THEN
                  DO 775 IJ  = 1,IEJE
                  UFELD(IJ) = AZ0L  (IJ)
 775              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'AZ0W     ') THEN
                  DO 780 IJ  = 1,IEJE
                  UFELD(IJ) = AZ0W  (IJ)
 780              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'AZ0I     ') THEN
                  DO 785 IJ  = 1,IEJE
                  UFELD(IJ) = AZ0I  (IJ)
 785              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'ALSOL    ') THEN
                  DO 790 IJ  = 1,IEJE
                  UFELD(IJ) = ALSOL (IJ)
 790              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'ALSOW    ') THEN
                  DO 795 IJ  = 1,IEJE
                  UFELD(IJ) = ALSOW (IJ)
 795              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'ALSOI    ') THEN
                  DO 800 IJ  = 1,IEJE
                  UFELD(IJ) = ALSOI (IJ)
 800              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'AHFICE   ') THEN
                  DO 805 IJ  = 1,IEJE
                  UFELD(IJ) = AHFICE(IJ)
 805              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'QRES     ') THEN
                  DO 810 IJ  = 1,IEJE
                  UFELD(IJ) = QRES  (IJ)
 810              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TMCHL    ') THEN
                  DO 830 IJ  = 1,IEJE
                  UFELD(IJ) = TMCHL (IJ)
 830              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TMCHW    ') THEN
                  DO 835 IJ  = 1,IEJE
                  UFELD(IJ) = TMCHW (IJ)
 835              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TMCHI    ') THEN
                  DO 840 IJ  = 1,IEJE
                  UFELD(IJ) = TMCHI (IJ)
 840              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'BFLHSL   ') THEN
                  DO 845 IJ  = 1,IEJE
                  UFELD(IJ) = BFLHSL(IJ)
 845              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'BFLHSW   ') THEN
                  DO 850 IJ  = 1,IEJE
                  UFELD(IJ) = BFLHSW(IJ)
 850              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'BFLHSI   ') THEN
                  DO 855 IJ  = 1,IEJE
                  UFELD(IJ) = BFLHSI(IJ)
 855              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'BFLQDSL  ') THEN
                  DO 860 IJ  = 1,IEJE
                  UFELD(IJ) = BFLQDSL(IJ)
 860              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'BFLQDSW  ') THEN
                  DO 865 IJ  = 1,IEJE
                  UFELD(IJ) = BFLQDSW(IJ)
 865              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'BFLQDSI  ') THEN
                  DO 870 IJ  = 1,IEJE
                  UFELD(IJ) = BFLQDSI(IJ)
 870              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'QDBOXS   ') THEN
                  DO 875 IJ  = 1,IEJE
                  UFELD(IJ) = QDBOXS(IJ)
 875              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'QWBOXS   ') THEN
                  DO 880 IJ  = 1,IEJE
                  UFELD(IJ) = QWBOXS(IJ)
 880              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'EKBOXS   ') THEN
                  DO 885 IJ  = 1,IEJE
                  UFELD(IJ) = EKBOXS(IJ)
 885              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'FHBOXS   ') THEN
                  DO 890 IJ  = 1,IEJE
                  UFELD(IJ) = FHBOXS(IJ)
 890              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'FIBOXS   ') THEN
                  DO 895 IJ  = 1,IEJE
                  UFELD(IJ) = FIBOXS(IJ)
 895              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'GHPBL    ') THEN
                  DO 900 IJ  = 1,IEJE
                  UFELD(IJ) = GHPBL (IJ)
 900              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'BETA     ') THEN
                  DO 905 IJ  = 1,IEJE
                  UFELD(IJ) = BETA  (IJ)
 905              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'WMINLOK  ') THEN
                  DO 910 IJ  = 1,IEJE
                  UFELD(IJ) = WMINLOK(IJ)
 910              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'WMAXLOK  ') THEN
                  DO 915 IJ  = 1,IEJE
                  UFELD(IJ) = WMAXLOK(IJ)
 915              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'VBM10M   ') THEN
                  DO 920 IJ  = 1,IEJE
                  UFELD(IJ) = VBM10M(IJ)
 920              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'CAPE     ') THEN
                  DO 925 IJ  = 1,IEJE
                  UFELD(IJ) = CAPE  (IJ)
 925              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'WS1      ') THEN
                  DO 930 IJ  = 1,IEJE
                  UFELD(IJ) = WS1   (IJ)
 930              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'WS2      ') THEN
                  DO 932 IJ  = 1,IEJE
                  UFELD(IJ) = WS2  (IJ)
 932              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'WS3      ') THEN
                  DO 934 IJ  = 1,IEJE
                  UFELD(IJ) = WS3   (IJ)
 934              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'WS4      ') THEN
                  DO 936 IJ  = 1,IEJE
                  UFELD(IJ) = WS4  (IJ)
 936              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'WS5      ') THEN
                  DO 938 IJ  = 1,IEJE
                  UFELD(IJ) = WS5   (IJ)
 938              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'DZR      ') THEN
                  DO 940 IJ  = 1,IEJE
                  UFELD(IJ) = DZR  (IJ)
 940              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'DZS      ') THEN
                  DO 942 IJ  = 1,IEJE
                  UFELD(IJ) = DZS  (IJ)
 942              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'FKSAT    ') THEN
                  DO 944 IJ  = 1,IEJE
                  UFELD(IJ) = FKSAT(IJ)
 944              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'FMPOT    ') THEN
                  DO 946 IJ  = 1,IEJE
                  UFELD(IJ) = FMPOT(IJ)
 946              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'BCLAPP   ') THEN
                  DO 948 IJ  = 1,IEJE
                  UFELD(IJ) = BCLAPP(IJ)
 948              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'VPOR     ') THEN
                  DO 950 IJ  = 1,IEJE
                  UFELD(IJ) = VPOR (IJ)
 950              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'ETRANS   ') THEN
                  DO 952 IJ  = 1,IEJE
                  UFELD(IJ) = ETRANS(IJ)
 952              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'EBSOIL   ') THEN
                  DO 954 IJ  = 1,IEJE
                  UFELD(IJ) = EBSOIL(IJ)
 954              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'ESNOW    ') THEN
                  DO 956 IJ  = 1,IEJE
                  UFELD(IJ) = ESNOW(IJ)
 956              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'ESKIN    ') THEN
                  DO 958 IJ  = 1,IEJE
                  UFELD(IJ) = ESKIN(IJ)
 958              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'ERES     ') THEN
                  DO 960 IJ  = 1,IEJE
                  UFELD(IJ) = ERES(IJ)
 960              CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'QIVI     ') THEN
                  DO 962 IJ  = 1,IEJE
                  UFELD(IJ) = QIVI(IJ)
 962              CONTINUE
          ENDIF

          IF (YEMNAME(IVLOC).EQ.'QIBOXS   ') THEN
                  DO 964 IJ  = 1,IEJE
                  UFELD(IJ) = QIBOXS(IJ)
 964              CONTINUE
          ENDIF
C
          ENDIF
C
C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
          CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, 1)
C
C     3. MULTI LEVEL FELDER
      ELSE IF (KAKE(IVLOC,1).NE.KAKE(IVLOC,2)) THEN
C
C     3 ZEITEBENEN
          IF (NZEMTLV(IVLOC).EQ.3) THEN
C
          IF (YEMNAME(IVLOC).EQ.'U        ') THEN
              DO 540 K = KAKE(IVLOC,1),KAKE(IVLOC,2)
              IF(IPDB(8).EQ.109) THEN
                  IPDB(10) = K
              ELSE IF(IPDB(8).EQ.110) THEN
                  IPDB( 9) = K
                  IPDB(10) = K + 1
              ENDIF
C
                  DO 545 IJ  = 1,IEJE
                  UFELD (IJ) = U    (IJ,K,NX)
 545              CONTINUE
C
C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
              CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, K)
C
 540          CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'V        ') THEN
              DO 550 K = KAKE(IVLOC,1),KAKE(IVLOC,2)
              IF(IPDB(8).EQ.109) THEN
                  IPDB(10) = K
              ELSE IF(IPDB(8).EQ.110) THEN
                  IPDB( 9) = K
                  IPDB(10) = K + 1
              ENDIF
C
                  DO 555 IJ  = 1,IEJE
                  UFELD (IJ) = V    (IJ,K,NX)
 555              CONTINUE
C
C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
              CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, K)
C
 550          CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'T        ') THEN
              DO 560 K = KAKE(IVLOC,1),KAKE(IVLOC,2)
              IF(IPDB(8).EQ.109) THEN
                  IPDB(10) = K
              ELSE IF(IPDB(8).EQ.110) THEN
                  IPDB( 9) = K
                  IPDB(10) = K + 1
              ENDIF
C
                  DO 565 IJ  = 1,IEJE
                  UFELD (IJ) = T    (IJ,K,NX)
 565              CONTINUE
C
C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
              CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, K)
C
 560          CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'QD       ') THEN
              DO 570 K = KAKE(IVLOC,1),KAKE(IVLOC,2)
              IF(IPDB(8).EQ.109) THEN
                  IPDB(10) = K
              ELSE IF(IPDB(8).EQ.110) THEN
                  IPDB( 9) = K
                  IPDB(10) = K + 1
              ENDIF
C
                  DO 575 IJ  = 1,IEJE
                  UFELD (IJ) = QD   (IJ,K,NX)
 575              CONTINUE
C
C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
              CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, K)
C
 570          CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'QW       ') THEN
              DO 580 K = KAKE(IVLOC,1),KAKE(IVLOC,2)
              IF(IPDB(8).EQ.109) THEN
                  IPDB(10) = K
              ELSE IF(IPDB(8).EQ.110) THEN
                  IPDB( 9) = K
                  IPDB(10) = K + 1
              ENDIF
C
                  DO 585 IJ  = 1,IEJE
                  UFELD (IJ) = QW   (IJ,K,NX)
 585              CONTINUE
C
C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
              CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, K)
C
 580          CONTINUE
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TKE      ') THEN
              DO 590 K = KAKE(IVLOC,1),KAKE(IVLOC,2)
              IF(IPDB(8).EQ.109) THEN
                  IPDB(10) = K
              ELSE IF(IPDB(8).EQ.110) THEN
                  IPDB( 9) = K
                  IPDB(10) = K + 1
              ENDIF
C
                  DO 595 IJ  = 1,IEJE
                  UFELD (IJ) = TKE  (IJ,K,NX)
 595              CONTINUE
C
C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
              CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, K)
C
 590          CONTINUE
          ENDIF

          IF (YEMNAME(IVLOC).EQ.'PINT     ') THEN
            DO K = KAKE(IVLOC,1),KAKE(IVLOC,2)
              IF(IPDB(8).EQ.109) THEN
                IPDB(10) = K
              ELSE IF(IPDB(8).EQ.110) THEN
                IPDB( 9) = K
                IPDB(10) = K + 1
              ENDIF
C
              DO IJ  = 1,IEJE
                UFELD (IJ) = PINT (IJ,K,NX)
              END DO
C
C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
                CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, K)
C
            END DO
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'W        ') THEN
            DO K = KAKE(IVLOC,1),KAKE(IVLOC,2)
              IF(IPDB(8).EQ.109) THEN
                IPDB(10) = K
              ELSE IF(IPDB(8).EQ.110) THEN
                IPDB( 9) = K
                IPDB(10) = K + 1
              ENDIF
C
              DO IJ  = 1,IEJE
                UFELD (IJ) = W    (IJ,K,NX)
              END DO
C
C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
              CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, K)
C
            END DO
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'DWDT     ') THEN
            DO K = KAKE(IVLOC,1),KAKE(IVLOC,2)
              IF(IPDB(8).EQ.109) THEN
                  IPDB(10) = K
              ELSE IF(IPDB(8).EQ.110) THEN
                  IPDB( 9) = K
                  IPDB(10) = K + 1
              ENDIF
C
              DO IJ  = 1,IEJE
                UFELD (IJ) = DWDT (IJ,K,NX)
              END DO
C
C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
                CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, K)
C
            END DO
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'QI       ') THEN
              DO 1580 K = KAKE(IVLOC,1),KAKE(IVLOC,2)
              IF(IPDB(8).EQ.109) THEN
                  IPDB(10) = K
              ELSE IF(IPDB(8).EQ.110) THEN
                  IPDB( 9) = K
                  IPDB(10) = K + 1
              ENDIF
C
                  DO 1585 IJ  = 1,IEJE
                  UFELD (IJ) = QI   (IJ,K,NX)
 1585              CONTINUE
C
C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
              CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, K)
C
 1580          CONTINUE
          ENDIF
C
C     2 ZEITEBENEN
          ELSE IF (NZEMTLV(IVLOC).EQ.2) THEN
C
          IF (YEMNAME(IVLOC).EQ.'FI       ') THEN
              DO 600 K = KAKE(IVLOC,1),KAKE(IVLOC,2)
              IF(IPDB(8).EQ.109) THEN
                  IPDB(10) = K
              ELSE IF(IPDB(8).EQ.110) THEN
                  IPDB( 9) = K
                  IPDB(10) = K + 1
              ENDIF

                  DO 605 IJ  = 1,IEJE
                  UFELD (IJ) = FI   (IJ,K,NX2)
  605             CONTINUE

C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
              CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, K)

  600         CONTINUE
          ENDIF
C
C     1 ZEITEBENE
          ELSE IF (NZEMTLV(IVLOC).EQ.1) THEN
C
          IF (YEMNAME(IVLOC).EQ.'VERVEL   ') THEN
              DO 610 K = KAKE(IVLOC,1),KAKE(IVLOC,2)
              IF(IPDB(8).EQ.109) THEN
                  IPDB(10) = K
              ELSE IF(IPDB(8).EQ.110) THEN
                  IPDB( 9) = K
                  IPDB(10) = K + 1
              ENDIF
C
                  DO 615 IJ  = 1,IEJE
                  UFELD (IJ) = VERVEL(IJ,K)
  615             CONTINUE
C
C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
              CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, K)

  610         CONTINUE
C
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'FTKVM    ') THEN
              DO 620 K = KAKE(IVLOC,1),KAKE(IVLOC,2)
              IF(IPDB(8).EQ.109) THEN
                  IPDB(10) = K
              ELSE IF(IPDB(8).EQ.110) THEN
                  IPDB( 9) = K
                  IPDB(10) = K + 1
              ENDIF
C
                  DO 625 IJ  = 1,IEJE
                  UFELD (IJ) = FTKVM(IJ,K)
  625             CONTINUE
C
C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
              CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, K)

  620         CONTINUE
C
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'FTKVH    ') THEN
              DO 630 K = KAKE(IVLOC,1),KAKE(IVLOC,2)
              IF(IPDB(8).EQ.109) THEN
                  IPDB(10) = K
              ELSE IF(IPDB(8).EQ.110) THEN
                  IPDB( 9) = K
                  IPDB(10) = K + 1
              ENDIF
C
                  DO 635 IJ  = 1,IEJE
                  UFELD (IJ) = FTKVH(IJ,K)
  635             CONTINUE
C
C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
              CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, K)

  630         CONTINUE
C
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'ACLC     ') THEN
              DO 640 K = KAKE(IVLOC,1),KAKE(IVLOC,2)
              IF(IPDB(8).EQ.109) THEN
                  IPDB(10) = K
              ELSE IF(IPDB(8).EQ.110) THEN
                  IPDB( 9) = K
                  IPDB(10) = K + 1
              ENDIF
C
                  DO 645 IJ  = 1,IEJE
                  UFELD (IJ) = ACLC (IJ,K)
  645             CONTINUE
C
C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
              CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, K)

  640         CONTINUE
C
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'ACLCAC   ') THEN
              DO 650 K = KAKE(IVLOC,1),KAKE(IVLOC,2)
              IF(IPDB(8).EQ.109) THEN
                  IPDB(10) = K
              ELSE IF(IPDB(8).EQ.110) THEN
                  IPDB( 9) = K
                  IPDB(10) = K + 1
              ENDIF
C
                  DO 655 IJ  = 1,IEJE
                  UFELD (IJ) = ACLCAC(IJ,K)
  655             CONTINUE
C
C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
              CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, K)

  650         CONTINUE
C
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'EMTER    ') THEN
              DO 660 K = KAKE(IVLOC,1),KAKE(IVLOC,2)
              IF(IPDB(8).EQ.109) THEN
                  IPDB(10) = K
              ELSE IF(IPDB(8).EQ.110) THEN
                  IPDB( 9) = K
                  IPDB(10) = K + 1
              ENDIF
C
                  DO 665 IJ  = 1,IEJE
                  UFELD (IJ) = EMTER(IJ,K)
  665             CONTINUE
C
C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
              CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, K)

  660         CONTINUE
C
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TRSOL    ') THEN
              DO 670 K = KAKE(IVLOC,1),KAKE(IVLOC,2)
              IF(IPDB(8).EQ.109) THEN
                  IPDB(10) = K
              ELSE IF(IPDB(8).EQ.110) THEN
                  IPDB( 9) = K
                  IPDB(10) = K + 1
              ENDIF
C
                  DO 675 IJ  = 1,IEJE
                  UFELD (IJ) = TRSOL(IJ,K)
  675             CONTINUE
C
C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
              CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, K)

  670         CONTINUE
C
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'EMTEF    ') THEN
              DO 680 K = KAKE(IVLOC,1),KAKE(IVLOC,2)
              IF(IPDB(8).EQ.109) THEN
                  IPDB(10) = K
              ELSE IF(IPDB(8).EQ.110) THEN
                  IPDB( 9) = K
                  IPDB(10) = K + 1
              ENDIF
C
                  DO 685 IJ  = 1,IEJE
                  UFELD (IJ) = EMTEF(IJ,K)
  685             CONTINUE
C
C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
              CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, K)

  680         CONTINUE
C
          ENDIF
C
          IF (YEMNAME(IVLOC).EQ.'TRSOF    ') THEN
              DO 690 K = KAKE(IVLOC,1),KAKE(IVLOC,2)
              IF(IPDB(8).EQ.109) THEN
                  IPDB(10) = K
              ELSE IF(IPDB(8).EQ.110) THEN
                  IPDB( 9) = K
                  IPDB(10) = K + 1
              ENDIF
C
                  DO 695 IJ  = 1,IEJE
                  UFELD (IJ) = TRSOF(IJ,K)
  695             CONTINUE
C
C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
              CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, K)

  690         CONTINUE
C
          ENDIF
CCTBEKS
          IF (YEMNAME(IVLOC).EQ.'RPRAC   ') THEN
              DO 1000 K = KAKE(IVLOC,1),KAKE(IVLOC,2)
              IF(IPDB(8).EQ.109) THEN
                  IPDB(10) = K
              ELSE IF(IPDB(8).EQ.110) THEN
                  IPDB( 9) = K
                  IPDB(10) = K + 1
              ENDIF
C
                  DO 1005 IJ  = 1,IEJE
                  UFELD (IJ) = RPRAC(IJ,K)
 1005             CONTINUE
C
C     FELD IN GRIB-CODE-PACKEN UND AUSSCHREIBEN AUF DATEI 'NUNIT'
              CALL WRITEGB(NUNIT, YTYP, IVLOC, UFELD, AK, BK, K)
 1000         CONTINUE
C
          ENDIF
C
          ENDIF
      ENDIF
C
      RETURN
      END
