      SUBROUTINE PUTECA
     &        (YTYP  , UFELD  , IEJE  , KO    , KU    , KE, KE1,
     &        IVLOC  , YVARN  , LVARDA, NZV   , NZFELD, NZT    ,
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
C**** PUTECA   -   UP:BESETZEN DER FELDER IM LANGZEIT-SPEICHER
C**   AUFRUF   :   CALL PUTECA
C**
C**   ENTRIES  :   KEINE
C**   ZWECK    :   BESETZEN DER FELDER IM LANGZEIT-SPEICHER
C**   VERSIONS-
C**   DATUM    :   05.10.04
C**
C**   EXTERNALS:   KEINE
C**
C**   EINGABE-
C**   PARAMETER:   YTYP:   TYP DER DATEN:
C**                        A1: 1.ANFANGSDATENSATZ A2: 2.ANFANGSDATENSATZ
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
      INCLUDE "corg.h"
      INCLUDE "emgbch.h"
      INCLUDE "emgbri.h"
      INCLUDE "grib.h"
      INCLUDE "higkon.h"
C     ------------------------------------------------------------------
C     Dummy Arguments
C
      INTEGER,   INTENT(IN)    :: IEJE,KO,KU,KE,KE1,NZT,NZV,IVLOC
      INTEGER,   INTENT(INOUT) :: NZFELD
      REAL,      INTENT(INOUT) :: UFELD(IEJE)
      LOGICAL,   INTENT(INOUT) :: LVARDA(NZV)
      CHARACTER, INTENT(IN)    :: YVARN (NZV)*8, YTYP*(*)
C
      REAL, INTENT(INOUT) ::
     & U      (IEJE,KE ,3), V      (IEJE,KE ,3), T      (IEJE,KE ,3),
     & QD     (IEJE,KE ,3), QW     (IEJE,KE ,3), PS     (IEJE,    3),
     & SN     (IEJE,    3), WSECH  (IEJE,    3), TSECH  (IEJE,    3),
     & TSLECH (IEJE,    3), TSWECH (IEJE,    3), TSIECH (IEJE,    3),
     & TD     (IEJE,    3), TDCL   (IEJE,    3), WL     (IEJE,    3),
     & TD3    (IEJE,    3), TD4    (IEJE,    3), TD5    (IEJE,    3),
     & TSN    (IEJE,    3), QDB    (IEJE,    3), TKE    (IEJE,KE ,3),
     & QDBL   (IEJE,    3), QDBW   (IEJE,    3), QDBI   (IEJE,    3),
     & USTRL  (IEJE      ), USTRW  (IEJE      ), USTRI  (IEJE      ),
     & VSTRL  (IEJE      ), VSTRW  (IEJE      ), VSTRI  (IEJE      ),
     & EVAPL  (IEJE      ), EVAPW  (IEJE      ), EVAPI  (IEJE      ),
     & AHFSL  (IEJE      ), AHFSW  (IEJE      ), AHFSI  (IEJE      ),
     & AZ0L   (IEJE      ), AZ0W   (IEJE      ), AZ0I   (IEJE      ),
     & ALSOL  (IEJE      ), ALSOW  (IEJE      ), ALSOI  (IEJE      ),
     & AHFICE (IEJE      ), QRES   (IEJE      ), SRFL   (IEJE      ),
     & TMCHL  (IEJE      ), TMCHW  (IEJE      ), TMCHI  (IEJE      ),
     & BFLHSL (IEJE      ), BFLHSW (IEJE      ), BFLHSI (IEJE      ),
     & BFLQDSL(IEJE      ), BFLQDSW(IEJE      ), BFLQDSI(IEJE      ),
     & BLA    (IEJE      ), AZ0    (IEJE      ), FIB    (IEJE      ),
     & GLAC   (IEJE      ), VGRAT  (IEJE      ), FOREST (IEJE      ),
     & ALBECH (IEJE      ), WSMX   (IEJE      ), VLT    (IEJE      ),
     & FAO    (IEJE      ), BFLHS  (IEJE      ), BFLQDS (IEJE      ),
     & RUNOFF (IEJE      ), TMCM   (IEJE      ), TMCH   (IEJE      ),
     & VERVEL (IEJE,KE   ), EMTEF  (IEJE,2    ), TRSOF  (IEJE,2    ),
     & APRL   (IEJE      ), APRC   (IEJE      ), APRS   (IEJE      ),
     & ACLC   (IEJE,KE   ), ACLCV  (IEJE      ), ALBEDO (IEJE      ),
     & EVAP   (IEJE      ), VAROR  (IEJE      ), T2MAX  (IEJE      ),
     & T2MIN  (IEJE      ), TSMAX  (IEJE      ), TSMIN  (IEJE      ),
     & WIMAX  (IEJE      ), SEAICE (IEJE      ), SICED  (IEJE      ),
     & TOPMAX (IEJE      ), RGCGN  (IEJE      ),
     & EVAPM  (IEJE      ), DSNAC  (IEJE      ), DRAIN  (IEJE      ),
     & RLA    (IEJE      ), PHI    (IEJE      ), AHFS   (IEJE      ),
     & AHFL   (IEJE      ), VDIS   (IEJE      ), USTAR3 (IEJE      ),
     & ACLCOV (IEJE      ), U10    (IEJE      ), V10    (IEJE      ),
     & TEMP2  (IEJE      ), DEW2   (IEJE      ), TSURF  (IEJE      ),
     & WIND10 (IEJE      ), SRADS  (IEJE      ), TRADS  (IEJE      ),
     & USTR   (IEJE      ), SRAD0  (IEJE      ), TRAD0  (IEJE      ),
     & VSTR   (IEJE      ), SRAFS  (IEJE      ), TRAFS  (IEJE      ),
     & USTRGW (IEJE      ), SRAF0  (IEJE      ), TRAF0  (IEJE      ),
     & VSTRGW (IEJE      ), SCLFS  (IEJE      ), TCLFS  (IEJE      ),
     & VDISGW (IEJE      ), SCLF0  (IEJE      ), TCLF0  (IEJE      ),
     & SRAD0U (IEJE      ), SRADSU (IEJE      ), TRADSU (IEJE      ),
     & SNMEL  (IEJE      ), TSLIN  (IEJE      ), ACLCAC (IEJE,KE   ),
     & QVI    (IEJE      ), ALWCVI (IEJE      ), FI     (IEJE,KE ,2),
     & TEFF   (IEJE      ), FTKVM  (IEJE,KE   ), FTKVH  (IEJE,KE   ),
     & EMTER  (IEJE,KE1  ), TRSOL  (IEJE,KE1  ), QDBOXS (IEJE      ),
     & QWBOXS (IEJE      ), EKBOXS (IEJE      ), FHBOXS (IEJE      ),
     & FIBOXS (IEJE      ),TLAMBDA (IEJE      ),DLAMBDA (IEJE      ),
     & PORVOL (IEJE      ), FCAP   (IEJE      ), WI3    (IEJE,    3),
     & WI4    (IEJE,    3), WI5    (IEJE,    3), WI     (IEJE,    3),
     & WICL   (IEJE,    3), BETA   (IEJE      ), WMINLOK(IEJE      ),
     & WMAXLOK(IEJE      ), VBM10M (IEJE      ), CAPE   (IEJE      )
      REAL, INTENT(INOUT) ::
     & WS1   (IEJE       ), WS2   (IEJE       ), WS3   (IEJE       ),
     & WS4   (IEJE       ), WS5   (IEJE       ), DZR   (IEJE       ),
     & DZS   (IEJE       ), FKSAT (IEJE       ), FMPOT (IEJE       ),
     & BCLAPP(IEJE       ), VPOR  (IEJE       ), ETRANS(IEJE       ),
     & EBSOIL(IEJE       ), ESNOW (IEJE       ), ESKIN (IEJE       ),
     & ERES  (IEJE       ), QI    (IEJE, KE ,3), QIVI  (IEJE       ),
     & QIBOXS(IEJE       ), RPRAC (IEJE, KE   ), PINT  (IEJE,KE1, 3),
     & DWDT  (IEJE, KE ,3), W     (IEJE, KE1,3)
C
C     ------------------------------------------------------------------
C     Local Declarations
C
      LOGICAL :: LLANF
      REAL    :: ZFAK,ZBIAS
      INTEGER :: NTLEV,NTAB,K,IJ
C     ------------------------------------------------------------------
C     ANFANGSFELDER ABSPEICHERN
          READ(YTYP(2:2),'(I1)')  NTLEV
          IF(NZEMTLV(IVLOC).EQ.1) NTLEV = 1

C     ENTSPRECHENDES FELD IM EM-LANGZEITSPEICHER SUCHEN
          DO 10 NTAB = 1,NZV
          IF(YEMNAME(IVLOC).EQ.YVARN(NTAB)) THEN

C     BEI ZEITLICH INTEGRIERTEN GROESSEN (TIME RANGE INDIKATOR = 3)
C     DEN SKALIERUNGSFAKTOR SO AENDERN, DASS MIT VV(H) MULTIPLIZIERT
C     WIRD
              IF(NGBTRI(IVLOC).EQ.3) THEN
                  IF(IPDB(16).EQ.0) THEN
                      IF(IPDB(18).EQ.0) ZFAK = DTDEH/EMGBFK(IVLOC)
                      IF(IPDB(18).GT.0) ZFAK = FLOAT(IPDB(18))/
     &                                        (EMGBFK(IVLOC)*60.0)
                  ELSE IF(IPDB(16).EQ.1) THEN
                      IF(IPDB(18).EQ.0) ZFAK = DTDEH/EMGBFK(IVLOC)
                      IF(IPDB(18).GT.0) ZFAK = FLOAT(IPDB(18))/
     &                                         EMGBFK(IVLOC)
                  ENDIF
              ELSE
                  ZFAK = 1.0/EMGBFK(IVLOC)
              ENDIF

C     FELDER, VON DENEN NUR EIN NIVEAU VORLIEGT (SINGLE LEVEL FIELDS)
              IF(KO.EQ.KU) THEN
                  IF(IPDB(8).NE.  1 .AND. IPDB(8).NE.105 .AND.
     &               IPDB(8).NE.109) RETURN

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
              LLANF=NTLEV.EQ.1.AND.NZT.EQ.0.AND.NZEMTLV(IVLOC).GT.1
              ZBIAS=-EMGBBS(IVLOC)*EMGBFK(IVLOC)
C
C     BEI FORTSETZUNGSLAUF KEINE SKALIERUNG DER FELDER
C
              IF (NZT.EQ.0) THEN
                 DO 20 IJ=1,IEJE
                 UFELD(IJ)=(UFELD(IJ)+ZBIAS)*ZFAK
 20              CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'FI      ') THEN
                  IF (LLANF) THEN
                     DO 90 IJ=1,IEJE
                     FI(IJ,K,NTLEV)=UFELD(IJ)
                     FI(IJ,K,2    )=FI(IJ,K,NTLEV)
 90                  CONTINUE
                  ELSE
                     DO 95 IJ=1,IEJE
                     FI(IJ,K,NTLEV)=UFELD(IJ)
 95                  CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'U       ') THEN
                  IF (LLANF) THEN
                     DO 100 IJ=1,IEJE
                     U(IJ,K,NTLEV)=UFELD(IJ)
                     U(IJ,K,2    )=U(IJ,K,NTLEV)
 100                 CONTINUE
                  ELSE
                     DO 105 IJ=1,IEJE
                     U(IJ,K,NTLEV)=UFELD(IJ)
 105                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'V       ') THEN
                  IF (LLANF) THEN
                     DO 110 IJ=1,IEJE
                     V(IJ,K,NTLEV)=UFELD(IJ)
                     V(IJ,K,2    )=V(IJ,K,NTLEV)
 110                 CONTINUE
                  ELSE
                     DO 115 IJ=1,IEJE
                     V(IJ,K,NTLEV)=UFELD(IJ)
 115                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'T       ') THEN
                  IF (LLANF) THEN
                     DO 120 IJ=1,IEJE
                     T(IJ,K,NTLEV)=UFELD(IJ)
                     T(IJ,K,2    )=T(IJ,K,NTLEV)
 120                 CONTINUE
                  ELSE
                     DO 125 IJ=1,IEJE
                     T(IJ,K,NTLEV)=UFELD(IJ)
 125                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'QD      ') THEN
                  IF (LLANF) THEN
                     DO 130 IJ=1,IEJE
                     QD(IJ,K,NTLEV)=UFELD(IJ)
                     QD(IJ,K,2    )=QD(IJ,K,NTLEV)
 130                 CONTINUE
                  ELSE
                     DO 135 IJ=1,IEJE
                     QD(IJ,K,NTLEV)=UFELD(IJ)
 135                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'QW      ') THEN
                  IF (LLANF) THEN
                     DO 140 IJ=1,IEJE
                     QW(IJ,K,NTLEV)=UFELD(IJ)
                     QW(IJ,K,2    )=QW(IJ,K,NTLEV)
 140                 CONTINUE
                  ELSE
                     DO 145 IJ=1,IEJE
                     QW(IJ,K,NTLEV)=UFELD(IJ)
 145                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'PS      ') THEN
                  IF (LLANF) THEN
                     DO 150 IJ=1,IEJE
                     PS(IJ,NTLEV)=UFELD(IJ)
                     PS(IJ,2    )=PS(IJ,NTLEV)
 150                 CONTINUE
                  ELSE
                     DO 155 IJ=1,IEJE
                     PS(IJ,NTLEV)=UFELD(IJ)
 155                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'SN      ') THEN
                  IF (LLANF) THEN
                     DO 160 IJ=1,IEJE
                     SN(IJ,NTLEV)=UFELD(IJ)
                     SN(IJ,2    )=SN(IJ,NTLEV)
 160                 CONTINUE
                  ELSE
                     DO 165 IJ=1,IEJE
                     SN(IJ,NTLEV)=UFELD(IJ)
 165                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'WI3     ') THEN
                  IF (LLANF) THEN
                     DO 161 IJ=1,IEJE
                     WI3(IJ,NTLEV)=UFELD(IJ)
                     WI3(IJ,2    )=WI3(IJ,NTLEV)
 161                 CONTINUE
                  ELSE
                     DO 162 IJ=1,IEJE
                     WI3(IJ,NTLEV)=UFELD(IJ)
 162                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'WI4     ') THEN
                  IF (LLANF) THEN
                     DO 163 IJ=1,IEJE
                     WI4(IJ,NTLEV)=UFELD(IJ)
                     WI4(IJ,2    )=WI4(IJ,NTLEV)
 163                 CONTINUE
                  ELSE
                     DO 164 IJ=1,IEJE
                     WI4(IJ,NTLEV)=UFELD(IJ)
 164                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'WI5     ') THEN
                  IF (LLANF) THEN
                     DO 166 IJ=1,IEJE
                     WI5(IJ,NTLEV)=UFELD(IJ)
                     WI5(IJ,2    )=WI5(IJ,NTLEV)
 166                 CONTINUE
                  ELSE
                     DO 167 IJ=1,IEJE
                     WI5(IJ,NTLEV)=UFELD(IJ)
 167                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'WI      ') THEN
                  IF (LLANF) THEN
                     DO 168 IJ=1,IEJE
                     WI(IJ,NTLEV)=UFELD(IJ)
                     WI(IJ,2    )=WI(IJ,NTLEV)
 168                 CONTINUE
                  ELSE
                     DO 169 IJ=1,IEJE
                     WI(IJ,NTLEV)=UFELD(IJ)
 169                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'WICL    ') THEN
                  IF (LLANF) THEN
                     DO 171 IJ=1,IEJE
                     WICL(IJ,NTLEV)=UFELD(IJ)
                     WICL(IJ,2    )=WICL(IJ,NTLEV)
 171                 CONTINUE
                  ELSE
                     DO 172 IJ=1,IEJE
                     WICL(IJ,NTLEV)=UFELD(IJ)
 172                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'WSECH   ') THEN
                  IF (LLANF) THEN
                     DO 170 IJ=1,IEJE
                     WSECH(IJ,NTLEV)=UFELD(IJ)
                     WSECH(IJ,2    )=WSECH(IJ,NTLEV)
 170                 CONTINUE
                  ELSE
                     DO 175 IJ=1,IEJE
                     WSECH(IJ,NTLEV)=UFELD(IJ)
 175                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TSECH   ') THEN
                  IF (LLANF) THEN
                     DO 180 IJ=1,IEJE
                     TSECH(IJ,NTLEV)=UFELD(IJ)
                     TSECH(IJ,2    )=TSECH(IJ,NTLEV)
 180              CONTINUE
                  ELSE
                     DO 185 IJ=1,IEJE
                     TSECH(IJ,NTLEV)=UFELD(IJ)
 185                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TSLECH  ') THEN
                  IF (LLANF) THEN
                     DO 183 IJ=1,IEJE
                     TSLECH(IJ,NTLEV)=UFELD(IJ)
                     TSLECH(IJ,2    )=TSLECH(IJ,NTLEV)
 183                 CONTINUE
                  ELSE
                     DO 184 IJ=1,IEJE
                     TSLECH(IJ,NTLEV)=UFELD(IJ)
 184                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TSWECH  ') THEN
                  IF (LLANF) THEN
                     DO 186 IJ=1,IEJE
                     TSWECH(IJ,NTLEV)=UFELD(IJ)
                     TSWECH(IJ,2    )=TSWECH(IJ,NTLEV)
 186                 CONTINUE
                  ELSE
                     DO 187 IJ=1,IEJE
                     TSWECH(IJ,NTLEV)=UFELD(IJ)
 187                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TSIECH  ') THEN
                  IF (LLANF) THEN
                     DO 188 IJ=1,IEJE
                     TSIECH(IJ,NTLEV)=UFELD(IJ)
                     TSIECH(IJ,2    )=TSIECH(IJ,NTLEV)
 188                 CONTINUE
                  ELSE
                     DO 189 IJ=1,IEJE
                     TSIECH(IJ,NTLEV)=UFELD(IJ)
 189                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TD      ') THEN
                  IF (LLANF) THEN
                     DO 190 IJ=1,IEJE
                     TD(IJ,NTLEV)=UFELD(IJ)
                     TD(IJ,2    )=TD(IJ,NTLEV)
 190                 CONTINUE
                  ELSE
                     DO 195 IJ=1,IEJE
                     TD(IJ,NTLEV)=UFELD(IJ)
 195                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TDCL    ') THEN
                  IF (LLANF) THEN
                     DO 200 IJ=1,IEJE
                     TDCL(IJ,NTLEV)=UFELD(IJ)
                     TDCL(IJ,2    )=TDCL(IJ,NTLEV)
 200                 CONTINUE
                  ELSE
                     DO 205 IJ=1,IEJE
                     TDCL(IJ,NTLEV)=UFELD(IJ)
 205                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'WL      ') THEN
                  IF (LLANF) THEN
                     DO 210 IJ=1,IEJE
                     WL(IJ,NTLEV)=UFELD(IJ)
                     WL(IJ,2    )=WL(IJ,NTLEV)
 210                 CONTINUE
                  ELSE
                     DO 215 IJ=1,IEJE
                     WL(IJ,NTLEV)=UFELD(IJ)
 215                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TD3     ') THEN
                  IF (LLANF) THEN
                     DO 220 IJ=1,IEJE
                     TD3(IJ,NTLEV)=UFELD(IJ)
                     TD3(IJ,2    )=TD3(IJ,NTLEV)
 220                 CONTINUE
                  ELSE
                     DO 225 IJ=1,IEJE
                     TD3(IJ,NTLEV)=UFELD(IJ)
 225                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TD4     ') THEN
                  IF (LLANF) THEN
                     DO 230 IJ=1,IEJE
                     TD4(IJ,NTLEV)=UFELD(IJ)
                     TD4(IJ,2    )=TD4(IJ,NTLEV)
 230                 CONTINUE
                  ELSE
                     DO 235 IJ=1,IEJE
                     TD4(IJ,NTLEV)=UFELD(IJ)
 235                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TD5     ') THEN
                  IF (LLANF) THEN
                     DO 240 IJ=1,IEJE
                     TD5(IJ,NTLEV)=UFELD(IJ)
                     TD5(IJ,2    )=TD5(IJ,NTLEV)
 240                 CONTINUE
                  ELSE
                     DO 245 IJ=1,IEJE
                     TD5(IJ,NTLEV)=UFELD(IJ)
 245                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TSN     ') THEN
                  IF (LLANF) THEN
                     DO 250 IJ=1,IEJE
                     TSN(IJ,NTLEV)=UFELD(IJ)
                     TSN(IJ,2    )=TSN(IJ,NTLEV)
 250                 CONTINUE
                  ELSE
                     DO 255 IJ=1,IEJE
                     TSN(IJ,NTLEV)=UFELD(IJ)
 255                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'QDB     ') THEN
                  IF (LLANF) THEN
                     DO 260 IJ=1,IEJE
                     QDB(IJ,NTLEV)=UFELD(IJ)
                     QDB(IJ,2    )=QDB(IJ,NTLEV)
 260                 CONTINUE
                  ELSE
                     DO 265 IJ=1,IEJE
                     QDB(IJ,NTLEV)=UFELD(IJ)
 265                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'QDBL    ') THEN
                  IF (LLANF) THEN
                     DO 261 IJ=1,IEJE
                     QDBL(IJ,NTLEV)=UFELD(IJ)
                     QDBL(IJ,2    )=QDBL(IJ,NTLEV)
 261                 CONTINUE
                  ELSE
                     DO 262 IJ=1,IEJE
                     QDBL(IJ,NTLEV)=UFELD(IJ)
 262                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'QDBW    ') THEN
                  IF (LLANF) THEN
                     DO 263 IJ=1,IEJE
                     QDBW(IJ,NTLEV)=UFELD(IJ)
                     QDBW(IJ,2    )=QDBW(IJ,NTLEV)
 263                 CONTINUE
                  ELSE
                     DO 264 IJ=1,IEJE
                     QDBW(IJ,NTLEV)=UFELD(IJ)
 264                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'QDBI    ') THEN
                  IF (LLANF) THEN
                     DO 266 IJ=1,IEJE
                     QDBI(IJ,NTLEV)=UFELD(IJ)
                     QDBI(IJ,2    )=QDBI(IJ,NTLEV)
 266                 CONTINUE
                  ELSE
                     DO 267 IJ=1,IEJE
                     QDBI(IJ,NTLEV)=UFELD(IJ)
 267                 CONTINUE
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TKE     ') THEN
                  IF (LLANF) THEN
                     DO 270 IJ=1,IEJE
                     TKE(IJ,K,NTLEV)=UFELD(IJ)
                     TKE(IJ,K,2    )=TKE(IJ,K,NTLEV)
 270                 CONTINUE
                  ELSE
                     DO 275 IJ=1,IEJE
                     TKE(IJ,K,NTLEV)=UFELD(IJ)
 275                 CONTINUE
                  ENDIF
              ENDIF
CSP
              IF (YEMNAME(IVLOC).EQ.'QI      ') THEN
                  IF (LLANF) THEN
                     DO 1140 IJ=1,IEJE
                     QI(IJ,K,NTLEV)=UFELD(IJ)
                     QI(IJ,K,2    )=QI(IJ,K,NTLEV)
 1140                 CONTINUE
                  ELSE
                     DO 1145 IJ=1,IEJE
                     QI(IJ,K,NTLEV)=UFELD(IJ)
 1145                 CONTINUE
                  ENDIF
              ENDIF
CSP
              IF (YEMNAME(IVLOC).EQ.'VERVEL  ') THEN
                 DO 280 IJ=1,IEJE
                 VERVEL(IJ,K)=UFELD(IJ)
 280             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'ACLC    ') THEN
                 DO 285 IJ=1,IEJE
                 ACLC(IJ,K)=UFELD(IJ)
 285             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'EMTEF   ') THEN
                 DO 290 IJ=1,IEJE
                 EMTEF(IJ,K)=UFELD(IJ)
 290             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TRSOF   ') THEN
                 DO 295 IJ=1,IEJE
                 TRSOF(IJ,K)=UFELD(IJ)
 295             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'BLA     ') THEN
                 DO 300 IJ=1,IEJE
                 BLA(IJ)=UFELD(IJ)
 300             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'AZ0     ') THEN
                 DO 305 IJ=1,IEJE
                 AZ0(IJ)=UFELD(IJ)
 305             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'FIB     ') THEN
                 DO 310 IJ=1,IEJE
                 FIB(IJ)=UFELD(IJ)
 310             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'GLAC    ') THEN
                 DO 315 IJ=1,IEJE
                 GLAC(IJ)=UFELD(IJ)
 315             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'VGRAT   ') THEN
                 DO 320 IJ=1,IEJE
                 VGRAT(IJ)=UFELD(IJ)
 320             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'FOREST  ') THEN
                 DO 325 IJ=1,IEJE
                 FOREST(IJ)=UFELD(IJ)
 325             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'ALBECH  ') THEN
                 DO 330 IJ=1,IEJE
                 ALBECH(IJ)=UFELD(IJ)
 330             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'WSMX    ') THEN
                 DO 335 IJ=1,IEJE
                 WSMX(IJ)=UFELD(IJ)
 335             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'VLT     ') THEN
                 DO 340 IJ=1,IEJE
                 VLT(IJ)=UFELD(IJ)
 340             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'FAO     ') THEN
                 DO 345 IJ=1,IEJE
                 FAO(IJ)=UFELD(IJ)
 345             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'BFLHS   ') THEN
                 DO 350 IJ=1,IEJE
                 BFLHS(IJ)=UFELD(IJ)
 350             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'BFLQDS  ') THEN
                 DO 355 IJ=1,IEJE
                 BFLQDS(IJ)=UFELD(IJ)
 355             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TMCM    ') THEN
                 DO 360 IJ=1,IEJE
                 TMCM(IJ)=UFELD(IJ)
 360             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'RUNOFF  ') THEN
                 DO 365 IJ=1,IEJE
                 RUNOFF(IJ)=UFELD(IJ)
 365             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'DRAIN   ') THEN
                 DO 367 IJ=1,IEJE
                 DRAIN(IJ)=UFELD(IJ)
 367             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'SRFL    ') THEN
                 DO 368 IJ=1,IEJE
                 SRFL(IJ)=UFELD(IJ)
 368             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TMCH    ') THEN
                 DO 370 IJ=1,IEJE
                 TMCH(IJ)=UFELD(IJ)
 370             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'APRL    ') THEN
                 DO 375 IJ=1,IEJE
                 APRL(IJ)=UFELD(IJ)
 375             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'APRC    ') THEN
                 DO 380 IJ=1,IEJE
                 APRC(IJ)=UFELD(IJ)
 380             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'APRS    ') THEN
                 DO 385 IJ=1,IEJE
                 APRS(IJ)=UFELD(IJ)
 385             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'ACLCV   ') THEN
                 DO 390 IJ=1,IEJE
                 ACLCV(IJ)=UFELD(IJ)
 390             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'ALBEDO  ') THEN
                 DO 395 IJ=1,IEJE
                 ALBEDO(IJ)=UFELD(IJ)
 395             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'EVAP    ') THEN
                 DO 400 IJ=1,IEJE
                 EVAP(IJ)=UFELD(IJ)
 400             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'VAROR   ') THEN
                 DO 405 IJ=1,IEJE
                 VAROR(IJ)=UFELD(IJ)
 405             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'T2MAX   ') THEN
                 DO 410 IJ=1,IEJE
                 T2MAX(IJ)=UFELD(IJ)
 410             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'T2MIN   ') THEN
                 DO 415 IJ=1,IEJE
                 T2MIN(IJ)=UFELD(IJ)
 415             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TSMAX   ') THEN
                 DO 420 IJ=1,IEJE
                 TSMAX(IJ)=UFELD(IJ)
 420             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TSMIN   ') THEN
                 DO 425 IJ=1,IEJE
                 TSMIN(IJ)=UFELD(IJ)
 425             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'WIMAX   ') THEN
                 DO 430 IJ=1,IEJE
                 WIMAX(IJ)=UFELD(IJ)
 430             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'SEAICE  ') THEN
                 DO 435 IJ=1,IEJE
                 SEAICE(IJ)=UFELD(IJ)
 435             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'SICED   ') THEN
                 DO 440 IJ=1,IEJE
                 SICED(IJ)=UFELD(IJ)
 440             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TOPMAX  ') THEN
                 DO 445 IJ=1,IEJE
                 TOPMAX(IJ)=UFELD(IJ)
 445             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'RGCGN   ') THEN
                 DO 455 IJ=1,IEJE
                 RGCGN(IJ)=UFELD(IJ)
 455             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TLAMBDA ') THEN
                 DO 456 IJ=1,IEJE
                 TLAMBDA(IJ)=UFELD(IJ)
 456             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'DLAMBDA ') THEN
                 DO 457 IJ=1,IEJE
                 DLAMBDA(IJ)=UFELD(IJ)
 457             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'PORVOL  ') THEN
                 DO 458 IJ=1,IEJE
                 PORVOL (IJ)=UFELD(IJ)
 458             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'FCAP  ') THEN
                 DO 459 IJ=1,IEJE
                 FCAP   (IJ)=UFELD(IJ)
 459             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'EVAPM   ') THEN
                 DO 460 IJ=1,IEJE
                 EVAPM(IJ)=UFELD(IJ)
 460             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'DSNAC   ') THEN
                 DO 465 IJ=1,IEJE
                 DSNAC(IJ)=UFELD(IJ)
 465             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'RLA     ') THEN
                 DO 475 IJ=1,IEJE
                 RLA(IJ)=UFELD(IJ)
 475             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'PHI     ') THEN
                 DO 480 IJ=1,IEJE
                 PHI(IJ)=UFELD(IJ)
 480             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'VDIS    ') THEN
                 DO 485 IJ=1,IEJE
                 VDIS(IJ)=UFELD(IJ)
 485             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'AHFS    ') THEN
                 DO 490 IJ=1,IEJE
                 AHFS(IJ)=UFELD(IJ)
 490             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'AHFL    ') THEN
                 DO 495 IJ=1,IEJE
                 AHFL(IJ)=UFELD(IJ)
 495             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'USTAR3  ') THEN
                 DO 500 IJ=1,IEJE
                 USTAR3(IJ)=UFELD(IJ)
 500             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'ACLCOV  ') THEN
                 DO 505 IJ=1,IEJE
                 ACLCOV(IJ)=UFELD(IJ)
 505             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'U10     ') THEN
                 DO 510 IJ=1,IEJE
                 U10(IJ)=UFELD(IJ)
 510             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'V10     ') THEN
                 DO 515 IJ=1,IEJE
                 V10(IJ)=UFELD(IJ)
 515             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TEMP2   ') THEN
                 DO 520 IJ=1,IEJE
                 TEMP2(IJ)=UFELD(IJ)
 520             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'DEW2    ') THEN
                 DO 525 IJ=1,IEJE
                 DEW2(IJ)=UFELD(IJ)
 525             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TSURF   ') THEN
                 DO 530 IJ=1,IEJE
                 TSURF(IJ)=UFELD(IJ)
 530             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'WIND10  ') THEN
                 DO 535 IJ=1,IEJE
                 WIND10(IJ)=UFELD(IJ)
 535             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'SRADS   ') THEN
                 DO 540 IJ=1,IEJE
                 SRADS(IJ)=UFELD(IJ)
 540             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TRADS   ') THEN
                 DO 545 IJ=1,IEJE
                 TRADS(IJ)=UFELD(IJ)
 545             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'SRAD0   ') THEN
                 DO 550 IJ=1,IEJE
                 SRAD0(IJ)=UFELD(IJ)
 550             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TRAD0   ') THEN
                 DO 555 IJ=1,IEJE
                 TRAD0(IJ)=UFELD(IJ)
 555             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'USTR    ') THEN
                 DO 560 IJ=1,IEJE
                 USTR(IJ)=UFELD(IJ)
 560             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'VSTR    ') THEN
                 DO 565 IJ=1,IEJE
                 VSTR(IJ)=UFELD(IJ)
 565             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'SRAFS   ') THEN
                 DO 570 IJ=1,IEJE
                 SRAFS(IJ)=UFELD(IJ)
 570             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TRAFS   ') THEN
                 DO 575 IJ=1,IEJE
                 TRAFS(IJ)=UFELD(IJ)
 575             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'SRAF0   ') THEN
                 DO 580 IJ=1,IEJE
                 SRAF0(IJ)=UFELD(IJ)
 580             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TRAF0   ') THEN
                 DO 585 IJ=1,IEJE
                 TRAF0(IJ)=UFELD(IJ)
 585             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'SCLFS   ') THEN
                 DO 590 IJ=1,IEJE
                 SCLFS(IJ)=UFELD(IJ)
 590             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TCLFS   ') THEN
                 DO 595 IJ=1,IEJE
                 TCLFS(IJ)=UFELD(IJ)
 595             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'SCLF0   ') THEN
                 DO 600 IJ=1,IEJE
                 SCLF0(IJ)=UFELD(IJ)
 600             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TCLF0   ') THEN
                 DO 605 IJ=1,IEJE
                 TCLF0(IJ)=UFELD(IJ)
 605             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'USTRGW  ') THEN
                 DO 610 IJ=1,IEJE
                 USTRGW(IJ)=UFELD(IJ)
 610             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'VSTRGW  ') THEN
                 DO 615 IJ=1,IEJE
                 VSTRGW(IJ)=UFELD(IJ)
 615             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'VDISGW  ') THEN
                 DO 620 IJ=1,IEJE
                 VDISGW(IJ)=UFELD(IJ)
 620             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'SRAD0U  ') THEN
                 DO 625 IJ=1,IEJE
                 SRAD0U(IJ)=UFELD(IJ)
 625             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'SRADSU  ') THEN
                 DO 630 IJ=1,IEJE
                 SRADSU(IJ)=UFELD(IJ)
 630             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TRADSU  ') THEN
                 DO 635 IJ=1,IEJE
                 TRADSU(IJ)=UFELD(IJ)
 635             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TEFF    ') THEN
                 DO 640 IJ=1,IEJE
                 TEFF(IJ)=UFELD(IJ)
 640             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'SNMEL   ') THEN
                 DO 645 IJ=1,IEJE
                 SNMEL(IJ)=UFELD(IJ)
 645             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TSLIN   ') THEN
                 DO 650 IJ=1,IEJE
                 TSLIN(IJ)=UFELD(IJ)
 650             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'QVI     ') THEN
                 DO 655 IJ=1,IEJE
                 QVI(IJ)=UFELD(IJ)
 655             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'ALWCVI  ') THEN
                 DO 660 IJ=1,IEJE
                 ALWCVI(IJ)=UFELD(IJ)
 660             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'USTRL    ') THEN
                 DO 700 IJ  = 1,IEJE
                 USTRL (IJ) = UFELD(IJ)
 700             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'USTRW    ') THEN
                 DO 705 IJ  = 1,IEJE
                 USTRW (IJ) = UFELD(IJ)
 705             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'USTRI    ') THEN
                 DO 710 IJ  = 1,IEJE
                 USTRI (IJ) = UFELD(IJ)
 710             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'VSTRL    ') THEN
                 DO 715 IJ  = 1,IEJE
                 VSTRL (IJ) = UFELD(IJ)
 715             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'VSTRW    ') THEN
                 DO 720 IJ  = 1,IEJE
                 VSTRW (IJ) = UFELD(IJ)
 720             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'VSTRI    ') THEN
                 DO 725 IJ  = 1,IEJE
                 VSTRI (IJ) = UFELD(IJ)
 725             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'EVAPL    ') THEN
                 DO 730 IJ  = 1,IEJE
                 EVAPL (IJ) = UFELD(IJ)
 730             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'EVAPW    ') THEN
                 DO 735 IJ  = 1,IEJE
                 EVAPW (IJ) = UFELD(IJ)
 735             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'EVAPI    ') THEN
                 DO 740 IJ  = 1,IEJE
                 EVAPI (IJ) = UFELD(IJ)
 740             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'AHFSL    ') THEN
                 DO 760 IJ  = 1,IEJE
                 AHFSL (IJ) = UFELD(IJ)
 760             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'AHFSW    ') THEN
                 DO 765 IJ  = 1,IEJE
                 AHFSW (IJ) = UFELD(IJ)
 765             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'AHFSI    ') THEN
                 DO 770 IJ  = 1,IEJE
                 AHFSI (IJ) = UFELD(IJ)
 770             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'AZ0L     ') THEN
                 DO 775 IJ  = 1,IEJE
                 AZ0L  (IJ) = UFELD(IJ)
 775             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'AZ0W     ') THEN
                 DO 780 IJ  = 1,IEJE
                 AZ0W  (IJ) = UFELD(IJ)
 780             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'AZ0I     ') THEN
                 DO 785 IJ  = 1,IEJE
                 AZ0I  (IJ) = UFELD(IJ)
 785             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'ALSOL    ') THEN
                 DO 790 IJ  = 1,IEJE
                 ALSOL (IJ) = UFELD(IJ)
 790             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'ALSOW    ') THEN
                 DO 795 IJ  = 1,IEJE
                 ALSOW (IJ) = UFELD(IJ)
 795             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'ALSOI    ') THEN
                 DO 800 IJ  = 1,IEJE
                 ALSOI (IJ) = UFELD(IJ)
 800             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'AHFICE   ') THEN
                 DO 805 IJ  = 1,IEJE
                 AHFICE(IJ) = UFELD(IJ)
 805             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'QRES     ') THEN
                 DO 810 IJ  = 1,IEJE
                 QRES  (IJ) = UFELD(IJ)
 810             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TMCHL    ') THEN
                 DO 830 IJ  = 1,IEJE
                 TMCHL (IJ) = UFELD(IJ)
 830             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TMCHW    ') THEN
                 DO 835 IJ  = 1,IEJE
                 TMCHW (IJ) = UFELD(IJ)
 835             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TMCHI    ') THEN
                 DO 840 IJ  = 1,IEJE
                 TMCHI (IJ) = UFELD(IJ)
 840             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'BFLHSL   ') THEN
                 DO 845 IJ  = 1,IEJE
                 BFLHSL(IJ) = UFELD(IJ)
 845             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'BFLHSW   ') THEN
                 DO 850 IJ  = 1,IEJE
                 BFLHSW(IJ) = UFELD(IJ)
 850             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'BFLHSI   ') THEN
                 DO 855 IJ  = 1,IEJE
                 BFLHSI(IJ) = UFELD(IJ)
 855             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'BFLQDSL  ') THEN
                 DO 860 IJ  = 1,IEJE
                 BFLQDSL(IJ)= UFELD(IJ)
 860             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'BFLQDSW  ') THEN
                 DO 865 IJ  = 1,IEJE
                 BFLQDSW(IJ)= UFELD(IJ)
 865             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'BFLQDSI  ') THEN
                 DO 870 IJ  = 1,IEJE
                 BFLQDSI(IJ)= UFELD(IJ)
 870             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'ACLCAC  ') THEN
                 DO 670 IJ=1,IEJE
                 ACLCAC(IJ,K)=UFELD(IJ)
 670             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'FTKVM   ') THEN
                 DO 675 IJ=1,IEJE
                 FTKVM(IJ,K)=UFELD(IJ)
 675             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'FTKVH   ') THEN
                 DO 680 IJ=1,IEJE
                 FTKVH(IJ,K)=UFELD(IJ)
 680             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'EMTER   ') THEN
                 DO 685 IJ=1,IEJE
                 EMTER(IJ,K)=UFELD(IJ)
 685             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'TRSOL   ') THEN
                 DO 690 IJ=1,IEJE
                 TRSOL(IJ,K)=UFELD(IJ)
 690             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'QDBOXS  ') THEN
                 DO 875 IJ=1,IEJE
                 QDBOXS(IJ)=UFELD(IJ)
 875             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'QWBOXS  ') THEN
                 DO 880 IJ=1,IEJE
                 QWBOXS(IJ)=UFELD(IJ)
 880             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'EKBOXS  ') THEN
                 DO 885 IJ=1,IEJE
                 EKBOXS(IJ)=UFELD(IJ)
 885             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'FHBOXS  ') THEN
                 DO 890 IJ=1,IEJE
                 FHBOXS(IJ)=UFELD(IJ)
 890             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'FIBOXS  ') THEN
                 DO 895 IJ=1,IEJE
                 FIBOXS(IJ)=UFELD(IJ)
 895             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'BETA    ') THEN
                 DO 900 IJ=1,IEJE
                 BETA  (IJ)=UFELD(IJ)
 900             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'WMINLOK ') THEN
                 DO 905 IJ=1,IEJE
                 WMINLOK(IJ)=UFELD(IJ)
 905             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'WMAXLOK ') THEN
                 DO 910 IJ=1,IEJE
                 WMAXLOK(IJ)=UFELD(IJ)
 910             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'VBM10M  ') THEN
                 DO 915 IJ=1,IEJE
                 VBM10M(IJ)=UFELD(IJ)
 915             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'CAPE    ') THEN
                 DO 920 IJ=1,IEJE
                 CAPE  (IJ)=UFELD(IJ)
 920             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'WS1     ') THEN
                 DO 922 IJ=1,IEJE
                 WS1  (IJ)=UFELD(IJ)
 922             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'WS2     ') THEN
                 DO 924 IJ=1,IEJE
                 WS2  (IJ)=UFELD(IJ)
 924             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'WS3     ') THEN
                 DO 926 IJ=1,IEJE
                 WS3  (IJ)=UFELD(IJ)
 926             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'WS4     ') THEN
                 DO 928 IJ=1,IEJE
                 WS4  (IJ)=UFELD(IJ)
 928             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'WS5     ') THEN
                 DO 930 IJ=1,IEJE
                 WS5  (IJ)=UFELD(IJ)
 930             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'DZR     ') THEN
                 DO 932 IJ=1,IEJE
                 DZR  (IJ)=UFELD(IJ)
 932             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'DZS     ') THEN
                 DO 934 IJ=1,IEJE
                 DZS  (IJ)=UFELD(IJ)
 934             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'FKSAT   ') THEN
                 DO 936 IJ=1,IEJE
                 FKSAT(IJ)=UFELD(IJ)
 936             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'FMPOT   ') THEN
                 DO 938 IJ=1,IEJE
                 FMPOT(IJ)=UFELD(IJ)
 938             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'BCLAPP  ') THEN
                 DO 940 IJ=1,IEJE
                 BCLAPP(IJ)=UFELD(IJ)
 940             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'VPOR    ') THEN
                 DO 942 IJ=1,IEJE
                 VPOR (IJ)=UFELD(IJ)
 942             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'ETRANS  ') THEN
                 DO 944 IJ=1,IEJE
                 ETRANS(IJ)=UFELD(IJ)
 944             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'EBSOIL  ') THEN
                 DO 946 IJ=1,IEJE
                 EBSOIL(IJ)=UFELD(IJ)
 946             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'ESNOW   ') THEN
                 DO 948 IJ=1,IEJE
                 ESNOW(IJ)=UFELD(IJ)
 948             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'ESKIN   ') THEN
                 DO 950 IJ=1,IEJE
                 ESKIN(IJ)=UFELD(IJ)
 950             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'ERES    ') THEN
                 DO 952 IJ=1,IEJE
                 ERES (IJ)=UFELD(IJ)
 952             CONTINUE
              ENDIF
CSP
              IF (YEMNAME(IVLOC).EQ.'QIBOXS  ') THEN
                 DO 954 IJ=1,IEJE
                 QIBOXS(IJ)=UFELD(IJ)
 954             CONTINUE
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'QIVI  ') THEN
                 DO 956 IJ=1,IEJE
                 QIVI(IJ)=UFELD(IJ)
 956             CONTINUE
              ENDIF

              IF (YEMNAME(IVLOC).EQ.'W       ') THEN
                  W(1:IEJE,K,NTLEV) = UFELD(1:IEJE)
                  IF (LLANF) THEN
                    W(1:IEJE,K,2) = UFELD(1:IEJE)
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'DWDT    ') THEN
                  DWDT(1:IEJE,K,NTLEV) = UFELD(1:IEJE)
                  IF (LLANF) THEN
                    DWDT(1:IEJE,K,2) = UFELD(1:IEJE)
                  ENDIF
              ENDIF
C
              IF (YEMNAME(IVLOC).EQ.'PINT    ') THEN
                  PINT(1:IEJE,K,NTLEV) = UFELD(1:IEJE)
                  IF (LLANF) THEN
                    PINT(1:IEJE,K,2) = UFELD(1:IEJE)
                  ENDIF
              ENDIF
CSP
CCTBEKS
              IF (YEMNAME(IVLOC).EQ.'RPRAC   ') THEN
                 DO 957 IJ=1,IEJE
                 RPRAC(IJ,K)=UFELD(IJ)
 957             CONTINUE
              ENDIF
C
              RETURN
          ENDIF
 10       CONTINUE

      END SUBROUTINE PUTECA
