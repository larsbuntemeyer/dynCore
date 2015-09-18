      SUBROUTINE ECEPREP(YTYP  ,
     &   AK    , BK    , TMKVMH, U      , V      , T      , QD    ,
     &   QW    , FI    , VERVEL, ACLC   , ACLCAC , TKE    , EMTER ,
     &   TRSOL , EMTEF , TRSOF , PS     , QDB    , TSECH  , WSECH ,
     &   TSLECH, TSWECH, TSIECH, USTRL  , USTRW  , USTRI  , VSTRL ,
     &   VSTRW , VSTRI , EVAPL , EVAPW  , EVAPI  ,
     &   AHFSL , AHFSW , AHFSI , AZ0L   , AZ0W   , AZ0I   ,
     &   ALSOL , ALSOW , ALSOI , AHFICE , QRES   ,
     &   TMCHL , TMCHW , TMCHI , QDBL   , QDBW   , QDBI   ,
     &   BFLHSL, BFLHSW, BFLHSI, BFLQDSL, BFLQDSW, BFLQDSI, SN    ,
     &   TD    , TDCL  , WL    , TSN    , TD3    , TD4    , TD5   ,
     &   SRADS , TRADS , SRAD0 , TRAD0  , APRL   , APRC   , APRS  ,
     &   VDIS  , AHFS  , FIB   , BLA    , AHFL   , USTAR3 , RUNOFF,
     &   ACLCV , ACLCOV, TMCM  , TMCH   , PHI    , RLA    , BFLHS ,
     &   BFLQDS, U10   , V10   , TEMP2  , DEW2   , TSURF  , WIND10,
     &   AZ0   , ALBECH, ALBEDO, USTR   , VSTR   , EVAP   , EVAPM ,
     &   SRAFS , TRAFS , SRAF0 , TRAF0  , SCLFS  , TCLFS  , SCLF0 ,
     &   TCLF0 , USTRGW, VSTRGW, VDISGW , VGRAT  , VAROR  , VLT   ,
     &   T2MAX , SRAD0U, SRADSU, TRADSU , T2MIN  , SEAICE , SICED ,
     &   FOREST, TEFF  , TSMAX , TSMIN  , WIMAX  , TOPMAX , SNMEL ,
     &   TSLIN , DSNAC , FAO   , RGCGN  , WSMX   , QVI    ,
     &   ALWCVI, GLAC  , DRAIN , SRFL   , QDBOXS , QWBOXS , EKBOXS,
     &   FHBOXS, FIBOXS,TLAMBDA, DLAMBDA, PORVOL , FCAP   , WI3   ,
     &   WI4   , WI5   , WI    , WICL   , GHPBL  , BETA   , WMINLOK,
     &   WMAXLOK,VBM10M, CAPE  , WS1    , WS2    , WS3    , WS4   ,
     &   WS5   , DZR   , DZS   , FKSAT  , FMPOT  , BCLAPP , VPOR  ,
     &   ETRANS, EBSOIL, ESNOW , ESKIN  , ERES   , QI     , QIVI  ,
     &   QIBOXS, PINT  , DWDT  , W      , RPRAC)
C
C     SUBROUTINE ECEPREP
C
C**** ECEPREP  -   UP:ERSTELLUNG VON ERGEBNISDATEIEN DES REMO
C**
C**   AUFRUF   :   CALL ECEPREP(YTYP) IN UP *EC4ORG*
C**   ENTRIES  :   KEINE
C**   ZWECK    :   ERSTELLUNG VON ERGEBNISDATEIEN DES EM UND CATALOGI-
C**                SIEREN DER FELDER IN EINER GRIB-CODE-DATEI ODER IN
C**                EINER DATENBANK.
C**   VERSIONS-
C**   DATUM    :   07.03.05
C**
C**   EXTERNALS:   MAKEPN, WRITEC4, SEND, PRIV, PRCV
C**
C**   EINGABE-
C**   PARAMETER:   YTYP : TYP DER ERGEBNISDATEN:
C**                      'E': NORMALE ERGEBNISDATEN  FUER GESAMTGEBIET
C**                      'D':         ERGEBNISDATEN  FUER TEILGEBIET
C**                      'F': FORTSETZUNGSDATEN      FUER GESAMTGEBIET
C**                      'T': TRAJEKTORIEN-DATEN     FUER GESAMTGEBIET
C**   AUSGABE-
C**   PARAMETER:   KEINE
C**
C**   COMMON-
C**   BLOECKE  :   ORG, CORG, COMPHY, UNITCH, UNITNR
C**
C**   METHODE  :   SEQUENTIELLER AUFRUF DER PROGRAMME ZUM ERSTELLEN
C**                UND CATALOGISIEREN DER ERGEBNISDATEN
C**   FEHLERBE-
C**   HANDLUNG :   STOP/ABORT IM FEHLERFALL
C**   VERFASSER:   R.PODZUN
C
C
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
      INCLUDE "org.h"
      INCLUDE "corg.h"
      INCLUDE "unitch.h"
      INCLUDE "unitnr.h"
C
      CHARACTER, INTENT(IN) :: YTYP*(*)
C
      REAL, INTENT(IN) :: AK(KE1), BK(KE1)
C
      REAL, INTENT(IN) :: TMKVMH(IE*(KE-1),JE,2)
C
      REAL, INTENT(IN) ::
     & U     (IEJE,KE,3), V     (IEJE,KE,3), T     (IEJE,KE,3),
     & QD    (IEJE,KE,3), QW    (IEJE,KE,3), FI    (IEJE,KE,2),
     & VERVEL(IEJE,KE  ), ACLC  (IEJE,KE  ), ACLCAC(IEJE,KE  ),
     & TKE   (IEJE,KE,3), EMTER (IEJE,KE1 ), TRSOL (IEJE,KE1 ),
     & EMTEF (IEJE,2   ), TRSOF (IEJE,2   )
C
      REAL, INTENT(IN) ::
     & PS     (IEJE,3), QDB   (IEJE,3), TSECH (IEJE,3), WSECH (IEJE,3),
     & TSLECH (IEJE,3), TSWECH(IEJE,3), TSIECH(IEJE,3), TD5   (IEJE,3),
     & QDBL   (IEJE,3), QDBW  (IEJE,3), QDBI  (IEJE,3), TD4   (IEJE,3),
     & SN     (IEJE,3), TD    (IEJE,3), TDCL  (IEJE,3), WL    (IEJE,3),
     & TSN    (IEJE,3), TD3   (IEJE,3)
C
      REAL, INTENT(IN) ::
     & USTRL  (IEJE  ), USTRW (IEJE  ), USTRI (IEJE  ), VSTRL (IEJE  ),
     & VSTRW  (IEJE  ), VSTRI (IEJE  ), EVAPL (IEJE  ), EVAPW (IEJE  ),
     & EVAPI  (IEJE  ), AHFSL (IEJE  ), AHFSW (IEJE  ), AHFSI (IEJE  ),
     & AZ0L   (IEJE  ), AZ0W  (IEJE  ), AZ0I  (IEJE  ), ALSOL (IEJE  ),
     & ALSOW  (IEJE  ), ALSOI (IEJE  ), AHFICE(IEJE  ), QRES  (IEJE  ),
     & TMCHL  (IEJE  ), TMCHW (IEJE  ), U10   (IEJE  ), V10   (IEJE  ),
     & TMCHI  (IEJE  ), BFLHSL(IEJE  ), BFLHSW(IEJE  ), BFLHSI(IEJE  ),
     & BFLQDSL(IEJE  ), SRADS (IEJE  ), SRFL  (IEJE  ), VDISGW(IEJE  ),
     & BFLQDSW(IEJE  ), TRADS (IEJE  ), SRAD0 (IEJE  ), TRAD0 (IEJE  ),
     & APRL   (IEJE  ), APRC  (IEJE  ), APRS  (IEJE  ), VDIS  (IEJE  ),
     & AHFS   (IEJE  ), FIB   (IEJE  ), BLA   (IEJE  ), AHFL  (IEJE  ),
     & USTAR3 (IEJE  ), RUNOFF(IEJE  ), ACLCV (IEJE  ), ACLCOV(IEJE  ),
     & BFLQDSI(IEJE  ), TMCM  (IEJE  ), TMCH  (IEJE  ), VSTRGW(IEJE  ),
     & PHI    (IEJE  ), RLA   (IEJE  ), BFLHS (IEJE  ), BFLQDS(IEJE  ),
     & TEMP2  (IEJE  ), DEW2  (IEJE  ), TSURF (IEJE  ), WIND10(IEJE  ),
     & AZ0    (IEJE  ), ALBECH(IEJE  ), ALBEDO(IEJE  ), USTR  (IEJE  ),
     & VSTR   (IEJE  ), EVAP  (IEJE  ), EVAPM (IEJE  ), SRAFS (IEJE  ),
     & TRAFS  (IEJE  ), SRAF0 (IEJE  ), TRAF0 (IEJE  ), SCLFS (IEJE  ),
     & TCLFS  (IEJE  ), SCLF0 (IEJE  ), TCLF0 (IEJE  ), USTRGW(IEJE  ),
     & QDBOXS (IEJE  ), QWBOXS(IEJE  ), EKBOXS(IEJE  ), FHBOXS(IEJE  ),
     & FIBOXS (IEJE  ),TLAMBDA(IEJE  ),DLAMBDA(IEJE  ), PORVOL(IEJE  ),
     & FCAP   (IEJE  ), WI3   (IEJE,3), WI4   (IEJE,3), WI5   (IEJE,3),
     & WI     (IEJE,3), WICL  (IEJE,3), GHPBL (IEJE)  , BETA  (IEJE)  ,
     & WMINLOK(IEJE)  , WMAXLOK(IEJE) , VBM10M(IEJE)  , CAPE  (IEJE)
C
      REAL, INTENT(IN) ::
     & VGRAT  (IEJE  ), VAROR (IEJE  ), VLT   (IEJE  ), T2MAX (IEJE  ),
     & SRAD0U (IEJE  ), SRADSU(IEJE  ), TRADSU(IEJE  ), T2MIN (IEJE  ),
     & SEAICE (IEJE  ), SICED (IEJE  ), FOREST(IEJE  ), TEFF  (IEJE  ),
     & TSMAX  (IEJE  ), TSMIN (IEJE  ), WIMAX (IEJE  ), TOPMAX(IEJE  ),
     & SNMEL  (IEJE  ), TSLIN (IEJE  ), DSNAC (IEJE  ), FAO   (IEJE  ),
     & RGCGN  (IEJE  ), WSMX  (IEJE  ), QVI   (IEJE  ),
     & ALWCVI (IEJE  ), GLAC  (IEJE  ), DRAIN (IEJE  )
C
      REAL, INTENT(IN) ::
     & WS1    (IEJE  ), WS2   (IEJE  ), WS3   (IEJE  ), WS4   (IEJE  ),
     & WS5    (IEJE  ), DZR   (IEJE  ), DZS   (IEJE  ), FKSAT (IEJE  ),
     & FMPOT  (IEJE  ), BCLAPP(IEJE  ), VPOR  (IEJE  ), ETRANS(IEJE  ),
     & EBSOIL (IEJE  ), ESNOW (IEJE  ), ESKIN (IEJE  ), ERES  (IEJE  )

      REAL ::   PINT(IEJE,KE1,3), DWDT(IEJE,KE ,3), 
     1          W   (IEJE,KE1,3)
CSP
      REAL, INTENT(IN) ::
     & QI  (IEJE,KE,3), QIVI (IEJE)   , QIBOXS(IEJE  ),
     & RPRAC (IEJE,KE)
C     ANZAHL DER VARIABLEN FUER EINEN FORTSETZUNGSLAUF
      NZVF=NZMXVF
C     DATEINAMEN ERZEUGEN; DATEINAME WIRD IN 'YEDNAM' BZW.
C     'YDDNAM', 'YFDNAM', 'YTDNAM' GESCHRIEBEN
CCM
C      WRITE (*,*) "eceprep: E Datei schreiben:"
C      PRINT *, NZT
C      PRINT *, myid
C      PRINT *, "THIS IS ECEPREC IF STARTS"
      IF (YTYP.NE.'E'.AND.YTYP.NE.'N') CALL MAKEPN(YTYP, NZT + 1)
C
CBE
C      PRINT *, 'end of if makepn call done'
C     SCHREIBEN DER ERGEBNISDATEI
C     ERGEBNISDATEN IN BANK EINBRINGEN ODER IN GRIB-CODE-DATEI
C     CATALOGISIEREN
      IF (YTYP.EQ.'E') THEN
         CALL WRITEC4('E'   , YMVARN, NZVM   , NZT+1  , NE     ,
     &        AK    , BK    , TMKVMH, U      , V      , T      , QD    ,
     &        QW    , FI    , VERVEL, ACLC   , ACLCAC , TKE    , EMTER ,
     &        TRSOL , EMTEF , TRSOF , PS     , QDB    , TSECH  , WSECH ,
     &        TSLECH, TSWECH, TSIECH, USTRL  , USTRW  , USTRI  , VSTRL ,
     &        VSTRW , VSTRI , EVAPL , EVAPW  , EVAPI  ,
     &        AHFSL , AHFSW , AHFSI , AZ0L   , AZ0W   , AZ0I   ,
     &        ALSOL , ALSOW , ALSOI , AHFICE , QRES   ,
     &        TMCHL , TMCHW , TMCHI , QDBL   , QDBW   , QDBI   ,
     &        BFLHSL, BFLHSW, BFLHSI, BFLQDSL, BFLQDSW, BFLQDSI, SN    ,
     &        TD    , TDCL  , WL    , TSN    , TD3    , TD4    , TD5   ,
     &        SRADS , TRADS , SRAD0 , TRAD0  , APRL   , APRC   , APRS  ,
     &        VDIS  , AHFS  , FIB   , BLA    , AHFL   , USTAR3 , RUNOFF,
     &        ACLCV , ACLCOV, TMCM  , TMCH   , PHI    , RLA    , BFLHS ,
     &        BFLQDS, U10   , V10   , TEMP2  , DEW2   , TSURF  , WIND10,
     &        AZ0   , ALBECH, ALBEDO, USTR   , VSTR   , EVAP   , EVAPM ,
     &        SRAFS , TRAFS , SRAF0 , TRAF0  , SCLFS  , TCLFS  , SCLF0 ,
     &        TCLF0 , USTRGW, VSTRGW, VDISGW , VGRAT  , VAROR  , VLT   ,
     &        T2MAX , SRAD0U, SRADSU, TRADSU , T2MIN  , SEAICE , SICED ,
     &        FOREST, TEFF  , TSMAX , TSMIN  , WIMAX  , TOPMAX , SNMEL ,
     &        TSLIN , DSNAC , FAO   , RGCGN  , WSMX   , QVI    ,
     &        ALWCVI, GLAC  , DRAIN , SRFL   , QDBOXS , QWBOXS , EKBOXS,
     &        FHBOXS, FIBOXS,TLAMBDA, DLAMBDA, PORVOL , FCAP   , WI3   ,
     &        WI4   , WI5   , WI    , WICL   , GHPBL  , BETA   ,WMINLOK,
     &        WMAXLOK,VBM10M, CAPE  , WS1    , WS2    , WS3    , WS4   ,
     &        WS5   , DZR   , DZS   , FKSAT  , FMPOT  , BCLAPP , VPOR  ,
     &        ETRANS, EBSOIL, ESNOW , ESKIN  , ERES   , QI     , QIVI  ,
     &        QIBOXS, PINT  , DWDT   , W      , RPRAC)
      ENDIF
      IF (YTYP.EQ.'N') THEN
         CALL WRITEC4('N'   , YNVARN, NZVN   , NZT+1  , NE     ,
     &        AK    , BK    , TMKVMH, U      , V      , T      , QD    ,
     &        QW    , FI    , VERVEL, ACLC   , ACLCAC , TKE    , EMTER ,
     &        TRSOL , EMTEF , TRSOF , PS     , QDB    , TSECH  , WSECH ,
     &        TSLECH, TSWECH, TSIECH, USTRL  , USTRW  , USTRI  , VSTRL ,
     &        VSTRW , VSTRI , EVAPL , EVAPW  , EVAPI  ,
     &        AHFSL , AHFSW , AHFSI , AZ0L   , AZ0W   , AZ0I   ,
     &        ALSOL , ALSOW , ALSOI , AHFICE , QRES   ,
     &        TMCHL , TMCHW , TMCHI , QDBL   , QDBW   , QDBI   ,
     &        BFLHSL, BFLHSW, BFLHSI, BFLQDSL, BFLQDSW, BFLQDSI, SN    ,
     &        TD    , TDCL  , WL    , TSN    , TD3    , TD4    , TD5   ,
     &        SRADS , TRADS , SRAD0 , TRAD0  , APRL   , APRC   , APRS  ,
     &        VDIS  , AHFS  , FIB   , BLA    , AHFL   , USTAR3 , RUNOFF,
     &        ACLCV , ACLCOV, TMCM  , TMCH   , PHI    , RLA    , BFLHS ,
     &        BFLQDS, U10   , V10   , TEMP2  , DEW2   , TSURF  , WIND10,
     &        AZ0   , ALBECH, ALBEDO, USTR   , VSTR   , EVAP   , EVAPM ,
     &        SRAFS , TRAFS , SRAF0 , TRAF0  , SCLFS  , TCLFS  , SCLF0 ,
     &        TCLF0 , USTRGW, VSTRGW, VDISGW , VGRAT  , VAROR  , VLT   ,
     &        T2MAX , SRAD0U, SRADSU, TRADSU , T2MIN  , SEAICE , SICED ,
     &        FOREST, TEFF  , TSMAX , TSMIN  , WIMAX  , TOPMAX , SNMEL ,
     &        TSLIN , DSNAC , FAO   , RGCGN  , WSMX   , QVI    ,
     &        ALWCVI, GLAC  , DRAIN , SRFL   , QDBOXS , QWBOXS , EKBOXS,
     &        FHBOXS, FIBOXS,TLAMBDA, DLAMBDA, PORVOL , FCAP   , WI3   ,
     &        WI4   , WI5   , WI    , WICL   , GHPBL  , BETA   ,WMINLOK,
     &        WMAXLOK,VBM10M, CAPE  , WS1    , WS2    , WS3    , WS4   ,
     &        WS5   , DZR   , DZS   , FKSAT  , FMPOT  , BCLAPP , VPOR  ,
     &        ETRANS, EBSOIL, ESNOW , ESKIN  , ERES   , QI     , QIVI  ,
     &        QIBOXS, PINT  , DWDT   , W      , RPRAC)
      ENDIF
      IF (YTYP.EQ.'M') THEN
         CALL WRITEC4('M'   , YMVARN, NZVM   , NZT+1  , NE     ,
     &        AK    , BK    , TMKVMH, U      , V      , T      , QD    ,
     &        QW    , FI    , VERVEL, ACLC   , ACLCAC , TKE    , EMTER ,
     &        TRSOL , EMTEF , TRSOF , PS     , QDB    , TSECH  , WSECH ,
     &        TSLECH, TSWECH, TSIECH, USTRL  , USTRW  , USTRI  , VSTRL ,
     &        VSTRW , VSTRI , EVAPL , EVAPW  , EVAPI  ,
     &        AHFSL , AHFSW , AHFSI , AZ0L   , AZ0W   , AZ0I   ,
     &        ALSOL , ALSOW , ALSOI , AHFICE , QRES   ,
     &        TMCHL , TMCHW , TMCHI , QDBL   , QDBW   , QDBI   ,
     &        BFLHSL, BFLHSW, BFLHSI, BFLQDSL, BFLQDSW, BFLQDSI, SN    ,
     &        TD    , TDCL  , WL    , TSN    , TD3    , TD4    , TD5   ,
     &        SRADS , TRADS , SRAD0 , TRAD0  , APRL   , APRC   , APRS  ,
     &        VDIS  , AHFS  , FIB   , BLA    , AHFL   , USTAR3 , RUNOFF,
     &        ACLCV , ACLCOV, TMCM  , TMCH   , PHI    , RLA    , BFLHS ,
     &        BFLQDS, U10   , V10   , TEMP2  , DEW2   , TSURF  , WIND10,
     &        AZ0   , ALBECH, ALBEDO, USTR   , VSTR   , EVAP   , EVAPM ,
     &        SRAFS , TRAFS , SRAF0 , TRAF0  , SCLFS  , TCLFS  , SCLF0 ,
     &        TCLF0 , USTRGW, VSTRGW, VDISGW , VGRAT  , VAROR  , VLT   ,
     &        T2MAX , SRAD0U, SRADSU, TRADSU , T2MIN  , SEAICE , SICED ,
     &        FOREST, TEFF  , TSMAX , TSMIN  , WIMAX  , TOPMAX , SNMEL ,
     &        TSLIN , DSNAC , FAO   , RGCGN  , WSMX   , QVI    ,
     &        ALWCVI, GLAC  , DRAIN , SRFL   , QDBOXS , QWBOXS , EKBOXS,
     &        FHBOXS, FIBOXS,TLAMBDA, DLAMBDA, PORVOL , FCAP   , WI3   ,
     &        WI4   , WI5   , WI    , WICL   , GHPBL  , BETA   ,WMINLOK,
     &        WMAXLOK,VBM10M, CAPE  , WS1    , WS2    , WS3    , WS4   ,
     &        WS5   , DZR   , DZS   , FKSAT  , FMPOT  , BCLAPP , VPOR  ,
     &        ETRANS, EBSOIL, ESNOW , ESKIN  , ERES   , QI     , QIVI  ,
     &        QIBOXS, PINT  , DWDT   , W      , RPRAC)
      ENDIF
      IF (YTYP.EQ.'T') THEN
         IF (MYID .EQ. 0) THEN
            CALL SEND2(NUTDAT, YTDNAM, YTDCAT)
         ENDIF
         CALL WRITEC4('T'   , YTVARN, NZVT   , NZT+1  , NE     ,
     &        AK    , BK    , TMKVMH, U      , V      , T      , QD    ,
     &        QW    , FI    , VERVEL, ACLC   , ACLCAC , TKE    , EMTER ,
     &        TRSOL , EMTEF , TRSOF , PS     , QDB    , TSECH  , WSECH ,
     &        TSLECH, TSWECH, TSIECH, USTRL  , USTRW  , USTRI  , VSTRL ,
     &        VSTRW , VSTRI , EVAPL , EVAPW  , EVAPI  ,
     &        AHFSL , AHFSW , AHFSI , AZ0L   , AZ0W   , AZ0I   ,
     &        ALSOL , ALSOW , ALSOI , AHFICE , QRES   ,
     &        TMCHL , TMCHW , TMCHI , QDBL   , QDBW   , QDBI   ,
     &        BFLHSL, BFLHSW, BFLHSI, BFLQDSL, BFLQDSW, BFLQDSI, SN    ,
     &        TD    , TDCL  , WL    , TSN    , TD3    , TD4    , TD5   ,
     &        SRADS , TRADS , SRAD0 , TRAD0  , APRL   , APRC   , APRS  ,
     &        VDIS  , AHFS  , FIB   , BLA    , AHFL   , USTAR3 , RUNOFF,
     &        ACLCV , ACLCOV, TMCM  , TMCH   , PHI    , RLA    , BFLHS ,
     &        BFLQDS, U10   , V10   , TEMP2  , DEW2   , TSURF  , WIND10,
     &        AZ0   , ALBECH, ALBEDO, USTR   , VSTR   , EVAP   , EVAPM ,
     &        SRAFS , TRAFS , SRAF0 , TRAF0  , SCLFS  , TCLFS  , SCLF0 ,
     &        TCLF0 , USTRGW, VSTRGW, VDISGW , VGRAT  , VAROR  , VLT   ,
     &        T2MAX , SRAD0U, SRADSU, TRADSU , T2MIN  , SEAICE , SICED ,
     &        FOREST, TEFF  , TSMAX , TSMIN  , WIMAX  , TOPMAX , SNMEL ,
     &        TSLIN , DSNAC , FAO   , RGCGN  , WSMX   , QVI    ,
     &        ALWCVI, GLAC  , DRAIN , SRFL   , QDBOXS , QWBOXS , EKBOXS,
     &        FHBOXS, FIBOXS,TLAMBDA, DLAMBDA, PORVOL , FCAP   , WI3   ,
     &        WI4   , WI5   , WI    , WICL   , GHPBL  , BETA   ,WMINLOK,
     &        WMAXLOK,VBM10M, CAPE  , WS1    , WS2    , WS3    , WS4   ,
     &        WS5   , DZR   , DZS   , FKSAT  , FMPOT  , BCLAPP , VPOR  ,
     &        ETRANS, EBSOIL, ESNOW , ESKIN  , ERES   , QI     , QIVI  ,
     &        QIBOXS, PINT  , DWDT   , W      , RPRAC)
      ENDIF
      IF (YTYP.EQ.'D') THEN
         IF (MYID .EQ. 0) THEN
            CALL SEND2(NUDDAT, YDDNAM, YDDCAT)
         ENDIF
         CALL WRITEC4('D'   , YDVARN, NZVD   , NZT+1  , NE     ,
     &        AK    , BK    , TMKVMH, U      , V      , T      , QD    ,
     &        QW    , FI    , VERVEL, ACLC   , ACLCAC , TKE    , EMTER ,
     &        TRSOL , EMTEF , TRSOF , PS     , QDB    , TSECH  , WSECH ,
     &        TSLECH, TSWECH, TSIECH, USTRL  , USTRW  , USTRI  , VSTRL ,
     &        VSTRW , VSTRI , EVAPL , EVAPW  , EVAPI  ,
     &        AHFSL , AHFSW , AHFSI , AZ0L   , AZ0W   , AZ0I   ,
     &        ALSOL , ALSOW , ALSOI , AHFICE , QRES   ,
     &        TMCHL , TMCHW , TMCHI , QDBL   , QDBW   , QDBI   ,
     &        BFLHSL, BFLHSW, BFLHSI, BFLQDSL, BFLQDSW, BFLQDSI, SN    ,
     &        TD    , TDCL  , WL    , TSN    , TD3    , TD4    , TD5   ,
     &        SRADS , TRADS , SRAD0 , TRAD0  , APRL   , APRC   , APRS  ,
     &        VDIS  , AHFS  , FIB   , BLA    , AHFL   , USTAR3 , RUNOFF,
     &        ACLCV , ACLCOV, TMCM  , TMCH   , PHI    , RLA    , BFLHS ,
     &        BFLQDS, U10   , V10   , TEMP2  , DEW2   , TSURF  , WIND10,
     &        AZ0   , ALBECH, ALBEDO, USTR   , VSTR   , EVAP   , EVAPM ,
     &        SRAFS , TRAFS , SRAF0 , TRAF0  , SCLFS  , TCLFS  , SCLF0 ,
     &        TCLF0 , USTRGW, VSTRGW, VDISGW , VGRAT  , VAROR  , VLT   ,
     &        T2MAX , SRAD0U, SRADSU, TRADSU , T2MIN  , SEAICE , SICED ,
     &        FOREST, TEFF  , TSMAX , TSMIN  , WIMAX  , TOPMAX , SNMEL ,
     &        TSLIN , DSNAC , FAO   , RGCGN  , WSMX   , QVI    ,
     &        ALWCVI, GLAC  , DRAIN , SRFL   , QDBOXS , QWBOXS , EKBOXS,
     &        FHBOXS, FIBOXS,TLAMBDA, DLAMBDA, PORVOL , FCAP   , WI3   ,
     &        WI4   , WI5   , WI    , WICL   , GHPBL  , BETA   ,WMINLOK,
     &        WMAXLOK,VBM10M, CAPE  , WS1    , WS2    , WS3    , WS4   ,
     &        WS5   , DZR   , DZS   , FKSAT  , FMPOT  , BCLAPP , VPOR  ,
     &        ETRANS, EBSOIL, ESNOW , ESKIN  , ERES   , QI     , QIVI  ,
     &        QIBOXS, PINT  , DWDT   , W      , RPRAC)
      ENDIF
      IF (YTYP.EQ.'F') THEN
         IF (MYID .EQ. 0) THEN
            CALL SEND2(NUFDAT, YFDNAM, YFDCAT)
         ENDIF
         CALL WRITEC4('F'   , YFVARN, NZVF   , NZT    , NJ     ,
     &        AK    , BK    , TMKVMH, U      , V      , T      , QD    ,
     &        QW    , FI    , VERVEL, ACLC   , ACLCAC , TKE    , EMTER ,
     &        TRSOL , EMTEF , TRSOF , PS     , QDB    , TSECH  , WSECH ,
     &        TSLECH, TSWECH, TSIECH, USTRL  , USTRW  , USTRI  , VSTRL ,
     &        VSTRW , VSTRI , EVAPL , EVAPW  , EVAPI  ,
     &        AHFSL , AHFSW , AHFSI , AZ0L   , AZ0W   , AZ0I   ,
     &        ALSOL , ALSOW , ALSOI , AHFICE , QRES   ,
     &        TMCHL , TMCHW , TMCHI , QDBL   , QDBW   , QDBI   ,
     &        BFLHSL, BFLHSW, BFLHSI, BFLQDSL, BFLQDSW, BFLQDSI, SN    ,
     &        TD    , TDCL  , WL    , TSN    , TD3    , TD4    , TD5   ,
     &        SRADS , TRADS , SRAD0 , TRAD0  , APRL   , APRC   , APRS  ,
     &        VDIS  , AHFS  , FIB   , BLA    , AHFL   , USTAR3 , RUNOFF,
     &        ACLCV , ACLCOV, TMCM  , TMCH   , PHI    , RLA    , BFLHS ,
     &        BFLQDS, U10   , V10   , TEMP2  , DEW2   , TSURF  , WIND10,
     &        AZ0   , ALBECH, ALBEDO, USTR   , VSTR   , EVAP   , EVAPM ,
     &        SRAFS , TRAFS , SRAF0 , TRAF0  , SCLFS  , TCLFS  , SCLF0 ,
     &        TCLF0 , USTRGW, VSTRGW, VDISGW , VGRAT  , VAROR  , VLT   ,
     &        T2MAX , SRAD0U, SRADSU, TRADSU , T2MIN  , SEAICE , SICED ,
     &        FOREST, TEFF  , TSMAX , TSMIN  , WIMAX  , TOPMAX , SNMEL ,
     &        TSLIN , DSNAC , FAO   , RGCGN  , WSMX   , QVI    ,
     &        ALWCVI, GLAC  , DRAIN , SRFL   , QDBOXS , QWBOXS , EKBOXS,
     &        FHBOXS, FIBOXS,TLAMBDA, DLAMBDA, PORVOL , FCAP   , WI3   ,
     &        WI4   , WI5   , WI    , WICL   , GHPBL  , BETA   ,WMINLOK,
     &        WMAXLOK,VBM10M, CAPE  , WS1    , WS2    , WS3    , WS4   ,
     &        WS5   , DZR   , DZS   , FKSAT  , FMPOT  , BCLAPP , VPOR  ,
     &        ETRANS, EBSOIL, ESNOW , ESKIN  , ERES   , QI     , QIVI  ,
     &        QIBOXS, PINT  , DWDT   , W      , RPRAC)
      ENDIF
C     WENN FORTSETZUNGSDATEN (YTYP='F'), 2.ZEITPUNKT AUSGEBEN
      IF (YTYP.EQ.'F') THEN
         YFDNAM(8:8) = 'g'
C
         IF (MYID .EQ. 0) THEN
            CALL SEND2(NUFDAT, YFDNAM, YFDCAT)
         ENDIF
C
         CALL WRITEC4('F'   , YFVARN, NZVF   , NZT+1  , NE     ,
     &        AK    , BK    , TMKVMH, U      , V      , T      , QD    ,
     &        QW    , FI    , VERVEL, ACLC   , ACLCAC , TKE    , EMTER ,
     &        TRSOL , EMTEF , TRSOF , PS     , QDB    , TSECH  , WSECH ,
     &        TSLECH, TSWECH, TSIECH, USTRL  , USTRW  , USTRI  , VSTRL ,
     &        VSTRW , VSTRI , EVAPL , EVAPW  , EVAPI  ,
     &        AHFSL , AHFSW , AHFSI , AZ0L   , AZ0W   , AZ0I   ,
     &        ALSOL , ALSOW , ALSOI , AHFICE , QRES   ,
     &        TMCHL , TMCHW , TMCHI , QDBL   , QDBW   , QDBI   ,
     &        BFLHSL, BFLHSW, BFLHSI, BFLQDSL, BFLQDSW, BFLQDSI, SN    ,
     &        TD    , TDCL  , WL    , TSN    , TD3    , TD4    , TD5   ,
     &        SRADS , TRADS , SRAD0 , TRAD0  , APRL   , APRC   , APRS  ,
     &        VDIS  , AHFS  , FIB   , BLA    , AHFL   , USTAR3 , RUNOFF,
     &        ACLCV , ACLCOV, TMCM  , TMCH   , PHI    , RLA    , BFLHS ,
     &        BFLQDS, U10   , V10   , TEMP2  , DEW2   , TSURF  , WIND10,
     &        AZ0   , ALBECH, ALBEDO, USTR   , VSTR   , EVAP   , EVAPM ,
     &        SRAFS , TRAFS , SRAF0 , TRAF0  , SCLFS  , TCLFS  , SCLF0 ,
     &        TCLF0 , USTRGW, VSTRGW, VDISGW , VGRAT  , VAROR  , VLT   ,
     &        T2MAX , SRAD0U, SRADSU, TRADSU , T2MIN  , SEAICE , SICED ,
     &        FOREST, TEFF  , TSMAX , TSMIN  , WIMAX  , TOPMAX , SNMEL ,
     &        TSLIN , DSNAC , FAO   , RGCGN  , WSMX   , QVI    ,
     &        ALWCVI, GLAC  , DRAIN , SRFL   , QDBOXS , QWBOXS , EKBOXS,
     &        FHBOXS, FIBOXS,TLAMBDA, DLAMBDA, PORVOL , FCAP   , WI3   ,
     &        WI4   , WI5   , WI    , WICL   , GHPBL  , BETA   ,WMINLOK,
     &        WMAXLOK,VBM10M, CAPE  , WS1    , WS2    , WS3    , WS4   ,
     &        WS5   , DZR   , DZS   , FKSAT  , FMPOT  , BCLAPP , VPOR  ,
     &        ETRANS, EBSOIL, ESNOW , ESKIN  , ERES   , QI     , QIVI  ,
     &        QIBOXS, PINT  , DWDT   , W      , RPRAC)
C
      ENDIF
C
      END SUBROUTINE ECEPREP
