      SUBROUTINE EC4ORG(NOZ)
      !
      IMPLICIT NONE
      !
      !**** EC4ORG   -   UP:ORGANISATION DER REMO-PROGNOSE FUER EC4-PHYSIK
      !**   AUFRUF   :    --
      !**   ENTRIES  :   KEINE
      !**   ZWECK    :   ORGANISATION DER GESAMTEN REMO-PROGNOSE VOM HOLEN DER
      !**                ANFANGS- UND RANDDATEN BIS ZUM SCHREIBEN DER ERGEB-
      !**                NISSE. ANLEGEN ALLER FELDER.
      !**   VERSIONS-
      !**   DATUM    :   03.03.05
      !**                2007
      !**
      !**   EXTERNALS:   ECAPREP, ECRPREP, EC4INMI, PROGEC4, INBOXS
      !**                GUST   , NEAREC4, ECEPREP, ECACCU
      !**
      !**   EINGABE-
      !**   PARAMETER:   INPUT FUER REMORG UEBER NAMELIST
      !**   AUSGABE-
      !**   PARAMETER:   KEINE
      !**
      !**   COMMON-
      !**   BLOECKE  :   ORG, CORG, COMDYN, COMNMI, COMDIA, COMPHY
      !**
      !**   METHODE  :   SETZEN VON STEUERPARAMETERN ABHAENGIG VOM ZEITSCHRITT
      !**                SEQUENTIELLER AUFRUF VON UP'S ZUR PROGNOSE-DURCH-
      !**                FUEHRUNG
      !**   FEHLERBE-
      !**   HANDLUNG :   STOP IM FEHLERFALL
      !**   VERFASSER:   R.PODZUN
      !
      INCLUDE "parorg.h"
      INCLUDE "org.h"
      INCLUDE "corg.h"
      INCLUDE "comdyn.h"
      INCLUDE "comnmi.h"
      INCLUDE "comdia.h"
      INCLUDE "comphy.h"
      INCLUDE "comecphy.h"
      INCLUDE "faktinf.h"
      INCLUDE "parkon.h"
      INCLUDE "higkon.h"
      INCLUDE "unitnr.h"
      INCLUDE "comhyd.h"
      !
      INTEGER, INTENT(IN) :: NOZ
      !
      ! Local Variables
      !
      INTEGER             :: IJ,NSP
      INTEGER             :: I,J,K
      INTEGER (KIND=8), PARAMETER :: SECONDS_PER_HOUR = 3600
      INTEGER             :: NZTSTD
      CHARACTER(44)       :: YNOTIZ
      !
      REAL    ::
     &           T      (IE,JE,KE,3)       ,
     &           QD     (IE,JE,KE,3)       ,
     &           QW     (IE,JE,KE,3)       ,
     &           U      (IE,JE,KE,3)       ,
     &           V      (IE,JE,KE,3)       ,
     &           FI     (IE,JE,KE,2)       ,
     &           VERVEL (IE,JE,KE)         ,
     &           PS     (IEJE,3)         ,
     &           SEAICE (IEJE)           ,
     &           BLA    (IEJE)           ,
     &           OZONPL (IEJE,NOZ)       ,
     &           SO4ALL (IEJE,KE)        ,
     &           SO4NAT (IEJE,KE)        ,
     &           ZSO4ALL(IEJE,KE,14)     ,
     &           ZSO4NAT(IEJE,KE,14)     ,
     &           ZOZACT (IEJE,NOZ,14)
      INTEGER ::
     &           INFRL   (IEJE)           ,
     &           INFRW   (IEJE)           ,
     &           INFRI   (IEJE)
      REAL    ::
     &           FIB    (IEJE)           ,
     &           AK     (KE1)            ,
     &           BK     (KE1)            ,
     &           AKH    (KE)             ,
     &           BKH    (KE)             ,
     &           DAK    (KE)             ,
     &           DBK    (KE)             ,
     &           PHI    (IEJE)           ,
     &           RLA    (IEJE)           ,
     &           COSLAT (IEJE)
      REAL    ::
     &           GCPHI   (MOJE,2)        ,
     &           GACPHIR (MOJE,2)
      REAL    ::
     &           SINLAT (IEJE)           ,
     &           COSLON (IEJE)           ,
     &           SINLON (IEJE)           ,
     &           RMY    (IEJE,3)         ,
     &           A1T    (KE1)            ,
     &           A2T    (KE1)            ,
     &           ACPHIR (JE,2)           ,
     &           CPHI   (JE,2)           ,
     &           QDB    (IEJE,3)         ,
     &           QDBL   (IEJE,3)         ,
     &           QDBW   (IEJE,3)         ,
     &           QDBI   (IEJE,3)         ,
     &           TS     (IEJE,3)
      REAL    ::
     &           TB     (IEJE,3)         ,
     &           TG     (IEJE,3)         ,
     &           TGL    (IEJE,3)         ,
     &           TGW    (IEJE,3)         ,
     &           TGI    (IEJE,3)         ,
     &           SOTHDT (IEKE,JE,2)      ,
     &           UVTK   (IEKE,JE,2)      ,
     &           TMKVMH (IE*(KE-1),JE,2) ,
     &           TTK    (IEKE,JE)        ,
     &           QDTK   (IEKE,JE)        ,
     &           TTS    (IEKE,JE)        ,
     &           QDTS   (IEKE,JE)        ,
     &           QWTS   (IEKE,JE)
      REAL    ::
     &           TMCM   (IEJE)           ,
     &           TMCH   (IEJE)           ,
     &           TMCHL  (IEJE)           ,
     &           TMCHW  (IEJE)           ,
     &           TMCHI  (IEJE)           ,
     &           GLAC   (IEJE)           ,
     &           SICED  (IEJE)           ,
     &           TEFF   (IEJE)           ,
     &           TSECH  (IEJE,3)         ,
     &           TSLECH (IEJE,3)         ,
     &           TSWECH (IEJE,3)         ,
     &           TSIECH (IEJE,3)         ,
     &           WSECH  (IEJE,3)         ,
     &           SN     (IEJE,3)
      REAL    ::
     &           WL     (IEJE,3)         ,
     &           TD     (IEJE,3)         ,
     &           TDCL   (IEJE,3)         ,
     &           TD3    (IEJE,3)         ,
     &           TD4    (IEJE,3)         ,
     &           TD5    (IEJE,3)         ,
     &           TSN    (IEJE,3)         ,
     &           TSURF  (IEJE)           ,
     &           TSMAX  (IEJE)           ,
     &           TSMIN  (IEJE)
      REAL    ::
     &           TSLIN  (IEJE)           ,
     &           QRES   (IEJE)           ,
     &           DSNAC  (IEJE)           ,
     &           SNMEL  (IEJE)           ,
     &           RUNOFF (IEJE)           ,
     &           DRAIN  (IEJE)           ,
     &           VAROR  (IEJE)           ,
     &           SRFL   (IEJE)           ,
     &           THFL   (IEJE)           ,
     &           QHFL   (IEJE)           ,
     &           XHFL   (IEJE)           ,
     &           RSFC   (IEJE)           ,
     &           SSFC   (IEJE)
      REAL    ::
     &           RSFL   (IEJE)           ,
     &           SSFL   (IEJE)           ,
     &           AHFL   (IEJE)           ,
     &           AHFS   (IEJE)           ,
     &           AHFSL  (IEJE)           ,
     &           AHFSW  (IEJE)           ,
     &           AHFSI  (IEJE)           ,
     &           AHFICE (IEJE)           ,
     &           DHFT   (IEJE)           ,
     &           DHFQW  (IEJE)           ,
     &           DHFQS  (IEJE)           ,
     &           TOPMAX (IEJE)           ,
     &           AZ0    (IEJE)           ,
     &           AZ0L   (IEJE)           ,
     &           AZ0W   (IEJE)           ,
     &           AZ0I   (IEJE)
      REAL    ::
     &           APRC   (IEJE)           ,
     &           APRL   (IEJE)           ,
     &           APRLSAV(IEJE)           ,
     &           APRS   (IEJE)           ,
     &           EVAP   (IEJE)           ,
     &           EVAPM  (IEJE)           ,
     &           EVAPL  (IEJE)           ,
     &           EVAPW  (IEJE)           ,
     &           EVAPI  (IEJE)           ,
     &           ACLC   (IEJEKE)         ,
     &           ACLCAC (IEJEKE)         ,
     &           VGRAT  (IEJE)           ,
     &           FOREST (IEJE)           ,
     &           GHPBL  (IEJE)           ,
     &           BETA   (IEJE)           ,
     &           WMINLOK(IEJE)           ,
     &           WMAXLOK(IEJE)           ,
     &           CAPE   (IEJE)
      REAL    ::
     &           ALBECH (IEJE)           ,
     &           ALBEDO (IEJE)           ,
     &           ALSOL  (IEJE)           ,
     &           ALSOW  (IEJE)           ,
     &           ALSOI  (IEJE)           ,
     &           TKE    (IEJEKE,3)       ,
     &           DEW2   (IEJE)           ,
     &           WSMX   (IEJE)           ,
     &           VLT    (IEJE)           ,
     &           FAO    (IEJE)           ,
     &           RGCGN  (IEJE)           ,
     &           TLAMBDA(IEJE)           ,
     &           DLAMBDA(IEJE)           ,
     &           PORVOL (IEJE)           ,
     &           FCAP   (IEJE)           ,
     &           WI3    (IEJE,3)         ,
     &           WI4    (IEJE,3)         ,
     &           WI5    (IEJE,3)         ,
     &           WI     (IEJE,3)         ,
     &           WICL   (IEJE,3)         ,
     &           TEMP2  (IEJE)
      REAL    ::
     &           T2MAX  (IEJE)           ,
     &           T2MIN  (IEJE)           ,
     &           USTAR3 (IEJE)           ,
     &           USTR   (IEJE)           ,
     &           USTRL  (IEJE)           ,
     &           USTRW  (IEJE)           ,
     &           USTRI  (IEJE)           ,
     &           U10    (IEJE)           ,
     &           VDIS   (IEJE)           ,
     &           VSTR   (IEJE)           ,
     &           VSTRL  (IEJE)           ,
     &           VSTRW  (IEJE)           ,
     &           VSTRI  (IEJE)           ,
     &           V10    (IEJE)           ,
     &           WIND10 (IEJE)           ,
     &           WIMAX  (IEJE)           ,
     &           VBM10M (IEJE)
      REAL    ::
     &           WS1    (IEJE)           ,
     &           WS2    (IEJE)           ,
     &           WS3    (IEJE)           ,
     &           WS4    (IEJE)           ,
     &           WS5    (IEJE)           ,
     &           DZR    (IEJE)           ,
     &           DZS    (IEJE)           ,
     &           FKSAT  (IEJE)           ,
     &           FMPOT  (IEJE)           ,
     &           BCLAPP (IEJE)           ,
     &           VPOR   (IEJE)           ,
     &           ETRANS (IEJE)           ,
     &           EBSOIL (IEJE)           ,
     &           ESNOW  (IEJE)           ,
     &           ESKIN  (IEJE)           ,
     &           ERES   (IEJE)
      REAL    ::
     &           ACLCOV (IEJE)           ,
     &           ALWCVI (IEJE)           ,
     &           QVI    (IEJE)           ,
     &           EMTER  (IEJE*KE1)       ,
     &           TRSOL  (IEJE*KE1)       ,
     &           EMTEF  (IEJE*2)         ,
     &           TRSOF  (IEJE*2)         ,
     &           SCLF0  (IEJE)           ,
     &           SCLFS  (IEJE)           ,
     &           SRAF0  (IEJE)
      REAL    ::
     &           SRAFS  (IEJE)           ,
     &           TCLF0  (IEJE)           ,
     &           TCLFS  (IEJE)           ,
     &           TRAF0  (IEJE)           ,
     &           TRAFS  (IEJE)           ,
     &           ACLCV  (IEJE)           ,
     &           SRADS  (IEJE)           ,
     &           SRADSU (IEJE)           ,
     &           SRAD0  (IEJE)           ,
     &           SRAD0U (IEJE)
      REAL    ::
     &           TRADS  (IEJE)           ,
     &           TRADSU (IEJE)           ,
     &           TRAD0  (IEJE)           ,
     &           USTRGW (IEJE)           ,
     &           VDISGW (IEJE)           ,
     &           VSTRGW (IEJE)           ,
     &           VAR    (IEJE*4)         ,
     &           FC     (IEJE)           ,
     &           BFLHS  (IEJE)           ,
     &           BFLHSL (IEJE)           ,
     &           BFLHSW (IEJE)           ,
     &           BFLHSI (IEJE)
      REAL    ::
     &           BFLQDS (IEJE)           ,
     &           BFLQDSL(IEJE)           ,
     &           BFLQDSW(IEJE)           ,
     &           BFLQDSI(IEJE)           ,
     &           BFLUS  (IEJE)           ,
     &           BFLVS  (IEJE)           ,
     &           VVFH   (KE)             ,
     &           CEVAPCU(KE)             ,
     &           CVDAES (KE1)            ,
     &           CVDAEL (KE1)            ,
     &           CVDAEU (KE1)            ,
     &           CVDAED (KE1)            ,
     &           SISTM  (KE,KE)
      REAL    ::
     &           SIGAM  (KE,KE)          ,
     &           SITAU  (KE,KE)          ,
     &           SINUE  (KE   )          ,
     &           SIVMT  (KE,KE)          ,
     &           SIVMTI (KE,KE)          ,
     &           SICQ   (   KE)          ,
     &           TRIGSI (5*(MOIE-1)/2)     ,
     &           RZ1I   (6)              ,
     &           RZ2I   (6)               
      INTEGER :: IFAXI  (11)
      REAL    ::
     &           TRIGSJ (5*(MOJE-1)/2)     ,
     &           RZ1J   (6)              ,
     &           RZ2J   (6)              ,
     &           IFAXJ  (11)             ,
     &           UR     (IEJEKE,2)       ,
     &           VR     (IEJEKE,2)       ,
     &           TR     (IEJEKE,2)       ,
     &           QDR    (IEJEKE,2)       ,
     &           PSR    (IEJE,2)         ,
     &           QWR    (IEJEKE,2)       ,
     &           ZT2MAX (IEJE)           ,
     &           ZT2MIN (IEJE)           ,
     &           ZWIMAX (IEJE)           ,
     &           ZALB   (IEJE,12)        ,
     &           ZVGR   (IEJE,12)        ,
     &           ZVLT   (IEJE,12)
      REAL    ::
     &           TSWECHR(IEJE,2)         ,
     &           TSIECHR(IEJE,2)         ,
     &           SEAICER(IEJE,2)         ,
     &           QDBLR  (IEJE,2)         ,
     &           SICEDR (IEJE,2)
      REAL    ::
     &           QDBOXS (IEJE)            ,
     &           QWBOXS (IEJE)            ,
     &           EKBOXS (IEJE)            ,
     &           FHBOXS (IEJE)            ,
     &           FIBOXS (IEJE)
      REAL    ::
     &           QIBOXS (IEJE)            ,
     &           QI     (IEJEKE,3)        ,
     &           QIVI   (IEJE)            ,
     &           QITS   (IEKE,JE)         ,
     &           RPRAC  (IEJEKE)          ,
     &           ALPHABOUND(IE,JE,3)
      LOGICAL ::
     &           LOLAND (IEJE)            ,
     &           LOSEA  (IEJE)            ,
     &           LOICE  (IEJE)            ,
     &           LOGLAC (IEJE)            ,
     &           LALAND (IEJE)

C
C     ADDITIONAL FIELDS FOR NON-HYDROSTATIC
C     PINT          - PRESSURE
C     DWDT          - CONVERSION BETWEEN HYDROSTATIC AND NON-HYDROSTATIC
C     ETAS          - DIAGNOSTIC VERTICAL VELOCITY
C     W             - VERTICAL VELOCITY
      REAL ::   PINT(IE,JE,KE1,3), DWDT(IE,JE,KE ,3), 
     &          ETAS(IE,JE,KE1  ), W   (IE,JE,KE1,3)
C
C     STATEMENT - FUNCTIONS ZUR BERECHNUNG DER SAEETIGUNGSFEUCHTE
C     MAGNUS-FORMEL FUER WASSER
C
!      ZGEW(TT)        = B1 * EXP  ( B2W*(TT - B3)/(TT - B4W) )
!      ZGEE(TT)        = B1 * EXP  ( B2E*(TT - B3)/(TT - B4E) )
C
C     SPEZIFISCHE FEUCHTE AUS DAMPFDRUCK UND LUFTDRUCK
C
!      ZGQD(GE,PP)     = RDRD*GE/(PP - EMRDRD*GE)
C
C
      
      YNOTIZ='..... FORECASTTIME  = XXXXXXXX HOURS   .....'
C
C     ZUM STARTZEITPUNKT TENDENZEN INITIALISIEREN
C
      IF (NZT.EQ.0) THEN
         CALL SETRA(QDBOXS,IEJE          ,0.0)
         CALL SETRA(QWBOXS,IEJE          ,0.0)
         CALL SETRA(QIBOXS,IEJE          ,0.0)
         CALL SETRA(EKBOXS,IEJE          ,0.0)
         CALL SETRA(FHBOXS,IEJE          ,0.0)
         CALL SETRA(FIBOXS,IEJE          ,0.0)
         CALL SETRA(WI3   ,IEJE*3        ,0.0)
         CALL SETRA(WI4   ,IEJE*3        ,0.0)
         CALL SETRA(WI5   ,IEJE*3        ,0.0)
         CALL SETRA(WI    ,IEJE*3        ,0.0)
         CALL SETRA(WICL  ,IEJE*3        ,0.0)
      ENDIF
      IF (NZT.EQ.0.OR.NZT.EQ.NANF) THEN
         CALL SETRA(QW    ,IEJEKE*3      ,0.0)
         CALL SETRA(TMKVMH,IE*(KE-1)*JE*2,0.0)
         CALL SETRA(SOTHDT,IEJEKE*2      ,0.0)
         CALL SETRA(UVTK  ,IEJEKE*2      ,0.0)
         CALL SETRA(TTK   ,IEJEKE        ,0.0)
         CALL SETRA(QDTK  ,IEJEKE        ,0.0)
         CALL SETRA(TTS   ,IEJEKE        ,0.0)
         CALL SETRA(QDTS  ,IEJEKE        ,0.0)
         CALL SETRA(QWTS  ,IEJEKE        ,0.0)
         CALL SETRA(OZONPL,IEJE*NOZ      ,0.0)
         CALL SETRA(SO4ALL,IEJEKE        ,0.0)
         CALL SETRA(SO4NAT,IEJEKE        ,0.0)
         CALL SETRA(TMCH  ,IEJE          ,0.0)
         CALL SETRA(TMCHL ,IEJE          ,0.0)
         CALL SETRA(TMCHW ,IEJE          ,0.0)
         CALL SETRA(TMCHI ,IEJE          ,0.0)
         CALL SETRA(TMCM  ,IEJE          ,0.0)
         CALL SETRA(BFLUS ,IEJE          ,0.0)
         CALL SETRA(BFLVS ,IEJE          ,0.0)
         CALL SETRA(BFLQDS,IEJE          ,0.0)
         CALL SETRA(BFLQDSL,IEJE         ,0.0)
         CALL SETRA(BFLQDSW,IEJE         ,0.0)
         CALL SETRA(BFLQDSI,IEJE         ,0.0)
         CALL SETRA(BFLHS ,IEJE          ,0.0)
         CALL SETRA(BFLHSL,IEJE          ,0.0)
         CALL SETRA(BFLHSW,IEJE          ,0.0)
         CALL SETRA(BFLHSI,IEJE          ,0.0)
         CALL SETRA(APRL  ,IEJE          ,0.0)
         CALL SETRA(APRC  ,IEJE          ,0.0)
         CALL SETRA(APRS  ,IEJE          ,0.0)
         CALL SETRA(VDIS  ,IEJE          ,0.0)
         CALL SETRA(AHFS  ,IEJE          ,0.0)
         CALL SETRA(AHFSL ,IEJE          ,0.0)
         CALL SETRA(AHFSW ,IEJE          ,0.0)
         CALL SETRA(AHFSI ,IEJE          ,0.0)
         CALL SETRA(AHFL  ,IEJE          ,0.0)
         CALL SETRA(AHFICE,IEJE          ,0.0)
         CALL SETRA(USTAR3,IEJE          ,0.0)
         CALL SETRA(RUNOFF,IEJE          ,0.0)
         CALL SETRA(DRAIN ,IEJE          ,0.0)
         CALL SETRA(ACLCOV,IEJE          ,0.0)
         CALL SETRA(ACLCV ,IEJE          ,0.0)
         CALL SETRA(U10   ,IEJE          ,0.0)
         CALL SETRA(V10   ,IEJE          ,0.0)
         CALL SETRA(TEMP2 ,IEJE          ,0.0)
         CALL SETRA(DEW2  ,IEJE          ,0.0)
         CALL SETRA(TSURF ,IEJE          ,0.0)
         CALL SETRA(WIND10,IEJE          ,0.0)
         CALL SETRA(ALBEDO,IEJE          ,0.0)
         CALL SETRA(ALSOL ,IEJE          ,0.0)
         CALL SETRA(ALSOW ,IEJE          ,0.0)
         CALL SETRA(ALSOI ,IEJE          ,0.0)
         CALL SETRA(AZ0L  ,IEJE          ,0.0)
         CALL SETRA(AZ0W  ,IEJE          ,0.0)
         CALL SETRA(AZ0I  ,IEJE          ,0.0)
         CALL SETRA(SRADS ,IEJE          ,0.0)
         CALL SETRA(TRADS ,IEJE          ,0.0)
         CALL SETRA(SRAD0 ,IEJE          ,0.0)
         CALL SETRA(TRAD0 ,IEJE          ,0.0)
         CALL SETRA(EVAP  ,IEJE          ,0.0)
         CALL SETRA(EVAPM ,IEJE          ,0.0)
         CALL SETRA(EVAPL ,IEJE          ,0.0)
         CALL SETRA(EVAPW ,IEJE          ,0.0)
         CALL SETRA(EVAPI ,IEJE          ,0.0)
         CALL SETRA(USTR  ,IEJE          ,0.0)
         CALL SETRA(USTRL ,IEJE          ,0.0)
         CALL SETRA(USTRW ,IEJE          ,0.0)
         CALL SETRA(USTRI ,IEJE          ,0.0)
         CALL SETRA(VSTR  ,IEJE          ,0.0)
         CALL SETRA(VSTRL ,IEJE          ,0.0)
         CALL SETRA(VSTRW ,IEJE          ,0.0)
         CALL SETRA(VSTRI ,IEJE          ,0.0)
         CALL SETRA(SRAFS ,IEJE          ,0.0)
         CALL SETRA(TRAFS ,IEJE          ,0.0)
         CALL SETRA(SRAF0 ,IEJE          ,0.0)
         CALL SETRA(TRAF0 ,IEJE          ,0.0)
         CALL SETRA(SCLFS ,IEJE          ,0.0)
         CALL SETRA(TCLFS ,IEJE          ,0.0)
         CALL SETRA(SCLF0 ,IEJE          ,0.0)
         CALL SETRA(TCLF0 ,IEJE          ,0.0)
         CALL SETRA(USTRGW,IEJE          ,0.0)
         CALL SETRA(VSTRGW,IEJE          ,0.0)
         CALL SETRA(VDISGW,IEJE          ,0.0)
         CALL SETRA(SRAD0U,IEJE          ,0.0)
         CALL SETRA(SRADSU,IEJE          ,0.0)
         CALL SETRA(TRADSU,IEJE          ,0.0)
         CALL SETRA(SNMEL ,IEJE          ,0.0)
         CALL SETRA(TSLIN ,IEJE          ,0.0)
         CALL SETRA(QRES  ,IEJE          ,0.0)
         CALL SETRA(DSNAC ,IEJE          ,0.0)
         CALL SETRA(ACLCAC,IEJEKE        ,0.0)
         CALL SETRA(ACLC  ,IEJEKE        ,0.0)
         CALL SETRA(VERVEL,IEJEKE        ,0.0)
         CALL SETRA(GHPBL ,IEJE          ,0.0)
         CALL SETRA(CAPE  ,IEJE          ,0.0)
         CALL SETRA(QVI   ,IEJE          ,0.0)
         CALL SETRA(ALWCVI,IEJE          ,0.0)
         CALL SETRA(RSFL  ,IEJE          ,0.0)
         CALL SETRA(SSFL  ,IEJE          ,0.0)
         CALL SETRA(RSFC  ,IEJE          ,0.0)
         CALL SETRA(SSFC  ,IEJE          ,0.0)
         CALL SETRA(SRFL  ,IEJE          ,0.0)
         CALL SETRA(THFL  ,IEJE          ,0.0)
         CALL SETRA(QHFL  ,IEJE          ,0.0)
         CALL SETRA(XHFL  ,IEJE          ,0.0)
         CALL SETRA(TEFF  ,IEJE          ,0.0)
         CALL SETRA(RGCGN ,IEJE          ,0.0)
         CALL SETRA(TLAMBDA,IEJE         ,0.0)
         CALL SETRA(DLAMBDA,IEJE         ,0.0)
         CALL SETRA(PORVOL,IEJE          ,0.0)
         CALL SETRA(FCAP  ,IEJE          ,0.0)
         CALL SETRA(SICED ,IEJE          ,0.0)
         CALL SETRA(DHFT  ,IEJE          ,0.0)
         CALL SETRA(DHFQW ,IEJE          ,0.0)
         CALL SETRA(DHFQS ,IEJE          ,0.0)
         CALL SETRA(TKE   ,IEJEKE*3      ,0.0)
         CALL SETRA(EMTER ,IEJE*KE1      ,0.0)
         CALL SETRA(TRSOL ,IEJE*KE1      ,0.0)
         CALL SETRA(EMTEF ,IEJE*2        ,0.0)
         CALL SETRA(TRSOF ,IEJE*2        ,0.0)
         CALL SETRA(WIMAX ,IEJE          ,0.0)
         CALL SETRA(VBM10M,IEJE          ,0.0)
         CALL SETRA(TSMAX ,IEJE          ,-99999.)
         CALL SETRA(TSMIN ,IEJE          , 99999.)
         CALL SETRA(T2MAX ,IEJE          ,-99999.)
         CALL SETRA(T2MIN ,IEJE          , 99999.)
         CALL SETRA(TOPMAX,IEJE          , 99999.)
         CALL SETRA(ZT2MAX,IEJE          ,-99999.)
         CALL SETRA(ZT2MIN,IEJE          , 99999.)
         CALL SETRA(ZWIMAX,IEJE          ,0.0)
         CALL SETRA(QDB   ,IEJE*3        ,0.0)
         CALL SETRA(QDBL  ,IEJE*3        ,0.0)
         CALL SETRA(QDBW  ,IEJE*3        ,0.0)
         CALL SETRA(QDBI  ,IEJE*3        ,0.0)
         CALL SETRA(WS1   ,IEJE          ,0.0)
         CALL SETRA(WS2   ,IEJE          ,0.0)
         CALL SETRA(WS3   ,IEJE          ,0.0)
         CALL SETRA(WS4   ,IEJE          ,0.0)
         CALL SETRA(WS5   ,IEJE          ,0.0)
         CALL SETRA(DZR   ,IEJE          ,0.0)
         CALL SETRA(DZS   ,IEJE          ,0.0)
         CALL SETRA(FKSAT ,IEJE          ,0.0)
         CALL SETRA(FMPOT ,IEJE          ,0.0)
         CALL SETRA(BCLAPP,IEJE          ,0.0)
         CALL SETRA(VPOR  ,IEJE          ,0.0)
         CALL SETRA(ETRANS,IEJE          ,0.0)
         CALL SETRA(EBSOIL,IEJE          ,0.0)
         CALL SETRA(ESNOW ,IEJE          ,0.0)
         CALL SETRA(ESKIN ,IEJE          ,0.0)
         CALL SETRA(ERES  ,IEJE          ,0.0)
         CALL SETRA(QI    ,IEJEKE*3      ,0.0)
         CALL SETRA(QIVI  ,IEJE          ,0.0)
         CALL SETRA(QITS  ,IEJEKE        ,0.0)
         CALL SETRA(DWDT  ,IEJEKE*3      ,1.0)
         CALL SETRA(W     ,IEJE*(KE+1)*3 ,0.0)
         CALL SETRA(ETAS  ,IEJE*(KE+1)   ,0.0)
         CALL SETRA(RPRAC ,IEJEKE        ,0.0)
         CALL SETRA(WL    ,IEJE*3        ,0.0)
         IF (.NOT.LPHY) THEN
            CALL SETRA(TSECH  ,IEJE*3        ,0.0)
            CALL SETRA(TSLECH ,IEJE*3        ,0.0)
            CALL SETRA(TSWECH ,IEJE*3        ,0.0)
            CALL SETRA(TSIECH ,IEJE*3        ,0.0)
            CALL SETRA(WSECH  ,IEJE*3        ,0.0)
            CALL SETRA(SN     ,IEJE*3        ,0.0)
            CALL SETRA(TSN    ,IEJE*3        ,0.0)
            CALL SETRA(TD3    ,IEJE*3        ,0.0)
            CALL SETRA(TD4    ,IEJE*3        ,0.0)
            CALL SETRA(TD5    ,IEJE*3        ,0.0)
            CALL SETRA(TD     ,IEJE*3        ,0.0)
            CALL SETRA(TDCL   ,IEJE*3        ,0.0)
            CALL SETRA(TG     ,IEJE*3        ,0.0)
            CALL SETRA(TGL    ,IEJE*3        ,0.0)
            CALL SETRA(TGW    ,IEJE*3        ,0.0)
            CALL SETRA(TGI    ,IEJE*3        ,0.0)
         ENDIF
      ENDIF
      IF (NZT.GT.0) THEN
         CALL SETRA(TKE   ,IEJEKE*3      ,0.0)
         IF (.NOT.LPHY) THEN
            CALL SETRA(SOTHDT,IEJEKE*2      ,0.0)
            CALL SETRA(UVTK  ,IEJEKE*2      ,0.0)
            CALL SETRA(TTK   ,IEJEKE        ,0.0)
            CALL SETRA(QDTK  ,IEJEKE        ,0.0)
            CALL SETRA(TTS   ,IEJEKE        ,0.0)
            CALL SETRA(QDTS  ,IEJEKE        ,0.0)
            CALL SETRA(QWTS  ,IEJEKE        ,0.0)
            CALL SETRA(QITS  ,IEJEKE        ,0.0)
         ENDIF
      ENDIF
C
C     REMO-ANFANGSDATEN HOLEN UND ABSPEICHERN IN DIE BETREFFENDEN
C     FELDER
C
      CALL ECAPREP
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
CHG   ADDITIONAL FIELDS FOR NONHYDROSTATIC
     &     PINT   , DWDT   , W      , RPRAC)
C
C     WENN NZT=0, NUR KOPIE DER ANFANGSDATEN IN DIE ERSTEN RANDDATEN,
C     WENN MIT ANALYSEN ALS RANDWERTEN GERECHNET WIRD
C     WENN NZT>0, 1.RANDDATEI HOLEN UND ABSPEICHERN
C
      NZTRD = (NZT/NDR)*NDR
      CALL ECRPREP(NZTRD , 1     ,
     &     UR   , VR     , TR     , QDR    , QWR   , PSR   ,
     &     QDBLR, TSWECHR, TSIECHR, SEAICER, U     , V     , T     ,
     &     QD   , QW     , PS     , QDBL   , TSWECH, TSIECH, SEAICE,
     &     SICED, SICEDR)
C
C     2.RANDDATEI HOLEN UND ABSPEICHERN
C
      CALL ECRPREP(NZTRD+NDR, 2  ,
     &     UR   , VR     , TR     , QDR    , QWR   , PSR   ,
     &     QDBLR, TSWECHR, TSIECHR, SEAICER, U     , V     , T     ,
     &     QD   , QW     , PS     , QDBL   , TSWECH, TSIECH, SEAICE,
     &     SICED, SICEDR)

C
      IF (LVEG) THEN
C
C     EINLESEN DER ZEITLICH VARIABLEN BODENFELDER
C
      CALL READBOD(ZALB, ZVGR, ZVLT)
      ENDIF
C
      IF (LAEROZ) THEN
C
C     EINLESEN DER ZEITLICH VARIABLEN AEROSOLFELDER
C
      CALL READAER(ZSO4ALL, ZSO4NAT, YAKDAT1)
C
C     EINLESEN DER ZEITLICH VARIABLEN OZONFELDER
C
      CALL READO3(ZOZACT, NOZ, YAKDAT1)
      ENDIF
C
C     NORMAL-MODE-INITIALISIERUNG, FALLS GEWUENSCHT
C
      IF (LANMI)  THEN
          LAISTEP  = .TRUE.
          LDISTEP  = .FALSE.
      ELSE IF (LDNMI) THEN
          LAISTEP  = .FALSE.
          LDISTEP  = .TRUE.
      ENDIF
C
C     INFRL, INFRW UND INFRI BELEGEN,
C     TSECH AUS TSLECH, TSWECH UND TSIECH BERECHNEN
C     WI3, WI4, WI5, WI, WICL INITIALISIEREN
C
      IF (NZT.EQ.0) THEN
C
C     INFRL, INFRW UND INFRI BELEGEN
C
      DO  K = 1, KE1
        DO  J = 1, JE
          DO  I = 1, IE
            IJ = I + ( ( J - 1 ) * IE )
            PINT(I,J,K,1) = AK(K) + BK(K) * PS(IJ,1)
            PINT(I,J,K,2) = AK(K) + BK(K) * PS(IJ,2)
            PINT(I,J,K,3) = PINT(I,J,K,2)
          ENDDO
        ENDDO
      ENDDO
      DO IJ = 1,IEJE
         LOLAND(IJ)=.FALSE.
         LOSEA(IJ) =.FALSE.
         LOICE(IJ) =.FALSE.
         LOGLAC(IJ)=.FALSE.
         LALAND(IJ)=.FALSE.
      ENDDO
      DO IJ = 1,IEJE
         INFRL(IJ) = NINT(BLA(IJ)*FAKINF)
         INFRW(IJ) = NINT((1.-BLA(IJ))*(1.-SEAICE(IJ))*FAKINF)
         INFRI(IJ) = NINT(FAKINF) - INFRL(IJ) - INFRW(IJ)
         IF (INFRL(IJ).GT.0) LOLAND(IJ)=.TRUE.
         IF (INFRW(IJ).GT.0) LOSEA(IJ) =.TRUE.
         IF (INFRI(IJ).GT.0) LOICE(IJ) =.TRUE.
         IF (GLAC(IJ).GT.0.5) LOGLAC(IJ) =.TRUE.
         IF (BLA(IJ).GE.0.5) LALAND(IJ) =.TRUE.
         TSECH(IJ,1) = (FLOAT(INFRL(IJ)) * TSLECH(IJ,1)
     &               +  FLOAT(INFRW(IJ)) * TSWECH(IJ,1)
     &               +  FLOAT(INFRI(IJ)) * TSIECH(IJ,1))*EDFAKINF
         TSECH(IJ,2) = (FLOAT(INFRL(IJ)) * TSLECH(IJ,2)
     &               +  FLOAT(INFRW(IJ)) * TSWECH(IJ,2)
     &               +  FLOAT(INFRI(IJ)) * TSIECH(IJ,2))*EDFAKINF
         QDBW(IJ,1) = ZGQD( ZGEW( TSWECH(IJ,1) ), PS(IJ,1) )
         QDBW(IJ,2) = ZGQD( ZGEW( TSWECH(IJ,2) ), PS(IJ,2) )
         QDBI(IJ,1) = ZGQD( ZGEE( TSIECH(IJ,1) ), PS(IJ,1) )
         QDBI(IJ,2) = ZGQD( ZGEE( TSIECH(IJ,2) ), PS(IJ,2) )
         QDB (IJ,1) = (FLOAT(INFRL(IJ))*QDBL(IJ,1)
     &              +  FLOAT(INFRW(IJ))*QDBW(IJ,1)
     &              +  FLOAT(INFRI(IJ))*QDBI(IJ,1))*EDFAKINF
         QDB (IJ,2) = (FLOAT(INFRL(IJ))*QDBL(IJ,2)
     &              +  FLOAT(INFRW(IJ))*QDBW(IJ,2)
     &              +  FLOAT(INFRI(IJ))*QDBI(IJ,2))*EDFAKINF
C
         IF ((TD3(IJ,NE).GT.TCRITH).OR.(INFRL(IJ).EQ.0)) THEN
            WI3(IJ,NJ) = 0.
            WI3(IJ,NE) = 0.
         ELSE IF (TD3(IJ,NE).LT.TCRITL) THEN
            WI3(IJ,NJ) = 1.
            WI3(IJ,NE) = 1.
         ELSE
            WI3(IJ,NJ) = (TCRITH-TD3(IJ,NJ))/(TCRITH-TCRITL)
            WI3(IJ,NE) = (TCRITH-TD3(IJ,NE))/(TCRITH-TCRITL)
         ENDIF
C
         IF ((TD4(IJ,NE).GT.TCRITH).OR.(INFRL(IJ).EQ.0)) THEN
            WI4(IJ,NJ) = 0.
            WI4(IJ,NE) = 0.
         ELSE IF (TD4(IJ,NE).LT.TCRITL) THEN
            WI4(IJ,NJ) = 1.
            WI4(IJ,NE) = 1.
         ELSE
            WI4(IJ,NJ) = (TCRITH-TD4(IJ,NJ))/(TCRITH-TCRITL)
            WI4(IJ,NE) = (TCRITH-TD4(IJ,NE))/(TCRITH-TCRITL)
         ENDIF
C
         IF ((TD5(IJ,NE).GT.TCRITH).OR.(INFRL(IJ).EQ.0)) THEN
            WI5(IJ,NJ) = 0.
            WI5(IJ,NE) = 0.
         ELSE IF (TD5(IJ,NE).LT.TCRITL) THEN
            WI5(IJ,NJ) = 1.
            WI5(IJ,NE) = 1.
         ELSE
            WI5(IJ,NJ) = (TCRITH-TD5(IJ,NJ))/(TCRITH-TCRITL)
            WI5(IJ,NE) = (TCRITH-TD5(IJ,NE))/(TCRITH-TCRITL)
         ENDIF
C
         IF ((TD(IJ,NE).GT.TCRITH).OR.(INFRL(IJ).EQ.0)) THEN
            WI(IJ,NJ) = 0.
            WI(IJ,NE) = 0.
         ELSE IF (TD(IJ,NE).LT.TCRITL) THEN
            WI(IJ,NJ) = 1.
            WI(IJ,NE) = 1.
         ELSE
            WI(IJ,NJ) = (TCRITH-TD(IJ,NJ))/(TCRITH-TCRITL)
            WI(IJ,NE) = (TCRITH-TD(IJ,NE))/(TCRITH-TCRITL)
         ENDIF
C
         IF ((TDCL(IJ,NE).GT.TCRITH).OR.(INFRL(IJ).EQ.0)) THEN
            WICL(IJ,NJ) = 0.
            WICL(IJ,NE) = 0.
         ELSE IF (TDCL(IJ,NE).LT.TCRITL) THEN
            WICL(IJ,NJ) = 1.
            WICL(IJ,NE) = 1.
         ELSE
            WICL(IJ,NJ) = (TCRITH-TDCL(IJ,NJ))/(TCRITH-TCRITL)
            WICL(IJ,NE) = (TCRITH-TDCL(IJ,NE))/(TCRITH-TCRITL)
         ENDIF
C
      ENDDO
      ENDIF !IF (NZT.EQ.0)

CKS   INITIALIZATIONS FOR NOHALO SWITCH
      IF (NZT.EQ.0.OR.(NZT-1).EQ.NANF) THEN
        TSECH(:,3) = TSECH(:,2)
      ENDIF
C
      IF (NZT.EQ.0 .AND. (LANMI .OR. LDNMI)) THEN

         CALL EC4INMI
     &  (
     &   T      , QD     , QW    , U     , V      , FI     , VERVEL ,
     &   PS     , SEAICE , FIB   , AK    , BK     , AKH    ,
     &   BKH    , DAK    , DBK   , PHI   , RLA    , COSLAT , SINLAT ,
     &   COSLON , SINLON , RMY   , A1T   , A2T    , ACPHIR , CPHI   ,
     &   TGL    , TGW    , TGI   , QDBL  , QDBW   , QDBI   , QDB    ,
     &   TS     , TB     , TG    , SOTHDT, UVTK   , TMKVMH , TTK    ,
     &   QDTK   , TTS    , QDTS  , QWTS  , TMCM   , TMCH   ,
     &   TMCHL  , TMCHW  , TMCHI , SICED ,
     &   TEFF   , TSECH  , TSLECH, TSWECH, TSIECH , WSECH  , SN     ,
     &   WL     , TD     , TDCL  , TD3   , TD4    , TD5    , TSN    ,
     &   TSURF  , TSMAX  , TSMIN , TSLIN , DSNAC  , SNMEL  , RUNOFF ,
     &   VAROR  , SRFL   , THFL  , QHFL  , XHFL   , RSFC   , SSFC   ,
     &   RSFL   , SSFL   , AHFL  , AHFS  , DHFT   , DHFQW  , AHFSL  ,
     &   AHFSW  , AHFSI  , AHFICE, QRES  ,
     &   AZ0L   , AZ0W   , AZ0I  , DHFQS , TOPMAX , AZ0    , APRC   ,
     &   APRL   , APRS   , EVAP  , EVAPM , ACLC   , ACLCAC , VGRAT  ,
     &   FOREST , ALBECH , EVAPL , EVAPW , EVAPI  , ALSOL  , ALSOW  ,
     &   ALSOI  , ALBEDO , TKE   , DEW2  , WSMX   , VLT    , FAO    ,
     &   RGCGN  , TEMP2  , T2MAX , T2MIN , USTAR3 , USTR   ,
     &   U10    , USTRL  , USTRW , USTRI , VSTRL  , VSTRW  , VSTRI  ,
     &   INFRL  , INFRW  , INFRI , VDIS  , VSTR   , V10    , WIND10 ,
     &   WIMAX  , ACLCOV , ALWCVI, QVI   , EMTER  , TRSOL  , EMTEF  ,
     &   TRSOF  , SCLF0  , SCLFS , SRAF0 , SRAFS  , TCLF0  , TCLFS  ,
     &   TRAF0  , TRAFS  , ACLCV , SRADS , SRADSU , SRAD0  , SRAD0U ,
     &   TRADS  , TRADSU , TRAD0 , USTRGW, VDISGW , VSTRGW , VAR    ,
     &   CVDAES , CVDAEL , CVDAEU, CVDAED, CEVAPCU, SISTM  , SIGAM  ,
     &   SITAU  , SINUE  , SIVMT , SIVMTI, SICQ   , TRIGSI , RZ1I   ,
     &   RZ2I   , IFAXI  , UR    , VR    , TR     , QDR    , QWR    ,
     &   VVFH   , FC     , BFLHSL, BFLHSW, BFLHSI , BFLQDSL, BFLQDSW,
     &   BFLQDSI, BFLHS  , BFLQDS, BFLUS , BFLVS  , PSR    , TSWECHR,
     &   TSIECHR, SEAICER, QDBLR , DRAIN , TLAMBDA, DLAMBDA, PORVOL ,
     &   FCAP   , WI3    , WI4   , WI5   , WI     , WICL   , LOGLAC ,
     &   LOLAND , LOSEA  , LOICE , LALAND, SICEDR , GHPBL  , BETA   ,
     &   WMINLOK, WMAXLOK, CAPE  , OZONPL, NOZ    , SO4ALL , SO4NAT ,
     &   WS1    , WS2    , WS3   , WS4   , WS5    , DZR    , DZS    ,
     &   FKSAT  , FMPOT  , BCLAPP, VPOR  , ETRANS , EBSOIL , ESNOW  ,
     &   ESKIN  , ERES   , QI    , QIVI  , QITS   , GCPHI  , GACPHIR,
     &   PINT   , DWDT   , ETAS  , W     , RPRAC  , ALPHABOUND)
      ENDIF
      LAISTEP  = .FALSE.
      LDISTEP  = .FALSE.
C
C     ZEITSCHRITT HALBIEREN, WENN NZT = 0 (VORWAERTSZEITSCHRITT)
C
      IF (NZT.EQ.0) THEN
          DT = DT*0.5
      ENDIF
C
C     BERECHNUNG VON VBMAX FUER PROGEXP (AUCH BEI FORTSETZUNGSLAUF)
C
      CALL VBMXBER(U, V)
C
C     WENN OUTPUT-INTERVALL DER E-DATEIEN NICHT MIT DEM INITIALISIER-
C     UNGSINTERVALL DER INTERVALLBEZOGENEN GROESSEN (WIE Z.B T2MAX)
C     UEBEREINSTIMMT, WARNUNG AUSGEBEN.
C
      IF (MYID .EQ. 0) THEN
         IF (NDEA.NE.NDMXN) PRINT *,'!!!! WARNUNG: NDEA= ',NDEA,
     &        ' UNGLEICH NDMXN= ',NDMXN,' !!!!'
      ENDIF
C
C     AKTUELLE SUMMEN-BZW. QUADRATSUMMEN FELDER FUER DIE VERSCHIEDENEN
C     MITTELWERT BERECHNUNGEN EINLESEN FALLS NOETIG
C
      IF (MYID .EQ. 0) THEN
        IF (NZT.NE.0) THEN
          IF (.NOT. LGMON) THEN
            CALL MAKEPN('H',NZT+1)
            CALL GETD(NUHDAT,YHDNAM,YHDCAT)
            CALL READMF(NUHDAT)
          ENDIF
        ENDIF
      ENDIF
C
      IF(.NOT.LPHY) THEN
         DO IJ=1,IEJE
            TG(IJ,1)    = TSECH(IJ,1)
            TG(IJ,2)    = TSECH(IJ,1)
            TG(IJ,3)    = TSECH(IJ,1)
            TGL(IJ,1)   = TSLECH(IJ,1)
            TGL(IJ,2)   = TSLECH(IJ,1)
            TGL(IJ,3)   = TSLECH(IJ,1)
            TGW(IJ,1)   = TSWECH(IJ,1)
            TGW(IJ,2)   = TSWECH(IJ,1)
            TGW(IJ,3)   = TSWECH(IJ,1)
            TGI(IJ,1)   = TSIECH(IJ,1)
            TGI(IJ,2)   = TSIECH(IJ,1)
            TGI(IJ,3)   = TSIECH(IJ,1)
         ENDDO
      ENDIF
!
!     SET ALPHA FOR THE LATERAL BOUNDARY TREATMENT
!
      CALL GETALPHA(RMY, FAKRMY, ALPHABOUND)
C
C     ZEITSCHLEIFE
C
CKS   10 CONTINUE
C     --------------------------------------------------------
C     --------------------------------------------------------
C     Here starts the time loop                              |
C                                                            |
      DO WHILE (NZT+1.LE.NENDE)
      IF (LVEG) THEN
C
C     BERECHNUNG DER ZEITABHAENGIGEN BODENFELDER
C
         CALL PREPBOD(ALBECH, VGRAT, VLT, ZALB, ZVGR, ZVLT)
      ENDIF
C
      IF (LAEROZ) THEN
C
C     BERECHNUNG DES ZEITABHAENGIGEN AEROSOLFELDES
C
         CALL PREPAER(SO4ALL, ZSO4ALL, SO4NAT, ZSO4NAT)
C
C     BERECHNUNG DES ZEITABHAENGIGEN OZONFELDES AUF DRUCKFLAECHEN
C
         CALL PREPO3(OZONPL, ZOZACT, NOZ)
      ENDIF
C
C     NEUE RANDDATEI HOLEN, FALLS BENOETIGT
C
      IF (NZT+1.GT.NZTRD + NDR) THEN
          NRD1  = 3 - NRD1
          NRD2  = 3 - NRD2
          NZTRD = (NZT/NDR)*NDR
          CALL ECRPREP(NZTRD+NDR, 2  ,
     &         UR   , VR     , TR     , QDR    , QWR   , PSR   ,
     &         QDBLR, TSWECHR, TSIECHR, SEAICER, U     , V     , T     ,
     &         QD   , QW     , PS     , QDBL   , TSWECH, TSIECH, SEAICE,
     &         SICED, SICEDR)
      ENDIF
C
C     STEUERPARAMETER SETZEN FUER:
C     ZEITREIHEN UND MONATSMITTEL (GESAMT-GEBIET)
C
      LWRITEE = NZT+1.GE.NEAA .AND. MOD((NZT+1-NEAA),NDEA).EQ.0
C
C     ERGEBNISAUSGABEN (AUSSCHNITTS-GEBIET)
C
      LWRITED = NZT+1.GE.NDAA .AND. MOD((NZT+1-NDAA),NDDA).EQ.0
C
C     FORTSETZUNGSAUSGABEN (GESAMT-GEBIET)
C
      LWRITEF = NZT+1.GE.NFORA .AND. MOD((NZT+1-NFORA),NDFOR).EQ.0
C
C     ERGEBNISAUSGABEN (GESAMT-GEBIET)
C
      LWRITET = NZT+1.GE.NTAA .AND. MOD((NZT+1-NTAA),NDTA).EQ.0
C
C     ZYKLISCHE VERTAUSCHUNG DER ZEIT-INDIZES
C
      NSP = NA
      NA  = NJ
      NJ  = NE
      NE  = NSP
C
      NA2 = 3 - NA2
      NJ2 = 3 - NJ2
C
C     INFRL, INFRW UND INFRI BELEGEN
C
      DO IJ = 1,IEJE
         LOLAND(IJ)=.FALSE.
         LOSEA(IJ) =.FALSE.
         LOICE(IJ) =.FALSE.
         LOGLAC(IJ)=.FALSE.
         LALAND(IJ)=.FALSE.
      ENDDO
      DO IJ = 1,IEJE
         INFRL(IJ) = NINT (BLA(IJ)*FAKINF)
         INFRW(IJ) = NINT ((1.-BLA(IJ))*(1.-SEAICE(IJ))*FAKINF)
         INFRI(IJ) = NINT(FAKINF) - INFRL(IJ) - INFRW(IJ)
         IF (INFRL(IJ).GT.0) LOLAND(IJ)=.TRUE.
         IF (INFRW(IJ).GT.0) LOSEA(IJ) =.TRUE.
         IF (INFRI(IJ).GT.0) LOICE(IJ) =.TRUE.
         IF (GLAC(IJ).GT.0.5) LOGLAC(IJ) =.TRUE.
         IF (BLA(IJ).GE.0.5) LALAND(IJ) =.TRUE.
      ENDDO
C
C     PROGNOSE DURCHFUEHREN FUER EINEN ZEITSCHRITT
C
      CALL PROGEC4
     &  (
     &   T      , QD     , QW    , U     , V      , FI     , VERVEL ,
     &   PS     , SEAICE , FIB   , AK    , BK     , AKH    ,
     &   BKH    , DAK    , DBK   , PHI   , RLA    , COSLAT , SINLAT ,
     &   COSLON , SINLON , ALPHABOUND,A1T, A2T    , ACPHIR , CPHI   ,
     &   TGL    , TGW    , TGI   , QDBL  , QDBW   , QDBI   , QDB    ,
     &   TS     , TB     , TG    , SOTHDT, UVTK   , TMKVMH , TTK    ,
     &   QDTK   , TTS    , QDTS  , QWTS  , TMCM   , TMCH   ,
     &   TMCHL  , TMCHW  , TMCHI , SICED ,
     &   TEFF   , TSECH  , TSLECH, TSWECH, TSIECH , WSECH  , SN     ,
     &   WL     , TD     , TDCL  , TD3   , TD4    , TD5    , TSN    ,
     &   TSURF  , TSMAX  , TSMIN , TSLIN , DSNAC  , SNMEL  , RUNOFF ,
     &   VAROR  , SRFL   , THFL  , QHFL  , XHFL   , RSFC   , SSFC   ,
     &   RSFL   , SSFL   , AHFL  , AHFS  , DHFT   , DHFQW  , AHFSL  ,
     &   AHFSW  , AHFSI  , AHFICE, QRES  ,
     &   AZ0L   , AZ0W   , AZ0I  , DHFQS , TOPMAX , AZ0    , APRC   ,
     &   APRL   , APRS   , EVAP  , EVAPM , ACLC   , ACLCAC , VGRAT  ,
     &   FOREST , ALBECH , EVAPL , EVAPW , EVAPI  , ALSOL  , ALSOW  ,
     &   ALSOI  , ALBEDO , TKE   , DEW2  , WSMX   , VLT    , FAO    ,
     &   RGCGN  , TEMP2  , T2MAX , T2MIN , USTAR3 , USTR   ,
     &   U10    , USTRL  , USTRW , USTRI , VSTRL  , VSTRW  , VSTRI  ,
     &   INFRL  , INFRW  , INFRI , VDIS  , VSTR   , V10    , WIND10 ,
     &   WIMAX  , ACLCOV , ALWCVI, QVI   , EMTER  , TRSOL  , EMTEF  ,
     &   TRSOF  , SCLF0  , SCLFS , SRAF0 , SRAFS  , TCLF0  , TCLFS  ,
     &   TRAF0  , TRAFS  , ACLCV , SRADS , SRADSU , SRAD0  , SRAD0U ,
     &   TRADS  , TRADSU , TRAD0 , USTRGW, VDISGW , VSTRGW , VAR    ,
     &   CVDAES , CVDAEL , CVDAEU, CVDAED, CEVAPCU, SISTM  , SIGAM  ,
     &   SITAU  , SINUE  , SIVMT , SIVMTI, SICQ   , TRIGSI , RZ1I   ,
     &   RZ2I   , IFAXI  , UR    , VR    , TR     , QDR    , QWR    ,
     &   VVFH   , FC     , BFLHSL, BFLHSW, BFLHSI , BFLQDSL, BFLQDSW,
     &   BFLQDSI, BFLHS  , BFLQDS, BFLUS , BFLVS  , PSR    , TSWECHR,
     &   TSIECHR, SEAICER, QDBLR , DRAIN , TLAMBDA, DLAMBDA, PORVOL ,
     &   FCAP   , WI3    , WI4   , WI5   , WI     , WICL   , LOGLAC ,
     &   LOLAND , LOSEA  , LOICE , LALAND, SICEDR , GHPBL  , BETA   ,
     &   WMINLOK, WMAXLOK, CAPE  , OZONPL, NOZ    , SO4ALL , SO4NAT ,
     &   WS1    , WS2    , WS3   , WS4   , WS5    , DZR    , DZS    ,
     &   FKSAT  , FMPOT  , BCLAPP, VPOR  , ETRANS , EBSOIL , ESNOW  ,
     &   ESKIN  , ERES   , QI    , QIVI  , QITS   , GACPHIR, GCPHI  ,
     &   PINT   , DWDT   , ETAS  , W     , RPRAC)
C
C     EVT.ERGEBNISDATEI (GESAMT-GEBIET) AUSGEBEN
C
C     FLUESSE DURCH DIE RAENDER DES MODELLGEBIETES BERECHNEN,
C     WENN GEWUENSCHT
C
      IF (LINBOX) THEN
C
         CALL INBOXS(
     &        U     , V     , T     , QD    , QW  , PS    ,
     &        FI    , FIB   , DAK   , DBK   , CPHI, QDBOXS,
     &        QWBOXS, EKBOXS, FHBOXS, FIBOXS, QI  , QIBOXS)
C
      ENDIF
C
C     MAXIMALE WINDBOE BERECHNEN
C
      CALL GUST( NE    ,
C     &   IE    , JE    , KE    , IEJE  , IEJEKE, IEKE  , KE1   ,
     &     T     , QD    , U     , V     , PS    , TB    , TG    ,
     &     QDB   , TMCM  , FIB   , AKH   , BKH   , VBM10M)
C
C     DIAGNOSTISCHE GROESSEN BERECHNEN
C
      IF (LNEAR) THEN
C
         CALL NEAREC4(NE,
C     &   IE    , JE    , KE    , IEJE  , IEJEKE, IEKE  , KE1   ,
     &        T     , QD     , U     , V     , PS    , TB    , TG    ,
     &        QDB   , TMCM   , TMCH  , AZ0   , FIB   ,
     &        AKH   , BKH    , U10   , V10   , ZWIMAX, TEMP2 , DEW2  ,
     &        ZT2MIN, ZT2MAX)

         DO IJ=1,IEJE
            WIND10(IJ)=SQRT(U10(IJ)**2+V10(IJ)**2)
            T2MIN(IJ)=ZT2MIN(IJ)
            T2MAX(IJ)=ZT2MAX(IJ)
            WIMAX(IJ)=ZWIMAX(IJ)
            IF (LWRITEE) THEN
               ZWIMAX(IJ)=0.0
               ZT2MIN(IJ)=99999.
               ZT2MAX(IJ)=-99999.
            ENDIF
         ENDDO
      ENDIF

      APRLSAV(:) = APRL(:)

C
      IF (LWRITEE) CALL ECACCU(1,
     &    KE    , IEJE   ,
     &    ACLCAC, TEFF  , BFLHS , BFLQDS , AHFL   , AHFS   , APRC   ,
     &    APRL  , APRS  , EVAP  , DEW2   , TEMP2  , USTAR3 , USTR   ,
     &    U10   , VDIS  , VSTR  , V10    , WIND10 , ACLCOV , ALWCVI ,
     &    QVI   , SCLF0 , SCLFS , SRAF0  , SRAFS  , TCLF0  , TCLFS  ,
     &    TRAF0 , TRAFS , SRADS , SRADSU , SRAD0  , SRAD0U , TRADS  ,
     &    TRADSU, TRAD0 , DSNAC , RUNOFF , SNMEL  , TSLIN  , TSURF  ,
     &    USTRGW, VDISGW, VSTRGW, T2MAX  , T2MIN  , TSMAX  , TSMIN  ,
     &    WIMAX , TOPMAX, DRAIN , USTRL  , USTRW  , USTRI  , VSTRL  ,
     &    VSTRW , VSTRI , EVAPL , EVAPW  , EVAPI  ,
     &    AHFSL , AHFSW , AHFSI , BFLQDSL, BFLQDSW, BFLQDSI,
     &    BFLHSL, BFLHSW, BFLHSI, AHFICE , QRES   , QDBOXS , QWBOXS ,
     &    EKBOXS, FHBOXS, FIBOXS, VBM10M , CAPE   , ETRANS , EBSOIL ,
     &    ESNOW , ESKIN , ERES  , QIVI   , QIBOXS , RPRAC)
      IF (LWRITEE) CALL ECEPREP('E'    ,
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
     &   WMAXLOK, VBM10M,CAPE  , WS1    , WS2    , WS3    , WS4   ,
     &   WS5   , DZR   , DZS   , FKSAT  , FMPOT  , BCLAPP , VPOR  ,
     &   ETRANS, EBSOIL, ESNOW , ESKIN  , ERES   , QI     , QIVI  ,
     &   QIBOXS, PINT  , DWDT  , W      , RPRAC)
C
C     AUFSUMMIEREN FUER TAGESMITTEL (GESAMT-GEBIET), WENN
C     TAGESMITTEL GEWUENSCHT WERDEN
C
      IF (LTAMIT) THEN
C
         IF (LWRITEE) CALL ECEPREP('N'    ,
     &       AK    , BK    , TMKVMH, U      , V      , T      , QD    ,
     &       QW    , FI    , VERVEL, ACLC   , ACLCAC , TKE    , EMTER ,
     &       TRSOL , EMTEF , TRSOF , PS     , QDB    , TSECH  , WSECH ,
     &       TSLECH, TSWECH, TSIECH, USTRL  , USTRW  , USTRI  , VSTRL ,
     &       VSTRW , VSTRI , EVAPL , EVAPW  , EVAPI  ,
     &       AHFSL , AHFSW , AHFSI , AZ0L   , AZ0W   , AZ0I   ,
     &       ALSOL , ALSOW , ALSOI , AHFICE , QRES   ,
     &       TMCHL , TMCHW , TMCHI , QDBL   , QDBW   , QDBI   ,
     &       BFLHSL, BFLHSW, BFLHSI, BFLQDSL, BFLQDSW, BFLQDSI, SN    ,
     &       TD    , TDCL  , WL    , TSN    , TD3    , TD4    , TD5   ,
     &       SRADS , TRADS , SRAD0 , TRAD0  , APRL   , APRC   , APRS  ,
     &       VDIS  , AHFS  , FIB   , BLA    , AHFL   , USTAR3 , RUNOFF,
     &       ACLCV , ACLCOV, TMCM  , TMCH   , PHI    , RLA    , BFLHS ,
     &       BFLQDS, U10   , V10   , TEMP2  , DEW2   , TSURF  , WIND10,
     &       AZ0   , ALBECH, ALBEDO, USTR   , VSTR   , EVAP   , EVAPM ,
     &       SRAFS , TRAFS , SRAF0 , TRAF0  , SCLFS  , TCLFS  , SCLF0 ,
     &       TCLF0 , USTRGW, VSTRGW, VDISGW , VGRAT  , VAROR  , VLT   ,
     &       T2MAX , SRAD0U, SRADSU, TRADSU , T2MIN  , SEAICE , SICED ,
     &       FOREST, TEFF  , TSMAX , TSMIN  , WIMAX  , TOPMAX , SNMEL ,
     &       TSLIN , DSNAC , FAO   , RGCGN  , WSMX   , QVI    ,
     &       ALWCVI, GLAC  , DRAIN , SRFL   , QDBOXS , QWBOXS , EKBOXS,
     &       FHBOXS, FIBOXS,TLAMBDA, DLAMBDA, PORVOL , FCAP   , WI3   ,
     &       WI4   , WI5   , WI    , WICL   , GHPBL  , BETA   , WMINLOK,
     &       WMAXLOK, VBM10M,CAPE  , WS1    , WS2    , WS3    , WS4   ,
     &       WS5   , DZR   , DZS   , FKSAT  , FMPOT  , BCLAPP , VPOR  ,
     &       ETRANS, EBSOIL, ESNOW , ESKIN  , ERES   , QI     , QIVI  ,
     &       QIBOXS, PINT  , DWDT  , W      , RPRAC)
C
      ENDIF
C
C     AUFSUMMIEREN FUER MONATSMITTEL (GESAMT-GEBIET), WENN
C     MONATSMITTEL GEWUENSCHT WERDEN
C
      IF (LMOMIT) THEN
C
         IF (LWRITEE) CALL ECEPREP('M'    ,
     &       AK    , BK    , TMKVMH, U      , V      , T      , QD    ,
     &       QW    , FI    , VERVEL, ACLC   , ACLCAC , TKE    , EMTER ,
     &       TRSOL , EMTEF , TRSOF , PS     , QDB    , TSECH  , WSECH ,
     &       TSLECH, TSWECH, TSIECH, USTRL  , USTRW  , USTRI  , VSTRL ,
     &       VSTRW , VSTRI , EVAPL , EVAPW  , EVAPI  ,
     &       AHFSL , AHFSW , AHFSI , AZ0L   , AZ0W   , AZ0I   ,
     &       ALSOL , ALSOW , ALSOI , AHFICE , QRES   ,
     &       TMCHL , TMCHW , TMCHI , QDBL   , QDBW   , QDBI   ,
     &       BFLHSL, BFLHSW, BFLHSI, BFLQDSL, BFLQDSW, BFLQDSI, SN    ,
     &       TD    , TDCL  , WL    , TSN    , TD3    , TD4    , TD5   ,
     &       SRADS , TRADS , SRAD0 , TRAD0  , APRL   , APRC   , APRS  ,
     &       VDIS  , AHFS  , FIB   , BLA    , AHFL   , USTAR3 , RUNOFF,
     &       ACLCV , ACLCOV, TMCM  , TMCH   , PHI    , RLA    , BFLHS ,
     &       BFLQDS, U10   , V10   , TEMP2  , DEW2   , TSURF  , WIND10,
     &       AZ0   , ALBECH, ALBEDO, USTR   , VSTR   , EVAP   , EVAPM ,
     &       SRAFS , TRAFS , SRAF0 , TRAF0  , SCLFS  , TCLFS  , SCLF0 ,
     &       TCLF0 , USTRGW, VSTRGW, VDISGW , VGRAT  , VAROR  , VLT   ,
     &       T2MAX , SRAD0U, SRADSU, TRADSU , T2MIN  , SEAICE , SICED ,
     &       FOREST, TEFF  , TSMAX , TSMIN  , WIMAX  , TOPMAX , SNMEL ,
     &       TSLIN , DSNAC , FAO   , RGCGN  , WSMX   , QVI    ,
     &       ALWCVI, GLAC  , DRAIN , SRFL   , QDBOXS , QWBOXS , EKBOXS,
     &       FHBOXS, FIBOXS,TLAMBDA, DLAMBDA, PORVOL , FCAP   , WI3   ,
     &       WI4   , WI5   , WI    , WICL   , GHPBL  , BETA   , WMINLOK,
     &       WMAXLOK, VBM10M,CAPE  , WS1    , WS2    , WS3    , WS4   ,
     &       WS5   , DZR   , DZS   , FKSAT  , FMPOT  , BCLAPP , VPOR  ,
     &       ETRANS, EBSOIL, ESNOW , ESKIN  , ERES   , QI     , QIVI  ,
     &       QIBOXS, PINT  , DWDT  , W      , RPRAC)
C
      ENDIF
C
C     AKKUMULIERTE GROESSEN WIEDER ZURUECKKORRIGIEREN
C
      IF (LWRITEE) CALL ECACCU(2,
     &    KE    , IEJE   ,
     &    ACLCAC, TEFF  , BFLHS , BFLQDS , AHFL   , AHFS   , APRC   ,
     &    APRL  , APRS  , EVAP  , DEW2   , TEMP2  , USTAR3 , USTR   ,
     &    U10   , VDIS  , VSTR  , V10    , WIND10 , ACLCOV , ALWCVI ,
     &    QVI   , SCLF0 , SCLFS , SRAF0  , SRAFS  , TCLF0  , TCLFS  ,
     &    TRAF0 , TRAFS , SRADS , SRADSU , SRAD0  , SRAD0U , TRADS  ,
     &    TRADSU, TRAD0 , DSNAC , RUNOFF , SNMEL  , TSLIN  , TSURF  ,
     &    USTRGW, VDISGW, VSTRGW, T2MAX  , T2MIN  , TSMAX  , TSMIN  ,
     &    WIMAX , TOPMAX, DRAIN , USTRL  , USTRW  , USTRI  , VSTRL  ,
     &    VSTRW , VSTRI , EVAPL , EVAPW  , EVAPI  ,
     &    AHFSL , AHFSW , AHFSI , BFLQDSL, BFLQDSW, BFLQDSI,
     &    BFLHSL, BFLHSW, BFLHSI, AHFICE , QRES   , QDBOXS , QWBOXS ,
     &    EKBOXS, FHBOXS, FIBOXS, VBM10M , CAPE   , ETRANS , EBSOIL ,
     &    ESNOW , ESKIN , ERES  , QIVI   , QIBOXS , RPRAC)
C
C     EVT.ERGEBNISDATEI (AUSSCHNITTS-GEBIET) AUSGEBEN
C
      IF (LWRITED) CALL ECACCU(1,
     &    KE    , IEJE   ,
     &    ACLCAC, TEFF  , BFLHS , BFLQDS , AHFL   , AHFS   , APRC   ,
     &    APRL  , APRS  , EVAP  , DEW2   , TEMP2  , USTAR3 , USTR   ,
     &    U10   , VDIS  , VSTR  , V10    , WIND10 , ACLCOV , ALWCVI ,
     &    QVI   , SCLF0 , SCLFS , SRAF0  , SRAFS  , TCLF0  , TCLFS  ,
     &    TRAF0 , TRAFS , SRADS , SRADSU , SRAD0  , SRAD0U , TRADS  ,
     &    TRADSU, TRAD0 , DSNAC , RUNOFF , SNMEL  , TSLIN  , TSURF  ,
     &    USTRGW, VDISGW, VSTRGW, T2MAX  , T2MIN  , TSMAX  , TSMIN  ,
     &    WIMAX , TOPMAX, DRAIN , USTRL  , USTRW  , USTRI  , VSTRL  ,
     &    VSTRW , VSTRI , EVAPL , EVAPW  , EVAPI  ,
     &    AHFSL , AHFSW , AHFSI , BFLQDSL, BFLQDSW, BFLQDSI,
     &    BFLHSL, BFLHSW, BFLHSI, AHFICE , QRES   , QDBOXS , QWBOXS ,
     &    EKBOXS, FHBOXS, FIBOXS, VBM10M , CAPE   , ETRANS , EBSOIL ,
     &    ESNOW , ESKIN , ERES  , QIVI   , QIBOXS , RPRAC)
      IF (LWRITED) CALL ECEPREP('D'    ,
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
     &   WMAXLOK, VBM10M,CAPE  , WS1    , WS2    , WS3    , WS4   ,
     &   WS5   , DZR   , DZS   , FKSAT  , FMPOT  , BCLAPP , VPOR  ,
     &   ETRANS, EBSOIL, ESNOW , ESKIN  , ERES   , QI     , QIVI  ,
     &   QIBOXS, PINT  , DWDT  , W      , RPRAC)
C
C     AKKUMULIERTE GROESSEN WIEDER ZURUECKKORRIGIEREN
C
      IF (LWRITED) CALL ECACCU(2,
     &    KE    , IEJE   ,
     &    ACLCAC, TEFF  , BFLHS , BFLQDS , AHFL   , AHFS   , APRC   ,
     &    APRL  , APRS  , EVAP  , DEW2   , TEMP2  , USTAR3 , USTR   ,
     &    U10   , VDIS  , VSTR  , V10    , WIND10 , ACLCOV , ALWCVI ,
     &    QVI   , SCLF0 , SCLFS , SRAF0  , SRAFS  , TCLF0  , TCLFS  ,
     &    TRAF0 , TRAFS , SRADS , SRADSU , SRAD0  , SRAD0U , TRADS  ,
     &    TRADSU, TRAD0 , DSNAC , RUNOFF , SNMEL  , TSLIN  , TSURF  ,
     &    USTRGW, VDISGW, VSTRGW, T2MAX  , T2MIN  , TSMAX  , TSMIN  ,
     &    WIMAX , TOPMAX, DRAIN , USTRL  , USTRW  , USTRI  , VSTRL  ,
     &    VSTRW , VSTRI , EVAPL , EVAPW  , EVAPI  ,
     &    AHFSL , AHFSW , AHFSI , BFLQDSL, BFLQDSW, BFLQDSI,
     &    BFLHSL, BFLHSW, BFLHSI, AHFICE , QRES   , QDBOXS , QWBOXS ,
     &    EKBOXS, FHBOXS, FIBOXS, VBM10M , CAPE   , ETRANS , EBSOIL ,
     &    ESNOW , ESKIN , ERES  , QIVI   , QIBOXS , RPRAC)
C
C     EVT.ERGEBNISDATEI (GESAMT-GEBIET) FUER TRAJEKTORIEN AUSGEBEN
C
      IF (LWRITET) CALL ECACCU(1,
     &    KE    , IEJE   ,
     &    ACLCAC, TEFF  , BFLHS , BFLQDS , AHFL   , AHFS   , APRC   ,
     &    APRL  , APRS  , EVAP  , DEW2   , TEMP2  , USTAR3 , USTR   ,
     &    U10   , VDIS  , VSTR  , V10    , WIND10 , ACLCOV , ALWCVI ,
     &    QVI   , SCLF0 , SCLFS , SRAF0  , SRAFS  , TCLF0  , TCLFS  ,
     &    TRAF0 , TRAFS , SRADS , SRADSU , SRAD0  , SRAD0U , TRADS  ,
     &    TRADSU, TRAD0 , DSNAC , RUNOFF , SNMEL  , TSLIN  , TSURF  ,
     &    USTRGW, VDISGW, VSTRGW, T2MAX  , T2MIN  , TSMAX  , TSMIN  ,
     &    WIMAX , TOPMAX, DRAIN , USTRL  , USTRW  , USTRI  , VSTRL  ,
     &    VSTRW , VSTRI , EVAPL , EVAPW  , EVAPI  ,
     &    AHFSL , AHFSW , AHFSI , BFLQDSL, BFLQDSW, BFLQDSI,
     &    BFLHSL, BFLHSW, BFLHSI, AHFICE , QRES   , QDBOXS , QWBOXS ,
     &    EKBOXS, FHBOXS, FIBOXS, VBM10M , CAPE   , ETRANS , EBSOIL ,
     &    ESNOW , ESKIN , ERES  , QIVI   , QIBOXS , RPRAC)
      IF (LWRITET) CALL ECEPREP('T'    ,
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
     &   WMAXLOK, VBM10M,CAPE  , WS1    , WS2    , WS3    , WS4   ,
     &   WS5   , DZR   , DZS   , FKSAT  , FMPOT  , BCLAPP , VPOR  ,
     &   ETRANS, EBSOIL, ESNOW , ESKIN  , ERES   , QI     , QIVI  ,
     &   QIBOXS, PINT  , DWDT  , W      , RPRAC)
C
C     AKKUMULIERTE GROESSEN WIEDER ZURUECKKORRIGIEREN
C
      IF (LWRITET) CALL ECACCU(2,
     &    KE    , IEJE   ,
     &    ACLCAC, TEFF  , BFLHS , BFLQDS , AHFL   , AHFS   , APRC   ,
     &    APRL  , APRS  , EVAP  , DEW2   , TEMP2  , USTAR3 , USTR   ,
     &    U10   , VDIS  , VSTR  , V10    , WIND10 , ACLCOV , ALWCVI ,
     &    QVI   , SCLF0 , SCLFS , SRAF0  , SRAFS  , TCLF0  , TCLFS  ,
     &    TRAF0 , TRAFS , SRADS , SRADSU , SRAD0  , SRAD0U , TRADS  ,
     &    TRADSU, TRAD0 , DSNAC , RUNOFF , SNMEL  , TSLIN  , TSURF  ,
     &    USTRGW, VDISGW, VSTRGW, T2MAX  , T2MIN  , TSMAX  , TSMIN  ,
     &    WIMAX , TOPMAX, DRAIN , USTRL  , USTRW  , USTRI  , VSTRL  ,
     &    VSTRW , VSTRI , EVAPL , EVAPW  , EVAPI  ,
     &    AHFSL , AHFSW , AHFSI , BFLQDSL, BFLQDSW, BFLQDSI,
     &    BFLHSL, BFLHSW, BFLHSI, AHFICE , QRES   , QDBOXS , QWBOXS ,
     &    EKBOXS, FHBOXS, FIBOXS, VBM10M , CAPE   , ETRANS , EBSOIL ,
     &    ESNOW , ESKIN , ERES  , QIVI   , QIBOXS , RPRAC)
C
C     BERECHNUNG VON VBMAX FUER PROGEXP
C
      CALL VBMXBER(U, V)
CKS
C     GET MINIMUM TEMPERATURE (LOWER 150K WILL ISSUE A WARNING)
C
      CALL TEMPMIN(T)
CKS
C
C     WENN GEFORDERT, IN DIESEM ZEITSCHRITT DIE INTERVALLBEZOGENEN
C     GROESSEN NEU INITIALISIEREN
C
      IF (MOD(NZT+1,NDMXN).EQ.0) CALL ECACCU(3 ,
     &    KE    , IEJE   ,
     &    ACLCAC, TEFF  , BFLHS , BFLQDS , AHFL   , AHFS   , APRC   ,
     &    APRL  , APRS  , EVAP  , DEW2   , TEMP2  , USTAR3 , USTR   ,
     &    U10   , VDIS  , VSTR  , V10    , WIND10 , ACLCOV , ALWCVI ,
     &    QVI   , SCLF0 , SCLFS , SRAF0  , SRAFS  , TCLF0  , TCLFS  ,
     &    TRAF0 , TRAFS , SRADS , SRADSU , SRAD0  , SRAD0U , TRADS  ,
     &    TRADSU, TRAD0 , DSNAC , RUNOFF , SNMEL  , TSLIN  , TSURF  ,
     &    USTRGW, VDISGW, VSTRGW, T2MAX  , T2MIN  , TSMAX  , TSMIN  ,
     &    WIMAX , TOPMAX, DRAIN , USTRL  , USTRW  , USTRI  , VSTRL  ,
     &    VSTRW , VSTRI , EVAPL , EVAPW  , EVAPI  ,
     &    AHFSL , AHFSW , AHFSI , BFLQDSL, BFLQDSW, BFLQDSI,
     &    BFLHSL, BFLHSW, BFLHSI, AHFICE , QRES   , QDBOXS , QWBOXS ,
     &    EKBOXS, FHBOXS, FIBOXS, VBM10M , CAPE   , ETRANS , EBSOIL ,
     &    ESNOW , ESKIN , ERES  , QIVI   , QIBOXS , RPRAC)
C
C     EVT.FORTSETZUNGSDATEI (GESAMT-GEBIET) AUSGEBEN
C
      IF (LWRITEF) CALL ECEPREP('F'    ,
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
     &   WMAXLOK, VBM10M,CAPE  , WS1    , WS2    , WS3    , WS4   ,
     &   WS5   , DZR   , DZS   , FKSAT  , FMPOT  , BCLAPP , VPOR  ,
     &   ETRANS, EBSOIL, ESNOW , ESKIN  , ERES   , QI     , QIVI  ,
     &   QIBOXS, PINT  , DWDT  , W      , RPRAC)
C
C     AKTUELLE SUMMEN-BZW. QUADRATSUMMEN FELDER FUER DIE VERSCHIEDENEN
C     MITTELWERT BERECHNUNGEN SPEICHERN
C
      IF (MYID .EQ. 0) THEN
         IF (.NOT. LGMON .AND. LWRITEF) THEN
            CALL MAKEPN('H',NZT+1)
            CALL SEND2(NUHDAT,YHDNAM,YHDCAT)
            CALL WRITEMF(NUHDAT)
         ENDIF
      ENDIF
C
C     VORHERSAGESTUNDE IM DAYFILE PROTOKOLLIEREN
C
      IF (MYID .EQ. 0) THEN
         IF (MOD(NINT((NZT+1)*DT, 8),SECONDS_PER_HOUR).EQ.0) THEN
            NZTSTD = NINT((NZT+1)/3600.0*DT)
            WRITE(YNOTIZ(23:30),'(I8.8)') NZTSTD
            CALL REMARK(YNOTIZ)
         ENDIF
      ENDIF
C
C     NAECHSTEN ZEITSCHRITT BERECHNEN
C
      IF (NZT.EQ.0) THEN
        DT = DT*2.0
      ENDIF
      NZT    = NZT + 1

C
C     NEUES AKTUELLES DATUM BERECHNEN
C
      CALL DATUTC(NZT, YADAT, DT, YAKDAT1, YAKDAT2, NAKJATA, AKHH)
C
C     BERECHNUNG DES NAECHSTEN ZEITSCHRITTS ODER JOBENDE
C
CKS      IF (NZT+1.LE.NENDE) THEN
CKS         GOTO 10
CKS      ENDIF
      ENDDO ! DO WHILE (NZT+1.LE.NENDE)
C                                                            |
C     End of time loop                                       | 
C     --------------------------------------------------------
C     --------------------------------------------------------
C
      RETURN
         !
         !
         !
         CONTAINS
         !
         ! STATEMENT - FUNCTIONS ZUR BERECHNUNG DER SAEETIGUNGSFEUCHTE
         ! MAGNUS-FORMEL FUER WASSER
         !
         REAL FUNCTION ZGEW(TT)
            IMPLICIT NONE
            REAL, INTENT(IN) :: TT
            ZGEW = B1 * EXP  ( B2W*(TT - B3)/(TT - B4W) )
         END FUNCTION ZGEW
         !
         REAL FUNCTION ZGEE(TT)
            IMPLICIT NONE
            REAL, INTENT(IN) :: TT
            ZGEE = B1 * EXP  ( B2E*(TT - B3)/(TT - B4E) )
         END FUNCTION ZGEE
         !
         ! SPEZIFISCHE FEUCHTE AUS DAMPFDRUCK UND LUFTDRUCK
         !
         REAL FUNCTION ZGQD(GE,PP)
            IMPLICIT NONE
            REAL, INTENT(IN) :: GE,PP
            ZGQD = RDRD*GE/(PP - EMRDRD*GE)
         END FUNCTION ZGQD
         !
         !
         ! 
      END SUBROUTINE EC4ORG
