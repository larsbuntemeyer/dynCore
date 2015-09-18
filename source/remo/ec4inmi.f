C
C     SUBROUTINE EC4INMI
C
C**** EC4INMI  -   UP:IMPLIZITE NORMAL MODE INITIALISIERUNG 1. ORDNUNG
C**   AUFRUF   :   CALL EC4INMI
C**   ENTRIES  :   KEINE
C**   ZWECK    :   IMPLIZITE NORMAL MODE INITIALISIERUNG DER
C**                U,V,T,PS ANFANGSFELDER VON EM
C**   VERSIONS-
C**   DATUM    :   06.02.2007
C**
C**   EXTERNALS:   PROGEC4, GEOPOT, SETRA, GAELKO, SOLVE, CONGRA,
C**                MINV   , MXM   , INMI
C**   EINGABE-
C**   PARAMETER:   ---
C**   AUSGABE-
C**   PARAMETER:   ---
C**
C**   COMMON-
C**   BLOECKE  :   ORG   , PARAM ,  COMDYN,  COMNMI,  COMPHY,  HIGKON,
C**                COMDIA, PHYKON,  PARKON,  UNITNR,
C**                COMPCST, COMPMLF, COMPSLF
C**
C**   METHODE  :   NACH WERGEN
C**                RANDWERTE FUER SOLVE SIND IMMER NULL
C**
C**   FEHLERBE-
C**   HANDLUNG :   BERECHNUNG VON BAL FUER JEDE ITERATION UND
C**                JEDEN BEHANDELTEN MODE
C**
C**   VERFASSER:   P.PROHL
C
      SUBROUTINE EC4INMI
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
     &   PINT   , DWDT  , ETAS  , W      , RPRAC  , ALPHABOUND)
C
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
      INCLUDE "org.h"
      INCLUDE "comdyn.h"
      INCLUDE "comnmi.h"
      INCLUDE "comphy.h"
      INCLUDE "higkon.h"
      INCLUDE "comdia.h"
      INCLUDE "phykon.h"
      INCLUDE "parkon.h"
      INCLUDE "unitnr.h"
C
      INCLUDE "haloexchorg"
C
      LOGICAL      LZSITS,LZPHY,LZWRITD,LZWRITE,LZWRITF,LZWRITT

C-----------------------------------------------------------------------
C
C     PHYEC
C
      INTEGER, INTENT(IN) :: NOZ
C
      REAL, INTENT(INOUT) ::
     &          T  (IE,JE,KE,3) , QD(IE,JE,KE,3), QW(IE,JE,KE,3),
     &          TKE(IEJEKE,3)   , FI(IEJEKE,2)
C
      REAL, INTENT(INOUT) ::
     &          U(IE,JE,KE,3), V(IE,JE,KE,3)
C
      REAL, INTENT(IN) ::
     &          ACLC(IEJEKE), ACLCAC(IEJEKE), VERVEL(IEJEKE),
     &          OZONPL(IEJE*NOZ), SO4ALL(IEJEKE), SO4NAT(IEJEKE)
C
      REAL, INTENT(IN) ::
     &          EMTER(IEJE*KE1),TRSOL(IEJE*KE1),
     &          EMTEF(IEJE*2)  ,TRSOF(IEJE*2)  , VAR(IEJE*4)
C
      REAL, INTENT(IN) ::
     &          AK     (KE1), BK    (KE1), AKH   (KE)  , BKH   (KE)  ,
     &          DAK    (KE) , DBK   (KE) , PHI   (IEJE), RLA   (IEJE),
     &          A1T    (KE1), A2T   (KE1), ACPHIR(JE,2), CPHI  (JE,2),
     &          CVDAES (KE1), CVDAEL(KE1), CVDAEU(KE1) , CVDAED(KE1) ,
     &          CEVAPCU(KE)
C
      REAL, INTENT(IN) ::   GACPHIR(MOJE,2), GCPHI(MOJE,2)

      REAL, INTENT(IN) ::
     &          COSLAT(IEJE), SINLAT(IEJE), COSLON(IEJE),
     &          SINLON(IEJE), RMY(IE,JE,3)
C
      REAL, INTENT(IN) ::
     &          TMKVMH(IE*(KE-1),JE,2), SOTHDT(IEKE,JE,2),
     &          UVTK  (IEKE,JE,2)     , TTK   (IEKE,JE)  ,
     &          QDTK  (IEKE,JE)       , TTS   (IEKE,JE)  ,
     &          QDTS  (IEKE,JE)       , QWTS  (IEKE,JE)  ,
     &          TMCM  (IEJE)          , TMCH  (IEJE)     ,
     &          TMCHL (IEJE)          , TMCHW (IEJE)     ,
     &          TMCHI (IEJE)
C
      REAL, INTENT(INOUT) ::
     &          QDB   (IEJE,3), TSN   (IEJE,3), PS   (IE,JE,3),
     &          QDBL  (IEJE,3), QDBW  (IEJE,3), QDBI  (IEJE,3),
     &          TSECH (IEJE,3), WSECH (IEJE,3), SN    (IEJE,3),
     &          TSLECH(IEJE,3), TSWECH(IEJE,3), TSIECH(IEJE,3),
     &          WL    (IEJE,3), TD    (IEJE,3), TDCL  (IEJE,3),
     &          TD3   (IEJE,3), TD4   (IEJE,3), TD5   (IEJE,3),
     &          TS    (IEJE,3), TB    (IEJE,3), TG    (IEJE,3),
     &          TGL   (IEJE,3), TGW   (IEJE,3), TGI   (IEJE,3)
C
      REAL, INTENT(IN) ::
     &          FIB   (IEJE), SEAICE(IEJE),
     &          SICED (IEJE), TEFF  (IEJE), CAPE  (IEJE),
     &          TSMAX (IEJE), TSMIN (IEJE), TSLIN (IEJE),
     &          DSNAC (IEJE), SNMEL (IEJE), RUNOFF(IEJE),
     &          VAROR (IEJE), SRFL  (IEJE), THFL  (IEJE),
     &          QHFL  (IEJE), XHFL  (IEJE), RSFC  (IEJE),
     &          SSFC  (IEJE), RSFL  (IEJE), SSFL  (IEJE),
     &          AHFL  (IEJE), AHFS  (IEJE), DHFT  (IEJE),
     &          AHFSL (IEJE), AHFSW (IEJE), AHFSI (IEJE),
     &          AZ0L  (IEJE), AZ0W  (IEJE), AZ0I  (IEJE),
     &          AHFICE(IEJE), QRES  (IEJE), TSURF (IEJE),
     &          DHFQW (IEJE), DHFQS (IEJE), TOPMAX(IEJE),
     &          AZ0   (IEJE), APRC  (IEJE),
     &          APRL  (IEJE), APRS  (IEJE), EVAP  (IEJE),
     &          EVAPM (IEJE), VGRAT (IEJE), FOREST(IEJE),
     &          EVAPL (IEJE), EVAPW (IEJE), EVAPI (IEJE),
     &          ALSOL (IEJE), ALSOW (IEJE), ALSOI (IEJE),
     &          ALBECH(IEJE), ALBEDO(IEJE), DEW2  (IEJE),
     &          WSMX  (IEJE), VLT   (IEJE), FAO   (IEJE),
     &          RGCGN (IEJE), TEMP2 (IEJE),
     &          T2MAX (IEJE), T2MIN (IEJE), USTAR3(IEJE),
     &          USTR  (IEJE), U10   (IEJE), VDIS  (IEJE),
     &          VSTR  (IEJE), V10   (IEJE), WIND10(IEJE),
     &          USTRL (IEJE), USTRW (IEJE), USTRI (IEJE),
     &          VSTRL (IEJE), VSTRW (IEJE), VSTRI (IEJE),
     &          WIMAX (IEJE), ACLCOV(IEJE), ALWCVI(IEJE),
     &          QVI   (IEJE), SCLF0 (IEJE), SCLFS (IEJE),
     &          SRAF0 (IEJE), SRAFS (IEJE), TCLF0 (IEJE),
     &          TCLFS (IEJE), TRAF0 (IEJE), TRAFS (IEJE),
     &          ACLCV (IEJE), SRADS (IEJE), SRADSU(IEJE),
     &          SRAD0 (IEJE), SRAD0U(IEJE), TRADS (IEJE),
     &          TRADSU(IEJE), TRAD0 (IEJE), USTRGW(IEJE),
     &          VDISGW(IEJE), VSTRGW(IEJE), DRAIN (IEJE),
     &          TLAMBDA(IEJE),DLAMBDA(IEJE),PORVOL(IEJE),
     &          FCAP  (IEJE), WI3 (IEJE,3), WI4 (IEJE,3),
     &          WI5 (IEJE,3), WI  (IEJE,3), WICL(IEJE,3),
     &          GHPBL (IEJE), BETA  (IEJE), WMINLOK(IEJE),
     &          WMAXLOK(IEJE)

      REAL  ::  PINT(IE,JE,KE1,3), DWDT(IE,JE,KE ,3), 
     &          ETAS(IE,JE,KE1  ), W   (IE,JE,KE1,3)
C
      REAL, INTENT(IN) ::
     &          WS1    (IEJE), WS2  (IEJE), WS3   (IEJE),
     &          WS4    (IEJE), WS5  (IEJE), DZR   (IEJE),
     &          DZS    (IEJE), FKSAT(IEJE), FMPOT (IEJE),
     &          BCLAPP (IEJE), VPOR (IEJE), ETRANS(IEJE),
     &          EBSOIL (IEJE), ESNOW(IEJE), ESKIN (IEJE),
     &          ERES   (IEJE)
C
      LOGICAL, INTENT(IN) ::
     &          LOLAND(IEJE), LOSEA (IEJE), LOICE (IEJE),
     &          LOGLAC(IEJE), LALAND(IEJE)
C
C     SIKOR
C
      REAL, INTENT(IN) ::
     &          SISTM(KE,KE), SIGAM (KE,KE), SITAU(KE,KE), SINUE(KE),
     &          SIVMT(KE,KE), SIVMTI(KE,KE), SICQ (KE)
      REAL, INTENT(IN) ::
     &          TRIGSI(5*(MOIE-1)/2), RZ1I(6), RZ2I(6)
      INTEGER, INTENT(IN) :: IFAXI(11)
C
C     ECRANDAS
C
      REAL, INTENT(IN) ::
     &          UR (IEJEKE,2), VR (IEJEKE,2), TR(IEJEKE,2),
     &          QDR(IEJEKE,2), QWR(IEJEKE,2)
C
C     PROGEXP
C
      REAL, INTENT(IN) ::
     &          VVFH   (KE)  , FC    (IE,JE), BFLHS  (IEJE),
     &          BFLQDS (IEJE), BFLUS  (IEJE), BFLVS  (IEJE),
     &          BFLHSL (IEJE), BFLHSW (IEJE), BFLHSI (IEJE),
     &          BFLQDSL(IEJE), BFLQDSW(IEJE), BFLQDSI(IEJE)
C
C     ECRANDUP
C
      REAL, INTENT(IN) ::
     &          PSR    (IEJE,2), TSWECHR(IEJE,2), TSIECHR(IEJE,2),
     &          SEAICER(IEJE,2), QDBLR  (IEJE,2), SICEDR (IEJE,2)
      INTEGER, INTENT(IN) ::
     &          INFRL  (IEJE)  , INFRW  (IEJE)  , INFRI  (IEJE)
C
CSP
C
      REAL, INTENT(INOUT) ::   
     &          QI(IE,JE,KE,3), QITS  (IEKE,JE), QIVI(IEJE),
     &          RPRAC(IEJEKE), ALPHABOUND(IE,JE,3)
C
C     LOKALE FELDER
C
      REAL ::   SIGAMI(KE,2*KE)
      REAL ::   SISTMI(KE,2*KE)
      REAL ::   ZDUM2 (IE,JE)
      REAL ::   ZDUM3 (IE,JE)
      REAL ::   ZDPS  (IE,JE)
      REAL ::   ZDU   (IE,JE,KE)
      REAL ::   ZDV   (IE,JE,KE)
      REAL ::   ZDT   (IE,JE,KE)
      REAL ::   ZDP   (IE,JE,KE)
      REAL ::   ZDUTR (IE,JE,NVM)
      REAL ::   ZDVTR (IE,JE,NVM)
      REAL ::   ZDPTR (IE,JE,NVM)
      REAL ::   ZDDIV (IE,JE)
      REAL ::   ZDZET (IE,JE)
      REAL ::   ZDPSI (IE,JE)
      REAL ::   ZDCHI (IE,JE)
      REAL ::   CPDLDDP(JE)
      REAL ::   ACPDLR (JE)
      REAL ::   ZACDLR (JE)
      REAL ::   ACPHDLQ(JE)
      REAL ::   ACLDCQ (JE,NVM)
      REAL ::   ZCOEFP (IE,MOJE)
      REAL ::   ZCOEFH (IE,MOJE,NVM)
      REAL ::   ZBPTP  (IE,MOJE)
      REAL ::   ZBPTH  (IE,MOJE,NVM)
      REAL ::   ZAP    (MOJE)
      REAL ::   ZBBP   (MOJE)
      REAL ::   ZBBH   (MOJE,NVM)
      REAL ::   ZCP    (MOJE)
      REAL ::   ZALF   (IE,JE)
      REAL ::   ZCLF   (IE,JE)
      REAL ::   ZAPF   (IE,JE)
      REAL ::   ZBBF   (IE,JE)
      REAL ::   ZCPF   (IE,JE)
      REAL ::   ZDI    (MOIE)
      REAL ::   ZFCI   (IE,JE)
      REAL ::   ZFCJ   (IE,JE)
      REAL ::   ZRES   (IE,JE)
      REAL ::   ZRS2   (IE,JE)
      REAL ::   ZXNUE  (IE,JE)
      REAL ::   ZDPG   (IE,JE)
      REAL ::   ZDELPSI(IE,JE)
      REAL ::   ZDELCHI(IE,JE)
      REAL ::   ZDPSIG (IE,JE)
      REAL ::   ZDUG   (IE,JE)
      REAL ::   ZDVG   (IE,JE)
      REAL ::   ZW1I   (MOIE,JE)
      REAL ::   ZW2I   (MOIE,JE)
      REAL ::   ZDUMAX (KE)
      REAL ::   ZDVMAX (KE)
      REAL ::   ZDTMAX (KE)
      REAL ::   ZDUMIT (KE)
      REAL ::   ZDVMIT (KE)
      REAL ::   ZDTMIT (KE)
      REAL ::   ZF_FFT(MOIE,JE)
      REAL ::   ZX_GAUSS(MOIE,JE)
      REAL ::   GACLDCQ(MOJE,NVM)
      REAL ::   GACPHDLQ(MOJE)
C
      REAL ::   EDDT,ZFCQQ,ZFKPS,ZFMAX,ZFMIN,ZZADPHR,ZZDT,ZEPSASS
      REAL ::   ZDPSMIT,ZDPSMAX,ZDLDDPQ,ZDLAMRQ,ZCPLIJ,ZBBPLIJ
      REAL ::   ZBAL,ZAPLIJ,ZADPHIR
C
      INTEGER :: I,J,K,IAU1,IAV1,IEJEKV
      INTEGER :: MGAUSS,NZDR,NIT0,NIT,NFFT,MT,K2,K1,JEV1,JEU1
      INTEGER :: JAV1,JAU1,IJ,IEV1,IEU1
CKS
C     FOR HALO EXCHANGE
      INTEGER :: KSIZE(10), ITYPE
      
C-----------------------------------------------------------------------

C     VARIABLE FUER SPEZIELLEN FALL SETZEN, ORIGINALWERTE AUFHEBEN

      ZZDT     = DT
      ZEPSASS  = EPSASS
      LZSITS   = LSITS
      LZPHY    = LPHY
      LZWRITD  = LWRITED
      LZWRITE  = LWRITEE
      LZWRITF  = LWRITEF
      LZWRITT  = LWRITET

      NZDR     = NDR
      NDR      = NINT (NDR * DT / DTNMI)
      FAKRMY   = DTNMI / DT
      IF (NZT .EQ. 0)  FAKRMY = 0.5 * FAKRMY

      CALL GETALPHA(RMY, FAKRMY, ALPHABOUND)

      NA  = 1
      NJ  = 2
      NE  = 3
      NA2 = 1
      NJ2 = 2

      DT       = DTNMI
      IF (NZT .EQ. 0)  DT =0.5 * DT
      EPSASS   = 0.
      LSITS    = .FALSE.
      IF (LANMI)  THEN
          LPHY     = .FALSE.
      ELSE
          LPHY     = .TRUE.
      ENDIF
      LWRITED  = .FALSE.
      LWRITEE  = .FALSE.
      LWRITEF  = .FALSE.
      LWRITET  = .FALSE.

C     KONTROLLVARIABLE VORBESETZEN
      ZDPSMAX  = 0.0
      ZDPSMIT  = 0.0

C-----------------------------------------------------------------------

      CALL COPYRE(PS(1,1,NJ)  ,ZDPS,IEJE)
      CALL COPYRE(U (1,1,1,NJ),ZDU ,IEJEKE)
      CALL COPYRE(V (1,1,1,NJ),ZDV ,IEJEKE)
      CALL COPYRE(T (1,1,1,NJ),ZDT ,IEJEKE)
      CALL COPYRE(QD(1,1,1,NJ),ZDP ,IEJEKE)

C
C     KONTROLLVARIABLE VORBESETZEN
      ZDPSMAX   = 0.0
      ZDPSMIT   = 0.0
      DO K  = 1,KE
         ZDUMAX(K) = 0.0
         ZDVMAX(K) = 0.0
         ZDTMAX(K) = 0.0
         ZDUMIT(K) = 0.0
         ZDVMIT(K) = 0.0
         ZDTMIT(K) = 0.0
      ENDDO

C-----------------------------------------------------------------------

      IEJEKV = IE * JE * NVM
      MT     = 5 * (MOIE-1) / 2
      NFFT   = MYFFT_JUP - MYFFT_JLO + 3
      MGAUSS = MYGAUSS_IUP - MYGAUSS_ILO + 3

      IF (NEIGHBOR(1) .EQ. -1) THEN
         IAU1 = IAU + 1
         IAV1 = IAV + 1
      ELSE
         IAU1 = IAU
         IAV1 = IAV
      ENDIF

      IF (NEIGHBOR(3) .EQ. -1) THEN
         IEU1 = IEU - 1
         IEV1 = IEV - 1
      ELSE
         IEU1 = IEU
         IEV1 = IEV
      ENDIF

      IF (NEIGHBOR(4) .EQ. -1) THEN
         JAU1 = JAU + 1
         JAV1 = JAV + 1
      ELSE
         JAU1 = JAU
         JAV1 = JAV
      ENDIF

      IF (NEIGHBOR(2) .EQ. -1) THEN
         JEU1 = JEU - 1
         JEV1 = JEV - 1
      ELSE
         JEU1 = JEU
         JEV1 = JEV
      ENDIF

C     BERECHNUNG DER INVERSEN MATRIZEN VON SIGAM UND SISTM

      DO K2 = 1 , KE
!DIR$ SHORTLOOP
         DO K1 = 1 , KE
            SIGAMI(K1,K2) = SIGAM(K1,K2)
            SISTMI(K1,K2) = SISTM(K1,K2)
            SIGAMI(K1,K2+KE) = 0.
            SISTMI(K1,K2+KE) = 0.
         ENDDO
         SIGAMI(K2,K2+KE) = 1.
         SISTMI(K2,K2+KE) = 1.
      ENDDO

      CALL MINVLAP(SIGAMI,KE)
      CALL MINVLAP(SISTMI,KE)

C    VORBEREITUNG DER SCHLEIFE UEBER DIE VERTIKALEN MODES

      ZDLDDPQ = DLADDPH ** 2
      ZDLAMRQ = (RERD / EDDLAM) **2
      ZADPHIR = EDDPHI / RERD
      ZZADPHR = 0.5 * ZADPHIR

      DO J = JAHGG , JEHGG
         ZAP(J)     = ZDLDDPQ * GCPHI(J,1) * GCPHI(J-1,2)
         ZCP(J)     = ZDLDDPQ * GCPHI(J,1) * GCPHI(J  ,2)
         ZBBP(J)    = -(2. + ZAP(J) + ZCP(J))
         GACPHDLQ(J)= ZDLAMRQ * (GCPHI(J,1) **2)
      ENDDO
C
C     ** ACPHDLQ ** WURDE VORHER IN LOOP 50 BELEGT
C
      DO J = JAA, JEA
         ACPDLR(J)  = EDDLAM * ACPHIR(J,1)
         ZACDLR(J)  = 0.5 * ACPDLR(J)
         ACPHDLQ(J) = ZDLAMRQ * (CPHI(J,1) **2)
         CPDLDDP(J) = DLADDPH * CPHI(J,2)
      ENDDO

      CALL GAELKO (ZCOEFP,ZBPTP,MYGAUSS_ILO,MYGAUSS_IUP,
     &             MYGAUSS_JLO,MYGAUSS_JUP,MOIE,MOJE,
     &             ZAP,ZBBP,ZCP,ZDI,PI)

      ZFMIN =  1.0E+10
      ZFMAX = -1.0E+10
      DO J = JAH , JEH
         DO I = IAH , IEH
            ZFMIN = AMIN1 (ZFMIN, FC(I,J))
            ZFMAX = AMAX1 (ZFMAX, FC(I,J))
         ENDDO
      ENDDO

      COUNT = 1
      CALL PALLREDUCER(ZFMIN,MPI_MIN)
      CALL PALLREDUCER(ZFMAX,MPI_MAX)

      ZFCQQ = 0.5 * (ZFMIN **2 + ZFMAX **2)

      DO K = 1 , NVM
         DO J = JAHGG , JEHGG
            GACLDCQ(J,K) = GACPHDLQ(J) / SICQ(K)
            ZBBH(J,K)   = ZBBP(J) - ZFCQQ * GACLDCQ(J,K)
         ENDDO
         CALL GAELKO (ZCOEFH(1,1,K),ZBPTH(1,1,K),MYGAUSS_ILO,
     &        MYGAUSS_IUP,MYGAUSS_JLO,MYGAUSS_JUP,MOIE,MOJE,
     &        ZAP,ZBBH(1,K),ZCP,ZDI,PI)
      ENDDO

      DO K = 1 , NVM
         DO J = JAA , JEA
            ACLDCQ(J,K) = ACPHDLQ(J) / SICQ(K)
         ENDDO
      ENDDO

      DO J = JAA , JEA
         DO I = IAA , IEH
            ZFCI(I,J) = 0.5 * (FC(I,J) + FC(I+1,J))
         ENDDO
      ENDDO

      DO J = JAA , JEH
         DO I = IAA , IEA
            ZFCJ(I,J) = 0.5 * (FC(I,J) + FC(I,J+1))
         ENDDO
      ENDDO

      DO J = JAH , JEH
         ZAPLIJ = ZDLDDPQ * CPHI(J,1) * CPHI(J-1,2)
         ZCPLIJ = ZDLDDPQ * CPHI(J,1) * CPHI(J  ,2)
         DO I = IAH , IEH
            ZALF(I,J) = ZFCI(I-1,J)
            ZCLF(I,J) = ZFCI(I  ,J)
            ZAPF(I,J) = ZFCJ(I,J-1) * ZAPLIJ
            ZCPF(I,J) = ZFCJ(I,J  ) * ZCPLIJ
            ZBBF(I,J) = -(ZALF(I,J) + ZCLF(I,J) + ZAPF(I,J) + ZCPF(I,J))
         ENDDO
      ENDDO


C-----------------------------------------------------------------------

C     ITERATIVE WIEDERHOLUNG DER NORMAL MODE INITIALISIERUNG
C     ------------------------------------------------------
C     ------------------------------------------------------

C     UEBERSCHRIFT FUER KONTROLLAUSDRUCK
      IF (MYID .EQ. 0) THEN
         IF (LDIA) WRITE(NUAUFTR, 9900)
      ENDIF
 9900 FORMAT('1     KONTROLLAUSGABE DER NICHTLINEAREN NORMAL MODE',
     &       ' INITIALISIERUNG NACH TEMPERTON/WERGEN',/,
     &       '0     ITERATION   MODE-NR.    BAL',/)
C
      DO IJ = 1,IEJE
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
C
      DO NIT = 1 , NITNMI
C
C     PROGORG 1 ZEITSCHRITT LAUFEN LASSEN (NUR EXPLIZITE PROGNOSE)

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
     &   PINT   , DWDT  , ETAS  , W      , RPRAC)

      ITYPE = 20000
      KSIZE(1:10) = (/KE,KE,KE,1,KE,KE,KE,0,0,0/)
      CALL HALOEXCH(KSIZE, ITYPE, U(:,:,:,NE), V(:,:,:,NE), T(:,:,:,NE),
     &     PS(:,:,NE), QD(:,:,:,NE), QW(:,:,:,NE), QI(:,:,:,NE))

C
C     ZEITLICHE AENDERUNG VON U, V, PS, T
C     -----------------------------------

      EDDT = 1. / DTNMI
      DO J = JAA , JEA
         DO I = IAA , IEA
            ZDPS(I,J) = EDDT * (PS(I,J,NE) - PS(I,J,NJ))
         ENDDO
      ENDDO

      DO K = 1 , KE
         DO J = JAA , JEA
            DO I = IAA , IEA
               ZDU(I,J,K) = EDDT * (U(I,J,K,NE) - U(I,J,K,NJ))
               ZDV(I,J,K) = EDDT * (V(I,J,K,NE) - V(I,J,K,NJ))
               ZDT(I,J,K) = EDDT * (T(I,J,K,NE) - T(I,J,K,NJ))
            ENDDO
         ENDDO
      ENDDO

C     ZEITLICHE AENDERUNG VON P
C     -------------------------
C     MATRIZENMULTIPLIKATION ZDT * SIGAM UNTER BERUECKSICHTIGUNG
C     DER FORM VON SIGAM

      ZFKPS = R * SITR / SIPSR
      CALL SETRA (ZDUM2,IEJE,0.)

      DO J = JAA , JEA
         DO I = IAA , IEA
            ZDUM3(I,J) = ZDPS(I,J)  * ZFKPS
         ENDDO
      ENDDO

      DO K = KE , 1 , -1
         DO J = JAA , JEA
            DO I = IAA , IEA
               ZDP(I,J,K) = ZDUM2(I,J) + ZDT(I,J,K) * SIGAM(K,K)
               ZDP(I,J,K) = ZDP(I,J,K) + ZDUM3(I,J)
               ZDUM2(I,J) = ZDUM2(I,J) + ZDT(I,J,K) * SIGAM(K,1)
            ENDDO
         ENDDO
      ENDDO

C     TRANSFORMATION IN DEN NORMAL MODE RAUM
C     --------------------------------------

      CALL SETRA (ZDUTR,IEJEKV,0.)
      CALL SETRA (ZDVTR,IEJEKV,0.)
      CALL SETRA (ZDPTR,IEJEKV,0.)
      DO K2 = 1 , NVM
         DO K1 = 1 , KE
            DO J = JAA , JEA
               DO I = IAA , IEA
      ZDUTR(I,J,K2) = ZDUTR(I,J,K2) + ZDU(I,J,K1) * SIVMT(K1,K2)
      ZDVTR(I,J,K2) = ZDVTR(I,J,K2) + ZDV(I,J,K1) * SIVMT(K1,K2)
      ZDPTR(I,J,K2) = ZDPTR(I,J,K2) + ZDP(I,J,K1) * SIVMT(K1,K2)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C     AUSGABE VON BAL (ANFANGSFELDER) ALS KONTROLLGROESSE
C     ---------------------------------------------------

      IF (NIT .EQ. 1) THEN
         DO K = 1 , NVM
            ZBAL = 0.
            DO J = JAH , JEH
               DO I = IAH , IEH
                  ZBAL = ZBAL + ZDPTR(I,J,K) **2
     &                 + SICQ(K) * (ZDUTR(I,J,K) **2 + ZDVTR(I,J,K) **2)
               ENDDO
            ENDDO

            COUNT=1
            CALL PREDUCERR(ZBAL,MPI_SUM)

            NIT0 = 0
            IF (MYID .EQ. 0) THEN
               IF (LDIA) WRITE(NUAUFTR, 9901) NIT0, K, ZBAL
            ENDIF

 9901       FORMAT('0', I10, I11, E14.5)

         ENDDO

C     SCHLEIFE UEBER DIE VERTIKALEN MODES
C     -----------------------------------
C     -----------------------------------

      ENDIF

      CALL SETRA (ZDDIV  ,IEJE,0.)
      CALL SETRA (ZDZET  ,IEJE,0.)
      CALL SETRA (ZRES   ,IEJE,0.)
      CALL SETRA (ZRS2   ,IEJE,0.)
      CALL SETRA (ZDPSI  ,IEJE,0.)
      CALL SETRA (ZDCHI  ,IEJE,0.)
      CALL SETRA (ZXNUE  ,IEJE,0.)
      CALL SETRA (ZDELPSI,IEJE,0.)
      CALL SETRA (ZDELCHI,IEJE,0.)
      CALL SETRA (ZDPSIG ,IEJE,0.)
      CALL SETRA (ZDUG   ,IEJE,0.)
      CALL SETRA (ZDVG   ,IEJE,0.)

      DO K = 1 , NVM

C     BERECHNUNG DER ZEITLICHEN AENDERUNG VON DIVERGENZ UND VORTICITY
C     ---------------------------------------------------------------

         DO J = JAH , JEH
            DO I = IAH , IEH
         ZDDIV(I,J) = ACPHIR(J,1) *
     &   (EDDLAM *(          ZDUTR(I,J,K) -             ZDUTR(I-1,J,K))
     &  + EDDPHI *(CPHI(J,2)*ZDVTR(I,J,K) - CPHI(J-1,2)*ZDVTR(I,J-1,K)))
            ENDDO
         ENDDO

         DO J = JAA , JEH
            DO I = IAA , IEH
         ZDUM2(I,J) = ACPHIR(J,2) *
     &   (EDDLAM *(            ZDVTR(I+1,J,K) -           ZDVTR(I,J,K))
     &  - EDDPHI *(CPHI(J+1,1)*ZDUTR(I,J+1,K) - CPHI(J,1)*ZDUTR(I,J,K)))
            ENDDO
         ENDDO

         DO J = JAH , JEH
            DO I = IAH , IEH
               ZDZET(I,J) = 0.25 * (ZDUM2(I-1,J-1) + ZDUM2(I,J-1)
     &                            + ZDUM2(I-1,J)   + ZDUM2(I,J))
            ENDDO
         ENDDO

C     LOESUNG EINER POISSONGLEICHUNG ZUR BERECHNUNG DER ZEITLICHEN
C     AENDERUNG VON STROMFUNKTION UND GESCHWINDIGKEITSPOTENTIAL
C     ------------------------------------------------------------

         DO J = JAA , JEA
            DO I = IAA , IEA
               ZDUM2(I,J) = ACPHDLQ(J) * ZDZET(I,J)
               ZDUM3(I,J) = ACPHDLQ(J) * ZDDIV(I,J)
            ENDDO
         ENDDO

         CALL SOLVE (ZDPSI,ZDUM2,ZF_FFT,ZW1I,ZW2I,ZX_GAUSS,ZCOEFP,ZBPTP,
     &        ZAP,TRIGSI,RZ1I,RZ2I,IFAXI,NFFT,MGAUSS,MT)

         ITYPE = 10000
         KSIZE(1:10) = (/1,0,0,0,0,0,0,0,0,0/)
         CALL HALOEXCH(KSIZE, ITYPE, ZDPSI)

         CALL SOLVE (ZDCHI,ZDUM3,ZF_FFT,ZW1I,ZW2I,ZX_GAUSS,ZCOEFP,ZBPTP,
     &        ZAP,TRIGSI,RZ1I,RZ2I,IFAXI,NFFT,MGAUSS,MT)

         ITYPE = 10000
         KSIZE(1:10) = (/1,0,0,0,0,0,0,0,0,0/)
         CALL HALOEXCH(KSIZE, ITYPE, ZDCHI)

C     BERECHNUNG DER RECHTEN SEITE DER GLEICHUNG FUER
C     DIE ZEITLICHE AENDERUNG VON PG
C     -----------------------------------------------

         DO J = JAH , JEH
            ZAPLIJ    = ZDLDDPQ*CPHI(J,1)*CPHI(J-1,2)
            ZCPLIJ    = ZDLDDPQ*CPHI(J,1)*CPHI(J  ,2)
            ZBBPLIJ   = - (2.0 + ZAPLIJ + ZCPLIJ)
            DO I = IAH , IEH
      ZRES(I,J) = -(ZALF(I,J)*ZDPSI(I-1,J) + ZCLF(I,J)*ZDPSI(I+1,J)
     &            + ZAPF(I,J)*ZDPSI(I,J-1) + ZCPF(I,J)*ZDPSI(I,J+1)
     &            + ZBBF(I,J)*ZDPSI(I,J))
     &            +           ZDPTR(I-1,J,K) +           ZDPTR(I+1,J,K)
     &            + ZAPLIJ  * ZDPTR(I,J-1,K) + ZCPLIJ  * ZDPTR(I,J+1,K)
     &            + ZBBPLIJ * ZDPTR(I,J,K)
            ENDDO
         ENDDO

C     BERECHNUNG DER ZEITLICHEN AENDERUNG VON PG
C     ------------------------------------------

         CALL SETRA (ZDPG,IEJE,0.)

         CALL CONGRA(ZDPG,ZRES,ZRS2,ZXNUE,
     &        ZAP,ZBBP,ZCOEFP,ZBPTP,
     &        ZBBH(1:MOJE,K),ZCOEFH(1:IE,1:MOJE,K),ZBPTH(1:IE,1:MOJE,K),
     &        ACLDCQ(1:JE,K),
     &        ZALF,ZCLF,ZBBF,ZAPF,ZCPF,ZW1I,ZW2I,
     &        ZF_FFT,ZX_GAUSS,TRIGSI,RZ1I,RZ2I,IFAXI,
     &        NFFT,MGAUSS,MT,MYPOSGRD(2))


C     LOESUNG DER POISSONGLEICHUNG ZUR BERECHNUNG DER AENDERUNG VON CHI
C     ------------------------------------------------------------------

         DO J = JAA , JEA
            DO I = IAA , IEA
               ZDUM2(I,J) = ACLDCQ(J,K) * ZDPG(I,J)
            ENDDO
         ENDDO

         CALL SOLVE (ZDELCHI,ZDUM2,ZF_FFT,ZW1I,ZW2I,ZX_GAUSS,ZCOEFP,
     &        ZBPTP,ZAP,TRIGSI,RZ1I,RZ2I,IFAXI,NFFT,MGAUSS,MT)

         ITYPE = 10000
         KSIZE(1:10) = (/1,0,0,0,0,0,0,0,0,0/)
         CALL HALOEXCH(KSIZE, ITYPE, ZDELCHI)

C     BERECHNUNG DER RECHTEN SEITE DER GLEICHUNG FUER AENDERUNG VON PTR
C     -----------------------------------------------------------------

         DO J = JAA , JEH
            DO I = IAA , IEH
      ZDUM2(I,J) = CPDLDDP(J) *
     &            (ZDELCHI(I  ,J  ) * (-ZFCI(I,J  ) + ZFCJ(I  ,J))
     &           + ZDELCHI(I+1,J+1) * (-ZFCI(I,J+1) + ZFCJ(I+1,J))
     &           + ZDELCHI(I+1,J  ) * ( ZFCI(I,J  ) - ZFCJ(I+1,J))
     &           + ZDELCHI(I  ,J+1) * ( ZFCI(I,J+1) - ZFCJ(I  ,J)))
            ENDDO
         ENDDO

         DO J = JAH , JEH
            DO I = IAH , IEH
               ZRES(I,J) = 0.25 * (ZDUM2(I-1,J-1) + ZDUM2(I,J-1)
     &                           + ZDUM2(I-1,J)   + ZDUM2(I,J))
     &                           + ACPHDLQ(J) * ZDDIV(I,J)
            ENDDO
         ENDDO

C     BERECHNUNG DER AENDERUNG VON PTR
C     --------------------------------

         CALL SETRA (ZDPTR(1,1,K),IEJE,0.)
         CALL CONGRA(ZDPTR(1,1,K),ZRES,ZRS2,ZXNUE,
     &        ZAP,ZBBP,ZCOEFP,ZBPTP,
     &        ZBBH(1,K),ZCOEFH(1,1,K),ZBPTH(1,1,K),ACLDCQ(1,K),
     &        ZALF,ZCLF,ZBBF,ZAPF,ZCPF,ZW1I,ZW2I,
     &        ZF_FFT,ZX_GAUSS,TRIGSI,RZ1I,RZ2I,IFAXI,
     &        NFFT,MGAUSS,MT,MYPOSGRD(2))


C     LOESUNG VON 2 POISSONGL. ZUR BERECHNUNG DER AENDERUNG VON PSI
C     -------------------------------------------------------------

         CALL SETRA (ZDUM2,IEJE,0.)
         CALL SOLVE (ZXNUE,ZDPTR(1,1,K),ZF_FFT,ZW1I,ZW2I,ZX_GAUSS,
     &        ZCOEFP,ZBPTP,ZAP,TRIGSI,RZ1I,RZ2I,IFAXI,NFFT,MGAUSS,MT)

         ITYPE = 10000
         KSIZE(1:10) = (/1,0,0,0,0,0,0,0,0,0/)
         CALL HALOEXCH(KSIZE, ITYPE, ZXNUE)

         DO J = JAH , JEH
            DO I = IAH , IEH
      ZDUM2(I,J) = ACLDCQ(J,K) *
     &           (ZALF(I,J) * ZXNUE(I-1,J) + ZCLF(I,J) * ZXNUE(I+1,J)
     &          + ZAPF(I,J) * ZXNUE(I,J-1) + ZCPF(I,J) * ZXNUE(I,J+1)
     &          + ZBBF(I,J) * ZXNUE(I,J))
            ENDDO
         ENDDO

         CALL SOLVE (ZDELPSI,ZDUM2,ZF_FFT,ZW1I,ZW2I,ZX_GAUSS,ZCOEFP,
     &        ZBPTP,ZAP,TRIGSI,RZ1I,RZ2I,IFAXI,NFFT,MGAUSS,MT)

         ITYPE = 10000
         KSIZE(1:10) = (/1,0,0,0,0,0,0,0,0,0/)
         CALL HALOEXCH(KSIZE, ITYPE, ZDELPSI)

C     AUSGABE VON BAL ALS KONTROLLGROESSE
C     -----------------------------------

         DO J = JAH , JEH
            DO I = IAH , IEH
      ZDUM2(I,J) = ZALF(I,J)*ZDELCHI(I-1,J) + ZCLF(I,J)*ZDELCHI(I+1,J)
     &           + ZAPF(I,J)*ZDELCHI(I,J-1) + ZCPF(I,J)*ZDELCHI(I,J+1)
     &           + ZBBF(I,J)*ZDELCHI(I,J)
            ENDDO
         ENDDO

         CALL SOLVE (ZDPSIG,ZDUM2,ZF_FFT,ZW1I,ZW2I,ZX_GAUSS,ZCOEFP,
     &        ZBPTP,ZAP,TRIGSI,RZ1I,RZ2I,IFAXI,NFFT,MGAUSS,MT)

         ITYPE = 10000
         KSIZE(1:10) = (/1,0,0,0,0,0,0,0,0,0/)
         CALL HALOEXCH(KSIZE, ITYPE, ZDPSIG)

         DO J = JAU1 , JEU1
            DO I = IAU1 , IEU1
               ZDUG(I,J) =  -ZZADPHR * (ZDPSIG(I,J+1) - ZDPSIG(I,J-1))
     &                   + ZACDLR(J) * (ZDCHI(I+1,J)  - ZDCHI(I-1,J))
            ENDDO
         ENDDO

         DO J = JAV1 , JEV1
            DO I = IAV1 , IEV1
               ZDVG(I,J) = ZACDLR(J) * (ZDPSIG(I+1,J) - ZDPSIG(I-1,J))
     &                     + ZZADPHR * (ZDCHI(I,J+1)  - ZDCHI(I,J-1))
            ENDDO
         ENDDO

         ZBAL = 0.
         DO J = JAH , JEH
            DO I = IAH , IEH
               ZBAL = ZBAL + ZDPG(I,J) **2
     &              + SICQ(K) * (ZDUG(I,J) **2 + ZDVG(I,J) **2)
            ENDDO
         ENDDO

         COUNT=1
         CALL PREDUCERR(ZBAL,MPI_SUM)

         IF (MYID .EQ. 0) THEN
            IF (LDIA) WRITE(NUAUFTR, 9901) NIT, K, ZBAL
         ENDIF


C     UMRECHNUNG VON DELPSI UND DELCHI IN DELUTR UND DELVTR
C     -----------------------------------------------------

         CALL SETRA (ZDUTR(1,1,K),IEJE,0.)
         CALL SETRA (ZDVTR(1,1,K),IEJE,0.)
         DO J = JAA , JEH
            DO I = IAH , IEH + 1
               ZDUM2(I,J) = ZADPHIR * (ZDELPSI(I,J+1) - ZDELPSI(I,J))
            ENDDO
         ENDDO

         DO J = JAU , JEU
            DO I = IAU , IEU
               ZDUTR(I,J,K) = -0.25 * (ZDUM2(I,J-1) + ZDUM2(I+1,J-1)
     &                               + ZDUM2(I,J)   + ZDUM2(I+1,J))
     &                  + ACPDLR(J) * (ZDELCHI(I+1,J) - ZDELCHI(I,J))
            ENDDO
         ENDDO

         DO J = JAH , JEH + 1
            DO I = IAA , IEH
               ZDUM2(I,J) = ACPDLR(J) * (ZDELPSI(I+1,J) - ZDELPSI(I,J))
            ENDDO
         ENDDO

         DO J = JAV , JEV
            DO I = IAV , IEV
               ZDVTR(I,J,K) =  0.25 * (ZDUM2(I-1,J  )   + ZDUM2(I,J  )
     &                               + ZDUM2(I-1,J+1)   + ZDUM2(I,J+1))
     &                    + ZADPHIR * (ZDELCHI(I,J+1) - ZDELCHI(I,J))
            ENDDO
         ENDDO
      ENDDO !K = 1 , NVM

C     RUECKTRANSFORMATION IN DEN PHYSIKALISCHEN RAUM
C     ----------------------------------------------

      CALL SETRA (ZDU,IEJEKE,0.)
      CALL SETRA (ZDV,IEJEKE,0.)
      CALL SETRA (ZDP,IEJEKE,0.)
      DO K2 = 1 , KE
         DO K1 = 1 , NVM
            DO J = JAH , JEH
               DO I = IAH , IEH
      ZDU(I,J,K2) = ZDU(I,J,K2) + ZDUTR(I,J,K1) * SIVMTI(K1,K2)
      ZDV(I,J,K2) = ZDV(I,J,K2) + ZDVTR(I,J,K1) * SIVMTI(K1,K2)
      ZDP(I,J,K2) = ZDP(I,J,K2) + ZDPTR(I,J,K1) * SIVMTI(K1,K2)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

C     AENDERUNG DER TEMPERATUR UND DES BODENDRUCKS
C     --------------------------------------------

      CALL SETRA (ZDPS,IEJE,0.)
      DO K2 = 1 , KE
         CALL SETRA (ZDUM2,IEJE,0.)
         DO K1 = 1 , KE
            DO J = JAH , JEH
               DO I = IAH , IEH
                  ZDUM2(I,J) = ZDUM2(I,J) + ZDP(I,J,K1) * SISTMI(K1,K2)
               ENDDO
            ENDDO
         ENDDO

         DO J = JAH , JEH
            DO I = IAH , IEH
               ZDPS(I,J) = ZDPS(I,J) + SINUE(K2) * ZDUM2(I,J)
            ENDDO
         ENDDO
      ENDDO

C     MATRIZENMULTIPLIKATION (ZDP - ..) * SIGAMI UNTER BERUECKSICHTIGUNG
C     DER FORM VON SIGAMI (OBERE DREIECKSMATRIX)

      DO K = 1 , KE
         DO J = JAH , JEH
            DO I = IAH , IEH
               ZDP(I,J,K) = ZDP(I,J,K) - ZFKPS * ZDPS(I,J)
            ENDDO
         ENDDO
      ENDDO

      CALL SETRA (ZDT,IEJEKE,0.)
      DO K2 = 1 , KE
      DO K1 = K2, KE
      DO J = JAH , JEH
      DO I = IAH , IEH
         ZDT(I,J,K2) = ZDT(I,J,K2) + ZDP(I,J,K1) * SIGAMI(K1,K2)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C     AENDERUNGEN VON U, V, PS, T SUMMIEREN FUER KONTROLLAUSGABE
C     ----------------------------------------------------------

      IF(NIT.EQ.NITNMI) THEN
          DO J = JAH , JEH
             DO I = IAH , IEH
                ZDPSMAX   = MAX(ZDPSMAX, ABS(ZDPS(I,J)))
                ZDPSMIT   =     ZDPSMIT  +   ZDPS(I,J)
             ENDDO
          ENDDO

          COUNT=1
          CALL PREDUCERR(ZDPSMIT,MPI_SUM)
          CALL PREDUCERR(ZDPSMAX,MPI_MAX)
          IF (MYID .EQ. 0) THEN
             ZDPSMIT   = ZDPSMIT/FLOAT((JEHGG-JAHGG+1)*(IEHGG-IAHGG+1))
          ENDIF

      ENDIF
      DO K =   1 , KE
         DO J = JAH , JEH
            DO I = IAH , IEH
               ZDUMAX(K) = MAX(ZDUMAX(K), ABS(ZDU(I,J,K)))
               ZDUMIT(K) =     ZDUMIT(K)  +   ZDU(I,J,K)
               ZDVMAX(K) = MAX(ZDVMAX(K), ABS(ZDV(I,J,K)))
               ZDVMIT(K) =     ZDVMIT(K)  +   ZDV(I,J,K)
               ZDTMAX(K) = MAX(ZDTMAX(K), ABS(ZDT(I,J,K)))
               ZDTMIT(K) =     ZDTMIT(K)  +   ZDT(I,J,K)
            ENDDO
         ENDDO
      ENDDO

      COUNT=KE
      CALL PREDUCER(ZDUMIT(1),MPI_SUM)
      CALL PREDUCER(ZDVMIT(1),MPI_SUM)
      CALL PREDUCER(ZDTMIT(1),MPI_SUM)
      CALL PREDUCER(ZDUMAX(1),MPI_MAX)
      CALL PREDUCER(ZDVMAX(1),MPI_MAX)
      CALL PREDUCER(ZDTMAX(1),MPI_MAX)

      IF (MYID .EQ. 0) THEN
      DO K=1,KE
         ZDUMIT(K) = ZDUMIT(K)/FLOAT((JEHGG-JAHGG+1)*(IEHGG-IAHGG+1))
         ZDVMIT(K) = ZDVMIT(K)/FLOAT((JEHGG-JAHGG+1)*(IEHGG-IAHGG+1))
         ZDTMIT(K) = ZDTMIT(K)/FLOAT((JEHGG-JAHGG+1)*(IEHGG-IAHGG+1))
      ENDDO

C     AENDERUNGEN VON U, V, PS, T AUSDRUCKEN IN KONTROLLAUSGABE
C     ---------------------------------------------------------

      IF (LDIA) THEN
      IF(NIT.EQ.NITNMI) THEN
          WRITE(NUAUFTR, 9902)
 9902     FORMAT(//,'0     AENDERUNG DER ANFANGSFELDER DURCH DIE INMI:',
     &           /)
          WRITE(NUAUFTR, 9903) ZDPSMIT*0.01, ZDPSMAX*0.01
 9903     FORMAT('0     PS(HPA):  MITTEL = ',F9.6,'  ABS(MAXIMUM) = ',
     &           F9.5,/,
     &           '0     K     U(M/S):    MITTEL       ABS(MAXIMUM)',
     &                 '      V(M/S):    MITTEL       ABS(MAXIMUM)',
     &                 '      T(K):      MITTEL       ABS(MAXIMUM)',/)
          DO K = 1,KE
             WRITE(NUAUFTR, 9904) K, ZDUMIT(K), ZDUMAX(K),
     &            ZDVMIT(K), ZDVMAX(K),
     &            ZDTMIT(K), ZDTMAX(K)
 9904        FORMAT(' ',I6,F22.5,F16.5,F26.5,F16.5,F26.5,F16.5)
          ENDDO
      ENDIF
      ENDIF
      ENDIF

C     AENDERUNGEN VON U, V, PS, T ADDIEREN (FUER ZEITPUNKTE -1 UND 0)
C     ------------------------------------

      IF (NIT .EQ. NITNMI) THEN
        DO J = JAH , JEH
          DO I = IAH , IEH
            PS(I,J,NA) = PS(I,J,NA) + ZDPS(I,J)
          ENDDO
        ENDDO
      ENDIF

      DO K = 1 , KE
        DO J = JAH , JEH
          DO I = IAH , IEH
            U (I,J,K,NA) = U (I,J,K,NA) + ZDU(I,J,K)
            V (I,J,K,NA) = V (I,J,K,NA) + ZDV(I,J,K)
            T (I,J,K,NA) = T (I,J,K,NA) + ZDT(I,J,K)
          ENDDO
        ENDDO
      ENDDO

      ITYPE = 20000
      KSIZE(1:10) = (/KE,KE,KE,1,KE,KE,KE,0,0,0/)
      CALL HALOEXCH(KSIZE, ITYPE, U(:,:,:,NA), V(:,:,:,NA), T(:,:,:,NA),
     &     PS(:,:,NA), QD(:,:,:,NA), QW(:,:,:,NA), QI(:,:,:,NA))

C     NJ-WERT DURCH NA-WERT IM GESAMTGEBIET ERSETZEN
      DO J = JAA , JEA
         DO I = IAA , IEA
            PS(I,J,NJ) = PS(I,J,NA)
         ENDDO
      ENDDO

      DO K = 1 , KE
        DO J = JAA , JEA
          DO I = IAA , IEA
            U (I,J,K,NJ) = U (I,J,K,NA)
            V (I,J,K,NJ) = V (I,J,K,NA)
            T (I,J,K,NJ) = T (I,J,K,NA)
            QD(I,J,K,NJ) = QD(I,J,K,NA)
            QW(I,J,K,NJ) = QW(I,J,K,NA)
            QI(I,J,K,NJ) = QI(I,J,K,NA)
          ENDDO
        ENDDO
      ENDDO

      DO K = 1 , KE+1
        DO J = JAA , JEA
          DO I = IAA , IEA
            PINT(I,J,K,NA) = AK(K) + BK(K)*PS(I,J,NA)
            PINT(I,J,K,NJ) = PINT(I,J,K,NA)
            W   (I,J,K,NA) = 0.
            W   (I,J,K,NJ) = 0.
          ENDDO
        ENDDO
      ENDDO

C     NEUES GEOPOTENTIAL BERECHNEN
C     ----------------------------

      CALL GEOPOT(NA, NA2, 1, IE, 1, JE,
     &            T  , QD, QW  , FI    , FIB ,
     &            QI , PINT    , DWDT)
      CALL GEOPOT(NJ, NJ2, 1, IE, 1, JE,
     &            T  , QD, QW  , FI    , FIB ,
     &            QI , PINT    , DWDT)

      ENDDO !NIT = 1 , NITNMI

C     VARIABLE AUF ORIGINALWERT ZURUECKSETZEN

      DT       = ZZDT
      EPSASS   = ZEPSASS
      LSITS    = LZSITS
      LPHY     = LZPHY
      LWRITED  = LZWRITD
      LWRITEE  = LZWRITE
      LWRITEF  = LZWRITF
      LWRITET  = LZWRITT
      LAISTEP  = .FALSE.
      LDISTEP  = .FALSE.

      NDR      = NZDR
      FAKRMY   = 1.

      NA  = 3
      NJ  = 1
      NE  = 2
      NA2 = 2
      NJ2 = 1

      END SUBROUTINE EC4INMI
