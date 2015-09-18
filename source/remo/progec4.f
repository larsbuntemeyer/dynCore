      SUBROUTINE PROGEC4
     &  (
     &   T      , QD     , QW    , U     , V      , FI     , VERVEL ,
     &   PS     , SEAICE , FIB   , AK    , BK     , AKH    ,
     &   BKH    , DAK    , DBK   , PHI   , RLA    , COSLAT , SINLAT ,
     &   COSLON , SINLON , ALPHABOUND, A1T, A2T   , ACPHIR , CPHI   ,
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
     &   PINT   , DWDT  , ETAS  , W      , RPRAC )
      !
      IMPLICIT NONE
      !
      !C
      !C**** PROGEC4  -   DRIVING ROUTINE FOR ONE TIMESTEP
      !C**   CALLS    :   CALL PROGEC4 IN EC4ORG, EC4INMI
      !C**   ENTRIES  :      ---
      !C**   PURPOSE  :   COMPUTES ALL PROGNOSTIC VARIABLES FOR TIME NZT+1
      !C**
      !C**   DATE     :   05.10.04
      !C**
      !C**   EXTERNALS:   ECRANDUP, PHYEC , PROGEXP, SIKOR, ECRANDAS, SETRA
      !C**
      !C**   INPUT-
      !C**   PARAMETER:      ---
      !C**
      !C**   OUTPUT-
      !C**   PARAMETER:      ---
      !C**
      !C**   COMMON   :   PARAM  , ORG, COMDYN, COMPYH, PHYKON, HIGKON,
      !C**                PROGCHK, COMDIA
      !C**   METHOD   :   THE FORECAST AREA FROM JAH TO JEH FOR AN EXPLICIT
      !C**                FORECAST WILL BE DIVIDED INTO NTASKS SUBAREAS WHICH
      !C**                ARE USED PARALLEL. IN EACH SUBAREA THE EXPLICIT
      !C**                CALCULATION IS DONE IN SLABS. ALLOCATION OF MEMORY
      !C**                FOR ALL SUBAREAS (NTASKS) FOLLOWED BY THE SEMI-IMPLICIT
      !C**                CORRECTION, THE BOUNDARY RELAXATYION AND THE ASSELIN
      !C**                FILTER
      !C**   ERROR--
      !C**   HANDLING :      ---
      !C**   AUTHOR   :   R.PODZUN
      !C
      INCLUDE "parorg.h"
      INCLUDE "org.h"
      INCLUDE "comdyn.h"
      INCLUDE "comphy.h"
      INCLUDE "phykon.h"
      INCLUDE "higkon.h"
      INCLUDE "progchk.h"
      INCLUDE "comdia.h"
      INCLUDE "comhyd.h"
C
      INCLUDE "haloexchorg"
      !
      !----------------------------------------------------------------
      ! Dummy Arguments
      !
      INTEGER, INTENT(IN)  :: NOZ
      !
      !     PHYEC
      !
      REAL, INTENT(INOUT) ::
     &          T  (IEJEKE,3) , QD(IEJEKE,3), QW(IEJEKE,3),
     &          FI(IEJEKE,2)
      !
      REAL, INTENT(IN) ::
     &          TKE(IEJEKE,3) 
      !
      REAL, INTENT(INOUT) ::
     &          U(IE,JE,KE,3), V(IE,JE,KE,3)
      !
      REAL, INTENT(IN) ::
     &          ACLC(IEJEKE), ACLCAC(IEJEKE), VERVEL(IEJEKE),
     &          OZONPL(IEJE*NOZ), SO4ALL(IEJEKE), SO4NAT(IEJEKE)
      !
      REAL, INTENT(IN) ::
     &          EMTER(IEJE*KE1),TRSOL(IEJE*KE1),
     &          EMTEF(IEJE*2)  ,TRSOF(IEJE*2)  , VAR(IEJE*4)
      !
      REAL, INTENT(IN) ::
     &          AK     (KE1), BK    (KE1), AKH   (KE)  , BKH   (KE)  ,
     &          DAK    (KE) , DBK   (KE) , PHI   (IEJE), RLA   (IEJE),
     &          A1T    (KE1), A2T   (KE1), ACPHIR(JE,2), CPHI  (JE,2),
     &          CVDAES (KE1), CVDAEL(KE1), CVDAEU(KE1) , CVDAED(KE1) ,
     &          CEVAPCU(KE)
      !
      REAL, INTENT(IN) ::   GACPHIR(MOJE,2), GCPHI(MOJE,2)
      !
      REAL, INTENT(IN) ::
     &          COSLAT(IEJE), SINLAT(IEJE), COSLON(IEJE),
     &          SINLON(IEJE), ALPHABOUND(IE,JE,3)
      !
      REAL, INTENT(IN) ::
     &          TMKVMH(IE*(KE-1),JE,2), SOTHDT(IEKE,JE,2),
     &          UVTK  (IEKE,JE,2)     , TTK   (IEKE,JE)  ,
     &          QDTK  (IEKE,JE)       , TTS   (IEKE,JE)  ,
     &          QDTS  (IEKE,JE)       , QWTS  (IEKE,JE)  ,
     &          TMCM  (IEJE)          , TMCH  (IEJE)     ,
     &          TMCHL (IEJE)          , TMCHW (IEJE)     ,
     &          TMCHI (IEJE)
      !
      REAL, INTENT(IN) ::
     &          TSN   (IEJE,3), 
     &          TSECH (IEJE,3), WSECH (IEJE,3), SN    (IEJE,3),
     &          TSLECH(IEJE,3), 
     &          WL    (IEJE,3), TD    (IEJE,3), TDCL  (IEJE,3),
     &          TD3   (IEJE,3), TD4   (IEJE,3), TD5   (IEJE,3),
     &          TS    (IEJE,3), TB    (IEJE,3), TG    (IEJE,3),
     &          TGL   (IEJE,3), TGW   (IEJE,3), TGI   (IEJE,3)
      REAL, INTENT(INOUT) :: 
     &          PS    (IEJE,3), TSWECH(IEJE,3), TSIECH(IEJE,3),
     &          QDB   (IEJE,3), 
     &          QDBL  (IEJE,3), QDBW  (IEJE,3), QDBI  (IEJE,3)
      !
      REAL, INTENT(IN) ::
     &          FIB   (IEJE), 
     &          TEFF  (IEJE), CAPE  (IEJE),
     &          TSMAX (IEJE), TSMIN (IEJE), TSLIN (IEJE),
     &          DSNAC (IEJE), SNMEL (IEJE), RUNOFF(IEJE),
     &          VAROR (IEJE), SRFL  (IEJE), THFL  (IEJE),
     &          QHFL  (IEJE), XHFL  (IEJE), RSFC  (IEJE),
     &          SSFC  (IEJE), RSFL  (IEJE), SSFL  (IEJE),
     &          DHFT  (IEJE),
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
      REAL, INTENT(INOUT) ::
     &          SEAICE(IEJE), SICED (IEJE)
      !
      REAL, INTENT(INOUT) ::
     &          AHFL  (IEJE), AHFS  (IEJE),
     &          AHFSL (IEJE), AHFSW (IEJE), AHFSI (IEJE)
      !
      REAL, INTENT(IN) ::
     &          WS1    (IEJE), WS2  (IEJE), WS3   (IEJE),
     &          WS4    (IEJE), WS5  (IEJE), DZR   (IEJE),
     &          DZS    (IEJE), FKSAT(IEJE), FMPOT (IEJE),
     &          BCLAPP (IEJE), VPOR (IEJE), ETRANS(IEJE),
     &          EBSOIL (IEJE), ESNOW(IEJE), ESKIN (IEJE),
     &          ERES   (IEJE)
      REAL ::   PINT(IE,JE,KE1,3), DWDT(IE,JE,KE ,3),
     &          ETAS(IE,JE,KE1  ), W   (IE,JE,KE1,3)

      LOGICAL, INTENT(IN) ::
     &          LOGLAC(IEJE), LOLAND(IEJE), LOSEA (IEJE),
     &          LOICE (IEJE), LALAND(IEJE)
      !
      !     SIKOR
      !
      REAL, INTENT(IN) ::
     &          SISTM(KE,KE), SIGAM (KE,KE), SITAU(KE,KE), SINUE(KE),
     &          SIVMT(KE,KE), SIVMTI(KE,KE), SICQ (KE)
      REAL, INTENT(IN) ::
     &          TRIGSI(5*(MOIE-1)/2), RZ1I(6), RZ2I(6), IFAXI(11)
      !
      !     ECRANDAS
      !
      REAL, INTENT(IN) ::
     &          UR (IEJEKE,2), VR (IEJEKE,2), TR (IEJEKE,2),
     &          QDR(IEJEKE,2), QWR(IEJEKE,2)
      !
      !     PROGEXP
      !
      REAL, INTENT(IN) ::
     &          VVFH   (KE)  , FC     (IEJE), BFLHS  (IEJE),
     &          BFLQDS (IEJE), BFLUS  (IEJE), BFLVS  (IEJE),
     &          BFLHSL (IEJE), BFLHSW (IEJE), BFLHSI (IEJE),
     &          BFLQDSL(IEJE), BFLQDSW(IEJE), BFLQDSI(IEJE)
      !
      !     ECRANDUP
      !
      REAL, INTENT(IN) ::
     &          PSR    (IEJE,2), TSWECHR(IEJE,2), TSIECHR(IEJE,2),
     &          SEAICER(IEJE,2), QDBLR  (IEJE,2), SICEDR (IEJE,2)
      INTEGER, INTENT(INOUT) :: INFRL(IEJE), INFRW(IEJE), INFRI(IEJE)
CSP
      REAL, INTENT(INOUT) :: QI(IEJEKE,3)
      REAL, INTENT(IN) :: QITS  (IEKE,JE), QIVI (IEJE),
     &                    RPRAC(IEJEKE)
      !
      ! End of Dummy Arguments
      !----------------------------------------------------------------
      ! Local Variables
      !
      REAL    :: QIR(IEJEKE,2)
      !
      INTEGER :: I,J,IJ,K
      REAL    :: ZDTIRD,ZEMZDR
      REAL    :: ZPW (IE,JE,KE   ),
     &           DIFW(IE,JE,KE   ), DIFD(IE,JE,KE   )
CKS
C     FOR HALO EXCHANGE
      INTEGER :: KSIZE(10), ITYPE
      !
      ! End of Local Variables
      !----------------------------------------------------------------
      ! 
      ! Allocate Dimensions
      !
      ! DIMENSION OM850MZ(16), OM500MZ(16), OM300MZ(16)
      ! These are declared in progchk.h
      DIMENSION OM850M(1), OM500M(1), OM300M(1)
C
C     BFLUS, BFLVS WERDEN AUF 0. GESETZT, DA SIE BEI ECHAM4-PHYSIK
C     NICHT AUSGEWERTET WERDEN
C
      IF (NZT.EQ.NANF+1) THEN
        CALL SETRA(BFLUS,IEJE,0.)
        CALL SETRA(BFLVS,IEJE,0.)
      ENDIF
      !
      !     DT2,ED2DT,DTDEH AND NEHDDT IN COMMON *HIGKON*
      DT2    = 2.0*DT
      ED2DT  = 1.0/DT2
      DTDEH  = DT/3600.
      NEHDDT = NINT( 3600./DT )
      !
      OM850M = 0.0
      OM500M = 0.0
      OM300M = 0.0
C
CHG   INIT OF SOME VALUES
      DIFW(:,:,:)   = 0.0
      DIFD(:,:,:)   = 0.0

C     COEFFICIENT FOR LINEAR TIME INTERPOLATION OF BOUNDARIES
      ZDTIRD = REAL( NZT + 1 - NZTRD ) / REAL( NDR )
      ZEMZDR = 1. - ZDTIRD
C
C     DEFINE ALL PROGNOSTIC VARIABLE FOR TIME NZT+1 WITH
C     CORRESPONDING BOUNDARY VALUES

      CALL ECRANDUP(ZEMZDR , ZDTIRD,
     &     U      , V     , T    , QD     , QW     , QI     ,
     &     UR     , VR    , TR   , QDR    , QWR    , PS     , TSWECH,
     &     TSIECH , SEAICE, INFRL, INFRW  , INFRI  , QDB    , QDBL  ,
     &     QDBW   , QDBI  , PSR  , TSWECHR, TSIECHR, SEAICER, QDBLR ,
     &     SICED  , SICEDR, QIR  , PINT   , DWDT   , W      , AK    ,
     &     BK)
      !
      !     CALL SUBROUTINES FOR PHYSICS
      IF ( LPHY ) THEN
        !
        CALL PHYEC
        !
        ! actual parameter names correspond to dummy argument names in:
        ! phyec
        !
     &       (
     &        T     , QD    , QW     , U     , V     , VERVEL,
     &        PS    , AK    , BK     , AKH   , BKH   , DAK   ,
     &        DBK   , PHI   , RLA    , COSLAT, SINLAT, COSLON, SINLON,
     &        ALPHABOUND, A1T, A2T   , ACPHIR, CPHI  , TGL   , TGW   ,
     &        TGI   , QDBL  , QDBW   , QDBI  , QDB   , TS    , TB    ,
     &        TG    , SOTHDT, UVTK   , TMKVMH, TTK   , QDTK  , TTS   ,
     &        QDTS  , QWTS  , TMCM   , TMCH  ,
     &        TMCHL , TMCHW , TMCHI  , SICED , TEFF  , TSECH ,
     &        TSLECH, TSWECH, TSIECH , WSECH , SN    , WL    , TD    ,
     &        TDCL  , TD3   , TD4    , TD5   , TSN   , TSURF , TSMAX ,
     &        TSMIN , TSLIN , DSNAC  , SNMEL , RUNOFF, VAROR , SRFL  ,
     &        THFL  , QHFL  , XHFL   , RSFC  , SSFC  , RSFL  , SSFL  ,
     &        AHFL  , AHFS  , DHFT   , DHFQW , AHFSL , AHFSW , AHFSI ,
     &        AHFICE, QRES  , AZ0L   , AZ0W  ,
     &        AZ0I  , DHFQS , TOPMAX , AZ0   , APRC  , APRL  , APRS  ,
     &        EVAP  , EVAPM , ACLC   , ACLCAC, VGRAT , FOREST, ALBECH,
     &        EVAPL , EVAPW , EVAPI  , ALSOL , ALSOW , ALSOI , ALBEDO,
     &        TKE   , DEW2  , WSMX   , VLT   , FAO   , RGCGN ,
     &        TEMP2 , T2MAX , T2MIN  , USTAR3, USTR  , U10   , USTRL ,
     &        USTRW , USTRI , VSTRL  , VSTRW , VSTRI , INFRL , INFRW ,
     &        INFRI , VDIS  , VSTR   , V10   , WIND10, WIMAX , ACLCOV,
     &        ALWCVI, QVI   , EMTER  , TRSOL , EMTEF , TRSOF , SCLF0 ,
     &        SCLFS , SRAF0 , SRAFS  , TCLF0 , TCLFS , TRAF0 , TRAFS ,
     &        ACLCV , SRADS , SRADSU , SRAD0 , SRAD0U, TRADS , TRADSU,
     &        TRAD0 , USTRGW, VDISGW , VSTRGW, VAR   , CVDAES, CVDAEL,
     &        CVDAEU, CVDAED, CEVAPCU, DRAIN ,TLAMBDA,DLAMBDA, PORVOL,
     &        FCAP  , WI3   , WI4    , WI5   , WI    , WICL  , LOGLAC,
     &        LOLAND, LOSEA , LOICE  , LALAND, GHPBL , BETA  , WMINLOK,
     &        WMAXLOK,CAPE  , OZONPL , NOZ   , SO4ALL, SO4NAT, WS1   ,
     &        WS2   , WS3   , WS4    , WS5   , DZR   , DZS   , FKSAT ,
     &        FMPOT , BCLAPP, VPOR   , ETRANS, EBSOIL, ESNOW , ESKIN ,
     &        ERES  , QI    , QITS   , QIVI  , RPRAC)

      END IF

      CALL PROGEXP(
        !
        ! WARNING! Here are 3 dummy argument names that
        ! differ from actual parameter names
        !
        !     |        |        |
        !     v        v        v
     &     OM850M(1), OM500M(1), OM300M(1),
        !     ^        ^        ^
        !     |        |        |
     &     AK    , BK    , AKH    , BKH    , DAK    , DBK  , A1T   ,
     &     A2T   , VVFH  , FIB    , FC     , ACPHIR , CPHI , PS    ,
     &     TGL   , TGW   , TGI    , QDBL   , QDBW   , QDBI , BFLHSL,
     &     BFLHSW, BFLHSI, BFLQDSL, BFLQDSW, BFLQDSI, TG   , QDB   ,
     &     BFLHS , BFLQDS, BFLUS  , BFLVS  , U      , V    , T     ,
     &     QD    , QW    , FI     , TMKVMH , TMCHL  , TMCHW,
     &     TMCHI , TMCM  , TMCH   , SOTHDT , TTK    , QDTK , UVTK  ,
     &     TTS   , QDTS  , QWTS   , INFRL  , INFRW  , INFRI, QI    ,
     &     QITS  , PINT  , DWDT   , ETAS )
C
C     COMPUTES THE NON-HYDROSTATIC CORRECTION AFTER JANJIC ET AL. (2001)
C
      IF (.NOT.LHYDRO) THEN
C
C     UPDATE HALOS FOR PROGNOSTIC VARIABLES
C
        ITYPE = 95
        KSIZE(1:10) = (/KE,KE,KE,1,KE,KE,KE,KE1,KE1,0/)
        CALL HALOEXCH(KSIZE, ITYPE, U(:,:,:,NE), V(:,:,:,NE), T(:,NE),
     &       PS(:,NE), QD(:,NE), QW(:,NE), QI(:,NE),
     &       PINT(:,:,:,NE), ETAS)
C
C     COMPUTE W=DZ/DT DIAGNOSTICALLY
C
        CALL CDZDT( 
     &       DAK    , DBK   , U      , V   ,
     &       ACPHIR , DT    , FIB    , QI  ,
     &       DWDT   , PINT  , W      , T   , QD  , QW,
     &       PS     , ETAS  , ZPW    )
C
C       UPDATE W IN HALOS
C
        ITYPE = 295
        KSIZE(1:10) = (/KE1,0,0,0,0,0,0,0,0,0/)
        CALL HALOEXCH(KSIZE, ITYPE, W(:,:,:,NE))
C
        IF( NZT > 3 ) THEN
C
C       COMPUTE DW/DT
C
         CALL CDWDT( 
     &         DAK    , DBK   ,
     &         ACPHIR , DT    , ETAS   ,
     &         DWDT   , W     , U   , V   , PS )
C
C        UPDATE DW/DT IN HALOS
C
         ITYPE = 195
         KSIZE(1:10) = (/KE,0,0,0,0,0,0,0,0,0/)
         CALL HALOEXCH(KSIZE, ITYPE, DWDT(:,:,:,NE))
C
C        CALCULATE DIFFUSION TERMS FOR DW/DT
C
         CALL CDWDT2( DIFW  , DIFD   , DWDT   , W  )
C
C        CALCULATE THE FINAL UPDATE OF THE PROGNOSTIC VARIABLES DUE TO
C        NON-HYDROSTATIC CORRECTION
C
         CALL FINAL( 
     &         AK     , BK    , DAK    , DBK  ,
     &         DWDT   , PINT  , W      , T      , QD  , QW, PS ,
     &         DIFW   , DIFD  , QI)
C
       ELSE
         DO K = 1 , KE1
           DO J = JAA, JEA
             IJ = (J - 1)*IE
             DO I = IAA , IEA
               PINT(I,J,K,NE) = AK(K) + BK(K) * PS(IJ+I,NE)
             ENDDO
           ENDDO
         ENDDO

         DO K = 1 , KE
           DO J = JAA , JEA
             DO I = IAA , IEA
               DWDT(I,J,K,NE) = 1.
             ENDDO
           ENDDO
         ENDDO
        ENDIF
C
      ELSE !LHYDRO=true
      !
      !     ANTEILE AN DER VERTIKALBEWEGUNG SUMMIEREN
C
        DO K = 1 , KE1
          DO J = JAA, JEA
            IJ = (J - 1)*IE
            DO I = IAA , IEA
              PINT(I,J,K,NE) = AK(K) + BK(K) * PS(IJ+I,NE)
              W   (I,J,K,NE) = 0.
            ENDDO
          ENDDO
        ENDDO

        DO K = 1 , KE
          DO J = JAA , JEA
            DO I = IAA , IEA
              DWDT(I,J,K,NE) = 1.
            ENDDO
          ENDDO
        ENDDO

      ENDIF !LHYDRO
C
      AHFS (1:IEJE) = BFLHS (1:IEJE)
      AHFSL(1:IEJE) = BFLHSL(1:IEJE)
      AHFSW(1:IEJE) = BFLHSW(1:IEJE)
      AHFSI(1:IEJE) = BFLHSI(1:IEJE)
      AHFL (1:IEJE) = BFLQDS(1:IEJE)

C     ANTEILE AN DER VERTIKALBEWEGUNG SUMMIEREN
      COUNT = 1
      CALL PREDUCER(OM850M,MPI_SUM)
      CALL PREDUCER(OM500M,MPI_SUM)
      CALL PREDUCER(OM300M,MPI_SUM)
      !
      !-----------------------------------------------------------------------
      !     PROGNOSTISCHER TEIL, SCHRITT 2: SEMI-IMPLIZITE KORREKTUREN
      !     ------------------------------------------------------------------
      !
      !     SEMI-IMPLIZITE KORREKTURTERME BERECHNEN UND DIE WERTE DER
      !     EXPLIZITEN PROGNOSE (NE) ENTSPRECHEND MODIFIZIEREN.
      IF ( LSITS )  THEN
        CALL SIKOR
     &       (
     &        U     , V    , T    , PS   , ACPHIR, CPHI  ,
     &        SISTM , SIGAM, SITAU, SINUE, SIVMT , SIVMTI, SICQ,
     &        TRIGSI, RZ1I , RZ2I , IFAXI, GACPHIR, GCPHI,
     &        PINT  , DWDT , AK   , BK    )
      ENDIF
      !
      !-----------------------------------------------------------------------
      !     PROGNOSTISCHER TEIL, SCHRITT 3: RANDRELAXATION UND ASSELIN-FILTER
      !     ------------------------------------------------------------------
      !     INITIALISIERUNGSZEITSCHRITTE UEBERSPRINGEN RANDRELAXATION UND
      !     ASSELIN-FILTERUNG
      !KS   UPDATE OF THE HALOS
      IF (.NOT.LAISTEP) THEN
        CALL ECRANDAS
     &       (
     &        AKH , BKH , ALPHABOUND, PS  ,
     &        U   , V   , T   , QD  , QW    , UR    ,
     &        VR  , TR  , QDR , QWR , PSR   , QI    ,
     &        QIR , PINT, DWDT, W   , AK    , BK    , DAK   , DBK)
      ELSE
C
C       1. LINKEN RAND BELEGEN
C
        IF (NEIGHBOR(1) .EQ. -1) THEN
          DWDT(1:2,:,:,NE)=1.
        ENDIF
C
C       2. OBEREN RAND BELEGEN
C
        IF (NEIGHBOR(2) .EQ. -1) THEN
          DWDT(:,JE-1:JE,:,NE)=1.
        ENDIF
C
C       3. RECHTEN RAND BELEGEN
C
        IF (NEIGHBOR(3) .EQ. -1) THEN
          DWDT(IE-1:IE,:,:,NE)=1.
        ENDIF
C
C       4. UNTEREN RAND BELEGEN
C
        IF (NEIGHBOR(4) .EQ. -1) THEN
          DWDT(:,1:2,:,NE)=1.
        ENDIF
      ENDIF !LAISTEP
C
C-----------------------------------------------------------------------
C     DIAGNOSTISCHER TEIL:    BERECHNUNG DES GEOPOTENTIALS FUER DEN
C     -------------------     ZEITPUNKT NE
      CALL GEOPOT(NE, NA2, 1, IE, 1, JE,
     &            T  , QD, QW  , FI    , FIB ,
     &            QI , PINT    , DWDT)
      !-----------------------------------------------------------------------
      !
      RETURN
      END SUBROUTINE PROGEC4
