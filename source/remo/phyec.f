      SUBROUTINE PHYEC
              !
              ! actual parameter names correspond to dummy argument names
              ! when called in :
              ! progec4
              !
     &        (
     &         T     , QD    , QW     , U     , V     , VERVEL,
     &         PS    , AK    , BK     , AKH   , BKH   , DAK   ,
     &         DBK   , PHI   , RLA    , COSLAT, SINLAT, COSLON, SINLON,
     &         ALPHABOUND, A1T, A2T   , ACPHIR, CPHI  , TGL   , TGW   ,
     &         TGI   , QDBL  , QDBW   , QDBI  , QDB   , TS    , TB    ,
     &         TG    , SOTHDT, UVTK   , TMKVMH, TTK   , QDTK  , TTS   ,
     &         QDTS  , QWTS  , TMCM   , TMCH  ,
     &         TMCHL , TMCHW , TMCHI  , SICED , TEFF  , TSECH ,
     &         TSLECH, TSWECH, TSIECH , WSECH , SN    , WL    , TD    ,
     &         TDCL  , TD3   , TD4    , TD5   , TSN   , TSURF , TSMAX ,
     &         TSMIN , TSLIN , DSNAC  , SNMEL , RUNOFF, VAROR , SRFL  ,
     &         THFL  , QHFL  , XHFL   , RSFC  , SSFC  , RSFL  , SSFL  ,
     &         AHFL  , AHFS  , DHFT   , DHFQW , AHFSL , AHFSW , AHFSI ,
     &         AHFICE, QRES  , AZ0L   , AZ0W  ,
     &         AZ0I  , DHFQS , TOPMAX , AZ0   , APRC  , APRL  , APRS  ,
     &         EVAP  , EVAPM , ACLC   , ACLCAC, VGRAT , FOREST, ALBECH,
     &         EVAPL , EVAPW , EVAPI  , ALSOL , ALSOW , ALSOI , ALBEDO,
     &         TKE   , DEW2  , WSMX   , VLT   , FAO   , RGCGN ,
     &         TEMP2 , T2MAX , T2MIN  , USTAR3, USTR  , U10   , USTRL ,
     &         USTRW , USTRI , VSTRL  , VSTRW , VSTRI , INFRL , INFRW ,
     &         INFRI , VDIS  , VSTR   , V10   , WIND10, WIMAX , ACLCOV,
     &         ALWCVI, QVI   , EMTER  , TRSOL , EMTEF , TRSOF , SCLF0 ,
     &         SCLFS , SRAF0 , SRAFS  , TCLF0 , TCLFS , TRAF0 , TRAFS ,
     &         ACLCV , SRADS , SRADSU , SRAD0 , SRAD0U, TRADS , TRADSU,
     &         TRAD0 , USTRGW, VDISGW , VSTRGW, VAR   , CVDAES, CVDAEL,
     &         CVDAEU, CVDAED, CEVAPCU, DRAIN ,TLAMBDA,DLAMBDA, PORVOL,
     &         FCAP  , WI3   , WI4    , WI5   , WI    , WICL  , LOGLAC,
     &         LOLAND, LOSEA , LOICE  , LALAND, GHPBL , BETA  , WMINLOK,
     &         WMAXLOK,CAPE  , OZONPL , NOZ   , SO4ALL, SO4NAT, WS1   ,
     &         WS2   , WS3   , WS4    , WS5   , DZR   , DZS   , FKSAT ,
     &         FMPOT , BCLAPP, VPOR   , ETRANS, EBSOIL, ESNOW , ESKIN ,
     &         ERES  , QI    , QITS   , QIVI  , RPRAC)
      !
      IMPLICIT NONE
      !
C
C     COMMON-BLOCKS FROM DWD-PHYSICS
C
C     ORG: NZT, NHORPH
C     HIGKON: RDDRM1
C     PHYKON: R
C     COMDYN: DT
C
      INCLUDE "parorg.h"
      INCLUDE "org.h"
      INCLUDE "higkon.h"
      INCLUDE "phykon.h"
      INCLUDE "comdia.h"
      INCLUDE "comdyn.h"
      INCLUDE "comecphy.h"
      INCLUDE "YOTLUC"
C
      INCLUDE "haloexchorg"
C
      INTEGER, INTENT(IN)  :: NOZ
      !
      REAL, INTENT(INOUT) ::
     &          T  (IEJE,KE,3) , QD(IEJE,KE,3), QW(IEJE,KE,3),
     &          TKE(IEJE,KE,3)
C
      REAL, INTENT(INOUT) ::
     &          U(IEJE,KE,3), V(IEJE,KE,3)
C
      REAL, INTENT(INOUT) ::
     &          ACLC(IEJE,KE), ACLCAC(IEJE,KE), VERVEL(IEJE,KE),
     &          OZONPL(IEJE,NOZ), SO4ALL(IEJE,KE), SO4NAT(IEJE,KE)
C
      REAL, INTENT(INOUT) ::
     &          EMTER(IEJE,KE1),TRSOL(IEJE,KE1),
     &          EMTEF(IEJE,2)  ,TRSOF(IEJE,2)  , VAR(IEJE,4)
C
      REAL, INTENT(INOUT) ::
     &          AK     (KE1), BK    (KE1), AKH   (KE)  , BKH   (KE)  ,
     &          DAK    (KE) , DBK   (KE) , PHI   (IEJE), RLA   (IEJE),
     &          A1T    (KE1), A2T   (KE1), ACPHIR(JE,2), CPHI  (JE,2),
     &          CVDAES (KE1), CVDAEL(KE1), CVDAEU(KE1) , CVDAED(KE1) ,
     &          CEVAPCU(KE)
C
      REAL, INTENT(OUT) ::
     &          COSLAT(IEJE), SINLAT(IEJE), COSLON(IEJE),
     &          SINLON(IEJE), TMKVMH(IE*(KE-1),JE,2)
      REAL, INTENT(IN) :: ALPHABOUND(IE,JE,3)
C
      REAL, INTENT(OUT) ::
     &          SOTHDT(IEKE,JE,2), UVTK(IEKE,JE,2), TTS(IEKE,JE)
      REAL, INTENT(OUT) ::
     &          TTK   (IEKE,JE)  ,
     &          QDTK  (IEKE,JE)       , 
     &          QDTS  (IEKE,JE)       , QWTS  (IEKE,JE)  
      REAL, INTENT(OUT) :: 
     &          TMCM  (IEJE)          , TMCH  (IEJE)     
      REAL, INTENT(OUT) ::
     &          TMCHL (IEJE)          , TMCHW (IEJE)     ,
     &          TMCHI (IEJE)
C
      REAL, INTENT(OUT) :: WSECH (IEJE,3), WL    (IEJE,3)
      REAL, INTENT(OUT) ::
     &          QDB   (IEJE,3), PS   (IEJE,3),
     &          QDBL  (IEJE,3), QDBW  (IEJE,3), QDBI (IEJE,3)
      REAL, INTENT(INOUT) ::
     &          TSN   (IEJE,3), 
     &          TSECH (IEJE,3), SN    (IEJE,3),
     &          TSLECH(IEJE,3), TSWECH(IEJE,3) ,TSIECH(IEJE,3),
     &          TD    (IEJE,3), TDCL (IEJE,3),
     &          TD3   (IEJE,3), TD4   (IEJE,3), TD5  (IEJE,3)
      REAL, INTENT(OUT) ::
     &          TS    (IEJE,3), TB    (IEJE,3), TG   (IEJE,3),
     &          TGL   (IEJE,3), TGW   (IEJE,3), TGI  (IEJE,3)
C
      REAL, INTENT(INOUT) ::
     &          SICED (IEJE), TEFF  (IEJE), CAPE  (IEJE),
     &          TSMAX (IEJE), TSMIN (IEJE), TSLIN (IEJE),
     &          DSNAC (IEJE), SNMEL (IEJE), RUNOFF(IEJE),
     &          VAROR (IEJE), SRFL  (IEJE), THFL  (IEJE),
     &          QHFL  (IEJE), XHFL  (IEJE)
      REAL, INTENT(OUT)::
     &          RSFC  (IEJE),
     &          SSFC  (IEJE), RSFL  (IEJE), SSFL  (IEJE)
      REAL, INTENT(INOUT) ::
     &          AHFL  (IEJE), AHFS  (IEJE), DHFT  (IEJE),
     &          AHFSL (IEJE), AHFSW (IEJE), AHFSI (IEJE)
      REAL, INTENT(OUT) ::
     &          AZ0L  (IEJE), AZ0W  (IEJE), AZ0I  (IEJE),
     &          AZ0   (IEJE), WSMX  (IEJE)
      REAL, INTENT(INOUT)  ::
     &          AHFICE(IEJE), QRES  (IEJE), TSURF (IEJE),
     &          DHFQW (IEJE), DHFQS (IEJE), TOPMAX(IEJE)
      REAL, INTENT(OUT) ::
     &          APRC  (IEJE),
     &          APRL  (IEJE), APRS  (IEJE)
      REAL, INTENT(INOUT)  ::
     &          EVAP  (IEJE),
     &          EVAPM (IEJE), VGRAT (IEJE), FOREST(IEJE),
     &          EVAPL (IEJE), EVAPW (IEJE), EVAPI (IEJE),
     &          ALSOL (IEJE), ALSOW (IEJE), ALSOI (IEJE),
     &          ALBECH(IEJE), ALBEDO(IEJE), DEW2  (IEJE),
     &          VLT   (IEJE), FAO   (IEJE),
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
     &          WS1    (IEJE), WS2  (IEJE), WS3   (IEJE),
     &          WS4    (IEJE), WS5  (IEJE), DZR   (IEJE),
     &          DZS    (IEJE), FKSAT(IEJE), FMPOT (IEJE),
     &          BCLAPP (IEJE), VPOR (IEJE), ETRANS(IEJE),
     &          EBSOIL (IEJE), ESNOW(IEJE), ESKIN (IEJE),
     &          ERES   (IEJE)
      LOGICAL, INTENT(INOUT) ::
     &          LOGLAC(IEJE), LOLAND(IEJE), LOSEA (IEJE),
     &          LOICE (IEJE), LALAND(IEJE)
      REAL, INTENT(INOUT)    ::   
     &          QI(IEJE,KE,3), QITS(IEKE,JE), QIVI(IEJE),
     &          RPRAC(IEJE,KE)
C
C     DECLARATION OF LOCAL WORK ARRAYS
C
      REAL   ::
     &          VCT(2*KE1)
C
      REAL   ::
     &          ZTRAD (IEJE,KE)   , ZTGWDR(IEJE,KE)   ,
     &          ZUGWDR(IEJE,KE)   , ZVGWDR(IEJE,KE)   ,
     &          ZTTS  (IEJE,KE),
     &          ZQDTS (IEJE,KE)   , ZQWTS (IEJE,KE)   , ZTTK  (IEJE,KE),
     &          ZUTK  (IEJE,KE)   , ZVTK  (IEJE,KE)   , ZQDTK (IEJE,KE),
     &          Z2ZM  (IEJE,KE-1), Z2ZH  (IEJE,KE-1), ZDQDT (IEJE,KE)
C
      REAL   ::
     &          ZW5    (IEKE)     ,
     &          ZW7    (IE*(KE-1)), ZW8  (IE*(KE-1)) ,
     &          ZW9    (IEKE)     , ZW10 (IEKE)      ,
     &          ZW11   (IEKE)     , ZW12 (IEKE)      ,
     &          ZW13   (IEKE)     , ZQDBNA (IEJE)    ,
     &          ZQDBLNA(IEJE)     , ZQDBWNA(IEJE)    ,
     &          ZQDBINA(IEJE)     , ZW3(IEKE)        ,
     &          ZW4(IEKE)
C
      REAL   ::
     &          ZU(IEJE,KE), ZV(IEJE,KE), ZTMCM(IEJE,KE), ZTMCH(IEJE,KE)
C
      LOGICAL ::   LDIUR  , LVDIFF, LSURF, LCONV, LCOND
      LOGICAL ::   LGWDRAG, LSOLC , LAER , LCFC , LGADSRH
      LOGICAL ::   LO1, LO2
      INTEGER ::   INFRL(IEJE), INFRW(IEJE), INFRI(IEJE)
CSP
      REAL    ::   ZQITS(IEJE,KE)
      REAL    ::   ZW14(IEKE)
      !
      REAL    :: EPS,ZCONACC,ZDPN2,ZDPO,ZDPU,ZGROP,ZRHO,ZTN,ZTVO
      REAL    :: ZTVU,ZDT
      INTEGER :: I,ICOUNT,IJ,IJJ,IPHYSTART,IPHYSTOP
      INTEGER :: J,JX,JY
      INTEGER :: K,KB,KLEV,KS,KTASK
      INTEGER :: NCBASE,NCDATA,NMONTH,NPOINT,NPTASK,NRADFR,NRADIA
      INTEGER :: NRADPFR,NTBASE,NTDATA,NTIMST,NVDIFF
CKS
C     FOR HALO EXCHANGE
      INTEGER :: KSIZE(10), ITYPE
C
C     INITIALISIERUNG VON FELDERN UND KONSTANTEN
C
      TMCHL   (:) = 0.0
      TMCHW   (:) = 0.0
      TMCHI   (:) = 0.0
C
      ZTRAD (:,:) = 0.0
      ZTGWDR(:,:) = 0.0
      ZUGWDR(:,:) = 0.0
      ZVGWDR(:,:) = 0.0
      ZDQDT (:,:) = 0.0
      ZU    (:,:) = 0.0
      ZV    (:,:) = 0.0
      ZTTS  (:,:) = 0.0
      ZQDTS (:,:) = 0.0
      ZQWTS (:,:) = 0.0
      ZTTK  (:,:) = 0.0
      ZUTK  (:,:) = 0.0
      ZVTK  (:,:) = 0.0
      ZQDTK (:,:) = 0.0
      ZTMCM (:,:) = 0.0
      ZTMCH (:,:) = 0.0
      ZQITS (:,:) = 0.0
      Z2ZM  (:,:) = 0.0
      Z2ZH  (:,:) = 0.0
C
      VCT(1:KE1) = AK(:)
      VCT(1+KE1:2*KE1) = BK(:)
C
      ZQDBNA (:) = QDB(:,NA)
      ZQDBLNA(:) = QDBL(:,NA)
      ZQDBWNA(:) = QDBW(:,NA)
      ZQDBINA(:) = QDBI(:,NA)
C
      LO1 = NZT.EQ.0
      LO2 = (NZT.EQ.NANF+1).AND.(NANF.NE.0)
C
      IF (LO1) THEN
C
         AZ0(:) = MAX(AZ0(:),1.5E-05)
         AZ0L(:) = AZ0(:)
         AZ0W(:) = AZ0(:)
         AZ0I(:) = AZ0(:)
         WSMX(:) = AMAX1(WSMX(:),1.E-10)
         WSECH(:,NA) = AMIN1(WSECH(:,NA),WSMX(:))
         WSECH(:,NJ) = AMIN1(WSECH(:,NJ),WSMX(:))
C
         ZCONACC=1.0
         ZDT=2.*DT
C
      ELSE
C
         ZDT=DT
         ZCONACC=0.5
C
      ENDIF
C
      IF (LO1.OR.LO2) THEN
C
         COSLAT(:) = COS(PHI(:))
         SINLAT(:) = SIN(PHI(:))
         COSLON(:) = COS(RLA(:))
         SINLON(:) = SIN(RLA(:))
C
      ENDIF
C
C     INITIALIZATION OF ....
C
      CALL INIT
     &     (KE     , IEJE  , KE1    ,
     &      NZT    , NANF  , NRADIA , NVDIFF, CEVAPCU, VCT   , FAO    ,
     &      RGCGN  , CVDAES, CVDAEL , CVDAEU, CVDAED , NRADFR, LGADSRH,
        !
        !   WARNING: LECRAD differs from dummy name 'LRAD'
        !                                                     ↓    
     &      NRADPFR, LSOLC , LAER   , LCFC  , EPS    , LDIUR , LECRAD ,
        !                                                       ↓    
        !                                                       LRAD              
     &      LVDIFF , LSURF , LGWDRAG, LCONV , LCOND  , NCBASE, NCDATA ,
     &      NTBASE , NTDATA, NMONTH , NTIMST, TLAMBDA,DLAMBDA, PORVOL ,
     &      FCAP)
        !   All other dummy names correspond to their actual parameters
C
C     PREPARING GREENHOUSE-GASES AND CFCS IF REQUIRED
C
      IF (LSCEN) THEN
         IF (LO1.OR.LO2) CALL READGTS
         CALL PREPGRG
      ENDIF
C
C     OMEGA BERECHNEN, WENN CUMULUSKONVEKTION AN
C
      IF (LCONV) THEN
         CALL ECKONT(
              !
              ! WARNING: Following actual parameters have different
              ! dummy names in ECKONT:
              ! VERVEL,ZDQDT
              !
     &        VERVEL, ZDQDT,
     &        AK , BK , AKH   , BKH   , DAK   ,
     &        DBK   , ALPHABOUND, A1T, A2T, ACPHIR, CPHI  , PS    ,
     &        QDB   , U    , V  , T  , QD    , TMKVMH, TMCH)
      END IF
C
CJ-PP IF HALOS IGNORED IN PHYSICS, FOLLOWING CALCULATIONS UNNECESSARY
      IF(.NOT. LNOHALO) THEN
         IF (NHORPH .NE. 0) THEN
            IF (MOD(IEJE,NHORPH).EQ.0) THEN
               NPTASK = IEJE/NHORPH
            ELSE
               NPTASK = IEJE/NHORPH+1
            ENDIF
         ELSE
            NPTASK = 1
         ENDIF
      ENDIF

C
C     FUER DIE PARALLELE VERSION:
C     DIE ROUTINE PHYEC WIRD VON JEDEM PE EINZELN AUFGERUFEN, ALSO FUER JEDES
C     TEILGEBIET.
C     DA ES PRO TEILGEBIET KEIN MULTITASKING GIBT, WIRD NPTASK=1 GESETZT.
C     (ALLERDINGS KANN ES AUS SPEICHERPLATZGRUENDEN DOCH ERFORDERLICH SEIN,
C     VON NHORPH GEBRAUCH ZU MACHEN, ALSO DAS BERECHNETE NPTASK ZU BENUTZEN)
C
C     BESTIMMUNG DER BERECHNUNGSINDICES (DIE AEUSSERSTE RANDZEILE,
C     DIE DURCH MPP DAZUGENOMMEN WURDE, WIRD WEGGELASSEN)
C
C!CDIR CONCUR
C!CDIR NODEP
      DO KLEV=1,KE
         DO K=1,IEJE
            JY = (K-1)/IE+1
            JX = K-(JY-1)*IE
            IF (((JX .EQ. 1) .OR. (JX .EQ. 2)) .AND.
     &           (NEIGHBOR(1) .EQ. -1)) THEN
               ZU(K,KLEV)  = U(K,KLEV,NA)
            ELSE IF ((JX .EQ. 1) .AND. (NEIGHBOR(1) .NE. -1)) THEN
               ZU(K,KLEV)  = U(K,KLEV,NA)
            ELSE
               ZU(K,KLEV)  = 0.5*(U(K,KLEV,NA)+U(K-1,KLEV,NA))
            ENDIF

            IF (((JY .EQ. 1) .OR. (JY .EQ. 2)) .AND.
     &           (NEIGHBOR(4) .EQ. -1)) THEN
               ZV(K,KLEV)  = V(K,KLEV,NA)
            ELSE IF ((JY .EQ. 1) .AND. (NEIGHBOR(4) .NE. -1)) THEN
               ZV(K,KLEV)  = V(K,KLEV,NA)
            ELSE
               ZV(K,KLEV)  = 0.5*(V(K,KLEV,NA)+V(K-IE,KLEV,NA))
            ENDIF
         ENDDO
      ENDDO
C
C
      IF (LDIA .AND. NZT.EQ.NANF+1 .AND. MYID .EQ. 0) THEN
         WRITE(*,'(12A8)') 'LRAD','LVDIFF','LSURF','LCONV','LCOND',
     &                     'LGWDRAG','LDIUR','LSOLC','LAER','LCFC',
     &                     'LGADSRH','LAEROZ'
         WRITE(*,'(12L8)') LECRAD,LVDIFF,LSURF,LCONV,LCOND,
     &                     LGWDRAG,LDIUR,LSOLC,LAER,LCFC,
     &                     LGADSRH,LAEROZ
      ENDIF
C
C LOOP OVER ALL TASKS (TO BE MULTI-TASKED)
C
C!CDIR CNCALL
C!CDIR CONCUR( BY=1 )
CJ-PP IF NO HALOS, LOOPING OVER HALO FREE REGIONS
C     IF HALOS INCLUDED, LOOPING OVER THE NHORP VECTOR LENGTHS
      IF(LNOHALO) THEN
         IPHYSTART = JAHCOMP
         IPHYSTOP = JEHCOMP
      ELSE
         IPHYSTART = 1
         IPHYSTOP = NPTASK
      ENDIF
C
C     PHYSICS LOOP
      DO KTASK=IPHYSTART,IPHYSTOP
C
CJ-PP IF NOT, VECTOR LENGTH ONE ROW MINUS HALOS
C     IF HALOS, NORMAL LOOPING
         IF(LNOHALO) THEN
            KB = (KTASK-1)*IE+IAHCOMP
            KS = (KTASK-1)*IE+IEHCOMP
         ELSE
            KB = (KTASK-1)*NHORPH+1
            IF(KTASK.EQ.NPTASK) THEN
               KS = IEJE
            ELSE
               KS = KTASK*NHORPH
            ENDIF
         ENDIF
C
C     LENGTH OF THE SUBVECTOR
         NPOINT = KS-KB+1
C
         CALL PHYORG(
              !
              ! WARNING!! Here the actual parameter names differ from
              ! the names of the dummy arguments in phyorg
              !
              ! Here, one should only pass the range of the subector (KB,KS)
              ! and not the whole vectors themselves. They should be moved
              ! to a module from which PHYORG can work.
              !
     &        NZT    , NPOINT  , KE     , KE1     , ZCONACC, VCT     ,
     &        NRADIA , CEVAPCU , CVDAES , CVDAEL  , CVDAEU ,
     &        CVDAED , ZDT     , 2.*ZDT , EPS     , NCBASE ,
     &        NTBASE , NTIMST  , NMONTH  , LVDIFF , LSURF   ,
     &        LGWDRAG, LCONV   , LCOND  , NRADFR  , LECRAD  ,
     &        LSOLC   , LAER   , LCFC    , LSICED , IEXC    ,
     &        NOZ    , LGADSRH , LAEROZ , L5LAY   , LWDIF  ,
C38
C             -- PROGNOSTIC FIELDS AND TENDENCIES --
     &        ZU    (KB:KS,:)   , ZV    (KB:KS,:) , T     (KB:KS,:,NA),
     &        QD    (KB:KS,:,NA), ZDQDT (KB:KS,:) , QW    (KB:KS,:,NA),
     &        ZTRAD (KB:KS,:)   , ZTGWDR(KB:KS,:) , ZUGWDR(KB:KS,:)   ,
     &        ZVGWDR(KB:KS,:)   , ZTTS  (KB:KS,:) , ZQDTS (KB:KS,:)   ,
     &        ZQWTS (KB:KS,:)   , ZTTK  (KB:KS,:) , ZUTK  (KB:KS,:)   ,
     &        ZVTK  (KB:KS,:)   , ZQDTK (KB:KS,:) ,
C             -- PRESSURE FIELDS AND COORDINATES --
     &        PS    (KB:KS,NA)  , SINLAT(KB:KS)   ,
     &        COSLAT(KB:KS)     , SINLON(KB:KS)   , COSLON(KB:KS)     ,
C             -- SKIN-TEMPERATURE RELATED FIELDS --
     &        SICED (KB:KS)     , TEFF  (KB:KS)   ,
C             -- SURF TEMP AND MOIST FIELDS --
     &        TSECH (KB:KS,NE)  , TSECH (KB:KS,NJ), TSECH (KB:KS,NA)  ,
     &        TSLECH(KB:KS,NE)  , TSLECH(KB:KS,NJ), TSLECH(KB:KS,NA)  ,
     &        TSWECH(KB:KS,NE)  , TSWECH(KB:KS,NJ), TSWECH(KB:KS,NA)  ,
     &        TSIECH(KB:KS,NE)  , TSIECH(KB:KS,NJ), TSIECH(KB:KS,NA)  ,
     &        WSECH (KB:KS,NE)  , WSECH (KB:KS,NJ), WSECH (KB:KS,NA)  ,
     &        SN    (KB:KS,NE)  , SN    (KB:KS,NJ), SN    (KB:KS,NA)  ,
     &        WL    (KB:KS,NE)  , WL    (KB:KS,NJ), WL    (KB:KS,NA)  ,
     &        TD    (KB:KS,NE)  , TD    (KB:KS,NJ), TD    (KB:KS,NA)  ,
C30
C             -- SURFACE FLUXES AND FLUX DERIVATIVES --
     &        SRFL  (KB:KS)   , THFL  (KB:KS)   , QHFL  (KB:KS)   ,
     &        XHFL  (KB:KS)   , RSFC  (KB:KS)   , SSFC  (KB:KS)   ,
     &        RSFL  (KB:KS)   , SSFL  (KB:KS)   , AHFS  (KB:KS)   ,
     &        AHFSL (KB:KS)   , AHFSW (KB:KS)   , AHFSI (KB:KS)   ,
     &        AHFL  (KB:KS)   , AHFICE(KB:KS)   , QDB   (KB:KS,NA),
     &        QDBL  (KB:KS,NA), QDBW  (KB:KS,NA), QDBI  (KB:KS,NA),
     &        DHFT  (KB:KS)   , DHFQW (KB:KS)   , DHFQS (KB:KS)   ,
C             -- FIELDS USED IN DEEP CONVECTION --
     &        TOPMAX(KB:KS)   , VERVEL(KB:KS,:) , AZ0   (KB:KS)   ,
     &        AZ0L  (KB:KS)   , AZ0W  (KB:KS)   , AZ0I  (KB:KS)   ,
C             -- RAIN AND SNOW FALL, EVAPORATION AND ALBEDO --
     &        APRC  (KB:KS)   , APRL  (KB:KS)   , APRS  (KB:KS)   ,
     &        EVAP  (KB:KS)   , EVAPM (KB:KS)   , EVAPL (KB:KS)   ,
     &        EVAPW (KB:KS)   , EVAPI (KB:KS)   , ACLC  (KB:KS,:) ,
     &        ACLCAC(KB:KS,:) , VGRAT (KB:KS)   , FOREST(KB:KS)   ,
     &        ALBECH(KB:KS)   , ALBEDO(KB:KS)   , ALSOL (KB:KS)   ,
     &        ALSOW (KB:KS)   , ALSOI (KB:KS)   ,
C             -- REMAINING ELEMENTS IN *VDIFF* --
     &        DEW2  (KB:KS)   ,
     &        TKE   (KB:KS,:,NE),TKE (KB:KS,:,NJ),TKE(KB:KS,:,NA) ,
     &        WSMX  (KB:KS)   , VLT   (KB:KS)   ,
C49
     &        RGCGN (KB:KS)   ,
     &        TEMP2 (KB:KS)   , T2MAX (KB:KS)   , T2MIN (KB:KS)   ,
     &        USTAR3(KB:KS)   , USTR  (KB:KS)   , U10   (KB:KS)   ,
     &        VDIS  (KB:KS)   , VSTR  (KB:KS)   , V10   (KB:KS)   ,
     &        WIND10(KB:KS)   , WIMAX (KB:KS)   , USTRL (KB:KS)   ,
     &        USTRW (KB:KS)   , USTRI (KB:KS)   , VSTRL (KB:KS)   ,
     &        VSTRW (KB:KS)   , VSTRI (KB:KS)   , INFRL (KB:KS)   ,
     &        INFRW (KB:KS)   , INFRI (KB:KS)   , ZTMCM (KB:KS,:) ,
     &        ZTMCH (KB:KS,:) , TMCHL (KB:KS)   , TMCHW (KB:KS)   ,
     &        TMCHI (KB:KS)   ,
C             -- REMAINING ELEMENTS --
     &        ACLCOV(KB:KS)   , ALWCVI(KB:KS)   ,
     &        QVI   (KB:KS)   , EMTER (KB:KS,:) , TRSOL (KB:KS,:) ,
     &        EMTEF (KB:KS,:) , TRSOF (KB:KS,:) , SCLF0 (KB:KS)   ,
     &        SCLFS (KB:KS)   , SRAF0 (KB:KS)   , SRAFS (KB:KS)   ,
     &        TCLF0 (KB:KS)   , TCLFS (KB:KS)   , TRAF0 (KB:KS)   ,
     &        TRAFS (KB:KS)   , ACLCV (KB:KS)   , SRADS (KB:KS)   ,
     &        SRADSU(KB:KS)   , SRAD0 (KB:KS)   , SRAD0U(KB:KS)   ,
     &        TRADS (KB:KS)   , TRADSU(KB:KS)   , TRAD0 (KB:KS)   ,
     &        DSNAC (KB:KS)   , RUNOFF(KB:KS)   , SNMEL (KB:KS)   ,
     &        TDCL  (KB:KS,NE), TDCL  (KB:KS,NJ), TDCL  (KB:KS,NA),
     &        TD3   (KB:KS,NE), TD3   (KB:KS,NJ), TD3   (KB:KS,NA),
     &        TD4   (KB:KS,NE), TD4   (KB:KS,NJ), TD4   (KB:KS,NA),
     &        TD5   (KB:KS,NE), TD5   (KB:KS,NJ), TD5   (KB:KS,NA),
     &        QRES  (KB:KS)   , TSLIN (KB:KS)   , TSMAX (KB:KS)   ,
     &        TSMIN (KB:KS)   , TSN   (KB:KS,NE), TSN   (KB:KS,NJ),
     &        TSN   (KB:KS,NA), TSURF (KB:KS)   , VAROR (KB:KS)   ,
     &        DRAIN (KB:KS)   ,
CC48
C             -- REMAINING ELEMENTS IN *GWDRAG* --
     &        USTRGW(KB:KS)   , VAR   (KB:KS,:) ,
     &        VDISGW(KB:KS)   , VSTRGW(KB:KS)   ,
C             -- SOIL PARAMETER --
     &        TLAMBDA(KB:KS)  ,
     &        DLAMBDA(KB:KS)  , PORVOL(KB:KS)   , FCAP  (KB:KS)   ,
     &        WI3   (KB:KS,NE), WI3   (KB:KS,NJ), WI3   (KB:KS,NA),
     &        WI4   (KB:KS,NE), WI4   (KB:KS,NJ), WI4   (KB:KS,NA),
     &        WI5   (KB:KS,NE), WI5   (KB:KS,NJ), WI5   (KB:KS,NA),
     &        WI    (KB:KS,NE), WI    (KB:KS,NJ), WI    (KB:KS,NA),
     &        WICL  (KB:KS,NE), WICL  (KB:KS,NJ), WICL  (KB:KS,NA),
     &        LOGLAC(KB:KS)   , LOLAND(KB:KS)   , LOSEA (KB:KS)   ,
     &        LOICE (KB:KS)   , LALAND(KB:KS)   , GHPBL (KB:KS)   ,
     &        BETA  (KB:KS)   , WMINLOK(KB:KS)  , WMAXLOK(KB:KS)  ,
     &        CAPE  (KB:KS)   , OZONPL(KB:KS,:) , SO4ALL(KB:KS,:) ,
     &        SO4NAT(KB:KS,:) ,
C             -- 5 LAYER SOIL SCHEME: MOISTURE LAYERS AND PARAMETERS --
     &        WS1   (KB:KS)   , WS2   (KB:KS)   ,
     &        WS3   (KB:KS)   , WS4   (KB:KS)   , WS5   (KB:KS)   ,
     &        DZR   (KB:KS)   , DZS   (KB:KS)   , FKSAT (KB:KS)   ,
     &        FMPOT (KB:KS)   , BCLAPP(KB:KS)   , VPOR  (KB:KS)   ,
C             -- EVAPORATION FLUXES --
     &        ETRANS(KB:KS)   , EBSOIL(KB:KS)   ,
     &        ESNOW (KB:KS)   , ESKIN (KB:KS)   , ERES  (KB:KS)   ,
C             -- CLOUD ICE --
     &        QI  (KB:KS,:,NA), ZQITS (KB:KS,:) , QIVI  (KB:KS)   ,
     &        RPRAC (KB:KS,:)  )
C
      ENDDO                     ! KTASK=IPHYSTART,IPHYSTOP
C
CKS
C     TURBULENT EXCHANGE COEFFICIENTS AND WIND TENDENCIES ARE USED IN THE
C     DYNAMICS (PROGEXP) AND ARE NEEDED ALSO FROM THE NEIGHBORING GRID BOXES.
C     THATS WHY THE HALOS NEED TO BE UPDATED (WHICH ARE NOT CONSISTENT ANYMORE
C     DUE TO ECKONT/ECTIED CALL).
C
      ITYPE = 60
      KSIZE(1:10) = (/KE,KE,KE,KE,0,0,0,0,0,0/)
      CALL HALOEXCH(KSIZE, ITYPE, ZTMCM, ZTMCH, ZUTK, ZVTK)
CKS
      Z2ZM(:,1:KE-1) = ZTMCM(:,1:KE-1)
      Z2ZH(:,1:KE-1) = ZTMCH(:,1:KE-1)
      TMCM(:)        = ZTMCM(:,KE)
      TMCH(:)        = ZTMCH(:,KE)
C
C!CDIR CONCUR
C!CDIR CNCALL
      APRL(:)=MAX(APRL(:),0.)
      APRC(:)=MAX(APRC(:),0.)
      APRS(:)=MAX(APRS(:),0.)
      WL(:,NE)=MAX(WL(:,NE),0.)
C
C
C
C
CDJ CHANGES FOR REMO1.0
C MODIFICATION OF THE DIFFUSION COEFFICIENTS IN ORDER TO ADJUST THE
C ECHAM 4 PHYSICS FORMULATION TO THE DWD DYNAMICS STRUCTURE (EG.PROGEXP)
C!CDIR CONCUR
      DO K = KE,2,-1
         DO IJ = 1,IEJE
            ZDPU      = DAK(K  ) + DBK(K  )*PS(IJ,NA)
            ZDPO      = DAK(K-1) + DBK(K-1)*PS(IJ,NA)
            ZDPN2     = ZDPU + ZDPO
            ZTVU      = T(IJ,K,NA) * ( 1.0 + RDDRM1*QD(IJ,K,NA) )
            ZTVO      = T(IJ,K-1,NA) * ( 1.0 + RDDRM1*QD(IJ,K-1,NA) )
            ZTN       = (ZDPU*ZTVO + ZDPO*ZTVU)/ZDPN2
            ZRHO      = (AK(K) + BK(K)*PS(IJ,NA))/(R*ZTN)
            ZGROP     = 2.0*G*G*ZRHO*ZRHO/ZDPN2
            Z2ZM(IJ,K-1)= Z2ZM(IJ,K-1)*ZGROP
            Z2ZH(IJ,K-1)= Z2ZH(IJ,K-1)*ZGROP
         ENDDO
      ENDDO
CDJ 11.10.94
C
C!CDIR CONCUR
      DO J=1,JE
         ICOUNT = 0
         DO K=1,KE-1
            DO I=1,IE
               ICOUNT = ICOUNT+1
               IJJ = (J-1)*IE+I
               ZW7(ICOUNT) =  Z2ZM(IJJ,K)
               ZW8(ICOUNT) = Z2ZH(IJJ,K)
            ENDDO
         ENDDO
         TMKVMH(:,J,1) = ZW7(:)
         TMKVMH(:,J,2) = ZW8(:)
      ENDDO
C!CDIR CONCUR
C!CDIR CNCALL
      DO J=1,JE
C
        ICOUNT = 0
        DO K=1,KE
          DO I=1,IE
            ICOUNT = ICOUNT+1
            IJJ = (J-1)*IE+I
            ZW3(ICOUNT) = ZUTK(IJJ,K)
            ZW4(ICOUNT) = ZVTK(IJJ,K)
            ZW5(ICOUNT) = ZTRAD(IJJ,K)
            ZW9(ICOUNT) = ZTTK(IJJ,K)
            ZW10(ICOUNT)= ZQDTK(IJJ,K)
            ZW11(ICOUNT)= ZTTS(IJJ,K)
            ZW12(ICOUNT)= ZQDTS(IJJ,K)
            ZW13(ICOUNT)= ZQWTS(IJJ,K)
            ZW14(ICOUNT)= ZQITS(IJJ,K)
          ENDDO
        ENDDO
C
        SOTHDT(:,J,1) = ZW5(:)
        SOTHDT(:,J,2) = 0.0     ! OLD ZTRAD0 VARIABLE
        UVTK(:,J,1)   = ZW3(:)
        UVTK(:,J,2)   = ZW4(:)
        TTK(:,J)      = ZW9(:)
        QDTK(:,J)     = ZW10(:)
        TTS(:,J)      = ZW11(:)
        QDTS(:,J)     = ZW12(:) 
        QWTS(:,J)     = ZW13(:)
        QITS(:,J)     = ZW14(:)
      ENDDO
C
CRP   BELEGEN VON QDB (FUERS TESTEN QDB=QD(KE))
C     UND TS BZW. TB AUF DEN ZEITEBENEN
C
      QDB(:,NE) = QDB(:,NA)
      QDB(:,NA) = ZQDBNA(:)
      QDBL(:,NE) = QDBL(:,NA)
      QDBL(:,NA) = ZQDBLNA(:)
      QDBW(:,NE) = QDBW(:,NA)
      QDBW(:,NA) = ZQDBWNA(:)
      QDBI(:,NE) = QDBI(:,NA)
      QDBI(:,NA) = ZQDBINA(:)
      TS(:,NA) = TSECH(:,NA)
      TS(:,NJ) = TSECH(:,NJ)
      TS(:,NE) = TSECH(:,NE)
      TB(:,NA) = TSECH(:,NA)
      TB(:,NJ) = TSECH(:,NJ)
      TB(:,NE) = TSECH(:,NE)
      TG(:,NA) = TSECH(:,NA)
      TG(:,NJ) = TSECH(:,NJ)
      TG(:,NE) = TSECH(:,NE)
      TGL(:,NA) = TSLECH(:,NA)
      TGL(:,NJ) = TSLECH(:,NJ)
      TGL(:,NE) = TSLECH(:,NE)
      TGW(:,NA) = TSWECH(:,NA)
      TGW(:,NJ) = TSWECH(:,NJ)
      TGW(:,NE) = TSWECH(:,NE)
      TGI(:,NA) = TSIECH(:,NA)
      TGI(:,NJ) = TSIECH(:,NJ)
      TGI(:,NE) = TSIECH(:,NE)
C
      RETURN
C
      END
