SUBROUTINE progexp                              &

! Code converted using TO_F90 by Alan Miller
! Date: 2015-09-07  Time: 11:10:31

!
! WARNING! Here are 3 dummy argument names that
! differ from actual parameter names.
! The rest if fine.
!
!  OM850M(1), OM500M(1), OM300M(1),
!     |        |        |
!     v        v        v  &
(zom850m, zom500m, zom300m,  &
    ak    , bk    , akh    , bkh    , dak    , dbk  , a1t   ,  &
    a2t   , vvfh  , fib    , fc     , acphir , cphi , ps    ,  &
    tgl   , tgw   , tgi    , qdbl   , qdbw   , qdbi , bflhsl,  &
    bflhsw, bflhsi, bflqdsl, bflqdsw, bflqdsi, tg   , qdb   ,  &
    bflhs , bflqds, bflus  , bflvs  , u      , v    , t     ,  &
    qd    , qw    , fi     , tmkvmh , tmchl  , tmchw,  &
    tmchi , tmcm  , tmch   , sothdt , ttk    , qdtk , uvtk  ,  &
    tts   , qdts  , qwts   , infrl  , infrw  , infri, qi    ,  &
    qits  , pint  , dwdt   , etas)

USE RuntimeParameters
USE Grid
USE Hydro

!**** PROGEXP  -   UP: PROGNOSE FUER EINEN ZEITSCHRITT
!**   AUFRUF   :   CALL PROGEXP ( JATPROG, JETPROG,
!**               1               ZOM850M, ZOM500M, ZOM300M )
!**   ENTRIES  :      ---
!**   ZWECK    :   BERECHNUNG ALLER PROGNOSEVARIABLEN FUER
!**                DEN ZEITPUNKT NZT+1 (EXPLIZITE PROGNOSE)
!**   VERSIONS-
!**   DATUM    :   22.7.1991
!**
!**   EXTERNALS:   GAUSS  , HQTOTQ , COPYRE
!**
!**   EINGABE-
!**   PARAMETER:   JATPROG: ANFANGS-J-INDEX DES PROGNOSEBEREICHS;
!**                JETPROG: END    -J-INDEX DES PROGNOSEBEREICHS;
!**                         DIE FESTLEGUNG ERFOLGT IN UP *PROGORG*.
!**   AUSGABE-
!**   PARAMETER    ZOM850M: OMEGA-MITTELWERT IN ETWA 850 HPA (KONTROLLE)
!**                ZOM500M: OMEGA-MITTELWERT IN ETWA 500 HPA (KONTROLLE)
!**                ZOM300M: OMEGA-MITTELWERT IN ETWA 300 HPA (KONTROLLE)
!**
!**   COMMON-
!**   BLOECKE  :   PARAM  , ORG, COMDYN, COMPYH, PHYKON, HIGKON, PARKON
!**                PROGCHK, COMDIA, COMPCST, COMPMLF, COMPSLF, COMPGP3
!**
!**   METHODE  :   ZEITLICH:  SEMI-IMPLIZIT/LEAP-FROG
!**                RAEUMLICH: FINITE DIFFERENZEN 2.ORDNUNG
!**                           ARAKAWA-C-GITTER (HORIZONTAL)
!**                           ETA-KOORDINATE (VERTIKAL)
!**   FEHLERBE-
!**   HANDLUNG :      ---
!**   VERFASSER:   G. DOMS UND D. MAJEWSKI


IMPLICIT NONE

!      INCLUDE "parorg.h"
!      INCLUDE "org.h"
!      INCLUDE "comdyn.h"
!      INCLUDE "phykon.h"
!      INCLUDE "higkon.h"
!      INCLUDE "parkon.h"
!      INCLUDE "progchk.h"
!      INCLUDE "comdia.h"
!      INCLUDE "faktinf.h"



!-----------------------------------------------------------------------

!     ERFORDERLICHE EM-FELDER AUS DEM LANGZEITSPEICHER DIMENSIONIEREN
!     ---------------------------------------------------------------
!     VERTIKAL-KOORDINATEN-PARAMETER
REAL :: ak(ke1), bk(ke1), akh(ke), bkh(ke), dak(ke), dbk(ke)

!     VERTIKAL VARIIERENDER IMPLIZITHEITSGRAD DER V-DIFFUSION
REAL :: a1t(ke1), a2t(ke1)

!     EVT. VERTIKAL VARIIERENDER VERSTAERKUNGSFAKTOR DER H-DIFFUSION
REAL :: vvfh(ke)

!     EXTERNE PARAMETER
REAL :: fib(ie,je), fc(ie,je), acphir(je,2), cphi(je,2)

!     PROGNOSTISCHE BODENFELDER

REAL :: ps(ie,je,3),  &
    tg (ie,je,3), tgl (ie,je,3), tgw (ie,je,3), tgi (ie,je,3),  &
    qdb(ie,je,3), qdbl(ie,je,3), qdbw(ie,je,3), qdbi(ie,je,3)
!     DIAGNOSTISCHE BODENFELDER
REAL ::   bflhs  (ie,je),  &
    bflhsl (ie,je), bflhsw (ie,je), bflhsi (ie,je), bflqds (ie,je),  &
    bflqdsl(ie,je), bflqdsw(ie,je), bflqdsi(ie,je), bflus  (ie,je), bflvs  (ie,je)

INTEGER, intent(in) :: infrl  (ie,je), infrw  (ie,je), infri  (ie,je)

!     ATMOSPHAEREN-FELDER
REAL ::   u    (ie,je,ke,3), v (ie,je,ke,3),  &
    t    (ie,je,ke,3), qd(ie,je,ke,3), qw   (ie,je,ke,3), fi(ie,je,ke,2)

!     FELDER DER KOEFFIZIENTEN UND FLUESSE (PHYSIKALISCHE UPS)
REAL :: tmkvmh(ie*(ke-1),je,2)

REAL ::   tmcm(ie,je),  &
    tmch(ie,je), tmchl(ie,je), tmchw(ie,je), tmchi(ie,je)

REAL ::   sothdt(ieke,je,2)

REAL ::   ttk(ieke,je), qdtk(ieke,je), uvtk(ieke,je,2)

REAL ::   tts(ieke,je), qdts(ieke,je), qwts(ieke,je)

REAL, intent(in) :: dwdt(ie,je,ke ,3)
REAL, intent(inout) :: pint(ie,je,ke1,3)
REAL, intent(inout) :: etas(ie,je,ke1)

!     LOKALE FELDER DIMENSIONIEREN
!     ----------------------------
REAL  :: psdt(ie)
REAL  :: ztpa(ie,ke,5), zgqd (ie,ke,3), ztv (ie,ke,5), zfihf(ie,ke,5)

!     NEUE FELDER FUER GEWICHTETE SIGMA-HORIZONTALDIFFUSION:
REAL  ::   zphf(ie,ke,5), zpnf(ie,ke1,5), zdp(ie,ke,3),  &
    zdp5(ie,ke,5), zdpu(ie,ke ,3), zdpv(ie,ke,4)

REAL  ::   ztmkvm(ie,2:ke,3), ztmkvh(ie,2:ke,3),  &
    zttk  (ie,  ke  ), zqdtk (ie,  ke  ),  &
    zutk  (ie,  ke,2), zvtk  (ie,  ke,2),  &
    ztts  (ie,  ke  ), zqdts (ie,  ke  ),  &
    zqwts (ie,  ke  ), zsodta(ie,  ke  ),  &
    zthdta(ie,  ke  )

REAL  ::   zalopn(ie,ke1,3), zaloph(ie,ke ,3),  &
    zbetak(ie,ke ,3), zfih  (ie,ke ,3), zpphi (ie,ke ,2), zplam (ie,ke   ),  &
    zgu   (ie,ke ,2), zgv   (ie,ke ,3), zekin (ie,ke ,2), zzeta (ie,ke ,2),  &
    zetas (ie,ke1,2), zsdiv (ie,ke1,2)

REAL  ::   zlapt(ie,ke,3), zlapqd(ie,ke,3), zlapqw(ie,ke,3),  &
    zlapu(ie,ke,3), zlapv (ie,ke,3)

REAL  ::   ztadv (ie,ke), ztdifh(ie,ke), zqddih(ie,ke),  &
    zqwdih(ie,ke), zalpom(ie,ke), zqdadv(ie,ke)

REAL  ::   zgrad (ie,ke), zvzeq (ie,ke), zeddpq(ie,ke), zudifh(ie,ke),  &
    zvdifh(ie,ke), zuzeq (ie,ke)

REAL  ::   zqkor(ie), zfczet(ie)

!     FELDER FUER SUBROUTINE GAUSS
REAL  ::   aga(ie,ke,4), agb(ie,ke,4), agc(ie,ke,4), agd(ie,ke,4),  &
    age(ie,ke,4)

!     FELDER FUER SUBROUTINE HQTOTQ
REAL  ::   he    (ie,ke), qdwe  (ie,ke), phfe  (ie,ke),  &
    tstart(ie,ke), gqdsta(ie,ke), phfsta(ie,ke),qdle  (ie,ke),  &
    qdie  (ie,ke)


CHARACTER (LEN=26) :: ytxt

REAL  ::   qi(ie,je,ke,3), qits  (ieke,je),  &
    zqits(ie,ke  ), zlapqi(ie,ke,3), zqidih(ie,ke)


LOGICAL  :: lmassf
REAL ::  z1, z2, z3, z4, z4drerd, za1, za1a, za2, za2a,  &
    za3, zagam, zagat, zagcm, zagct, zalog2, zdeltf, zepsray, zfadvx
REAL ::  zfadvy, zfdivo, zfdivu, zfdivx, zfkio, zfkiu, zfvzo,  &
    zgvco, zgvcu, zom300m, zom500m, zom850m, zqd1, zqd2, zqd3,  &
    zqd4, zqw1, zqw2, zt1, zfvzu
REAL ::  zt2, zt3, zt4, ztkvl, ztkvz, ztmkvhm, ztrc, zx1, zx2, zxo,  &
    zxu, zzgew
INTEGER  :: i, index1, index2, index3, index4, j, jm2, jm3, jm4,  &
    jn5, jo2, jo3, jo4, jo5, js4, js5, jsp, ju3, ju4, jm5
INTEGER  :: ju5, jzm, jzn, jzo, jzs, jzu, k
INTEGER  :: km1, kp1

!     STATEMENT - FUNCTIONS ZUR BERECHNUNG DER SAEETIGUNGSFEUCHTE
!     TODO: Find out, why this macro statement gives slightly different
!           results than the function statement below.

REAL ::  fgew,fgqd,tt,GE,pp

!     MAGNUS-FORMEL FUER WASSER
fgew(tt)        = b1 * EXP  ( b2w*(tt - b3)/(tt - b4w) )
!     SPEZIFISCHE FEUCHTE AUS DAMPFDRUCK UND LUFTDRUCK
fgqd(GE,pp)     = rdrd*GE/(pp - emrdrd*GE)


!DIR$ NOTASK
!-----------------------------------------------------------------------
!     VORBEREITENDE MASSNAHMEN
index1 = 1
index2 = 2
index3 = 3
index4 = 4
!     PARAMETER ZUR BERECHNUNG DER DIFFUSION
!     UND LN(2) SETZEN.
alcnva  = 0.5
za1a    = alcnva
IF(za1a < 0.5) za1a = 0.5
za2a    = (1.-za1a)
zalog2  = ALOG(2.)
z4drerd = 4.0/rerd

!     WENN DIE MAXIMALE WINDGESCHWINDIGKEIT 95 % DES STABILITAETS-
!     KRITISCHEN WERTES UEBERSCHREITET, WIRD IN DER U- UND V-GLEICHUNG
!     RAYLEIGH-REIBUNG EINGEFUEHRT, UM DEN WIND ABZUBREMSEN.
IF ( vbmxv > 0.95*vbcfl ) THEN
  zepsray = 0.0005*ed2dt*(vbmxv - 0.95*vbcfl)/(0.05*vbcfl)
  IF ( vbmxv > vbcfl*1.05 ) THEN
    WRITE (ytxt,'(A,F5.1,A)') 'WARNING: VBMAX=',vbmxv,' M/S'
    CALL remark ( ytxt )
    IF (vbmxv > 250.0) THEN
      PRINT *,'VBMAX EXCEEDS 250 M/S'
      STOP 1
    END IF
  END IF
ELSE
  zepsray = 0.0
END IF

!     ZTKVZ , ZTKVL : GEWICHTSFAKTOREN ZUR HORIZONTALEN MITTELUNG DER
!                     VERTIKALEN TRANSPORTKOEFFIZIENTEN
ztkvz = 0.9
ztkvl = (1.-ztkvz)*0.25

!     OMEGA FLAECHENMITTELWERTE VORBEREITEN
zom850m = 0.
zom500m = 0.
zom300m = 0.

!     BERECHNUNG MIT ODER OHNE MASSENFLUSS-KORREKTURSCHEMA
lmassf = .true.

!     SETZEN DES STEUERPARAMETERS FUER TROCKENE (QD=QW=0) BZW. FEUCHTE
!     (QD>0, QW>0) PROGNOSE; WENN EIN ADIABATISCHER INITIALISIERUNGS-
!     ZEITSCHRITT GERECHNET WIRD, IST ZTRC=0.0, SONST ZTRC=1.0
IF(laistep) THEN
  ztrc = 0.0
ELSE
  ztrc = 1.0
END IF

!-----------------------------------------------------------------------

!     EXPLIZITE PROGNOSE OHNE RANDRELAXATION
!     DIE PROGNOSE ERFOLGT 'SCHEIBENWEISE' VON J = JATPROG BIS JETPROG.
!     IN EINEM TASK WERDEN JEWEILS ZUSAMMENHAENGENDE BEREICHE VON
!     J = JATPROG BIS JETPROG BEHANDELT; DIE TASK-EINTEILUNG UND
!     STEUERUNG ERFOLGT IM UP *PROGORG*; MAXIMAL KOENNEN VIER TASKS
!     PARALLEL LAUFEN; D.H. DER BEREICH JAH BIS JEH WIRD DANN IN VIER
!     UNTERBEREICHE GETEILT.

!***********************************************************************
!*                                                                     *
!*    ZUORDNUNG VON SCHEIBENINDIZES  'SUED'--------->-----------'NORD' *
!*    -----------------------------                                    *
!*                                          J-2 J-1      J      J+1 J+2*
!*    U, T, QD, QW, PHY.UP'S                     +       +       +     *
!*    V                                              +       +         *
!*                                                                     *
!*    ZALOPN(K+1) ALOG( P(K+1/2) )              JU3     JM3     JO3    *
!*    ZPHF(K)                               JS5 JU5     JM5     JO5 JN5*
!*    ZDP (K)                                   JU3     JM3     JO3    *
!*    ZPNF(K)                               JS5 JU5     JM5     JO5 JN5*
!*    ZTV (K)                               JS5 JU5     JM5     JO5 JN5*
!*    ZGQD(K)                                   JU3     JM3     JO3    *
!*    ZALOPH(K)   ALOG( P(K) )                  JU3     JM3     JO3    *
!*    ZBETAK(K)   HILFSGROESSE F. FIH, ALOPN    JU3     JM3     JO3    *
!*    ZFIHF(K)    FI(K) AN HAUPTFLAECHEN    JS5 JU3     JM3     JO3 JN5*
!*                FUER PHYS. PARAMETRISIER.                            *
!*    ZFIH(K)     FI(K) AN HAUPTFLAECHEN        JU3     JM3     JO3    *
!*                FUER DRUCKGRADIENTTERM                               *
!*    ZGU(K)      .5*(DP(I+1)+DP(I))*U                  JM2     JO2    *
!*    ZGV(K)      .5*(DP(J+1)+DP(J))*V              JU3     JM3     JO3*
!*    ZZETA(K)    POTENTIELLE VORTICITY             JM2     JO2        *
!*    ZEKIN(K)    KINETISCHE ENERGIE                    JM2     JO2    *
!*    ZPPHI(K)    R*TV*GRADY(LN(P))                 JM2     JO2        *
!*    ZSDIV(K+1)  VERT. SUMME DER DIVERGENZEN           JM2     JO2    *
!*    ZETAS(K+1)  MODIF. VERTIKALGESCHWINDIGK.          JM2     JO2    *
!*    ZSKHH(K)    DP*SKH/DLAM**2 AM H-GP.               JM2     JO2    *
!*    ZSKHU(K)    DP*SKH/DLAM**2 AM U-GP.               JM2     JO2    *
!*    ZSKHV(K)    DP*SKH/DLAM**2 AM V-GP.           JM2     JO2        *
!*    ZSKHZ(K)    DP*SKH/DLAM**2 AM Z-GP.           JM2     JO2        *
!*    ZHJ(K)      H(NJ)                         JU3     JM3     JO3    *
!*    ZQDWJ(K)    QDW(NJ)                       JU3     JM3     JO3    *
!*    ZTMKVM(K)   TMKVM(K)                      JU3     JM3     JO3    *
!*    ZTMKVH(K)   TMKVH(K)                      JU3     JM3     JO3    *
!*    ZUTK(K)     U-TENDENZ (KONVEKTIV)                 JM2     JO2    *
!*    ZVTK(K)     V-TENDENZ (KONVEKTIV)                 JM2     JO2    *
!*                                                                     *
!*    ZPLAM(K)    R*TV*GRADX(LN(P))                      +             *
!*    ZALPOM(K)   ALPHA*OMEGA                            +             *
!*    ZTPA (K)    HP (NA)                  JS5  JU5     JM5     JO5 JN5*
!*    ZQDWA(K)    QDW(NA)                  JS5  JU5     JM5     JO5 JN5*
!*    ZPLAM(K)    DRUCKGRADIENT                          +             *
!*    SOWIE WEITERE LOKALE FELDER OHNE SCHEI-            +             *
!*    BENINDEX SIND IN DER J-FLAECHE DEFINIERT           +             *
!*                                                                     *
!*                                         J-2  J-1      J      J+1 J+2*
!*                                                                     *
!***********************************************************************

!     VORBESETZUNG LOKALER FELDER AM OBER- UND UNTERRAND
!     --------------------------------------------------

DO i  = iaa , iea
  
  zpnf  (i,1,1) = 0.
  zpnf  (i,1,2) = 0.
  zpnf  (i,1,3) = 0.
  zpnf  (i,1,4) = 0.
  zpnf  (i,1,5) = 0.
  zalopn(i,1,1) = 0.
  zalopn(i,1,2) = 0.
  zalopn(i,1,3) = 0.
  zbetak(i,1,1) = zalog2
  zbetak(i,1,2) = zalog2
  zbetak(i,1,3) = zalog2
  zsdiv (i,1,1) = 0.
  zsdiv (i,1,2) = 0.
  zetas (i,1,1) = 0.
  zetas (i,1,2) = 0.
  zetas(i,ke1,1)= 0.
  zetas(i,ke1,2)= 0.
  aga(i,1,1)    = 0.
  aga(i,1,2)    = 0.
  aga(i,1,3)    = 0.
  agc(i,ke,1)   = 0.
  agc(i,ke,2)   = 0.
  agc(i,ke,3)   = 0.
  aga(i,1,4)    = 0.
  agc(i,ke,4)   = 0.
END DO

etas(iaa:iea,:,1:ke1) = 0.

!     ANFANGSINDIZES FUER ZYKLISCHES UMSPEICHERN SETZEN

jm2 = 1
jo2 = 2

ju3 = 1
jm3 = 2
jo3 = 3

js4 = 1
ju4 = 2
jm4 = 3
jo4 = 4

js5 = 1
ju5 = 2
jm5 = 3
jo5 = 4
jn5 = 5

!     LOKALE FELDER AM 'SUEDRAND' (J=JATPROG, JATPROG-1) VORBESETZEN
!     --------------------------------------------------------------

!     ACHTUNG:
!     UM BEIM AUTOTASKING (PARALLER DURCHLAUF VON UP *PROGEXP*) IDENTI-
!     SCHE ERGEBNISSE ZU ERZIELEN, MUSS DARAUF GEACHTET WERDEN, DASS
!     DIE FORMULIERUNG DER TERME FUER DEN 'SUEDRAND' UND FUER DAS INNE-
!     RE DES PROGNOSEGEBIETES EXAKT IDENTISCH SIND. SONST KANN DURCH DIE
!     UNVERMEIDLICHEN RUNDUNGSFEHLER DAS ERGEBNIS JE NACH ANZAHL DER
!     TASKS VERSCHIEDEN SEIN.

jzs = jah - 2
jzu = jah - 1
jzm = jah
jzo = jah + 1

!     GEPACKTE FELDER DER PHYSIKALISCHEN TENDENZEN UND KOEFFIZIENTEN
!     AUSPACKEN
CALL copyre(tmkvmh(1,jzu,1),ztmkvm(1,2,ju3),ie*(ke-1))
CALL copyre(tmkvmh(1,jzu,2),ztmkvh(1,2,ju3),ie*(ke-1))
CALL copyre(tmkvmh(1,jzm,1),ztmkvm(1,2,jm3),ie*(ke-1))
CALL copyre(tmkvmh(1,jzm,2),ztmkvh(1,2,jm3),ie*(ke-1))
CALL copyre(uvtk  (1,jzm,1),zutk  (1,1,jm2),ieke)
CALL copyre(uvtk  (1,jzm,2),zvtk  (1,1,jm2),ieke)

zx2 = .5*rerd*(cphi(jzu,1)+cphi(jzm,1))
zxo = cphi(jzm,1)*eddphi
zxu = cphi(jzu,1)*eddphi

DO k = 1 , ke
  kp1 = k + 1
  DO i  = iaa , iea
    zdp   (i,k  ,ju3) = dak(k  ) + dbk(k  )*ps(i,jzu,nj)
    zdp   (i,k  ,jm3) = dak(k  ) + dbk(k  )*ps(i,jzm,nj)
!     DEFINITION VON ZDP5 FUER GEWICHTETE SIGMA HDIFFUSION:
!==============================================================
!    !!!! Should ZDP5 be based on PINT?
!==============================================================
    zdp5  (i,k  ,js5) = dak(k  ) + dbk(k  )*ps(i,jzs,na)
    zdp5  (i,k  ,ju5) = dak(k  ) + dbk(k  )*ps(i,jzu,na)
    zdp5  (i,k  ,jm5) = dak(k  ) + dbk(k  )*ps(i,jzm,na)
    zdp5  (i,k  ,jo5) = dak(k  ) + dbk(k  )*ps(i,jzo,na)
!     *************
    zphf  (i,k  ,js5) = 0.5 * ( pint(i,jzs,k,nj) + pint(i,jzs,kp1,nj))
    zphf  (i,k  ,ju5) = 0.5 * ( pint(i,jzu,k,nj) + pint(i,jzu,kp1,nj))
    zphf  (i,k  ,jm5) = 0.5 * ( pint(i,jzm,k,nj) + pint(i,jzm,kp1,nj))
    zphf  (i,k  ,jo5) = 0.5 * ( pint(i,jzo,k,nj) + pint(i,jzo,kp1,nj))
    zpnf  (i,kp1,js5) = pint(i,jzs,kp1,nj)
    zpnf  (i,kp1,ju5) = pint(i,jzu,kp1,nj)
    zpnf  (i,kp1,jm5) = pint(i,jzm,kp1,nj)
    zpnf  (i,kp1,jo5) = pint(i,jzo,kp1,nj)
    
    zalopn(i,kp1,ju3) = ALOG(zpnf(i,kp1,ju5))
    zalopn(i,kp1,jm3) = ALOG(zpnf(i,kp1,jm5))
    ztv   (i,k  ,js5) =     t(i,jzs,k,nj) * ( 1.0 + rddrm1*  &
        qd(i,jzs,k,nj) - (qw(i,jzs,k,nj)+qi(i,jzs,k,nj)))
    ztv   (i,k  ,ju5) =     t(i,jzu,k,nj) * ( 1.0 + rddrm1*  &
        qd(i,jzu,k,nj) - (qw(i,jzu,k,nj)+qi(i,jzu,k,nj)) )
    ztv   (i,k  ,jm5) =     t(i,jzm,k,nj) * ( 1.0 + rddrm1*  &
        qd(i,jzm,k,nj) - (qw(i,jzm,k,nj)+qi(i,jzm,k,nj)) )
    ztv   (i,k  ,jo5) =     t(i,jzo,k,nj) * ( 1.0 + rddrm1*  &
        qd(i,jzo,k,nj) - (qw(i,jzo,k,nj)+qi(i,jzo,k,nj)) )
    zzgew             = fgew ( t(i,jzu,k,nj) )
    zgqd (i,k,ju3)    = fgqd ( zzgew, zphf(i,k,ju5) )
    zzgew             = fgew ( t(i,jzm,k,nj) )
    zgqd (i,k,jm3)    = fgqd ( zzgew, zphf(i,k,jm5) )
  END DO
END DO

!     DEFINITION VON ZDPU UND ZDPV FUER GEWICHTETE SIGMA HDIFFUSION
DO k = 1 , ke
!     FEHLERKORREKTUR FUER ZDPU-DEFINITION
  DO i = iaa , ieh+1
    zdpu (i,k,ju3) = 0.5*(zdp5(i,k,ju5) + zdp5(i+1,k,ju5))
    zdpu (i,k,jm3) = 0.5*(zdp5(i,k,jm5) + zdp5(i+1,k,jm5))
  END DO
END DO

!KS   SPLITTED THE FOLLOWING LOOP INTO TWO TO REMOVE THE IF INSIDE A LOOP

!KS   FIRST PART WITH FOR ALL K<KE
DO k = 1 , ke-1
  DO i = iaa , iea
    zdpv (i,k,js4) = 0.5*(zdp5(i,k,js5) + zdp5(i,k,ju5))
    zdpv (i,k,ju4) = 0.5*(zdp5(i,k,ju5) + zdp5(i,k,jm5))
    zdpv (i,k,jm4) = 0.5*(zdp5(i,k,jm5) + zdp5(i,k,jo5))
!           ***************
    zalpom(i,k) = 0.0
    zfihf(i,k,js5) = fi(i,jzs,k+1,na2) + r*ztv(i,k,js5)*  &
        ALOG  ( zpnf(i,k+1,js5)/zphf(i,k,js5) ) / dwdt(i,jzs,k,nj)
    zfihf(i,k,ju5) = fi(i,jzu,k+1,na2) + r*ztv(i,k,ju5)*  &
        ALOG  ( zpnf(i,k+1,ju5)/zphf(i,k,ju5) ) / dwdt(i,jzu,k,nj)
    zfihf(i,k,jm5) = fi(i,jzm,k+1,na2) + r*ztv(i,k,jm5)*  &
        ALOG  ( zpnf(i,k+1,jm5)/zphf(i,k,jm5) ) / dwdt(i,jzm,k,nj)
    zfihf(i,k,jo5) = fi(i,jzo,k+1,na2) + r*ztv(i,k,jo5)*  &
        ALOG  ( zpnf(i,k+1,jo5)/zphf(i,k,jo5) ) / dwdt(i,jzo,k,nj)
    ztpa (i,k,js5) = t(i,jzs,k,na) + zfihf(i,k,js5)*wcpr
    ztpa (i,k,ju5) = t(i,jzu,k,na) + zfihf(i,k,ju5)*wcpr
    ztpa (i,k,jm5) = t(i,jzm,k,na) + zfihf(i,k,jm5)*wcpr
    ztpa (i,k,jo5) = t(i,jzo,k,na) + zfihf(i,k,jo5)*wcpr
  END DO
END DO

!KS   SECOND PART FOR K=KE

DO i = iaa , iea
  zdpv (i,ke,js4) = 0.5*(zdp5(i,ke,js5) + zdp5(i,ke,ju5))
  zdpv (i,ke,ju4) = 0.5*(zdp5(i,ke,ju5) + zdp5(i,ke,jm5))
  zdpv (i,ke,jm4) = 0.5*(zdp5(i,ke,jm5) + zdp5(i,ke,jo5))
!        ***************
  zalpom(i,ke) = 0.0
  zfihf(i,ke,js5) = fib(i,jzs)        + r*ztv(i,ke,js5)*  &
      ALOG  ( zpnf(i,ke+1,js5)/zphf(i,ke,js5) ) / dwdt(i,jzs,ke,nj)
  zfihf(i,ke,ju5) = fib(i,jzu)        + r*ztv(i,ke,ju5)*  &
      ALOG  ( zpnf(i,ke+1,ju5)/zphf(i,ke,ju5) ) / dwdt(i,jzu,ke,nj)
  zfihf(i,ke,jm5) = fib(i,jzm)        + r*ztv(i,ke,jm5)*  &
      ALOG  ( zpnf(i,ke+1,jm5)/zphf(i,ke,jm5) ) / dwdt(i,jzm,ke,nj)
  zfihf(i,ke,jo5) = fib(i,jzo)        + r*ztv(i,ke,jo5)*  &
      ALOG  ( zpnf(i,ke+1,jo5)/zphf(i,ke,jo5) ) / dwdt(i,jzo,ke,nj)
  ztpa (i,ke,js5) = t(i,jzs,ke,na) + zfihf(i,ke,js5)*wcpr
  ztpa (i,ke,ju5) = t(i,jzu,ke,na) + zfihf(i,ke,ju5)*wcpr
  ztpa (i,ke,jm5) = t(i,jzm,ke,na) + zfihf(i,ke,jm5)*wcpr
  ztpa (i,ke,jo5) = t(i,jzo,ke,na) + zfihf(i,ke,jo5)*wcpr
END DO
!KS


!KS   SPLITTED THE FOLLOWING LOOP INTO TWO TO REMOVE THE IF INSIDE A LOOP

!KS   FIRST PART FOR K=1
DO i = iaa, iea
  zaloph(i,1 ,ju3) = zalopn(i,2,ju3) - zbetak(i,1 ,ju3)
  zaloph(i,1 ,jm3) = zalopn(i,2,jm3) - zbetak(i,1 ,jm3)
  zfih(i,1,ju3) = fi(i,jzu,2,nj2) +  &
      r*ztv(i,1,ju5)*zbetak(i,1,ju3) / dwdt(i,jzu,1,nj)
  zfih(i,1,jm3) = fi(i,jzm,2,nj2) +  &
      r*ztv(i,1,jm5)*zbetak(i,1,jm3) / dwdt(i,jzm,1,nj)
END DO

!KS   SECOND PART FOR K>1
DO k = 2, ke-1
  DO i = iaa , iea
    zbetak(i,k,ju3) = 1. - zpnf(i,k  ,ju5)/zdp   (i,k,ju3)  &
        * (zalopn(i,k+1,ju3)-zalopn(i,k,ju3))
    zbetak(i,k,jm3) = 1. - zpnf(i,k  ,jm5)/zdp   (i,k,jm3)  &
        * (zalopn(i,k+1,jm3)-zalopn(i,k,jm3))
    zfih(i,k  ,ju3) = fi(i,jzu,k+1,nj2) +  r*ztv(i,k,ju5)  &
        * zbetak(i,k,ju3) / dwdt(i,jzu,k,nj)
    zfih(i,k  ,jm3) = fi(i,jzm,k+1,nj2) +  r*ztv(i,k,jm5)  &
        *zbetak(i,k,jm3) / dwdt(i,jzm,k,nj)
    zaloph(i,k,ju3) = zalopn(i,k+1,ju3) - zbetak(i,k,ju3)
    zaloph(i,k,jm3) = zalopn(i,k+1,jm3) - zbetak(i,k,jm3)
  END DO
END DO
!KS

!     CORIOLISPARAMETER AM ZETA-PUNKT BEREITSTELLEN
DO i  = iaa, ieh + 1
  zfczet(i) = 0.25*( fc(i,jzu) + fc(i+1,jzu) + fc(i,jzm) + fc(i+1,jzm) )
END DO

DO i = iaa , iea
  zbetak(i,ke,ju3) = 1. - zpnf(i,ke ,ju5) / zdp   (i,ke,ju3)  &
      * (zalopn(i,ke1,ju3) - zalopn(i,ke,ju3))
  zbetak(i,ke,jm3) = 1. - zpnf(i,ke ,jm5) / zdp   (i,ke,jm3)  &
      * (zalopn(i,ke1,jm3) - zalopn(i,ke,jm3))
  zaloph(i,ke,ju3) =    zalopn(i,ke1,ju3) - zbetak(i,ke,ju3)
  zaloph(i,ke,jm3) =    zalopn(i,ke1,jm3) - zbetak(i,ke,jm3)
  zfih  (i,ke,ju3) = fib(i,jzu) +  &
      r*ztv(i,ke,ju5)*zbetak(i,ke,ju3) / dwdt(i,jzu,ke,nj)
  zfih  (i,ke,jm3) = fib(i,jzm) +  &
      r*ztv(i,ke,jm5)*zbetak(i,ke,jm3) / dwdt(i,jzm,ke,nj)
END DO

!     WEITERE LOKALE FELDER AN U/ZETA-GITTERPUNKTEN (HF) VORBESETZEN

DO k = 1, ke
  DO i = iaa , ieh + 1
    zgu(i,k,jm2) = .5*(zdp(i,k,jm3)+ zdp(i+1,k,jm3))*u(i,jzm,k,nj)
    z1 = zx2*zfczet(i)
    z2 = eddlam*( v(i+1,jzu,k,nj) -     v(i,jzu,k,nj) )
    z3 = zxo   *  u(i  ,jzm,k,nj) - zxu*u(i,jzu,k,nj)
    z4 = (zdp(i,k,ju3)+zdp(i+1,k,ju3))*cphi(jzu,1) +  &
        (zdp(i,k,jm3)+zdp(i+1,k,jm3))*cphi(jzm,1)
    zzeta(i,k,jm2) = z4drerd * ( z1 + z2 - z3 )/ z4
  END DO
  
!     WEITERE LOKALE FELDER AN H/V-GITTERPUNKTEN VORBESETZEN
  zx1    = 0.5*r*edadphi
  zfkio  = cphi(jzm,2)/cphi(jzm,1)
  zfkiu  = cphi(jzu,2)/cphi(jzm,1)
  zfdivx = acphir(jzm,1)*eddlam
  zfdivo = acphir(jzm,1)*eddphi*cphi(jzm,2)
  zfdivu = acphir(jzm,1)*eddphi*cphi(jzu,2)
  
  DO i = iah - 1 , ieh + 1
    zgv(i,k,ju3) = .5*(zdp(i,k,ju3)+(dak(k)+dbk(k)*ps(i,jzm  ,nj)))*  &
        v(i,jzu,k,nj)
    zgv(i,k,jm3) = .5*(zdp(i,k,jm3)+(dak(k)+dbk(k)*ps(i,jzm+1,nj)))*  &
        v(i,jzm,k,nj)
    zpphi(i,k,jm2) = zx1*(ztv(i,k,jm5) + ztv   (i,k,ju5))*  &
        (zaloph(i,k,jm3) - zaloph(i,k,ju3))
    zekin(i,k,jm2) = .25*(  u(i-1,jzm,k,nj)**2 + u(i,jzm,k,nj)**2 +  &
        zfkiu*v(i,jzu,k,nj)**2 + zfkio*v(i,jzm,k,nj)**2)
  END DO
END DO

DO k = 1 , ke
  DO i = iah - 1, ieh + 1
    zsdiv(i,k+1,jm2) = zsdiv(i,k,jm2) +  &
        zfdivx*( zgu(i,k,jm2)-zgu(i-1,k,jm2) ) +  &
        zfdivo*zgv(i,k,jm3) - zfdivu*zgv(i,k,ju3)
  END DO
END DO

DO k = 2 , ke
  DO i = iah - 1, ieh + 1
    zetas(i,k,jm2) = bk(k)*zsdiv(i,ke1,jm2) - zsdiv(i,k,jm2)
  END DO
  etas(iah-1:ieh+1,jzm,k)  = zetas(iah-1:ieh+1,k,jm2)
END DO

!     LOKALE HILFSFELDER FUER HORIZONTALDIFFUSION
DO k = 1 , ke
  DO i = iah - 1, ieh + 1
    
!     GEWICHTUNG VON ZLAPT, ZLAPQD, ZLAPQW MIT GEWICHT DER SCHICHT
    zlapt(i,k,ju3) = ((ztpa(i+1,k,ju5)-ztpa(i  ,k,ju5))*zdpu(i  ,k,ju3)-  &
        (ztpa(i  ,k,ju5)-ztpa(i-1,k,ju5))*zdpu(i-1,k,ju3)+  &
        (cphi(jzu,2)*(ztpa(i  ,k,jm5)-ztpa(i  ,k,ju5))*zdpv(i  ,k,ju4)-  &
        cphi(jzs,2)*(ztpa(i  ,k,ju5)-ztpa(i  ,k,js5))*zdpv(i  ,k,js4))/  &
        cphi(jzu,1)) / zdp5(i,k,ju5)
    zlapqd(i,k,ju3) = ((qd(i+1,jzu,k,na)-qd(i  ,jzu,k,na))*zdpu(i  ,k,ju3)-  &
        (qd(i  ,jzu,k,na)-qd(i-1,jzu,k,na))*zdpu(i-1,k,ju3)+ (cphi(jzu,2)*  &
        (qd(i  ,jzm,k,na)-qd(i  ,jzu,k,na))*zdpv(i  ,k,ju4)- cphi(jzs,2)*  &
        (qd(i  ,jzu,k,na)-qd(i  ,jzs,k,na))*zdpv(i  ,k,js4))/  &
        cphi(jzu,1)) / zdp5(i,k,ju5)
    zlapqw(i,k,ju3) = ((qw(i+1,jzu,k,na)-qw(i  ,jzu,k,na))*zdpu(i  ,k,ju3)-  &
        (qw(i  ,jzu,k,na)-qw(i-1,jzu,k,na))*zdpu(i-1,k,ju3)+ (cphi(jzu,2)*  &
        (qw(i  ,jzm,k,na)-qw(i  ,jzu,k,na))*zdpv(i  ,k,ju4)- cphi(jzs,2)*  &
        (qw(i  ,jzu,k,na)-qw(i  ,jzs,k,na))*zdpv(i  ,k,js4))/  &
        cphi(jzu,1)) / zdp5(i,k,ju5)
    zlapqi(i,k,ju3) = ((qi(i+1,jzu,k,na)-qi(i  ,jzu,k,na))*zdpu(i  ,k,ju3)-  &
        (qi(i  ,jzu,k,na)-qi(i-1,jzu,k,na))*zdpu(i-1,k,ju3)+ (cphi(jzu,2)*  &
        (qi(i  ,jzm,k,na)-qi(i  ,jzu,k,na))*zdpv(i  ,k,ju4)- cphi(jzs,2)*  &
        (qi(i  ,jzu,k,na)-qi(i  ,jzs,k,na))*zdpv(i  ,k,js4))/  &
        cphi(jzu,1)) / zdp5(i,k,ju5)
!     ***************
    zlapu (i,k,ju3) = u (i+1,jzu,k,na) + u (i-1,jzu,k,na) -  &
        2.0*u (i  ,jzu,k,na) +  &
        ( cphi(jzu,2)*(u (i  ,jzm,k,na) - u (i  ,jzu,k,na)) -  &
        cphi(jzs,2)*(u (i  ,jzu,k,na) - u (i  ,jzs,k,na)) )/ cphi(jzu,1)
    zlapv (i,k,ju3) = v (i+1,jzu,k,na) + v (i-1,jzu,k,na) -  &
        2.0*v (i  ,jzu,k,na) +  &
        ( cphi(jzm,1)*(v (i  ,jzm,k,na) - v (i  ,jzu,k,na)) -  &
        cphi(jzu,1)*(v (i  ,jzu,k,na) - v (i  ,jzs,k,na)) )/ cphi(jzu,2)
    
    zlapt(i,k,jm3) = ((ztpa(i+1,k,jm5)-ztpa(i  ,k,jm5))*zdpu(i  ,k,jm3)-  &
        (ztpa(i  ,k,jm5)-ztpa(i-1,k,jm5))*zdpu(i-1,k,jm3)+  &
        (cphi(jzm,2)*(ztpa(i  ,k,jo5)-ztpa(i  ,k,jm5))*zdpv(i  ,k,jm4)-  &
        cphi(jzu,2)*(ztpa(i  ,k,jm5)-ztpa(i  ,k,ju5))*zdpv(i  ,k,ju4))/  &
        cphi(jzm,1)) / zdp5(i,k,jm5)
    zlapqd(i,k,jm3) = ((qd(i+1,jzm,k,na)-qd(i  ,jzm,k,na))*zdpu(i  ,k,jm3)-  &
        (qd(i  ,jzm,k,na)-qd(i-1,jzm,k,na))*zdpu(i-1,k,jm3)+ (cphi(jzm,2)*  &
        (qd(i  ,jzo,k,na)-qd(i  ,jzm,k,na))*zdpv(i  ,k,jm4)- cphi(jzu,2)*  &
        (qd(i  ,jzm,k,na)-qd(i  ,jzu,k,na))*zdpv(i  ,k,ju4))/  &
        cphi(jzm,1)) / zdp5(i,k,jm5)
    zlapqw(i,k,jm3) = ((qw(i+1,jzm,k,na)-qw(i  ,jzm,k,na))*zdpu(i  ,k,jm3)-  &
        (qw(i  ,jzm,k,na)-qw(i-1,jzm,k,na))*zdpu(i-1,k,jm3)+ (cphi(jzm,2)*  &
        (qw(i  ,jzo,k,na)-qw(i  ,jzm,k,na))*zdpv(i  ,k,jm4)- cphi(jzu,2)*  &
        (qw(i  ,jzm,k,na)-qw(i  ,jzu,k,na))*zdpv(i  ,k,ju4))/  &
        cphi(jzm,1)) / zdp5(i,k,jm5)
    zlapqi(i,k,jm3) = ((qi(i+1,jzm,k,na)-qi(i  ,jzm,k,na))*zdpu(i  ,k,jm3)-  &
        (qi(i  ,jzm,k,na)-qi(i-1,jzm,k,na))*zdpu(i-1,k,jm3)+ (cphi(jzm,2)*  &
        (qi(i  ,jzo,k,na)-qi(i  ,jzm,k,na))*zdpv(i  ,k,jm4)- cphi(jzu,2)*  &
        (qi(i  ,jzm,k,na)-qi(i  ,jzu,k,na))*zdpv(i  ,k,ju4))/  &
        cphi(jzm,1)) / zdp5(i,k,jm5)
    zlapu (i,k,jm3) = u (i+1,jzm,k,na) + u (i-1,jzm,k,na) -  &
        2.0*u (i  ,jzm,k,na) +  &
        ( cphi(jzm,2)*(u (i  ,jzo,k,na) - u (i  ,jzm,k,na)) -  &
        cphi(jzu,2)*(u (i  ,jzm,k,na) - u (i  ,jzu,k,na)) )/ cphi(jzm,1)
    zlapv (i,k,jm3) = v (i+1,jzm,k,na) + v (i-1,jzm,k,na) -  &
        2.0*v (i  ,jzm,k,na) +  &
        ( cphi(jzo,1)*(v (i  ,jzo,k,na) - v (i  ,jzm,k,na)) -  &
        cphi(jzm,1)*(v (i  ,jzm,k,na) - v (i  ,jzu,k,na)) )/ cphi(jzm,2)
  END DO
END DO

!     RAENDER DER SCHEIBENFELDER FUER DIE HQ-ZERLEGUNG VORBESETZEN;
!     (BELIEBIGE KONSISTENTE WERTE)
!VDIR NOVECTOR
DO i  = 1 , 2
  DO k  = 1 , ke
    he    (i,k) = wcp*t(1,jzu,k,nj) + wlk*qd(1,jzu,k,nj)
    qdwe  (i,k) =    qd(1,jzu,k,nj) +     qw(1,jzu,k,nj) + qi(1,jzu,k,nj)
    tstart(i,k) =     t(1,jzu,k,nj)
    gqdsta(i,k) = zgqd (1,k,ju3)
    phfsta(i,k) = zphf (1,k,ju5)
    phfe  (i,k) = zphf (1,k,ju5)
  END DO
END DO

!VDIR NOVECTOR
DO i = ie - 1 , ie
  DO k = 1, ke
    he    (i,k) = he    (1,k)
    qdwe  (i,k) = qdwe  (1,k)
    tstart(i,k) = tstart(1,k)
    gqdsta(i,k) = gqdsta(1,k)
    phfsta(i,k) = phfsta(1,k)
    phfe  (i,k) = phfe  (1,k)
  END DO
END DO


!     BEGINN DER SCHLEIFE IN J-RICHTUNG
!     ---------------------------------

!     ACHTUNG:
!     UM BEIM AUTOTASKING (PARALLER DURCHLAUF VON UP *PROGEXP*) IDENTI-
!     SCHE ERGEBNISSE ZU ERZIELEN, MUSS DARAUF GEACHTET WERDEN, DASS
!     DIE FORMULIERUNG DER TERME FUER DEN 'SUEDRAND' UND FUER DAS INNE-
!     RE DES PROGNOSEGEBIETES EXAKT IDENTISCH SIND. SONST KANN DURCH DIE
!     UNVERMEIDLICHEN RUNDUNGSFEHLER DAS ERGEBNIS JE NACH ANZAHL DER
!     TASKS VERSCHIEDEN SEIN.

DO j = jah, jeh
  
  jzn = j+2
  
!        LOKALE FELDER FUER J   BESETZEN
!        GEPACKTE FELDER DER PHYSIKALISCHEN TENDENZEN UND KOEFFIZIENTEN
!        AUSPACKEN
  CALL copyre(ttk(1,j)   ,zttk  ,ieke)
  CALL copyre(qdtk(1,j)  ,zqdtk ,ieke)
  CALL copyre(tts(1,j)   ,ztts  ,ieke)
  CALL copyre(qdts(1,j)  ,zqdts ,ieke)
  CALL copyre(qwts(1,j)  ,zqwts ,ieke)
  CALL copyre(sothdt(1,j,1), zsodta, ieke)
  CALL copyre(sothdt(1,j,2), zthdta, ieke)
  CALL copyre(qits(1,j)  ,zqits , ieke)
!        LOKALE FELDER FUER J+1 BESETZEN
!        GEPACKTE FELDER DER PHYSIKALISCHEN TENDENZEN UND KOEFFIZIENTEN
!        AUSPACKEN
  CALL copyre(tmkvmh(1,j+1,1),ztmkvm(1,2,jo3),ie*(ke-1))
  CALL copyre(tmkvmh(1,j+1,2),ztmkvh(1,2,jo3),ie*(ke-1))
  CALL copyre(uvtk  (1,j+1,1),zutk  (1,1,jo2),ieke)
  CALL copyre(uvtk  (1,j+1,2),zvtk  (1,1,jo2),ieke)
!        LOKALE FELDER FUER J+2 BESETZEN
  
  zx1 = .5*r*acphir(j,1)*eddlam
  zx2 = .5*rerd*(cphi(j,1)+cphi(j+1,1))
  zxo = cphi(j+1,1)*eddphi
  zxu = cphi(j  ,1)*eddphi
  
  DO k = 1 , ke
    kp1 = k + 1
    DO i = iaa , iea
      zdp  (i,k,jo3)    = dak(k  ) + dbk(k  )*ps(i,j+1,nj)
!     FUER GEWICHTETE HDIFF
!==============================================================
!    !!!! Should ZDP5 be based on PINT?
!==============================================================
      zdp5 (i,k,jn5)    = dak(k  ) + dbk(k  )*ps(i,jzn,na)
      zdpv (i,k,jo4)  = 0.5*(zdp5(i,k,jo5) + zdp5(i,k,jn5))
!     ************
      zphf (i,k,jn5)    = 0.5 * ( pint(i,jzn,k,nj) + pint(i,jzn,kp1,nj))
      zpnf (i,kp1,jn5)  = pint(i,jzn,kp1,nj)
      zalopn(i,kp1,jo3) = ALOG(zpnf(i,kp1,jo5))
      ztv   (i,k  ,jn5) =     t(i,jzn,k,nj) * ( 1.0 + rddrm1*  &
          qd(i,jzn,k,nj) - (qw(i,jzn,k,nj) + qi(i,jzn,k,nj)) )
      zzgew             = fgew ( t(i,j+1,k,nj) )
      zgqd (i,k,jo3)    = fgqd ( zzgew, zphf(i,k,jo5) )
    END DO
  END DO
  
!        FUER GEWICHTETE HDIFF
  DO k = 1 , ke
!           FEHLERKORREKTUR FUER ZDPU
    DO i = iaa , ieh+1
      zdpu (i,k,jo3) = 0.5*(zdp5(i,k,jo5) + zdp5(i+1,k,jo5))
    END DO
  END DO
!        ***************
  
!KS      SPLITTED LOOP TO REMOVE IF INSIDE LOOP
  
!KS      FIRST PART K<KE
  DO k = 1 , ke-1
    DO i = iaa , iea
      zfihf(i,k,jn5) = fi(i,jzn,k+1,na2) + r*ztv(i,k,jn5)*  &
          ALOG  ( zpnf(i,k+1,jn5)/zphf(i,k,jn5) ) / dwdt(i,jzn,k,nj)
      ztpa (i,k,jn5)    = t(i,jzn,k,na) + zfihf(i,k,jn5)*wcpr
    END DO
  END DO
!KS      SECOND PART K=KE
  DO i = iaa , iea
    zfihf(i,ke,jn5) = fib(i,jzn)        + r*ztv(i,ke,jn5)*  &
        ALOG  ( zpnf(i,ke+1,jn5)/zphf(i,ke,jn5) ) / dwdt(i,jzn,ke,nj)
    ztpa (i,ke,jn5)    = t(i,jzn,ke,na) + zfihf(i,ke,jn5)*wcpr
  END DO
!KS
  
!KS      SPLITTED LOOP TO REMOVE IF INSIDE LOOP
  
!KS      FIRST PART K=1
  DO i  = iaa , iea
    zaloph(i,1,jo3) = zalopn(i,2,jo3) -zbetak(i,1,jo3)
    zfih (i,1,jo3) = fi(i,j+1,2,nj2) + r*ztv(i,1,jo5)  &
        *zbetak(i,1,jo3) / dwdt(i,j+1,1,nj)
  END DO
!KS      SECOND PART K>1
  DO k = 2 , ke-1
    DO i  = iaa , iea
      zbetak(i,k,jo3) = 1. - zpnf(i,k  ,jo5)/zdp   (i,k,jo3)  &
          * (zalopn(i,k+1,jo3)-zalopn(i,k,jo3))
      zfih  (i,k,jo3) = fi(i,j+1,k+1,nj2) + r*ztv (i,k,jo5)  &
          *zbetak(i,k,jo3) / dwdt(i,j+1,k,nj)
      zaloph(i,k,jo3) = zalopn(i,k+1,jo3) - zbetak(i,k,jo3)
    END DO
  END DO
  
!        CORIOLISPARAMETER AM ZETA-PUNKT BEREITSTELLEN
  DO i = iaa, ieh + 1
    zfczet(i) = 0.25*( fc(i,j  ) + fc(i+1,j  ) + fc(i,j+1) + fc(i+1,j+1) )
  END DO
  
  DO i = iaa , iea
    zbetak(i,ke,jo3) = 1. - zpnf(i,ke ,jo5)  /zdp   (i,ke,jo3)  &
        * (zalopn(i,ke1,jo3) - zalopn(i,ke,jo3))
    zaloph(i,ke,jo3) =    zalopn(i,ke1,jo3) - zbetak(i,ke,jo3)
    zfih  (i,ke,jo3) = fib(i,j+1) + r*ztv(i,ke,jo5)  &
        *zbetak(i,ke,jo3) / dwdt(i,j+1,ke,nj)
  END DO
  
!        WEITERE LOKALE FELDER AN U/ZETA-GITTERPUNKTEN (HF) BESETZEN , JO
  
  DO k = 1 , ke
    DO i = iaa , ieh + 1
      zgu(i,k,jo2) = .5*(zdp(i,k,jo3)+zdp(i+1,k,jo3)) *u(i,j+1,k,nj)
      zplam(i,k)   = zx1*( ztv(i+1,k,jm5) +    ztv(i,k,jm5) )  &
          *(zaloph(i+1,k,jm3) - zaloph(i,k,jm3))
      z1          = zx2*zfczet(i)
      z2          = eddlam*( v(i+1,j,k,nj) -     v(i,j,k,nj) )
      z3          = zxo   *  u(i,j+1,k,nj) - zxu*u(i,j,k,nj)
      z4          = (zdp(i,k,jm3)+zdp(i+1,k,jm3))*cphi(j,  1) +  &
          (zdp(i,k,jo3)+zdp(i+1,k,jo3))*cphi(j+1,1)
      zzeta(i,k,jo2) = z4drerd * ( z1 + z2 - z3 )/ z4
    END DO
  END DO
  
!        WEITERE LOKALE FELDER AN H/V-GITTERPUNKTEN (HF) BESETZEN
  
  zx1 = 0.5*r*edadphi
  IF((j < jeh) .OR. (neighbor(2) /= -1)) THEN
    zfkio  = cphi(j+1,2)/cphi(j+1,1)
    zfkiu  = cphi(j  ,2)/cphi(j+1,1)
    zfdivx = acphir(j+1,1)*eddlam
    zfdivo = acphir(j+1,1)*eddphi*cphi(j+1,2)
    zfdivu = acphir(j+1,1)*eddphi*cphi(j  ,2)
    
    DO k = 1 , ke
      DO i = iah - 1, ieh + 1
        zgv(i,k,jo3) = .5*(zdp(i,k,jo3)+(dak(k)+dbk(k)*ps(i,j+2,nj)))*  &
            v(i,j+1,k,nj)
        zpphi(i,k,jo2) = zx1*(ztv(i,k,jo5) +    ztv(i,k,jm5))*  &
            (zaloph(i,k,jo3) - zaloph(i,k,jm3))
        zekin(i,k,jo2) = .25*(  u(i-1,j+1,k,nj)**2 + u(i,j+1,k,nj)**2+  &
            zfkiu*v(i,j  ,k,nj)**2 + zfkio*v(i,j+1,k,nj)**2)
      END DO
    END DO
    
    DO k = 1 , ke
      DO i = iah - 1, ieh + 1
        zsdiv(i,k+1,jo2) = zsdiv(i,k,jo2) +  &
            zfdivx*( zgu(i,k,jo2)-zgu(i-1,k,jo2) ) +  &
            zfdivo*zgv(i,k,jo3) - zfdivu*zgv(i,k,jm3)
      END DO
    END DO
    
    DO k = 2 , ke
      DO i = iah - 1, ieh + 1
        zetas(i,k,jo2) = bk(k)*zsdiv(i,ke1,jo2) - zsdiv(i,k,jo2)
      END DO
      etas(iah-1:ieh+1,j+1,k)  = zetas(iah-1:ieh+1,k,jo2)
    END DO
    
  ELSE
    DO k =   1 , ke
      DO i = iah , ieh
        zpphi(i,k,jo2) =  zx1*(ztv(i,k,jo5) +    ztv(i,k,jm5))  &
            *(zaloph(i,k,jo3) - zaloph(i,k,jm3))
      END DO
    END DO
  END IF
  
!        LOKALE HILFSFELDER FUER HORIZONTALDIFFUSION
  DO k = 1 , ke
    DO i = iah - 1 , ieh + 1
!              GEWICHTETE HDIFF
      zlapt(i,k,jo3) = ((ztpa(i+1,k,jo5)-ztpa(i  ,k,jo5))*zdpu(i  ,k,jo3)-  &
          (ztpa(i  ,k,jo5)-ztpa(i-1,k,jo5))*zdpu(i-1,k,jo3)+  &
          (cphi(j+1,2)*(ztpa(i  ,k,jn5)-ztpa(i  ,k,jo5))*zdpv(i  ,k,jo4)-  &
          cphi(j  ,2)*(ztpa(i  ,k,jo5)-ztpa(i  ,k,jm5))*zdpv(i  ,k,jm4))/  &
          cphi(j+1,1)) / zdp5(i,k,jo5)
      zlapqd(i,k,jo3) =  &
          ((qd(i+1,j+1,k,na)-qd(i  ,j+1,k,na))*zdpu(i  ,k,jo3)-  &
          (qd(i  ,j+1,k,na)-qd(i-1,j+1,k,na))*zdpu(i-1,k,jo3)+ (cphi(j+1,2)*  &
          (qd(i  ,jzn,k,na)-qd(i  ,j+1,k,na))*zdpv(i  ,k,jo4)- cphi(j  ,2)*  &
          (qd(i  ,j+1,k,na)-qd(i  ,j  ,k,na))*zdpv(i  ,k,jm4))/  &
          cphi(j+1,1)) / zdp5(i,k,jo5)
      zlapqw(i,k,jo3) =  &
          ((qw(i+1,j+1,k,na)-qw(i  ,j+1,k,na))*zdpu(i  ,k,jo3)-  &
          (qw(i  ,j+1,k,na)-qw(i-1,j+1,k,na))*zdpu(i-1,k,jo3)+ (cphi(j+1,2)*  &
          (qw(i  ,jzn,k,na)-qw(i  ,j+1,k,na))*zdpv(i  ,k,jo4)- cphi(j  ,2)*  &
          (qw(i  ,j+1,k,na)-qw(i  ,j  ,k,na))*zdpv(i  ,k,jm4))/  &
          cphi(j+1,1)) / zdp5(i,k,jo5)
      zlapqi(i,k,jo3) =  &
          ((qi(i+1,j+1,k,na)-qi(i  ,j+1,k,na))*zdpu(i  ,k,jo3)-  &
          (qi(i  ,j+1,k,na)-qi(i-1,j+1,k,na))*zdpu(i-1,k,jo3)+ (cphi(j+1,2)*  &
          (qi(i  ,jzn,k,na)-qi(i  ,j+1,k,na))*zdpv(i  ,k,jo4)- cphi(j  ,2)*  &
          (qi(i  ,j+1,k,na)-qi(i  ,j  ,k,na))*zdpv(i  ,k,jm4))/  &
          cphi(j+1,1)) / zdp5(i,k,jo5)
!     *********************
      zlapu (i,k,jo3) = u (i+1,j+1,k,na) + u (i-1,j+1,k,na) -  &
          2.0*u (i  ,j+1,k,na) +  &
          ( cphi(j+1,2)*(u (i  ,jzn,k,na) - u (i  ,j+1,k,na)) -  &
          cphi(j  ,2)*(u (i  ,j+1,k,na) - u (i  ,j  ,k,na)) )/ cphi(j+1,1)
      zlapv (i,k,jo3) = v (i+1,j+1,k,na) + v (i-1,j+1,k,na) -  &
          2.0*v (i  ,j+1,k,na) +  &
          ( cphi(jzn,1)*(v (i  ,jzn,k,na) - v (i  ,j+1,k,na)) -  &
          cphi(j+1,1)*(v (i  ,j+1,k,na) - v (i  ,j  ,k,na)) )/ cphi(j+1,2)
      
    END DO
  END DO
  
!        AM RECHTEN UND LINKEN RAND H-DIFFUSION 2. ORDNUNG
  DO k = 1 , ke
    
    IF( (j .EQ. jah) .AND. (neighbor(4) == -1)) THEN
!              SUEDLICHER RAND DES GESAMTGEBIETS
!              H-DIFFUSION 2. ORDNUNG FUER ALLE I
      DO i = iah , ieh
        ztdifh(i,k) = vvfh(k)*aks2*zlapt (i,k,jm3)
        zqddih(i,k) = vvfh(k)*aks2*zlapqd(i,k,jm3)
        zqwdih(i,k) = vvfh(k)*aks2*zlapqw(i,k,jm3)
        zudifh(i,k) = vvfh(k)*aks2*zlapu (i,k,jm3)
        zvdifh(i,k) = vvfh(k)*aks2*zlapv (i,k,jm3)
        zqidih(i,k) = vvfh(k)*aks2*zlapqi(i,k,jm3)
      END DO
    ELSE IF((j == jeh) .AND. (neighbor(2) == -1)) THEN
!              NOERDLICHER RAND DES GESAMTGEBIETS
!              H-DIFFUSION 2. ORDNUNG FUER ALLE I
      DO i = iah , ieh
        ztdifh(i,k) = vvfh(k)*aks2*zlapt (i,k,jm3)
        zqddih(i,k) = vvfh(k)*aks2*zlapqd(i,k,jm3)
        zqwdih(i,k) = vvfh(k)*aks2*zlapqw(i,k,jm3)
        zudifh(i,k) = vvfh(k)*aks2*zlapu (i,k,jm3)
        zvdifh(i,k) = vvfh(k)*aks2*zlapv (i,k,jm3)
        zqidih(i,k) = vvfh(k)*aks2*zlapqi(i,k,jm3)
      END DO
    ELSE
!              BREITENKREIS AUS DER GESAMTGEBIETSMITTE
      DO i = iah , ieh
        IF((i == iah) .AND. (neighbor(1) == -1)) THEN
!                    WESTLICHER RAND DES GESAMTGEBIETS
!                    H-DIFFUSION 2. ORDNUNG FUER DIESES I
          ztdifh(iah,k) = vvfh(k)*aks2*zlapt (iah,k,jm3)
          zqddih(iah,k) = vvfh(k)*aks2*zlapqd(iah,k,jm3)
          zqwdih(iah,k) = vvfh(k)*aks2*zlapqw(iah,k,jm3)
          zudifh(iah,k) = vvfh(k)*aks2*zlapu (iah,k,jm3)
          zvdifh(iah,k) = vvfh(k)*aks2*zlapv (iah,k,jm3)
          zqidih(iah,k) = vvfh(k)*aks2*zlapqi(iah,k,jm3)
        ELSE IF((i == ieh) .AND. (neighbor(3) == -1)) THEN
          ztdifh(ieh,k) = vvfh(k)*aks2*zlapt (ieh,k,jm3)
          zqddih(ieh,k) = vvfh(k)*aks2*zlapqd(ieh,k,jm3)
          zqwdih(ieh,k) = vvfh(k)*aks2*zlapqw(ieh,k,jm3)
          zudifh(ieh,k) = vvfh(k)*aks2*zlapu (ieh,k,jm3)
          zvdifh(ieh,k) = vvfh(k)*aks2*zlapv (ieh,k,jm3)
          zqidih(ieh,k) = vvfh(k)*aks2*zlapqi(ieh,k,jm3)
        ELSE
          
!                    IM INNEREN DES MODELLGEBIETES H-DIFFUSION 4. ORDNUNG
!                    GEWICHTETE HDIFF
          ztdifh(i,k) = -vvfh(k)*aks4*  &
              ((zlapt(i+1,k,jm3)-zlapt(i  ,k,jm3))*zdpu(i  ,k,jm3)-  &
              (zlapt(i  ,k,jm3)-zlapt(i-1,k,jm3))*zdpu(i-1,k,jm3)+  &
              ( cphi(j  ,2)*  &
              (zlapt(i  ,k,jo3)-zlapt(i  ,k,jm3))*zdpv(i  ,k,jm4)- cphi(j-1,2)*  &
              (zlapt(i  ,k,jm3)-zlapt(i  ,k,ju3))*zdpv(i  ,k,ju4))/  &
              cphi(j  ,1)) / zdp5(i,k,jm5)
          zqddih(i,k) = - vvfh(k)*aks4*  &
              ((zlapqd(i+1,k,jm3)-zlapqd(i  ,k,jm3))*zdpu(i  ,k,jm3)-  &
              (zlapqd(i  ,k,jm3)-zlapqd(i-1,k,jm3))*zdpu(i-1,k,jm3)+  &
              ( cphi(j  ,2)*  &
              (zlapqd(i  ,k,jo3)-zlapqd(i  ,k,jm3))*zdpv(i  ,k,jm4)-  &
              cphi(j-1,2)*  &
              (zlapqd(i  ,k,jm3)-zlapqd(i  ,k,ju3))*zdpv(i  ,k,ju4))/  &
              cphi(j  ,1)) / zdp5(i,k,jm5)
          zqwdih(i,k) = - vvfh(k)*aks4*  &
              ((zlapqw(i+1,k,jm3)-zlapqw(i  ,k,jm3))*zdpu(i  ,k,jm3)-  &
              (zlapqw(i  ,k,jm3)-zlapqw(i-1,k,jm3))*zdpu(i-1,k,jm3)+  &
              ( cphi(j  ,2)*  &
              (zlapqw(i  ,k,jo3)-zlapqw(i  ,k,jm3))*zdpv(i  ,k,jm4)-  &
              cphi(j-1,2)*  &
              (zlapqw(i  ,k,jm3)-zlapqw(i  ,k,ju3))*zdpv(i  ,k,ju4))/  &
              cphi(j  ,1)) / zdp5(i,k,jm5)
          zqidih(i,k) = - vvfh(k)*aks4*  &
              ((zlapqi(i+1,k,jm3)-zlapqi(i  ,k,jm3))*zdpu(i  ,k,jm3)-  &
              (zlapqi(i  ,k,jm3)-zlapqi(i-1,k,jm3))*zdpu(i-1,k,jm3)+  &
              ( cphi(j  ,2)*  &
              (zlapqi(i  ,k,jo3)-zlapqi(i  ,k,jm3))*zdpv(i  ,k,jm4)-  &
              cphi(j-1,2)*  &
              (zlapqi(i  ,k,jm3)-zlapqi(i  ,k,ju3))*zdpv(i  ,k,ju4))/  &
              cphi(j  ,1)) / zdp5(i,k,jm5)
!                    *****************
          zudifh(i,k) = - vvfh(k)*aks4*  &
              ( zlapu (i+1,k,jm3) + zlapu (i-1,k,jm3) - 2.0*zlapu (i  ,k,jm3) +  &
              ( cphi(j  ,2)*(zlapu (i  ,k,jo3) - zlapu (i  ,k,jm3)) -  &
              cphi(j-1,2)*(zlapu (i  ,k,jm3) - zlapu (i  ,k,ju3)) )/ cphi(j  ,1) )
          zvdifh(i,k) = - vvfh(k)*aks4*  &
              ( zlapv (i+1,k,jm3) + zlapv (i-1,k,jm3) - 2.0*zlapv (i  ,k,jm3) +  &
              ( cphi(j+1,1)*(zlapv (i  ,k,jo3) - zlapv (i  ,k,jm3)) -  &
              cphi(j  ,1)*(zlapv (i  ,k,jm3) - zlapv (i  ,k,ju3)) )/ cphi(j  ,2) )
        END IF
      END DO
    END IF
  END DO
  
!        EXPLIZITE PROGNOSE OHNE RANDRELAXATIONSTERME
!        --------------------------------------------
  
  zx2   =  rerd*acphir(j,1)
  zfadvx= .5*acphir(j,1)*eddlam
  zfadvy= .5*acphir(j,1)*eddphi
  
!        1. BODENDRUCK + DRUCK
!        ---------------------
  
  DO i = iah , ieh
    psdt(i) = zsdiv(i,ke1,jm2)
    ps(i,j,NE) = ps(i,j,na) - psdt(i)*dt2
    pint(i,j,1,NE) = ak(1) + bk(1)*ps(i,j,NE)
  END DO
  
  DO k = 2 , ke+1
    km1 = k - 1
    DO i = iah , ieh
      pint(i,j,k,NE) = -(bk(km1) + bk(k)) * dt2 * psdt(i) *  &
          dwdt(i,j,km1,nj) + pint(i,j,km1,na) + pint(i,j,k,na) -  &
          pint(i,j,km1,NE)
    END DO
  END DO
  
!        2. H - QDW - PROGNOSE
!        ---------------------
  
!        HORIZONTALADVEKTION UND QUELLTERME
  DO k = 1 , ke
    DO i = iah , ieh
      zeddpq(i,k) = 1./zdp(i,k,jm3)
      zgvco       = zgv(i,k,jm3)*cphi(j  ,2)
      zgvcu       = zgv(i,k,ju3)*cphi(j-1,2)
      zt1         = zgu(i  ,k,jm2)* ( t(i+1,j,k,nj) - t(i  ,j,k,nj) )
      zt2         = zgu(i-1,k,jm2)* ( t(i  ,j,k,nj) - t(i-1,j,k,nj) )
      zt3         = zgvco*( t(i,j+1,k,nj) - t(i,j,k,nj) )
      zt4         = zgvcu*( t(i,j,k,nj) - t(i,j-1,k,nj) )
      ztadv(i,k)  = -( zfadvx*(zt1+zt2) + zfadvy*(zt3+zt4) ) *zeddpq(i,k)
      zqd1        = zgu(i  ,k,jm2)* (qd(i+1,j,k,nj) - qd(i  ,j,k,nj))
      zqd2        = zgu(i-1,k,jm2)* (qd(i  ,j,k,nj) - qd(i-1,j,k,nj) )
      zqd3        = zgvco*( qd(i,j+1,k,nj)-qd(i,j,k,nj))
      zqd4        = zgvcu*( qd(i,j,k,nj)-qd(i,j-1,k,nj))
      zqdadv(i,k) = -( zfadvx*(zqd1+zqd2) + zfadvy*(zqd3+zqd4)) *zeddpq(i,k)
      
!              ALPHA*OMEGA
      za1 =  zgu(i  ,k,jm2)*zplam(i  ,k) + zgu(i-1,k,jm2)*zplam(i-1,k)
      za2 = ( zgvco*zpphi(i,k,jo2) + zgvcu*zpphi(i,k,jm2) )*zx2
      za3 = - r*ztv(i,k,jm5)*zeddpq(i,k)*  &
          (   (zalopn(i,k+1,jm3)-zalopn(i,k,jm3))*zsdiv (i,k,jm2)  &
          + (zsdiv (i,k+1,jm2)-zsdiv (i,k,jm2))*zbetak(i,k,jm3) )
      
      zalpom(i,k) = ( .5*( za1 + za2 )*zeddpq(i,k) + za3 )
      
      agb(i,k,1) = ed2dt
      agd(i,k,1) = ed2dt*t(i,j,k,na) + ztadv(i,k) +  &
          ztdifh(i,k) + zalpom(i,k) *wcpr+ zttk(i,k) +  &
          ztts(i,k) + zthdta(i,k) + zsodta(i,k)
      agd(i,k,2) = ed2dt*qd(i,j,k,na) + zqdadv(i,k) +  &
          zqddih(i,k) + zqdtk(i,k) + zqdts(i,k)
      agd(i,k,3) = ed2dt*qw(i,j,k,na) + zqwdih(i,k) + zqwts(i,k)
      agd(i,k,4) = ed2dt*qi(i,j,k,na) + zqidih(i,k) + zqits(i,k)
    END DO
  END DO
  
  
!        VERTIKALADVEKTION UND -DIFFUSION
  DO i = iah , ieh
    zagcm      = .5*zetas(i,2,jm2)*zeddpq(i,1)
    ztmkvhm    = ztkvz*ztmkvh(i  ,2,jm3) + ztkvl*  &
        ( ztmkvh(i+1,2,jm3) + ztmkvh(i,2,jo3) +  &
        ztmkvh(i-1,2,jm3) + ztmkvh(i,2,ju3) )
    zagct      = - ztmkvhm*zeddpq(i,1)
    agc(i,1,1) = zagcm*za1a + zagct*a1t(2)
    agb(i,1,1) = agb(i,1,1) - agc(i,1,1)
    agd(i,1,1) = agd(i,1,1) - ( za2a*zagcm + a2t(2)*zagct )*  &
        ( t(i,j,2,na) - t(i,j,1,na) )  &
        -wcpr* zagct*( zfihf(i,2,jm5) - zfihf(i,1,jm5) )
    agd(i,1,2) = agd(i,1,2) - ( za2a*zagcm + a2t(2)*zagct )*  &
        ( qd(i,j,2,na) - qd(i,j,1,na) )
    agd(i,1,3) = agd(i,1,3) - ( za2a*zagcm + a2t(2)*zagct )*  &
        ( qw(i,j,2,na) - qw(i,j,1,na) )
    agd(i,1,4) = agd(i,1,4) - ( za2a*zagcm + a2t(2)*zagct )*  &
        ( qi(i,j,2,na) - qi(i,j,1,na) )
    agb(i,1,2) = agb(i,1,1)
    agc(i,1,2) = agc(i,1,1)
    agb(i,1,3) = agb(i,1,1)
    agc(i,1,3) = agc(i,1,1)
    agb(i,1,4) = agb(i,1,1)
    agc(i,1,4) = agc(i,1,1)
  END DO
  
  DO k = 2 , ke-1
    DO i  = iah , ieh
      zagam      = -.5*zetas(i,k  ,jm2)*zeddpq(i,k)
      zagcm      =  .5*zetas(i,k+1,jm2)*zeddpq(i,k)
      ztmkvhm    = ztkvz*ztmkvh(i  ,k,jm3) + ztkvl*  &
          ( ztmkvh(i+1,k,jm3) + ztmkvh(i  ,k,jo3) +  &
          ztmkvh(i-1,k,jm3) + ztmkvh(i  ,k,ju3) )
      zagat      = - ztmkvhm*zeddpq(i,k)
      ztmkvhm    = ztkvz*ztmkvh(i  ,k+1,jm3) + ztkvl*  &
          ( ztmkvh(i+1,k+1,jm3) + ztmkvh(i  ,k+1,jo3) +  &
          ztmkvh(i-1,k+1,jm3) + ztmkvh(i  ,k+1,ju3) )
      zagct      = - ztmkvhm*zeddpq(i,k)
      aga(i,k,1) = zagam*za1a + zagat*a1t(k)
      agc(i,k,1) = zagcm*za1a + zagct*a1t(k+1)
      agb(i,k,1) = agb(i,k,1) - aga(i,k,1) - agc(i,k,1)
      agd(i,k,1) = agd(i,k,1) - ( za2a*zagam + a2t(k)*zagat )*  &
          ( t(i,j,k-1,na) - t(i,j,k,na) ) - ( za2a*zagcm + a2t(k+1)*zagct )*  &
          ( t(i,j,k+1,na) - t(i,j,k,na) )  &
          -wcpr* zagat*( zfihf(i,k-1,jm5) - zfihf(i,k,jm5))  &
          -wcpr* zagct*( zfihf(i,k+1,jm5) - zfihf(i,k,jm5))
      agd(i,k,2) = agd(i,k,2) - ( za2a*zagam + a2t(k)*zagat )*  &
          (qd(i,j,k-1,na) - qd(i,j,k,na) ) - ( za2a*zagcm + a2t(k+1)*zagct )*  &
          (qd(i,j,k+1,na) - qd(i,j,k,na) )
      agd(i,k,3) = agd(i,k,3) - ( za2a*zagam + a2t(k)*zagat )*  &
          (qw(i,j,k-1,na) - qw(i,j,k,na) ) - ( za2a*zagcm + a2t(k+1)*zagct )*  &
          (qw(i,j,k+1,na) - qw(i,j,k,na) )
      agd(i,k,4) = agd(i,k,4) - ( za2a*zagam + a2t(k)*zagat )*  &
          (qi(i,j,k-1,na) - qi(i,j,k,na) ) - ( za2a*zagcm + a2t(k+1)*zagct )*  &
          (qi(i,j,k+1,na) - qi(i,j,k,na) )
      aga(i,k,2) = aga(i,k,1)
      agb(i,k,2) = agb(i,k,1)
      agc(i,k,2) = agc(i,k,1)
      aga(i,k,3) = aga(i,k,1)
      agb(i,k,3) = agb(i,k,1)
      agc(i,k,3) = agc(i,k,1)
      aga(i,k,4) = aga(i,k,1)
      agb(i,k,4) = agb(i,k,1)
      agc(i,k,4) = agc(i,k,1)
    END DO
  END DO
  
  DO i = iah , ieh
    zagam       = -.5*zetas(i,ke,jm2)*zeddpq(i,ke)
    ztmkvhm     = ztkvz*ztmkvh(i  ,ke,jm3) + ztkvl*  &
        ( ztmkvh(i+1,ke,jm3) + ztmkvh(i  ,ke,jo3) +  &
        ztmkvh(i-1,ke,jm3) + ztmkvh(i  ,ke,ju3) )
    zagat       = -  ztmkvhm *zeddpq(i,ke)
    zagct       = - tmch(i,j)*zeddpq(i,ke)
    aga(i,ke,1) = zagam*za1a  + zagat*a1t(ke)
    agb(i,ke,1) = agb(i,ke,1) - aga(i,ke,1) - zagct*a1t(ke1)
    agd(i,ke,1) = agd(i,ke,1) - ( za2a*zagam + a2t(ke)*zagat )*  &
        ( t(i,j,ke-1,na) - t(i,j,ke,na) )  &
        - a2t(ke1)*zagct*( tg(i,j,na) - t(i,j,ke,na) )  &
        - a1t(ke1)*zagct*  tg(i,j,na)  &
        -wcpr* zagat*( zfihf(i,ke-1,jm5) - zfihf(i,ke,jm5) )  &
        -wcpr* zagct*( fib  (i,j)        - zfihf(i,ke,jm5) )
    agd(i,ke,2) = agd(i,ke,2) - ( za2a*zagam + a2t(ke)*zagat )*  &
        ( qd(i,j,ke-1,na) - qd(i,j,ke,na) )  &
        - a2t(ke1)*zagct*(qdb(i,j,na)  - qd(i,j,ke,na) )  &
        - a1t(ke1)*zagct* qdb(i,j,na)
    agd(i,ke,3) = agd(i,ke,3) - ( za2a*zagam + a2t(ke)*zagat )*  &
        ( qw(i,j,ke-1,na) - qw(i,j,ke,na) ) + a2t(ke1)*zagct* qw(i,j,ke,na)
    agd(i,ke,4) = agd(i,ke,4) - ( za2a*zagam + a2t(ke)*zagat )*  &
        ( qi(i,j,ke-1,na) - qi(i,j,ke,na) ) + a2t(ke1)*zagct* qi(i,j,ke,na)
    aga(i,ke,2) = aga(i,ke,1)
    agb(i,ke,2) = agb(i,ke,1)
    aga(i,ke,3) = aga(i,ke,1)
    agb(i,ke,3) = agb(i,ke,1)
    aga(i,ke,4) = aga(i,ke,1)
    agb(i,ke,4) = agb(i,ke,1)
  END DO
  
!        GAUSS - ELIMINATION UND H-QDW-ZERLEGUNG
  
  CALL gauss ( iah , ieh , index1, aga,agb,agc,agd,age)
  CALL gauss ( iah , ieh , index2, aga,agb,agc,agd,age)
  CALL gauss ( iah , ieh , index3, aga,agb,agc,agd,age)
  CALL gauss ( iah , ieh , index4, aga,agb,agc,agd,age)
  
  DO k = 1 , ke
    DO i = iah , ieh
      he(i,k)= age(i,k,1)*wcp + wlk * age(i,k,2)*ztrc
      qdwe(i,k)= (age(i,k,2) + age(i,k,3) + age(i,k,4))*ztrc
      qdle(i,k)= (age(i,k,2) + age(i,k,3))*ztrc
      qdie(i,k)= age(i,k,4)*ztrc
      tstart(i,k) = t(i,j,k,nj)
      gqdsta(i,k) = zgqd(i,k,jm3)
      phfsta(i,k) = zphf(i,k,jm5)
      phfe  (i,k) = akh(k) + bkh(k)*ps(i,j,NE)
    END DO
  END DO
  
!        BERECHNUNG DES FLUESSES LATENTER WAERME AM BODEN;
!        DIESER FLUESS WIRD AUFSUMMIERT
  
  DO i  = iah , ieh
    IF (infrl(i,j) > 0) THEN
      bflqdsl(i,j) = bflqdsl(i,j) - tmchl(i,j)*edg*(wlk*  &
          (a2t(ke1)*(qdbl(i,j,na) - (qd(i,j,ke,na)+ (qw(i,j,ke,na))) ) +  &
          a1t(ke1)*(qdbl(i,j,NE) - qdle (i,ke  )))-wls*(  &
          a2t(ke1)*qi(i,j,ke,na)+a1t(ke1)*qdie(i,ke)))*dtdeh
    END IF
    IF (infrw(i,j) > 0) THEN
      bflqdsw(i,j) = bflqdsw(i,j) - tmchw(i,j)*edg*(wlk*  &
          (a2t(ke1)*(qdbw(i,j,na) - (qd(i,j,ke,na)+ (qw(i,j,ke,na))) ) +  &
          a1t(ke1)*(qdbw(i,j,NE) - qdle (i,ke  )))-wls*(  &
          a2t(ke1)*qi(i,j,ke,na)+a1t(ke1)*qdie(i,ke)))*dtdeh
    END IF
    IF (infri(i,j) > 0) THEN
      bflqdsi(i,j) = bflqdsi(i,j) - tmchi(i,j)*edg*(wlk*  &
          (a2t(ke1)*(qdbi(i,j,na) - (qd(i,j,ke,na)+ (qw(i,j,ke,na))) ) +  &
          a1t(ke1)*(qdbi(i,j,NE) - qdle (i,ke  )))-wls*(  &
          a2t(ke1)*qi(i,j,ke,na)+a1t(ke1)*qdie(i,ke)))*dtdeh
    END IF
    bflqds(i,j) = (FLOAT(infrl(i,j))*bflqdsl(i,j)  &
        +  FLOAT(infrw(i,j))*bflqdsw(i,j)  &
        +  FLOAT(infri(i,j))*bflqdsi(i,j))*edfakinf
  END DO
  
!        MASSENFLUSS-KORREKTURSCHEMA  (FUER QDW<0)
  IF(lmassf) THEN
    DO i = iah , ieh
      zqkor(i) = 0.
    END DO
    DO k = 1 , ke
      DO i = iah , ieh
        zqkor(i)   = qdwe(i,k) + zqkor(i)*zeddpq(i,k)
        IF(zqkor(i) < 0.) THEN
          qdwe(i,k) = 0.
          zqkor(i)  = zqkor(i)*zdp(i,k,jm3)
        ELSE
          qdwe(i,k) = zqkor(i)
          zqkor(i)  = 0.
        END IF
      END DO
    END DO
  END IF
  
!        BEI ADIABATISCHER INITIALISIERUNG: T MIT UNVERAENDERTEM QD BERECH-
!        NEN; SONST IST DIE SKALIGE KONDENSATIONSRATE UNGLEICH NULL.
  IF (laistep) THEN
    DO k = 1 , ke
      DO i = iah , ieh
        qd(i,j,k,NE) = qd(i,j,k,na)
        qw(i,j,k,NE) = qw(i,j,k,na)
        t (i,j,k,NE) = wcpr*he(i,k)
        qi(i,j,k,NE) = qi(i,j,k,na)
      END DO
    END DO
    
  ELSE
    DO k = 1 , ke
      DO i = iah , ieh
        qd(i,j,k,NE) = age(i,k,2)
        qw(i,j,k,NE) = age(i,k,3)
        t (i,j,k,NE) = age(i,k,1)
        qi(i,j,k,NE) = age(i,k,4)
      END DO
    END DO
    DO i = iah , ieh
      zqkor(i) = 0.
    END DO
    DO k = 1 , ke
      DO i = iah , ieh
        zqw1 = qw(i,j,k,NE)
        zqw2 = qi(i,j,k,NE)
        IF (zqw1 < 0.) THEN
          t (i,j,k,NE) = t (i,j,k,NE) - wlk*qw(i,j,k,NE)*wcpr
          qd(i,j,k,NE) = qd(i,j,k,NE) + qw(i,j,k,NE)
          qw(i,j,k,NE) = 0.
        END IF
        IF (zqw2 < 0.) THEN
          t (i,j,k,NE) = t (i,j,k,NE) - wls*qi(i,j,k,NE)*wcpr
          qd(i,j,k,NE) = qd(i,j,k,NE) + qi(i,j,k,NE)
          qi(i,j,k,NE) = 0.
        END IF
        zqkor(i)   = qd(i,j,k,NE) + zqkor(i)*zeddpq(i,k)
        IF(zqkor(i) < 0.) THEN
          qd(i,j,k,NE) = 0.
          qw(i,j,k,NE) = 0.
          qi(i,j,k,NE) = 0.
          zqkor(i)  = zqkor(i)*zdp(i,k,jm3)
        ELSE
          qd(i,j,k,NE) = zqkor(i)
          zqkor(i)  = 0.
        END IF
      END DO
    END DO
    
!     BERECHNUNG DES FLUESSES SENSIBLER WAERME AM BODEN;
    
    DO i = iah , ieh
      IF (infrl(i,j) > 0) THEN
        bflhsl(i,j)= bflhsl(i,j)- tmchl(i,j)*edg*wcp*  &
            (a2t(ke1)*(tgl(i,j,na) - t    (i,j,ke,na)) +  &
            a1t(ke1)*(tgl(i,j,NE) - t    (i,j,ke,NE)) +  &
            wcpr*    (fib(i,j)    - zfihf(i,  ke,jm5)))*dtdeh
      END IF
      IF (infrw(i,j) > 0) THEN
        bflhsw(i,j)= bflhsw(i,j)- tmchw(i,j)*edg*wcp*  &
            (a2t(ke1)*(tgw(i,j,na) - t    (i,j,ke,na)) +  &
            a1t(ke1)*(tgw(i,j,NE) - t    (i,j,ke,NE)) +  &
            wcpr*    (fib(i,j)    - zfihf(i,  ke,jm5)))*dtdeh
      END IF
      IF (infri(i,j) > 0) THEN
        bflhsi(i,j)= bflhsi(i,j)- tmchi(i,j)*edg*wcp*  &
            (a2t(ke1)*(tgi(i,j,na) - t    (i,j,ke,na)) +  &
            a1t(ke1)*(tgi(i,j,NE) - t    (i,j,ke,NE)) +  &
            wcpr*    (fib(i,j)    - zfihf(i,  ke,jm5)))*dtdeh
      END IF
      bflhs(i,j) = (FLOAT(infrl(i,j))*bflhsl(i,j)  &
          +  FLOAT(infrw(i,j))*bflhsw(i,j)  &
          +  FLOAT(infri(i,j))*bflhsi(i,j))*edfakinf
    END DO
    
!         OMEGA-WERTE BESTIMMTER MODELLFLAECHEN SUMMIEREN (PROGCHK)
    DO i = iah , ieh
      zom850m   = zom850m +  &
          ABS(zalpom(i,kfl850)*zphf(i,kfl850,jm5)/(r*ztv(i,kfl850,jm5)))
      zom500m   = zom500m +  &
          ABS(zalpom(i,kfl500)*zphf(i,kfl500,jm5)/(r*ztv(i,kfl500,jm5)))
      zom300m   = zom300m +  &
          ABS(zalpom(i,kfl300)*zphf(i,kfl300,jm5)/(r*ztv(i,kfl300,jm5)))
    END DO
  END IF ! LAISTEP
  
  
!        3. U - PROGNOSE
!        ---------------
  
  zx1   = acphir(j,1)*eddlam
  zfvzo = 0.125*cphi(j  ,2)/cphi(j,1)
  zfvzu = 0.125*cphi(j-1,2)/cphi(j,1)
  
  DO k = 1 , ke
    DO i = iau , ieu
      zgrad(i,k) = -zplam(i,k) - zx1* ( (zfih(i+1,k,jm3) - zfih(i,k,jm3)) *  &
          0.5 * (dwdt(i,j,k,nj)+dwdt(i+1,j,k,nj))  &
          + zekin(i+1,k,jm2) - zekin(i,k,jm2) )
      zvzeq(i,k) = ( zzeta(i,k,jo2) + zzeta(i,k,jm2) )*  &
          (  zfvzo*( zgv(i+1,k,jm3) + zgv(i,k,jm3) )  &
          + zfvzu*( zgv(i+1,k,ju3) + zgv(i,k,ju3) ) )
      zeddpq(i,k) =  2./( zdp(i+1,k,jm3) + zdp(i,k,jm3) )
      agb(i,k,1) = ed2dt
      agd(i,k,1) =(ed2dt-zepsray)*u(i,j,k,na) + zgrad(i,k) +  &
          zvzeq(i,k) + zudifh(i,k) + 0.5*( zutk(i+1,k,jm2) + zutk(i,k,jm2) )
    END DO
  END DO
  
!        VERTIKALDIFFUSION UND VERTIKALADVEKTION
  DO i = iau , ieu
    zagcm = .25*( zetas (i,2,jm2) + zetas(i+1,2,jm2) )* zeddpq(i,1)
    zagct = -0.5*(ztmkvm(i,2,jm3)+ ztmkvm(i+1,2,jm3) )* zeddpq(i,1)
    agc(i,1,1) = zagcm*za1a + zagct*a1t(2)
    agb(i,1,1) = agb(i,1,1) - agc(i,1,1)
    agd(i,1,1) = agd(i,1,1) - ( za2a*zagcm + a2t(2)*zagct )*  &
        ( u(i,j,2,na) - u(i,j,1,na ) )
  END DO
  
  DO k = 2 , ke-1
    DO i = iau , ieu
      zagam = - 0.25*(zetas(i,k  ,jm2)+zetas(i+1,k  ,jm2))* zeddpq(i,k)
      zagcm =   0.25*(zetas(i,k+1,jm2)+zetas(i+1,k+1,jm2))* zeddpq(i,k)
      zagat = - 0.5 *(ztmkvm(i,k  ,jm3)+ztmkvm(i+1,k  ,jm3))* zeddpq(i,k)
      zagct = - 0.5 *(ztmkvm(i,k+1,jm3)+ztmkvm(i+1,k+1,jm3))* zeddpq(i,k)
      aga(i,k,1) = zagam*za1a + zagat*a1t(k)
      agc(i,k,1) = zagcm*za1a + zagct*a1t(k+1)
      agb(i,k,1) = agb(i,k,1) - aga(i,k,1) - agc(i,k,1)
      agd(i,k,1) = agd(i,k,1) - ( za2a*zagam + a2t(k)*zagat )*  &
          ( u(i,j,k-1,na) - u(i,j,k,na) ) - ( za2a*zagcm + a2t(k+1)*zagct )*  &
          ( u(i,j,k+1,na) - u(i,j,k,na) )
    END DO
  END DO
  
  DO i = iau , ieu
    zagam = - 0.25*(zetas  (i,ke,jm2)+zetas (i+1,ke,jm2))* zeddpq(i,ke)
    zagat = - 0.5 *( ztmkvm(i,ke,jm3)+ztmkvm(i+1,ke,jm3))* zeddpq(i,ke)
    zagct = - 0.5 *( tmcm (i,j   ) + tmcm (i+1,j  ) )    * zeddpq(i,ke)
    aga(i,ke,1) = zagam*za1a  + zagat*a1t(ke)
    agb(i,ke,1) = agb(i,ke,1) - aga(i,ke,1) - zagct*a1t(ke1)
    agd(i,ke,1) = agd(i,ke,1) - ( za2a*zagam + a2t(ke)*zagat )*  &
        ( u(i,j,ke-1,na) - u(i,j,ke,na) ) + a2t(ke1)*zagct*u(i,j,ke,na)
  END DO
  
!        GAUSS-ELIMINATION UND U-ENDE-WERTE AUF U(I,J,K,NE) ABLEGEN
  CALL gauss ( iau , ieu , index1, aga,agb,agc,agd,age)
  
  DO k = 1 , ke
    DO i = iau , ieu
      u(i,j,k,NE) = age(i,k,1)
    END DO
  END DO
  
!        BERECHNUNG DES U-IMPULSFLUSSES AM BODEN;
!        DIESER FLUESS WIRD AUFSUMMIERT
  DO i = iau , ieu
    bflus(i,j) = bflus(i,j) + 0.5*(tmcm(i,j) + tmcm(i+1,j))*edg*  &
        (a2t(ke1)*u(i,j,ke,na) + a1t(ke1)*u(i,j,ke,NE))*dtdeh
  END DO
  
!        4. V - PROGNOSE
!        ---------------
  
  IF(j <= jev) THEN
    
    DO k = 1 , ke
      DO i = iav , iev
        zgrad(i,k) = - zpphi(i,k,jo2) - edadphi*  &
            ( (zfih(i,k,jo3) - zfih(i,k,jm3) ) *  &
            0.5*(dwdt(i,j,k,nj)+dwdt(i,j+1,k,nj))  &
            + zekin(i,k,jo2) - zekin(i,k,jm2) )
        zuzeq(i,k) = - 0.125*( zzeta(i-1,k,jo2) + zzeta(i,k,jo2) )*  &
            (  zgu(i-1,k,jo2) + zgu  (i,k,jo2) + zgu(i-1,k,jm2) + zgu  (i,k,jm2) )
        zeddpq(i,k) =  2./( zdp(i,k,jo3) + zdp(i,k,jm3) )
        agb(i,k,2) = ed2dt
        agd(i,k,2) =(ed2dt-zepsray)*v(i,j,k,na) +  &
            zgrad(i,k) + zuzeq(i,k) + zvdifh(i,k) +  &
            0.5*( zvtk(i,k,jo2) + zvtk(i,k,jm2) )
      END DO
    END DO
    
!           VERTIKALDIFFUSION UND VERTIKALADVEKTION
    DO i = iav , iev
      zagcm =  0.25*( zetas(i,2,jm2) +  zetas(i,2,jo2) )* zeddpq(i,1)
      zagct = - 0.5*(ztmkvm(i,2,jm3) + ztmkvm(i,2,jo3) )* zeddpq(i,1)
      agc(i,1,2) = zagcm*za1a + zagct*a1t(2)
      agb(i,1,2) = agb(i,1,2) - agc(i,1,2)
      agd(i,1,2) = agd(i,1,2) - ( za2a*zagcm + a2t(2)*zagct )*  &
          ( v(i,j,2,na) - v(i,j,1,na ) )
    END DO
    
    DO k = 2 , ke-1
      DO i = iav , iev
        zagam = - 0.25*(zetas(i,k  ,jm2)+zetas(i,k  ,jo2))* zeddpq(i,k)
        zagcm =   0.25*(zetas(i,k+1,jm2)+zetas(i,k+1,jo2))* zeddpq(i,k)
        zagat = - 0.5 *( ztmkvm(i,k  ,jm3) + ztmkvm(i  ,k  ,jo3) ) *  &
            zeddpq(i,k)
        zagct = - 0.5 *( ztmkvm(i,k+1,jm3) + ztmkvm(i  ,k+1,jo3) ) *  &
            zeddpq(i,k)
        aga(i,k,2) = zagam*za1a + zagat*a1t(k)
        agc(i,k,2) = zagcm*za1a + zagct*a1t(k+1)
        agb(i,k,2) = agb(i,k,2) - aga(i,k,2) - agc(i,k,2)
        agd(i,k,2) = agd(i,k,2) - ( za2a*zagam + a2t(k)*zagat )*  &
            ( v(i,j,k-1,na) - v(i,j,k,na) )- ( za2a*zagcm + a2t(k+1)*zagct )*  &
            ( v(i,j,k+1,na) - v(i,j,k,na) )
      END DO
    END DO
    
    DO i = iav , iev
      zagam = - 0.25*(zetas (i,ke,jm2)+zetas (i,ke,jo2))* zeddpq(i,ke)
      zagat = - 0.5 *(ztmkvm(i,ke,jm3)+ztmkvm(i,ke,jo3))* zeddpq(i,ke)
      zagct = - 0.5 *( tmcm (i,j   ) + tmcm (i,j+1 ) )  * zeddpq(i,ke)
      aga(i,ke,2) = zagam*za1a  + zagat*a1t(ke)
      agb(i,ke,2) = agb(i,ke,2) - aga(i,ke,2) - zagct*a1t(ke1)
      agd(i,ke,2) = agd(i,ke,2) - ( za2a*zagam + a2t(ke)*zagat )*  &
          ( v(i,j,ke-1,na) - v(i,j,ke,na) ) + a2t(ke1)*zagct*v(i,j,ke,na)
    END DO
    
!           GAUSS ELIMINATION UND V-ENDE-WERTE AUF V(I,J,K,NE) ABLEGEN
    CALL gauss ( iav , iev , index2, aga,agb,agc,agd,age)
    
    DO k = 1 , ke
      DO i = iav , iev
        v(i,j,k,NE) = age(i,k,2)
      END DO
    END DO
    
!           BERECHNUNG DES V-IMPULSFLUSSES AM BODEN;
!           DIESER FLUESS WIRD AUFSUMMIERT
    DO i = iav , iev
      bflvs(i,j) = bflvs(i,j) + 0.5*(tmcm(i,j) + tmcm(i,j+1))*edg*  &
          (a2t(ke1)*v(i,j,ke,na) + a1t(ke1)*v(i,j,ke,NE))*dtdeh
    END DO
    
  END IF ! IF(J.LE.JEV)
  
!        ZYKLISCHE VERTAUSCHUNG DER SCHEIBENINDIZES
  
  jm2 = 3 - jm2
  jo2 = 3 - jo2
  
  jsp = ju3
  ju3 = jm3
  jm3 = jo3
  jo3 = jsp
  
  jsp = js4
  js4 = ju4
  ju4 = jm4
  jm4 = jo4
  jo4 = jsp
  
  jsp = js5
  js5 = ju5
  ju5 = jm5
  jm5 = jo5
  jo5 = jn5
  jn5 = jsp
  
END DO ! J = JAH, JEH

!     DIESER FLUESS WIRD AUFSUMMIERT
!     OMEGA-MITTELWERTE BESTIMMTER MODELLFLAECHEN BILDEN
zdeltf  = FLOAT( (iehgg-iahgg+1)*(jehgg-jahgg+1) )
zom850m = zom850m/zdeltf
zom500m = zom500m/zdeltf
zom300m = zom300m/zdeltf

!-----------------------------------------------------------------------
RETURN
!
!CONTAINS
!
! STATEMENT - FUNCTIONS ZUR BERECHNUNG DER SAEETIGUNGSFEUCHTE
! MAGNUS-FORMEL FUER WASSER
!
!REAL FUNCTION FGEW(TT)
!   IMPLICIT NONE
!   REAL, INTENT(IN) :: TT
!   FGEW = B1 * EXP( B2W*(TT - B3)/(TT - B4W) )
!END FUNCTION FGEW
!!
!! SPEZIFISCHE FEUCHTE AUS DAMPFDRUCK UND LUFTDRUCK
!!
!REAL FUNCTION FGQD(GE,PP)
!   IMPLICIT NONE
!   REAL, INTENT(IN) :: GE,PP
!   FGQD = RDRD*GE/(PP - EMRDRD*GE)
!END FUNCTION FGQD
!
END SUBROUTINE progexp
