      SUBROUTINE PROGEXP
        !
        ! WARNING! Here are 3 dummy argument names that
        ! differ from actual parameter names.
        ! The rest if fine.
        !
        !  OM850M(1), OM500M(1), OM300M(1),
        !     |        |        |
        !     v        v        v
     &    (ZOM850M, ZOM500M, ZOM300M,
     &     AK    , BK    , AKH    , BKH    , DAK    , DBK  , A1T   ,
     &     A2T   , VVFH  , FIB    , FC     , ACPHIR , CPHI , PS    ,
     &     TGL   , TGW   , TGI    , QDBL   , QDBW   , QDBI , BFLHSL,
     &     BFLHSW, BFLHSI, BFLQDSL, BFLQDSW, BFLQDSI, TG   , QDB   ,
     &     BFLHS , BFLQDS, BFLUS  , BFLVS  , U      , V    , T     ,
     &     QD    , QW    , FI     , TMKVMH , TMCHL  , TMCHW,
     &     TMCHI , TMCM  , TMCH   , SOTHDT , TTK    , QDTK , UVTK  ,
     &     TTS   , QDTS  , QWTS   , INFRL  , INFRW  , INFRI, QI    ,
     &     QITS  , PINT  , DWDT   , ETAS)
C
C
C**** PROGEXP  -   UP: PROGNOSE FUER EINEN ZEITSCHRITT
C**   AUFRUF   :   CALL PROGEXP ( JATPROG, JETPROG,
C**               1               ZOM850M, ZOM500M, ZOM300M )
C**   ENTRIES  :      ---
C**   ZWECK    :   BERECHNUNG ALLER PROGNOSEVARIABLEN FUER
C**                DEN ZEITPUNKT NZT+1 (EXPLIZITE PROGNOSE)
C**   VERSIONS-
C**   DATUM    :   22.7.1991
C**
C**   EXTERNALS:   GAUSS  , HQTOTQ , COPYRE
C**
C**   EINGABE-
C**   PARAMETER:   JATPROG: ANFANGS-J-INDEX DES PROGNOSEBEREICHS;
C**                JETPROG: END    -J-INDEX DES PROGNOSEBEREICHS;
C**                         DIE FESTLEGUNG ERFOLGT IN UP *PROGORG*.
C**   AUSGABE-
C**   PARAMETER    ZOM850M: OMEGA-MITTELWERT IN ETWA 850 HPA (KONTROLLE)
C**                ZOM500M: OMEGA-MITTELWERT IN ETWA 500 HPA (KONTROLLE)
C**                ZOM300M: OMEGA-MITTELWERT IN ETWA 300 HPA (KONTROLLE)
C**
C**   COMMON-
C**   BLOECKE  :   PARAM  , ORG, COMDYN, COMPYH, PHYKON, HIGKON, PARKON
C**                PROGCHK, COMDIA, COMPCST, COMPMLF, COMPSLF, COMPGP3
C**
C**   METHODE  :   ZEITLICH:  SEMI-IMPLIZIT/LEAP-FROG
C**                RAEUMLICH: FINITE DIFFERENZEN 2.ORDNUNG
C**                           ARAKAWA-C-GITTER (HORIZONTAL)
C**                           ETA-KOORDINATE (VERTIKAL)
C**   FEHLERBE-
C**   HANDLUNG :      ---
C**   VERFASSER:   G. DOMS UND D. MAJEWSKI
C

      IMPLICIT NONE

      INCLUDE "parorg.h"
      INCLUDE "org.h"
      INCLUDE "comdyn.h"
      INCLUDE "phykon.h"
      INCLUDE "higkon.h"
      INCLUDE "parkon.h"
      INCLUDE "progchk.h"
      INCLUDE "comdia.h"
      INCLUDE "faktinf.h"
C
C

C-----------------------------------------------------------------------

C     ERFORDERLICHE EM-FELDER AUS DEM LANGZEITSPEICHER DIMENSIONIEREN
C     ---------------------------------------------------------------
C     VERTIKAL-KOORDINATEN-PARAMETER
      REAL :: AK(KE1), BK(KE1), AKH(KE), BKH(KE), DAK(KE), DBK(KE)

C     VERTIKAL VARIIERENDER IMPLIZITHEITSGRAD DER V-DIFFUSION
      REAL :: A1T(KE1), A2T(KE1)

C     EVT. VERTIKAL VARIIERENDER VERSTAERKUNGSFAKTOR DER H-DIFFUSION
      REAL :: VVFH(KE)

C     EXTERNE PARAMETER
      REAL :: FIB(IE,JE), FC(IE,JE), ACPHIR(JE,2), CPHI(JE,2)

C     PROGNOSTISCHE BODENFELDER

      REAL :: PS(IE,JE,3),
     &        TG (IE,JE,3), TGL (IE,JE,3), TGW (IE,JE,3), TGI (IE,JE,3),
     &        QDB(IE,JE,3), QDBL(IE,JE,3), QDBW(IE,JE,3), QDBI(IE,JE,3)
C     DIAGNOSTISCHE BODENFELDER
      REAL ::   BFLHS  (IE,JE),
     &          BFLHSL (IE,JE), BFLHSW (IE,JE), BFLHSI (IE,JE),
     &          BFLQDS (IE,JE),
     &          BFLQDSL(IE,JE), BFLQDSW(IE,JE), BFLQDSI(IE,JE),
     &          BFLUS  (IE,JE), BFLVS  (IE,JE)         
C
      INTEGER, INTENT(IN) :: INFRL  (IE,JE), INFRW  (IE,JE), 
     &                       INFRI  (IE,JE)

C     ATMOSPHAEREN-FELDER
      REAL ::   U    (IE,JE,KE,3), V (IE,JE,KE,3),
     &          T    (IE,JE,KE,3), QD(IE,JE,KE,3),
     &          QW   (IE,JE,KE,3), FI(IE,JE,KE,2)

C     FELDER DER KOEFFIZIENTEN UND FLUESSE (PHYSIKALISCHE UPS)
      REAL :: TMKVMH(IE*(KE-1),JE,2)

      REAL ::   TMCM(IE,JE),
     &          TMCH(IE,JE), TMCHL(IE,JE), TMCHW(IE,JE), TMCHI(IE,JE)

      REAL ::   SOTHDT(IEKE,JE,2)

      REAL ::   TTK(IEKE,JE), QDTK(IEKE,JE), UVTK(IEKE,JE,2)

      REAL ::   TTS(IEKE,JE), QDTS(IEKE,JE), QWTS(IEKE,JE)

      REAL, INTENT(IN) :: DWDT(IE,JE,KE ,3)
      REAL, INTENT(INOUT) :: PINT(IE,JE,KE1,3)
      REAL, INTENT(INOUT) :: ETAS(IE,JE,KE1)

C     LOKALE FELDER DIMENSIONIEREN
C     ----------------------------
      REAL    :: PSDT(IE)
      REAL    :: ZTPA(IE,KE,5), ZGQD (IE,KE,3),
     &           ZTV (IE,KE,5), ZFIHF(IE,KE,5)

C     NEUE FELDER FUER GEWICHTETE SIGMA-HORIZONTALDIFFUSION:
      REAL ::   ZPHF(IE,KE,5), ZPNF(IE,KE1,5), ZDP(IE,KE,3),
     &          ZDP5(IE,KE,5), ZDPU(IE,KE ,3), ZDPV(IE,KE,4)

      REAL ::   ZTMKVM(IE,2:KE,3),
     &          ZTMKVH(IE,2:KE,3),
     &          ZTTK  (IE,  KE  ),
     &          ZQDTK (IE,  KE  ),
     &          ZUTK  (IE,  KE,2),
     &          ZVTK  (IE,  KE,2),
     &          ZTTS  (IE,  KE  ),
     &          ZQDTS (IE,  KE  ),
     &          ZQWTS (IE,  KE  ),
     &          ZSODTA(IE,  KE  ),
     &          ZTHDTA(IE,  KE  )

      REAL ::   ZALOPN(IE,KE1,3), ZALOPH(IE,KE ,3),
     &          ZBETAK(IE,KE ,3), ZFIH  (IE,KE ,3),
     &          ZPPHI (IE,KE ,2), ZPLAM (IE,KE   ),
     &          ZGU   (IE,KE ,2), ZGV   (IE,KE ,3),
     &          ZEKIN (IE,KE ,2), ZZETA (IE,KE ,2),
     &          ZETAS (IE,KE1,2), ZSDIV (IE,KE1,2)

      REAL ::   ZLAPT(IE,KE,3), ZLAPQD(IE,KE,3), ZLAPQW(IE,KE,3),
     &          ZLAPU(IE,KE,3), ZLAPV (IE,KE,3)

      REAL ::   ZTADV (IE,KE), ZTDIFH(IE,KE), ZQDDIH(IE,KE),
     &          ZQWDIH(IE,KE), ZALPOM(IE,KE), ZQDADV(IE,KE)

      REAL ::   ZGRAD (IE,KE), ZVZEQ (IE,KE),
     &          ZEDDPQ(IE,KE), ZUDIFH(IE,KE),
     &          ZVDIFH(IE,KE), ZUZEQ (IE,KE)

      REAL ::   ZQKOR(IE), ZFCZET(IE)

C     FELDER FUER SUBROUTINE GAUSS
      REAL ::   AGA(IE,KE,4), AGB(IE,KE,4),
     &          AGC(IE,KE,4), AGD(IE,KE,4),
     &          AGE(IE,KE,4)

C     FELDER FUER SUBROUTINE HQTOTQ
      REAL ::   HE    (IE,KE), QDWE  (IE,KE),
     &          PHFE  (IE,KE),
     &          TSTART(IE,KE), GQDSTA(IE,KE),
     &          PHFSTA(IE,KE),QDLE  (IE,KE),
     &          QDIE  (IE,KE)


      CHARACTER  YTXT*26
C
      REAL ::   QI(IE,JE,KE,3), QITS  (IEKE,JE),
     &          ZQITS(IE,KE  ), ZLAPQI(IE,KE,3), ZQIDIH(IE,KE)

C
      LOGICAL :: LMASSF
      REAL :: Z1, Z2, Z3, Z4, Z4DRERD, ZA1, ZA1A, ZA2, ZA2A, 
     &        ZA3, ZAGAM, ZAGAT, ZAGCM, ZAGCT, ZALOG2, ZDELTF, ZEPSRAY, 
     &        ZFADVX
      REAL :: ZFADVY, ZFDIVO, ZFDIVU, ZFDIVX, ZFKIO, ZFKIU, ZFVZO, 
     &        ZGVCO, ZGVCU, ZOM300M, ZOM500M, ZOM850M, ZQD1, ZQD2, ZQD3, 
     &        ZQD4, ZQW1, ZQW2, ZT1, ZFVZU
      REAL :: ZT2, ZT3, ZT4, ZTKVL, ZTKVZ, ZTMKVHM, ZTRC, ZX1, ZX2, ZXO, 
     &        ZXU, ZZGEW
      INTEGER :: I, INDEX1, INDEX2, INDEX3, INDEX4, J, JM2, JM3, JM4, 
     &           JN5, JO2, JO3, JO4, JO5, JS4, JS5, JSP, JU3, JU4, JM5
      INTEGER :: JU5, JZM, JZN, JZO, JZS, JZU, K
      INTEGER :: KM1, KP1
C
C     STATEMENT - FUNCTIONS ZUR BERECHNUNG DER SAEETIGUNGSFEUCHTE
C     TODO: Find out, why this macro statement gives slightly different
C           results than the function statement below.
C
      REAL :: FGEW,FGQD,TT,GE,PP
C
C     MAGNUS-FORMEL FUER WASSER
      FGEW(TT)        = B1 * EXP  ( B2W*(TT - B3)/(TT - B4W) )
C     SPEZIFISCHE FEUCHTE AUS DAMPFDRUCK UND LUFTDRUCK
      FGQD(GE,PP)     = RDRD*GE/(PP - EMRDRD*GE)


!DIR$ NOTASK
C-----------------------------------------------------------------------
C     VORBEREITENDE MASSNAHMEN
      INDEX1 = 1
      INDEX2 = 2
      INDEX3 = 3
      INDEX4 = 4
C     PARAMETER ZUR BERECHNUNG DER DIFFUSION
C     UND LN(2) SETZEN.
      ALCNVA  = 0.5
      ZA1A    = ALCNVA
      IF(ZA1A.LT.0.5) ZA1A = 0.5
      ZA2A    = (1.-ZA1A)
      ZALOG2  = ALOG(2.)
      Z4DRERD = 4.0/RERD

C     WENN DIE MAXIMALE WINDGESCHWINDIGKEIT 95 % DES STABILITAETS-
C     KRITISCHEN WERTES UEBERSCHREITET, WIRD IN DER U- UND V-GLEICHUNG
C     RAYLEIGH-REIBUNG EINGEFUEHRT, UM DEN WIND ABZUBREMSEN.
      IF ( VBMXV.GT.0.95*VBCFL ) THEN
         ZEPSRAY = 0.0005*ED2DT*(VBMXV - 0.95*VBCFL)/(0.05*VBCFL)
         IF ( VBMXV.GT.VBCFL*1.05 ) THEN
            WRITE (YTXT,'(A,F5.1,A)') 'WARNING: VBMAX=',VBMXV,' M/S'
            CALL REMARK ( YTXT )
            IF (VBMXV.GT.250.0) THEN
               PRINT *,'VBMAX EXCEEDS 250 M/S'
               STOP 1
            ENDIF
         ENDIF
      ELSE
         ZEPSRAY = 0.0
      ENDIF

C     ZTKVZ , ZTKVL : GEWICHTSFAKTOREN ZUR HORIZONTALEN MITTELUNG DER
C                     VERTIKALEN TRANSPORTKOEFFIZIENTEN
      ZTKVZ = 0.9
      ZTKVL = (1.-ZTKVZ)*0.25

C     OMEGA FLAECHENMITTELWERTE VORBEREITEN
      ZOM850M = 0.
      ZOM500M = 0.
      ZOM300M = 0.

C     BERECHNUNG MIT ODER OHNE MASSENFLUSS-KORREKTURSCHEMA
      LMASSF = .TRUE.

C     SETZEN DES STEUERPARAMETERS FUER TROCKENE (QD=QW=0) BZW. FEUCHTE
C     (QD>0, QW>0) PROGNOSE; WENN EIN ADIABATISCHER INITIALISIERUNGS-
C     ZEITSCHRITT GERECHNET WIRD, IST ZTRC=0.0, SONST ZTRC=1.0
      IF(LAISTEP) THEN
         ZTRC = 0.0
      ELSE
         ZTRC = 1.0
      ENDIF

C-----------------------------------------------------------------------

C     EXPLIZITE PROGNOSE OHNE RANDRELAXATION
C     DIE PROGNOSE ERFOLGT 'SCHEIBENWEISE' VON J = JATPROG BIS JETPROG.
C     IN EINEM TASK WERDEN JEWEILS ZUSAMMENHAENGENDE BEREICHE VON
C     J = JATPROG BIS JETPROG BEHANDELT; DIE TASK-EINTEILUNG UND
C     STEUERUNG ERFOLGT IM UP *PROGORG*; MAXIMAL KOENNEN VIER TASKS
C     PARALLEL LAUFEN; D.H. DER BEREICH JAH BIS JEH WIRD DANN IN VIER
C     UNTERBEREICHE GETEILT.

C***********************************************************************
C*                                                                     *
C*    ZUORDNUNG VON SCHEIBENINDIZES  'SUED'--------->-----------'NORD' *
C*    -----------------------------                                    *
C*                                          J-2 J-1      J      J+1 J+2*
C*    U, T, QD, QW, PHY.UP'S                     +       +       +     *
C*    V                                              +       +         *
C*                                                                     *
C*    ZALOPN(K+1) ALOG( P(K+1/2) )              JU3     JM3     JO3    *
C*    ZPHF(K)                               JS5 JU5     JM5     JO5 JN5*
C*    ZDP (K)                                   JU3     JM3     JO3    *
C*    ZPNF(K)                               JS5 JU5     JM5     JO5 JN5*
C*    ZTV (K)                               JS5 JU5     JM5     JO5 JN5*
C*    ZGQD(K)                                   JU3     JM3     JO3    *
C*    ZALOPH(K)   ALOG( P(K) )                  JU3     JM3     JO3    *
C*    ZBETAK(K)   HILFSGROESSE F. FIH, ALOPN    JU3     JM3     JO3    *
C*    ZFIHF(K)    FI(K) AN HAUPTFLAECHEN    JS5 JU3     JM3     JO3 JN5*
C*                FUER PHYS. PARAMETRISIER.                            *
C*    ZFIH(K)     FI(K) AN HAUPTFLAECHEN        JU3     JM3     JO3    *
C*                FUER DRUCKGRADIENTTERM                               *
C*    ZGU(K)      .5*(DP(I+1)+DP(I))*U                  JM2     JO2    *
C*    ZGV(K)      .5*(DP(J+1)+DP(J))*V              JU3     JM3     JO3*
C*    ZZETA(K)    POTENTIELLE VORTICITY             JM2     JO2        *
C*    ZEKIN(K)    KINETISCHE ENERGIE                    JM2     JO2    *
C*    ZPPHI(K)    R*TV*GRADY(LN(P))                 JM2     JO2        *
C*    ZSDIV(K+1)  VERT. SUMME DER DIVERGENZEN           JM2     JO2    *
C*    ZETAS(K+1)  MODIF. VERTIKALGESCHWINDIGK.          JM2     JO2    *
C*    ZSKHH(K)    DP*SKH/DLAM**2 AM H-GP.               JM2     JO2    *
C*    ZSKHU(K)    DP*SKH/DLAM**2 AM U-GP.               JM2     JO2    *
C*    ZSKHV(K)    DP*SKH/DLAM**2 AM V-GP.           JM2     JO2        *
C*    ZSKHZ(K)    DP*SKH/DLAM**2 AM Z-GP.           JM2     JO2        *
C*    ZHJ(K)      H(NJ)                         JU3     JM3     JO3    *
C*    ZQDWJ(K)    QDW(NJ)                       JU3     JM3     JO3    *
C*    ZTMKVM(K)   TMKVM(K)                      JU3     JM3     JO3    *
C*    ZTMKVH(K)   TMKVH(K)                      JU3     JM3     JO3    *
C*    ZUTK(K)     U-TENDENZ (KONVEKTIV)                 JM2     JO2    *
C*    ZVTK(K)     V-TENDENZ (KONVEKTIV)                 JM2     JO2    *
C*                                                                     *
C*    ZPLAM(K)    R*TV*GRADX(LN(P))                      +             *
C*    ZALPOM(K)   ALPHA*OMEGA                            +             *
C*    ZTPA (K)    HP (NA)                  JS5  JU5     JM5     JO5 JN5*
C*    ZQDWA(K)    QDW(NA)                  JS5  JU5     JM5     JO5 JN5*
C*    ZPLAM(K)    DRUCKGRADIENT                          +             *
C*    SOWIE WEITERE LOKALE FELDER OHNE SCHEI-            +             *
C*    BENINDEX SIND IN DER J-FLAECHE DEFINIERT           +             *
C*                                                                     *
C*                                         J-2  J-1      J      J+1 J+2*
C*                                                                     *
C***********************************************************************

C     VORBESETZUNG LOKALER FELDER AM OBER- UND UNTERRAND
C     --------------------------------------------------

      DO I  = IAA , IEA

         ZPNF  (I,1,1) = 0.
         ZPNF  (I,1,2) = 0.
         ZPNF  (I,1,3) = 0.
         ZPNF  (I,1,4) = 0.
         ZPNF  (I,1,5) = 0.
         ZALOPN(I,1,1) = 0.
         ZALOPN(I,1,2) = 0.
         ZALOPN(I,1,3) = 0.
         ZBETAK(I,1,1) = ZALOG2
         ZBETAK(I,1,2) = ZALOG2
         ZBETAK(I,1,3) = ZALOG2
         ZSDIV (I,1,1) = 0.
         ZSDIV (I,1,2) = 0.
         ZETAS (I,1,1) = 0.
         ZETAS (I,1,2) = 0.
         ZETAS(I,KE1,1)= 0.
         ZETAS(I,KE1,2)= 0.
         AGA(I,1,1)    = 0.
         AGA(I,1,2)    = 0.
         AGA(I,1,3)    = 0.
         AGC(I,KE,1)   = 0.
         AGC(I,KE,2)   = 0.
         AGC(I,KE,3)   = 0.
         AGA(I,1,4)    = 0.
         AGC(I,KE,4)   = 0.
      ENDDO

      ETAS(IAA:IEA,:,1:KE1) = 0.

C     ANFANGSINDIZES FUER ZYKLISCHES UMSPEICHERN SETZEN

      JM2 = 1
      JO2 = 2

      JU3 = 1
      JM3 = 2
      JO3 = 3
      
      JS4 = 1
      JU4 = 2
      JM4 = 3
      JO4 = 4

      JS5 = 1
      JU5 = 2
      JM5 = 3
      JO5 = 4
      JN5 = 5

C     LOKALE FELDER AM 'SUEDRAND' (J=JATPROG, JATPROG-1) VORBESETZEN
C     --------------------------------------------------------------

C     ACHTUNG:
C     UM BEIM AUTOTASKING (PARALLER DURCHLAUF VON UP *PROGEXP*) IDENTI-
C     SCHE ERGEBNISSE ZU ERZIELEN, MUSS DARAUF GEACHTET WERDEN, DASS
C     DIE FORMULIERUNG DER TERME FUER DEN 'SUEDRAND' UND FUER DAS INNE-
C     RE DES PROGNOSEGEBIETES EXAKT IDENTISCH SIND. SONST KANN DURCH DIE
C     UNVERMEIDLICHEN RUNDUNGSFEHLER DAS ERGEBNIS JE NACH ANZAHL DER
C     TASKS VERSCHIEDEN SEIN.

      JZS = JAH - 2
      JZU = JAH - 1
      JZM = JAH
      JZO = JAH + 1

C     GEPACKTE FELDER DER PHYSIKALISCHEN TENDENZEN UND KOEFFIZIENTEN
C     AUSPACKEN
      CALL COPYRE(TMKVMH(1,JZU,1),ZTMKVM(1,2,JU3),IE*(KE-1))
      CALL COPYRE(TMKVMH(1,JZU,2),ZTMKVH(1,2,JU3),IE*(KE-1))
      CALL COPYRE(TMKVMH(1,JZM,1),ZTMKVM(1,2,JM3),IE*(KE-1))
      CALL COPYRE(TMKVMH(1,JZM,2),ZTMKVH(1,2,JM3),IE*(KE-1))
      CALL COPYRE(UVTK  (1,JZM,1),ZUTK  (1,1,JM2),IEKE)
      CALL COPYRE(UVTK  (1,JZM,2),ZVTK  (1,1,JM2),IEKE)

      ZX2 = .5*RERD*(CPHI(JZU,1)+CPHI(JZM,1))
      ZXO = CPHI(JZM,1)*EDDPHI
      ZXU = CPHI(JZU,1)*EDDPHI

      DO K = 1 , KE
        KP1 = K + 1
        DO I  = IAA , IEA
          ZDP   (I,K  ,JU3) = DAK(K  ) + DBK(K  )*PS(I,JZU,NJ)
          ZDP   (I,K  ,JM3) = DAK(K  ) + DBK(K  )*PS(I,JZM,NJ)
C     DEFINITION VON ZDP5 FUER GEWICHTETE SIGMA HDIFFUSION:
C==============================================================
C    !!!! Should ZDP5 be based on PINT?
C==============================================================
          ZDP5  (I,K  ,JS5) = DAK(K  ) + DBK(K  )*PS(I,JZS,NA)
          ZDP5  (I,K  ,JU5) = DAK(K  ) + DBK(K  )*PS(I,JZU,NA)
          ZDP5  (I,K  ,JM5) = DAK(K  ) + DBK(K  )*PS(I,JZM,NA)
          ZDP5  (I,K  ,JO5) = DAK(K  ) + DBK(K  )*PS(I,JZO,NA)
C     *************
          ZPHF  (I,K  ,JS5) = 0.5 *
     &         ( PINT(I,JZS,K,NJ) + PINT(I,JZS,KP1,NJ))
          ZPHF  (I,K  ,JU5) = 0.5 *
     &         ( PINT(I,JZU,K,NJ) + PINT(I,JZU,KP1,NJ))
          ZPHF  (I,K  ,JM5) = 0.5 *
     &         ( PINT(I,JZM,K,NJ) + PINT(I,JZM,KP1,NJ))
          ZPHF  (I,K  ,JO5) = 0.5 *
     &         ( PINT(I,JZO,K,NJ) + PINT(I,JZO,KP1,NJ))
          ZPNF  (I,KP1,JS5) = PINT(I,JZS,KP1,NJ)
          ZPNF  (I,KP1,JU5) = PINT(I,JZU,KP1,NJ)
          ZPNF  (I,KP1,JM5) = PINT(I,JZM,KP1,NJ)
          ZPNF  (I,KP1,JO5) = PINT(I,JZO,KP1,NJ)

          ZALOPN(I,KP1,JU3) = ALOG(ZPNF(I,KP1,JU5))
          ZALOPN(I,KP1,JM3) = ALOG(ZPNF(I,KP1,JM5))
          ZTV   (I,K  ,JS5) =     T(I,JZS,K,NJ) * ( 1.0 + RDDRM1*
     &         QD(I,JZS,K,NJ) - (QW(I,JZS,K,NJ)+QI(I,JZS,K,NJ)))
          ZTV   (I,K  ,JU5) =     T(I,JZU,K,NJ) * ( 1.0 + RDDRM1*
     &         QD(I,JZU,K,NJ) - (QW(I,JZU,K,NJ)+QI(I,JZU,K,NJ)) )
          ZTV   (I,K  ,JM5) =     T(I,JZM,K,NJ) * ( 1.0 + RDDRM1*
     &         QD(I,JZM,K,NJ) - (QW(I,JZM,K,NJ)+QI(I,JZM,K,NJ)) )
          ZTV   (I,K  ,JO5) =     T(I,JZO,K,NJ) * ( 1.0 + RDDRM1*
     &         QD(I,JZO,K,NJ) - (QW(I,JZO,K,NJ)+QI(I,JZO,K,NJ)) )
          ZZGEW             = FGEW ( T(I,JZU,K,NJ) )
          ZGQD (I,K,JU3)    = FGQD ( ZZGEW, ZPHF(I,K,JU5) )
          ZZGEW             = FGEW ( T(I,JZM,K,NJ) )
          ZGQD (I,K,JM3)    = FGQD ( ZZGEW, ZPHF(I,K,JM5) )
        ENDDO
      ENDDO

C     DEFINITION VON ZDPU UND ZDPV FUER GEWICHTETE SIGMA HDIFFUSION
      DO K = 1 , KE
C     FEHLERKORREKTUR FUER ZDPU-DEFINITION
         DO I = IAA , IEH+1
            ZDPU (I,K,JU3) = 0.5*(ZDP5(I,K,JU5) + ZDP5(I+1,K,JU5))
            ZDPU (I,K,JM3) = 0.5*(ZDP5(I,K,JM5) + ZDP5(I+1,K,JM5))
         ENDDO
      ENDDO
C
CKS   SPLITTED THE FOLLOWING LOOP INTO TWO TO REMOVE THE IF INSIDE A LOOP
C
CKS   FIRST PART WITH FOR ALL K<KE
      DO K = 1 , KE-1
         DO I = IAA , IEA
            ZDPV (I,K,JS4) = 0.5*(ZDP5(I,K,JS5) + ZDP5(I,K,JU5))
            ZDPV (I,K,JU4) = 0.5*(ZDP5(I,K,JU5) + ZDP5(I,K,JM5))
            ZDPV (I,K,JM4) = 0.5*(ZDP5(I,K,JM5) + ZDP5(I,K,JO5))
C           ***************
            ZALPOM(I,K) = 0.0
            ZFIHF(I,K,JS5) = FI(I,JZS,K+1,NA2) + R*ZTV(I,K,JS5)*
     &           ALOG  ( ZPNF(I,K+1,JS5)/ZPHF(I,K,JS5) ) /
     &           DWDT(I,JZS,K,NJ)
            ZFIHF(I,K,JU5) = FI(I,JZU,K+1,NA2) + R*ZTV(I,K,JU5)*
     &           ALOG  ( ZPNF(I,K+1,JU5)/ZPHF(I,K,JU5) ) /
     &           DWDT(I,JZU,K,NJ)
            ZFIHF(I,K,JM5) = FI(I,JZM,K+1,NA2) + R*ZTV(I,K,JM5)*
     &           ALOG  ( ZPNF(I,K+1,JM5)/ZPHF(I,K,JM5) ) /
     &           DWDT(I,JZM,K,NJ)
            ZFIHF(I,K,JO5) = FI(I,JZO,K+1,NA2) + R*ZTV(I,K,JO5)*
     &           ALOG  ( ZPNF(I,K+1,JO5)/ZPHF(I,K,JO5) ) /
     &           DWDT(I,JZO,K,NJ)
            ZTPA (I,K,JS5) = T(I,JZS,K,NA) + ZFIHF(I,K,JS5)*WCPR
            ZTPA (I,K,JU5) = T(I,JZU,K,NA) + ZFIHF(I,K,JU5)*WCPR
            ZTPA (I,K,JM5) = T(I,JZM,K,NA) + ZFIHF(I,K,JM5)*WCPR
            ZTPA (I,K,JO5) = T(I,JZO,K,NA) + ZFIHF(I,K,JO5)*WCPR
         ENDDO
      ENDDO
C
CKS   SECOND PART FOR K=KE
C
      DO I = IAA , IEA
         ZDPV (I,KE,JS4) = 0.5*(ZDP5(I,KE,JS5) + ZDP5(I,KE,JU5))
         ZDPV (I,KE,JU4) = 0.5*(ZDP5(I,KE,JU5) + ZDP5(I,KE,JM5))
         ZDPV (I,KE,JM4) = 0.5*(ZDP5(I,KE,JM5) + ZDP5(I,KE,JO5))
C        ***************
         ZALPOM(I,KE) = 0.0
         ZFIHF(I,KE,JS5) = FIB(I,JZS)        + R*ZTV(I,KE,JS5)*
     &        ALOG  ( ZPNF(I,KE+1,JS5)/ZPHF(I,KE,JS5) ) /
     &        DWDT(I,JZS,KE,NJ)
         ZFIHF(I,KE,JU5) = FIB(I,JZU)        + R*ZTV(I,KE,JU5)*
     &        ALOG  ( ZPNF(I,KE+1,JU5)/ZPHF(I,KE,JU5) ) /
     &        DWDT(I,JZU,KE,NJ)
         ZFIHF(I,KE,JM5) = FIB(I,JZM)        + R*ZTV(I,KE,JM5)*
     &        ALOG  ( ZPNF(I,KE+1,JM5)/ZPHF(I,KE,JM5) ) /
     &        DWDT(I,JZM,KE,NJ)
         ZFIHF(I,KE,JO5) = FIB(I,JZO)        + R*ZTV(I,KE,JO5)*
     &        ALOG  ( ZPNF(I,KE+1,JO5)/ZPHF(I,KE,JO5) ) /
     &        DWDT(I,JZO,KE,NJ)
         ZTPA (I,KE,JS5) = T(I,JZS,KE,NA) + ZFIHF(I,KE,JS5)*WCPR
         ZTPA (I,KE,JU5) = T(I,JZU,KE,NA) + ZFIHF(I,KE,JU5)*WCPR
         ZTPA (I,KE,JM5) = T(I,JZM,KE,NA) + ZFIHF(I,KE,JM5)*WCPR
         ZTPA (I,KE,JO5) = T(I,JZO,KE,NA) + ZFIHF(I,KE,JO5)*WCPR
      ENDDO
CKS

C
CKS   SPLITTED THE FOLLOWING LOOP INTO TWO TO REMOVE THE IF INSIDE A LOOP
C
CKS   FIRST PART FOR K=1
      DO I = IAA, IEA
         ZALOPH(I,1 ,JU3) = ZALOPN(I,2,JU3) - ZBETAK(I,1 ,JU3)
         ZALOPH(I,1 ,JM3) = ZALOPN(I,2,JM3) - ZBETAK(I,1 ,JM3)
         ZFIH(I,1,JU3) = FI(I,JZU,2,NJ2) +
     &        R*ZTV(I,1,JU5)*ZBETAK(I,1,JU3) / DWDT(I,JZU,1,NJ)
         ZFIH(I,1,JM3) = FI(I,JZM,2,NJ2) +
     &        R*ZTV(I,1,JM5)*ZBETAK(I,1,JM3) / DWDT(I,JZM,1,NJ)
      ENDDO

CKS   SECOND PART FOR K>1
      DO K = 2, KE-1
         DO I = IAA , IEA
            ZBETAK(I,K,JU3) = 1. - ZPNF(I,K  ,JU5)/ZDP   (I,K,JU3)
     &           * (ZALOPN(I,K+1,JU3)-ZALOPN(I,K,JU3))
            ZBETAK(I,K,JM3) = 1. - ZPNF(I,K  ,JM5)/ZDP   (I,K,JM3)
     &           * (ZALOPN(I,K+1,JM3)-ZALOPN(I,K,JM3))
            ZFIH(I,K  ,JU3) = FI(I,JZU,K+1,NJ2) +  R*ZTV(I,K,JU5)
     &           * ZBETAK(I,K,JU3) / DWDT(I,JZU,K,NJ)
            ZFIH(I,K  ,JM3) = FI(I,JZM,K+1,NJ2) +  R*ZTV(I,K,JM5)
     &           *ZBETAK(I,K,JM3) / DWDT(I,JZM,K,NJ)
            ZALOPH(I,K,JU3) = ZALOPN(I,K+1,JU3) - ZBETAK(I,K,JU3)
            ZALOPH(I,K,JM3) = ZALOPN(I,K+1,JM3) - ZBETAK(I,K,JM3)
         ENDDO
      ENDDO
CKS

C     CORIOLISPARAMETER AM ZETA-PUNKT BEREITSTELLEN
      DO I  = IAA, IEH + 1
         ZFCZET(I) = 0.25*( FC(I,JZU) + FC(I+1,JZU) +
     &                      FC(I,JZM) + FC(I+1,JZM) )
      ENDDO

      DO I = IAA , IEA
         ZBETAK(I,KE,JU3) = 1. - ZPNF(I,KE ,JU5) / ZDP   (I,KE,JU3)
     &        * (ZALOPN(I,KE1,JU3) - ZALOPN(I,KE,JU3))
         ZBETAK(I,KE,JM3) = 1. - ZPNF(I,KE ,JM5) / ZDP   (I,KE,JM3)
     &        * (ZALOPN(I,KE1,JM3) - ZALOPN(I,KE,JM3))
         ZALOPH(I,KE,JU3) =    ZALOPN(I,KE1,JU3) - ZBETAK(I,KE,JU3)
         ZALOPH(I,KE,JM3) =    ZALOPN(I,KE1,JM3) - ZBETAK(I,KE,JM3)
         ZFIH  (I,KE,JU3) = FIB(I,JZU) +
     &        R*ZTV(I,KE,JU5)*ZBETAK(I,KE,JU3) / DWDT(I,JZU,KE,NJ)
         ZFIH  (I,KE,JM3) = FIB(I,JZM) +
     &        R*ZTV(I,KE,JM5)*ZBETAK(I,KE,JM3) / DWDT(I,JZM,KE,NJ)
      ENDDO

C     WEITERE LOKALE FELDER AN U/ZETA-GITTERPUNKTEN (HF) VORBESETZEN

      DO K = 1, KE
         DO I = IAA , IEH + 1
            ZGU(I,K,JM2) = .5*(ZDP(I,K,JM3)+
     &           ZDP(I+1,K,JM3))*U(I,JZM,K,NJ)
            Z1 = ZX2*ZFCZET(I)
            Z2 = EDDLAM*( V(I+1,JZU,K,NJ) -     V(I,JZU,K,NJ) )
            Z3 = ZXO   *  U(I  ,JZM,K,NJ) - ZXU*U(I,JZU,K,NJ)
            Z4 = (ZDP(I,K,JU3)+ZDP(I+1,K,JU3))*CPHI(JZU,1) +
     &           (ZDP(I,K,JM3)+ZDP(I+1,K,JM3))*CPHI(JZM,1)
            ZZETA(I,K,JM2) = Z4DRERD * ( Z1 + Z2 - Z3 )/ Z4
         ENDDO

C     WEITERE LOKALE FELDER AN H/V-GITTERPUNKTEN VORBESETZEN
         ZX1    = 0.5*R*EDADPHI
         ZFKIO  = CPHI(JZM,2)/CPHI(JZM,1)
         ZFKIU  = CPHI(JZU,2)/CPHI(JZM,1)
         ZFDIVX = ACPHIR(JZM,1)*EDDLAM
         ZFDIVO = ACPHIR(JZM,1)*EDDPHI*CPHI(JZM,2)
         ZFDIVU = ACPHIR(JZM,1)*EDDPHI*CPHI(JZU,2)

         DO I = IAH - 1 , IEH + 1
            ZGV(I,K,JU3) = 
     &           .5*(ZDP(I,K,JU3)+(DAK(K)+DBK(K)*PS(I,JZM  ,NJ)))*
     &                                            V(I,JZU,K,NJ)
            ZGV(I,K,JM3) =
     &           .5*(ZDP(I,K,JM3)+(DAK(K)+DBK(K)*PS(I,JZM+1,NJ)))*
     &                                            V(I,JZM,K,NJ)
            ZPPHI(I,K,JM2) = ZX1*(ZTV(I,K,JM5) + ZTV   (I,K,JU5))*
     &           (ZALOPH(I,K,JM3) - ZALOPH(I,K,JU3))
            ZEKIN(I,K,JM2) = 
     &           .25*(  U(I-1,JZM,K,NJ)**2 + U(I,JZM,K,NJ)**2 +
     &           ZFKIU*V(I,JZU,K,NJ)**2 + ZFKIO*V(I,JZM,K,NJ)**2)
         ENDDO
      ENDDO

      DO K = 1 , KE
         DO I = IAH - 1, IEH + 1
            ZSDIV(I,K+1,JM2) = ZSDIV(I,K,JM2) +
     &           ZFDIVX*( ZGU(I,K,JM2)-ZGU(I-1,K,JM2) ) +
     &           ZFDIVO*ZGV(I,K,JM3) - ZFDIVU*ZGV(I,K,JU3)
         ENDDO
      ENDDO

      DO K = 2 , KE
        DO I = IAH - 1, IEH + 1
          ZETAS(I,K,JM2) = BK(K)*ZSDIV(I,KE1,JM2) - ZSDIV(I,K,JM2)
        ENDDO
        ETAS(IAH-1:IEH+1,JZM,K)  = ZETAS(IAH-1:IEH+1,K,JM2)
      ENDDO

C     LOKALE HILFSFELDER FUER HORIZONTALDIFFUSION
      DO K = 1 , KE
         DO I = IAH - 1, IEH + 1

C     GEWICHTUNG VON ZLAPT, ZLAPQD, ZLAPQW MIT GEWICHT DER SCHICHT
            ZLAPT(I,K,JU3) = 
     &             ((ZTPA(I+1,K,JU5)-ZTPA(I  ,K,JU5))*ZDPU(I  ,K,JU3)-
     &              (ZTPA(I  ,K,JU5)-ZTPA(I-1,K,JU5))*ZDPU(I-1,K,JU3)+
     & (CPHI(JZU,2)*(ZTPA(I  ,K,JM5)-ZTPA(I  ,K,JU5))*ZDPV(I  ,K,JU4)-
     &  CPHI(JZS,2)*(ZTPA(I  ,K,JU5)-ZTPA(I  ,K,JS5))*ZDPV(I  ,K,JS4))/
     &  CPHI(JZU,1)) / ZDP5(I,K,JU5)
            ZLAPQD(I,K,JU3) = 
     &          ((QD(I+1,JZU,K,NA)-QD(I  ,JZU,K,NA))*ZDPU(I  ,K,JU3)-
     &           (QD(I  ,JZU,K,NA)-QD(I-1,JZU,K,NA))*ZDPU(I-1,K,JU3)+
     & (CPHI(JZU,2)*
     &           (QD(I  ,JZM,K,NA)-QD(I  ,JZU,K,NA))*ZDPV(I  ,K,JU4)-
     &  CPHI(JZS,2)*
     &           (QD(I  ,JZU,K,NA)-QD(I  ,JZS,K,NA))*ZDPV(I  ,K,JS4))/ 
     &  CPHI(JZU,1)) / ZDP5(I,K,JU5)
            ZLAPQW(I,K,JU3) =
     &          ((QW(I+1,JZU,K,NA)-QW(I  ,JZU,K,NA))*ZDPU(I  ,K,JU3)-
     &           (QW(I  ,JZU,K,NA)-QW(I-1,JZU,K,NA))*ZDPU(I-1,K,JU3)+
     & (CPHI(JZU,2)*
     &           (QW(I  ,JZM,K,NA)-QW(I  ,JZU,K,NA))*ZDPV(I  ,K,JU4)-
     &  CPHI(JZS,2)*
     &           (QW(I  ,JZU,K,NA)-QW(I  ,JZS,K,NA))*ZDPV(I  ,K,JS4))/ 
     &  CPHI(JZU,1)) / ZDP5(I,K,JU5)
            ZLAPQI(I,K,JU3) =
     &          ((QI(I+1,JZU,K,NA)-QI(I  ,JZU,K,NA))*ZDPU(I  ,K,JU3)-
     &           (QI(I  ,JZU,K,NA)-QI(I-1,JZU,K,NA))*ZDPU(I-1,K,JU3)+
     & (CPHI(JZU,2)*
     &           (QI(I  ,JZM,K,NA)-QI(I  ,JZU,K,NA))*ZDPV(I  ,K,JU4)-
     &  CPHI(JZS,2)*
     &           (QI(I  ,JZU,K,NA)-QI(I  ,JZS,K,NA))*ZDPV(I  ,K,JS4))/ 
     &  CPHI(JZU,1)) / ZDP5(I,K,JU5)
C     ***************
            ZLAPU (I,K,JU3) = U (I+1,JZU,K,NA) + U (I-1,JZU,K,NA) -
     &                    2.0*U (I  ,JZU,K,NA) +
     &         ( CPHI(JZU,2)*(U (I  ,JZM,K,NA) - U (I  ,JZU,K,NA)) -
     &           CPHI(JZS,2)*(U (I  ,JZU,K,NA) - U (I  ,JZS,K,NA)) )/
     &           CPHI(JZU,1)
            ZLAPV (I,K,JU3) = V (I+1,JZU,K,NA) + V (I-1,JZU,K,NA) -
     &                    2.0*V (I  ,JZU,K,NA) +
     &         ( CPHI(JZM,1)*(V (I  ,JZM,K,NA) - V (I  ,JZU,K,NA)) -
     &           CPHI(JZU,1)*(V (I  ,JZU,K,NA) - V (I  ,JZS,K,NA)) )/
     &           CPHI(JZU,2)

            ZLAPT(I,K,JM3) =
     &             ((ZTPA(I+1,K,JM5)-ZTPA(I  ,K,JM5))*ZDPU(I  ,K,JM3)-
     &              (ZTPA(I  ,K,JM5)-ZTPA(I-1,K,JM5))*ZDPU(I-1,K,JM3)+
     & (CPHI(JZM,2)*(ZTPA(I  ,K,JO5)-ZTPA(I  ,K,JM5))*ZDPV(I  ,K,JM4)-
     &  CPHI(JZU,2)*(ZTPA(I  ,K,JM5)-ZTPA(I  ,K,JU5))*ZDPV(I  ,K,JU4))/
     &  CPHI(JZM,1)) / ZDP5(I,K,JM5)
            ZLAPQD(I,K,JM3) = 
     &          ((QD(I+1,JZM,K,NA)-QD(I  ,JZM,K,NA))*ZDPU(I  ,K,JM3)-
     &           (QD(I  ,JZM,K,NA)-QD(I-1,JZM,K,NA))*ZDPU(I-1,K,JM3)+
     & (CPHI(JZM,2)*
     &           (QD(I  ,JZO,K,NA)-QD(I  ,JZM,K,NA))*ZDPV(I  ,K,JM4)-
     &  CPHI(JZU,2)*
     &           (QD(I  ,JZM,K,NA)-QD(I  ,JZU,K,NA))*ZDPV(I  ,K,JU4))/ 
     &  CPHI(JZM,1)) / ZDP5(I,K,JM5)
            ZLAPQW(I,K,JM3) =
     &          ((QW(I+1,JZM,K,NA)-QW(I  ,JZM,K,NA))*ZDPU(I  ,K,JM3)-
     &           (QW(I  ,JZM,K,NA)-QW(I-1,JZM,K,NA))*ZDPU(I-1,K,JM3)+
     & (CPHI(JZM,2)*
     &           (QW(I  ,JZO,K,NA)-QW(I  ,JZM,K,NA))*ZDPV(I  ,K,JM4)-
     &  CPHI(JZU,2)*
     &           (QW(I  ,JZM,K,NA)-QW(I  ,JZU,K,NA))*ZDPV(I  ,K,JU4))/ 
     &  CPHI(JZM,1)) / ZDP5(I,K,JM5)
            ZLAPQI(I,K,JM3) =
     &          ((QI(I+1,JZM,K,NA)-QI(I  ,JZM,K,NA))*ZDPU(I  ,K,JM3)-
     &           (QI(I  ,JZM,K,NA)-QI(I-1,JZM,K,NA))*ZDPU(I-1,K,JM3)+
     & (CPHI(JZM,2)*
     &           (QI(I  ,JZO,K,NA)-QI(I  ,JZM,K,NA))*ZDPV(I  ,K,JM4)-
     &  CPHI(JZU,2)*
     &           (QI(I  ,JZM,K,NA)-QI(I  ,JZU,K,NA))*ZDPV(I  ,K,JU4))/ 
     &  CPHI(JZM,1)) / ZDP5(I,K,JM5)
            ZLAPU (I,K,JM3) = U (I+1,JZM,K,NA) + U (I-1,JZM,K,NA) -
     &                    2.0*U (I  ,JZM,K,NA) +
     &         ( CPHI(JZM,2)*(U (I  ,JZO,K,NA) - U (I  ,JZM,K,NA)) -
     &           CPHI(JZU,2)*(U (I  ,JZM,K,NA) - U (I  ,JZU,K,NA)) )/
     &           CPHI(JZM,1)
            ZLAPV (I,K,JM3) = V (I+1,JZM,K,NA) + V (I-1,JZM,K,NA) -
     &                    2.0*V (I  ,JZM,K,NA) +
     &         ( CPHI(JZO,1)*(V (I  ,JZO,K,NA) - V (I  ,JZM,K,NA)) -
     &           CPHI(JZM,1)*(V (I  ,JZM,K,NA) - V (I  ,JZU,K,NA)) )/
     &           CPHI(JZM,2)
         ENDDO
      ENDDO

C     RAENDER DER SCHEIBENFELDER FUER DIE HQ-ZERLEGUNG VORBESETZEN;
C     (BELIEBIGE KONSISTENTE WERTE)
*VDIR NOVECTOR
      DO I  = 1 , 2
         DO K  = 1 , KE
            HE    (I,K) = WCP*T(1,JZU,K,NJ) + WLK*QD(1,JZU,K,NJ)
            QDWE  (I,K) =    QD(1,JZU,K,NJ) +     QW(1,JZU,K,NJ)
     &                     + QI(1,JZU,K,NJ)
            TSTART(I,K) =     T(1,JZU,K,NJ)
            GQDSTA(I,K) = ZGQD (1,K,JU3)
            PHFSTA(I,K) = ZPHF (1,K,JU5)
            PHFE  (I,K) = ZPHF (1,K,JU5)
         ENDDO
      ENDDO

*VDIR NOVECTOR
      DO I = IE - 1 , IE
         DO K = 1, KE
            HE    (I,K) = HE    (1,K)
            QDWE  (I,K) = QDWE  (1,K)
            TSTART(I,K) = TSTART(1,K)
            GQDSTA(I,K) = GQDSTA(1,K)
            PHFSTA(I,K) = PHFSTA(1,K)
            PHFE  (I,K) = PHFE  (1,K)
         ENDDO
      ENDDO


C     BEGINN DER SCHLEIFE IN J-RICHTUNG
C     ---------------------------------

C     ACHTUNG:
C     UM BEIM AUTOTASKING (PARALLER DURCHLAUF VON UP *PROGEXP*) IDENTI-
C     SCHE ERGEBNISSE ZU ERZIELEN, MUSS DARAUF GEACHTET WERDEN, DASS
C     DIE FORMULIERUNG DER TERME FUER DEN 'SUEDRAND' UND FUER DAS INNE-
C     RE DES PROGNOSEGEBIETES EXAKT IDENTISCH SIND. SONST KANN DURCH DIE
C     UNVERMEIDLICHEN RUNDUNGSFEHLER DAS ERGEBNIS JE NACH ANZAHL DER
C     TASKS VERSCHIEDEN SEIN.

      DO J = JAH, JEH

         JZN = J+2

C        LOKALE FELDER FUER J   BESETZEN
C        GEPACKTE FELDER DER PHYSIKALISCHEN TENDENZEN UND KOEFFIZIENTEN
C        AUSPACKEN
         CALL COPYRE(TTK(1,J)   ,ZTTK  ,IEKE)
         CALL COPYRE(QDTK(1,J)  ,ZQDTK ,IEKE)
         CALL COPYRE(TTS(1,J)   ,ZTTS  ,IEKE)
         CALL COPYRE(QDTS(1,J)  ,ZQDTS ,IEKE)
         CALL COPYRE(QWTS(1,J)  ,ZQWTS ,IEKE)
         CALL COPYRE(SOTHDT(1,J,1), ZSODTA, IEKE)
         CALL COPYRE(SOTHDT(1,J,2), ZTHDTA, IEKE)
         CALL COPYRE(QITS(1,J)  ,ZQITS , IEKE)
C        LOKALE FELDER FUER J+1 BESETZEN
C        GEPACKTE FELDER DER PHYSIKALISCHEN TENDENZEN UND KOEFFIZIENTEN
C        AUSPACKEN
         CALL COPYRE(TMKVMH(1,J+1,1),ZTMKVM(1,2,JO3),IE*(KE-1))
         CALL COPYRE(TMKVMH(1,J+1,2),ZTMKVH(1,2,JO3),IE*(KE-1))
         CALL COPYRE(UVTK  (1,J+1,1),ZUTK  (1,1,JO2),IEKE)
         CALL COPYRE(UVTK  (1,J+1,2),ZVTK  (1,1,JO2),IEKE)
C        LOKALE FELDER FUER J+2 BESETZEN

         ZX1 = .5*R*ACPHIR(J,1)*EDDLAM
         ZX2 = .5*RERD*(CPHI(J,1)+CPHI(J+1,1))
         ZXO = CPHI(J+1,1)*EDDPHI
         ZXU = CPHI(J  ,1)*EDDPHI

         DO K = 1 , KE
           KP1 = K + 1
           DO I = IAA , IEA
             ZDP  (I,K,JO3)    = DAK(K  ) + DBK(K  )*PS(I,J+1,NJ)
C     FUER GEWICHTETE HDIFF
C==============================================================
C    !!!! Should ZDP5 be based on PINT?
C==============================================================
             ZDP5 (I,K,JN5)    = DAK(K  ) + DBK(K  )*PS(I,JZN,NA)
             ZDPV (I,K,JO4)  = 0.5*(ZDP5(I,K,JO5) + ZDP5(I,K,JN5))
C     ************
             ZPHF (I,K,JN5)    = 0.5 *
     &            ( PINT(I,JZN,K,NJ) + PINT(I,JZN,KP1,NJ))
             ZPNF (I,KP1,JN5)  = PINT(I,JZN,KP1,NJ)
             ZALOPN(I,KP1,JO3) = ALOG(ZPNF(I,KP1,JO5))
             ZTV   (I,K  ,JN5) =     T(I,JZN,K,NJ) * ( 1.0 + RDDRM1*
     1            QD(I,JZN,K,NJ) - (QW(I,JZN,K,NJ) + QI(I,JZN,K,NJ)) )
             ZZGEW             = FGEW ( T(I,J+1,K,NJ) )
             ZGQD (I,K,JO3)    = FGQD ( ZZGEW, ZPHF(I,K,JO5) )
           ENDDO
         ENDDO

C        FUER GEWICHTETE HDIFF
         DO K = 1 , KE
C           FEHLERKORREKTUR FUER ZDPU
            DO I = IAA , IEH+1
               ZDPU (I,K,JO3) = 0.5*(ZDP5(I,K,JO5) + ZDP5(I+1,K,JO5))
            ENDDO
         ENDDO
C        ***************
C
CKS      SPLITTED LOOP TO REMOVE IF INSIDE LOOP
C
CKS      FIRST PART K<KE
         DO K = 1 , KE-1
            DO I = IAA , IEA
               ZFIHF(I,K,JN5) = FI(I,JZN,K+1,NA2) + R*ZTV(I,K,JN5)*
     &             ALOG  ( ZPNF(I,K+1,JN5)/ZPHF(I,K,JN5) ) /
     &             DWDT(I,JZN,K,NJ)
               ZTPA (I,K,JN5)    = T(I,JZN,K,NA) + ZFIHF(I,K,JN5)*WCPR
            ENDDO
         ENDDO
CKS      SECOND PART K=KE
         DO I = IAA , IEA
            ZFIHF(I,KE,JN5) = FIB(I,JZN)        + R*ZTV(I,KE,JN5)*
     &          ALOG  ( ZPNF(I,KE+1,JN5)/ZPHF(I,KE,JN5) ) /
     &          DWDT(I,JZN,KE,NJ)
            ZTPA (I,KE,JN5)    = T(I,JZN,KE,NA) + ZFIHF(I,KE,JN5)*WCPR
         ENDDO
CKS
C
CKS      SPLITTED LOOP TO REMOVE IF INSIDE LOOP
C
CKS      FIRST PART K=1
         DO I  = IAA , IEA
            ZALOPH(I,1,JO3) = ZALOPN(I,2,JO3) -ZBETAK(I,1,JO3)
            ZFIH (I,1,JO3) = FI(I,J+1,2,NJ2) + R*ZTV(I,1,JO5)
     &           *ZBETAK(I,1,JO3) / DWDT(I,J+1,1,NJ)
         ENDDO
CKS      SECOND PART K>1
         DO K = 2 , KE-1
            DO I  = IAA , IEA
               ZBETAK(I,K,JO3) = 1. - ZPNF(I,K  ,JO5)/ZDP   (I,K,JO3)
     &              * (ZALOPN(I,K+1,JO3)-ZALOPN(I,K,JO3))
               ZFIH  (I,K,JO3) = FI(I,J+1,K+1,NJ2) + R*ZTV (I,K,JO5)
     &              *ZBETAK(I,K,JO3) / DWDT(I,J+1,K,NJ)
               ZALOPH(I,K,JO3) = ZALOPN(I,K+1,JO3) - ZBETAK(I,K,JO3)
            ENDDO
         ENDDO

C        CORIOLISPARAMETER AM ZETA-PUNKT BEREITSTELLEN
         DO I = IAA, IEH + 1
            ZFCZET(I) = 0.25*( FC(I,J  ) + FC(I+1,J  ) +
     &                         FC(I,J+1) + FC(I+1,J+1) )
         ENDDO

         DO I = IAA , IEA
            ZBETAK(I,KE,JO3) = 1. - ZPNF(I,KE ,JO5)  /ZDP   (I,KE,JO3)
     &           * (ZALOPN(I,KE1,JO3) - ZALOPN(I,KE,JO3))
            ZALOPH(I,KE,JO3) =    ZALOPN(I,KE1,JO3) - ZBETAK(I,KE,JO3)
            ZFIH  (I,KE,JO3) = FIB(I,J+1) + R*ZTV(I,KE,JO5)
     &           *ZBETAK(I,KE,JO3) / DWDT(I,J+1,KE,NJ)
         ENDDO

C        WEITERE LOKALE FELDER AN U/ZETA-GITTERPUNKTEN (HF) BESETZEN , JO

         DO K = 1 , KE
            DO I = IAA , IEH + 1
               ZGU(I,K,JO2) = .5*(ZDP(I,K,JO3)+ZDP(I+1,K,JO3))
     &              *U(I,J+1,K,NJ)
               ZPLAM(I,K)   = ZX1*( ZTV(I+1,K,JM5) +    ZTV(I,K,JM5) )
     &              *(ZALOPH(I+1,K,JM3) - ZALOPH(I,K,JM3))
               Z1          = ZX2*ZFCZET(I)
               Z2          = EDDLAM*( V(I+1,J,K,NJ) -     V(I,J,K,NJ) )
               Z3          = ZXO   *  U(I,J+1,K,NJ) - ZXU*U(I,J,K,NJ)
               Z4          = (ZDP(I,K,JM3)+ZDP(I+1,K,JM3))*CPHI(J,  1) +
     &              (ZDP(I,K,JO3)+ZDP(I+1,K,JO3))*CPHI(J+1,1)
               ZZETA(I,K,JO2) = Z4DRERD * ( Z1 + Z2 - Z3 )/ Z4
            ENDDO
         ENDDO

C        WEITERE LOKALE FELDER AN H/V-GITTERPUNKTEN (HF) BESETZEN

         ZX1 = 0.5*R*EDADPHI
         IF((J.LT.JEH) .OR. (NEIGHBOR(2) .NE. -1)) THEN
            ZFKIO  = CPHI(J+1,2)/CPHI(J+1,1)
            ZFKIU  = CPHI(J  ,2)/CPHI(J+1,1)
            ZFDIVX = ACPHIR(J+1,1)*EDDLAM
            ZFDIVO = ACPHIR(J+1,1)*EDDPHI*CPHI(J+1,2)
            ZFDIVU = ACPHIR(J+1,1)*EDDPHI*CPHI(J  ,2)

            DO K = 1 , KE
               DO I = IAH - 1, IEH + 1
                  ZGV(I,K,JO3) =
     &                 .5*(ZDP(I,K,JO3)+(DAK(K)+DBK(K)*PS(I,J+2,NJ)))*
     &                                                  V(I,J+1,K,NJ)
                  ZPPHI(I,K,JO2) = ZX1*(ZTV(I,K,JO5) +    ZTV(I,K,JM5))*
     &                              (ZALOPH(I,K,JO3) - ZALOPH(I,K,JM3))
                  ZEKIN(I,K,JO2) =
     &                    .25*(  U(I-1,J+1,K,NJ)**2 + U(I,J+1,K,NJ)**2+
     &                 ZFKIU*V(I,J  ,K,NJ)**2 + ZFKIO*V(I,J+1,K,NJ)**2)
               ENDDO
            ENDDO

            DO K = 1 , KE
               DO I = IAH - 1, IEH + 1
                  ZSDIV(I,K+1,JO2) = ZSDIV(I,K,JO2) +
     &                 ZFDIVX*( ZGU(I,K,JO2)-ZGU(I-1,K,JO2) ) +
     &                 ZFDIVO*ZGV(I,K,JO3) - ZFDIVU*ZGV(I,K,JM3)
               ENDDO
            ENDDO

            DO K = 2 , KE
               DO I = IAH - 1, IEH + 1
                  ZETAS(I,K,JO2) = BK(K)*ZSDIV(I,KE1,JO2) -
     &                 ZSDIV(I,K,JO2)
               ENDDO
               ETAS(IAH-1:IEH+1,J+1,K)  = ZETAS(IAH-1:IEH+1,K,JO2)
            ENDDO

         ELSE
            DO K =   1 , KE
               DO I = IAH , IEH
                  ZPPHI(I,K,JO2) =  ZX1*(ZTV(I,K,JO5) +    ZTV(I,K,JM5))
     &                              *(ZALOPH(I,K,JO3) - ZALOPH(I,K,JM3))
               ENDDO
            ENDDO
         ENDIF

C        LOKALE HILFSFELDER FUER HORIZONTALDIFFUSION
         DO K = 1 , KE
            DO I = IAH - 1 , IEH + 1
C              GEWICHTETE HDIFF
               ZLAPT(I,K,JO3) =
     &             ((ZTPA(I+1,K,JO5)-ZTPA(I  ,K,JO5))*ZDPU(I  ,K,JO3)-
     &              (ZTPA(I  ,K,JO5)-ZTPA(I-1,K,JO5))*ZDPU(I-1,K,JO3)+
     & (CPHI(J+1,2)*(ZTPA(I  ,K,JN5)-ZTPA(I  ,K,JO5))*ZDPV(I  ,K,JO4)-
     &  CPHI(J  ,2)*(ZTPA(I  ,K,JO5)-ZTPA(I  ,K,JM5))*ZDPV(I  ,K,JM4))/
     &  CPHI(J+1,1)) / ZDP5(I,K,JO5)
               ZLAPQD(I,K,JO3) = 
     &          ((QD(I+1,J+1,K,NA)-QD(I  ,J+1,K,NA))*ZDPU(I  ,K,JO3)-
     &           (QD(I  ,J+1,K,NA)-QD(I-1,J+1,K,NA))*ZDPU(I-1,K,JO3)+
     & (CPHI(J+1,2)*
     &           (QD(I  ,JZN,K,NA)-QD(I  ,J+1,K,NA))*ZDPV(I  ,K,JO4)-
     &  CPHI(J  ,2)*
     &           (QD(I  ,J+1,K,NA)-QD(I  ,J  ,K,NA))*ZDPV(I  ,K,JM4))/ 
     &  CPHI(J+1,1)) / ZDP5(I,K,JO5)
               ZLAPQW(I,K,JO3) = 
     &          ((QW(I+1,J+1,K,NA)-QW(I  ,J+1,K,NA))*ZDPU(I  ,K,JO3)-
     &           (QW(I  ,J+1,K,NA)-QW(I-1,J+1,K,NA))*ZDPU(I-1,K,JO3)+
     & (CPHI(J+1,2)*
     &           (QW(I  ,JZN,K,NA)-QW(I  ,J+1,K,NA))*ZDPV(I  ,K,JO4)-
     &  CPHI(J  ,2)*
     &           (QW(I  ,J+1,K,NA)-QW(I  ,J  ,K,NA))*ZDPV(I  ,K,JM4))/ 
     &  CPHI(J+1,1)) / ZDP5(I,K,JO5)
               ZLAPQI(I,K,JO3) = 
     &          ((QI(I+1,J+1,K,NA)-QI(I  ,J+1,K,NA))*ZDPU(I  ,K,JO3)-
     &           (QI(I  ,J+1,K,NA)-QI(I-1,J+1,K,NA))*ZDPU(I-1,K,JO3)+
     & (CPHI(J+1,2)*
     &           (QI(I  ,JZN,K,NA)-QI(I  ,J+1,K,NA))*ZDPV(I  ,K,JO4)-
     &  CPHI(J  ,2)*
     &           (QI(I  ,J+1,K,NA)-QI(I  ,J  ,K,NA))*ZDPV(I  ,K,JM4))/ 
     &  CPHI(J+1,1)) / ZDP5(I,K,JO5)
C     *********************
               ZLAPU (I,K,JO3) = U (I+1,J+1,K,NA) + U (I-1,J+1,K,NA) -
     &                       2.0*U (I  ,J+1,K,NA) +
     &            ( CPHI(J+1,2)*(U (I  ,JZN,K,NA) - U (I  ,J+1,K,NA)) -
     &              CPHI(J  ,2)*(U (I  ,J+1,K,NA) - U (I  ,J  ,K,NA)) )/
     &              CPHI(J+1,1)
               ZLAPV (I,K,JO3) = V (I+1,J+1,K,NA) + V (I-1,J+1,K,NA) -
     &                       2.0*V (I  ,J+1,K,NA) +
     &            ( CPHI(JZN,1)*(V (I  ,JZN,K,NA) - V (I  ,J+1,K,NA)) -
     &              CPHI(J+1,1)*(V (I  ,J+1,K,NA) - V (I  ,J  ,K,NA)) )/
     &              CPHI(J+1,2)

            ENDDO
         ENDDO

C        AM RECHTEN UND LINKEN RAND H-DIFFUSION 2. ORDNUNG
         DO K = 1 , KE

            IF( (J. EQ. JAH) .AND. (NEIGHBOR(4) .EQ. -1)) THEN
C              SUEDLICHER RAND DES GESAMTGEBIETS
C              H-DIFFUSION 2. ORDNUNG FUER ALLE I
               DO I = IAH , IEH
                  ZTDIFH(I,K) = VVFH(K)*AKS2*ZLAPT (I,K,JM3)
                  ZQDDIH(I,K) = VVFH(K)*AKS2*ZLAPQD(I,K,JM3)
                  ZQWDIH(I,K) = VVFH(K)*AKS2*ZLAPQW(I,K,JM3)
                  ZUDIFH(I,K) = VVFH(K)*AKS2*ZLAPU (I,K,JM3)
                  ZVDIFH(I,K) = VVFH(K)*AKS2*ZLAPV (I,K,JM3)
                  ZQIDIH(I,K) = VVFH(K)*AKS2*ZLAPQI(I,K,JM3)
               ENDDO
            ELSE IF((J.EQ. JEH) .AND. (NEIGHBOR(2) .EQ. -1)) THEN
C              NOERDLICHER RAND DES GESAMTGEBIETS
C              H-DIFFUSION 2. ORDNUNG FUER ALLE I
               DO I = IAH , IEH
                  ZTDIFH(I,K) = VVFH(K)*AKS2*ZLAPT (I,K,JM3)
                  ZQDDIH(I,K) = VVFH(K)*AKS2*ZLAPQD(I,K,JM3)
                  ZQWDIH(I,K) = VVFH(K)*AKS2*ZLAPQW(I,K,JM3)
                  ZUDIFH(I,K) = VVFH(K)*AKS2*ZLAPU (I,K,JM3)
                  ZVDIFH(I,K) = VVFH(K)*AKS2*ZLAPV (I,K,JM3)
                  ZQIDIH(I,K) = VVFH(K)*AKS2*ZLAPQI(I,K,JM3)
               ENDDO
            ELSE
C              BREITENKREIS AUS DER GESAMTGEBIETSMITTE
               DO I = IAH , IEH
                  IF((I .EQ. IAH) .AND. (NEIGHBOR(1) .EQ. -1)) THEN
C                    WESTLICHER RAND DES GESAMTGEBIETS
C                    H-DIFFUSION 2. ORDNUNG FUER DIESES I
                     ZTDIFH(IAH,K) = VVFH(K)*AKS2*ZLAPT (IAH,K,JM3)
                     ZQDDIH(IAH,K) = VVFH(K)*AKS2*ZLAPQD(IAH,K,JM3)
                     ZQWDIH(IAH,K) = VVFH(K)*AKS2*ZLAPQW(IAH,K,JM3)
                     ZUDIFH(IAH,K) = VVFH(K)*AKS2*ZLAPU (IAH,K,JM3)
                     ZVDIFH(IAH,K) = VVFH(K)*AKS2*ZLAPV (IAH,K,JM3)
                     ZQIDIH(IAH,K) = VVFH(K)*AKS2*ZLAPQI(IAH,K,JM3)
                  ELSE IF((I .EQ. IEH) .AND. (NEIGHBOR(3) .EQ. -1)) THEN
                     ZTDIFH(IEH,K) = VVFH(K)*AKS2*ZLAPT (IEH,K,JM3)
                     ZQDDIH(IEH,K) = VVFH(K)*AKS2*ZLAPQD(IEH,K,JM3)
                     ZQWDIH(IEH,K) = VVFH(K)*AKS2*ZLAPQW(IEH,K,JM3)
                     ZUDIFH(IEH,K) = VVFH(K)*AKS2*ZLAPU (IEH,K,JM3)
                     ZVDIFH(IEH,K) = VVFH(K)*AKS2*ZLAPV (IEH,K,JM3)
                     ZQIDIH(IEH,K) = VVFH(K)*AKS2*ZLAPQI(IEH,K,JM3)
                  ELSE

C                    IM INNEREN DES MODELLGEBIETES H-DIFFUSION 4. ORDNUNG
C                    GEWICHTETE HDIFF
                     ZTDIFH(I,K) = -VVFH(K)*AKS4*
     &      ((ZLAPT(I+1,K,JM3)-ZLAPT(I  ,K,JM3))*ZDPU(I  ,K,JM3)-
     &       (ZLAPT(I  ,K,JM3)-ZLAPT(I-1,K,JM3))*ZDPU(I-1,K,JM3)+
     &    ( CPHI(J  ,2)*
     &       (ZLAPT(I  ,K,JO3)-ZLAPT(I  ,K,JM3))*ZDPV(I  ,K,JM4)-
     &      CPHI(J-1,2)*
     &       (ZLAPT(I  ,K,JM3)-ZLAPT(I  ,K,JU3))*ZDPV(I  ,K,JU4))/
     &      CPHI(J  ,1)) / ZDP5(I,K,JM5)
                     ZQDDIH(I,K) = - VVFH(K)*AKS4*
     &      ((ZLAPQD(I+1,K,JM3)-ZLAPQD(I  ,K,JM3))*ZDPU(I  ,K,JM3)-
     &       (ZLAPQD(I  ,K,JM3)-ZLAPQD(I-1,K,JM3))*ZDPU(I-1,K,JM3)+
     &    ( CPHI(J  ,2)*
     &       (ZLAPQD(I  ,K,JO3)-ZLAPQD(I  ,K,JM3))*ZDPV(I  ,K,JM4)-
     &      CPHI(J-1,2)*
     &       (ZLAPQD(I  ,K,JM3)-ZLAPQD(I  ,K,JU3))*ZDPV(I  ,K,JU4))/
     &      CPHI(J  ,1)) / ZDP5(I,K,JM5)
                     ZQWDIH(I,K) = - VVFH(K)*AKS4*
     &      ((ZLAPQW(I+1,K,JM3)-ZLAPQW(I  ,K,JM3))*ZDPU(I  ,K,JM3)-
     &       (ZLAPQW(I  ,K,JM3)-ZLAPQW(I-1,K,JM3))*ZDPU(I-1,K,JM3)+
     &    ( CPHI(J  ,2)*
     &       (ZLAPQW(I  ,K,JO3)-ZLAPQW(I  ,K,JM3))*ZDPV(I  ,K,JM4)-
     &      CPHI(J-1,2)*
     &       (ZLAPQW(I  ,K,JM3)-ZLAPQW(I  ,K,JU3))*ZDPV(I  ,K,JU4))/
     &      CPHI(J  ,1)) / ZDP5(I,K,JM5)
                     ZQIDIH(I,K) = - VVFH(K)*AKS4*
     &      ((ZLAPQI(I+1,K,JM3)-ZLAPQI(I  ,K,JM3))*ZDPU(I  ,K,JM3)-
     &       (ZLAPQI(I  ,K,JM3)-ZLAPQI(I-1,K,JM3))*ZDPU(I-1,K,JM3)+
     &    ( CPHI(J  ,2)*
     &       (ZLAPQI(I  ,K,JO3)-ZLAPQI(I  ,K,JM3))*ZDPV(I  ,K,JM4)-
     &      CPHI(J-1,2)*
     &       (ZLAPQI(I  ,K,JM3)-ZLAPQI(I  ,K,JU3))*ZDPV(I  ,K,JU4))/
     &      CPHI(J  ,1)) / ZDP5(I,K,JM5)
C                    *****************
                     ZUDIFH(I,K) = - VVFH(K)*AKS4*
     &               ( ZLAPU (I+1,K,JM3) + ZLAPU (I-1,K,JM3) -
     &             2.0*ZLAPU (I  ,K,JM3) +
     &  ( CPHI(J  ,2)*(ZLAPU (I  ,K,JO3) - ZLAPU (I  ,K,JM3)) -
     &    CPHI(J-1,2)*(ZLAPU (I  ,K,JM3) - ZLAPU (I  ,K,JU3)) )/
     &    CPHI(J  ,1) )
                     ZVDIFH(I,K) = - VVFH(K)*AKS4*
     &               ( ZLAPV (I+1,K,JM3) + ZLAPV (I-1,K,JM3) -
     &             2.0*ZLAPV (I  ,K,JM3) +
     &  ( CPHI(J+1,1)*(ZLAPV (I  ,K,JO3) - ZLAPV (I  ,K,JM3)) -
     &    CPHI(J  ,1)*(ZLAPV (I  ,K,JM3) - ZLAPV (I  ,K,JU3)) )/
     &    CPHI(J  ,2) )
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

C        EXPLIZITE PROGNOSE OHNE RANDRELAXATIONSTERME
C        --------------------------------------------

         ZX2   =  RERD*ACPHIR(J,1)
         ZFADVX= .5*ACPHIR(J,1)*EDDLAM
         ZFADVY= .5*ACPHIR(J,1)*EDDPHI

C        1. BODENDRUCK + DRUCK
C        ---------------------

         DO I = IAH , IEH
           PSDT(I) = ZSDIV(I,KE1,JM2)
           PS(I,J,NE) = PS(I,J,NA) - PSDT(I)*DT2
           PINT(I,J,1,NE) = AK(1) + BK(1)*PS(I,J,NE)
         ENDDO
  
         DO K = 2 , KE+1
           KM1 = K - 1
           DO I = IAH , IEH
             PINT(I,J,K,NE) = -(BK(KM1) + BK(K)) * DT2 * PSDT(I) *
     &            DWDT(I,J,KM1,NJ) +
     &            PINT(I,J,KM1,NA) + PINT(I,J,K,NA) - 
     &            PINT(I,J,KM1,NE)
           ENDDO
         ENDDO

C        2. H - QDW - PROGNOSE
C        ---------------------

C        HORIZONTALADVEKTION UND QUELLTERME
         DO K = 1 , KE
            DO I = IAH , IEH
               ZEDDPQ(I,K) = 1./ZDP(I,K,JM3)
               ZGVCO       = ZGV(I,K,JM3)*CPHI(J  ,2)
               ZGVCU       = ZGV(I,K,JU3)*CPHI(J-1,2)
               ZT1         = ZGU(I  ,K,JM2)*
     &              ( T(I+1,J,K,NJ) - T(I  ,J,K,NJ) )
               ZT2         = ZGU(I-1,K,JM2)*
     &              ( T(I  ,J,K,NJ) - T(I-1,J,K,NJ) )
               ZT3         = ZGVCO*( T(I,J+1,K,NJ) - T(I,J,K,NJ) )
               ZT4         = ZGVCU*( T(I,J,K,NJ) - T(I,J-1,K,NJ) )
               ZTADV(I,K)  = -( ZFADVX*(ZT1+ZT2) + ZFADVY*(ZT3+ZT4) )
     &              *ZEDDPQ(I,K)
               ZQD1        = ZGU(I  ,K,JM2)*
     &              (QD(I+1,J,K,NJ) - QD(I  ,J,K,NJ))
               ZQD2        = ZGU(I-1,K,JM2)*
     &              (QD(I  ,J,K,NJ) - QD(I-1,J,K,NJ) )
               ZQD3        = ZGVCO*( QD(I,J+1,K,NJ)-QD(I,J,K,NJ))
               ZQD4        = ZGVCU*( QD(I,J,K,NJ)-QD(I,J-1,K,NJ))
               ZQDADV(I,K) = -( ZFADVX*(ZQD1+ZQD2) + ZFADVY*(ZQD3+ZQD4))
     &              *ZEDDPQ(I,K)

C              ALPHA*OMEGA
               ZA1 =  ZGU(I  ,K,JM2)*ZPLAM(I  ,K) +
     &                ZGU(I-1,K,JM2)*ZPLAM(I-1,K)
               ZA2 = ( ZGVCO*ZPPHI(I,K,JO2) + ZGVCU*ZPPHI(I,K,JM2) )*ZX2
               ZA3 = - R*ZTV(I,K,JM5)*ZEDDPQ(I,K)*
     &         (   (ZALOPN(I,K+1,JM3)-ZALOPN(I,K,JM3))*ZSDIV (I,K,JM2)
     &           + (ZSDIV (I,K+1,JM2)-ZSDIV (I,K,JM2))*ZBETAK(I,K,JM3) )

               ZALPOM(I,K) = ( .5*( ZA1 + ZA2 )*ZEDDPQ(I,K) + ZA3 )

               AGB(I,K,1) = ED2DT
               AGD(I,K,1) = ED2DT*T(I,J,K,NA) + ZTADV(I,K) +
     &              ZTDIFH(I,K) + ZALPOM(I,K) *WCPR+ ZTTK(I,K) +
     &              ZTTS(I,K) + ZTHDTA(I,K) + ZSODTA(I,K)
               AGD(I,K,2) = ED2DT*QD(I,J,K,NA) + ZQDADV(I,K) +
     &              ZQDDIH(I,K) + ZQDTK(I,K) + ZQDTS(I,K)
               AGD(I,K,3) = ED2DT*QW(I,J,K,NA) + ZQWDIH(I,K)
     &              + ZQWTS(I,K)
               AGD(I,K,4) = ED2DT*QI(I,J,K,NA) + ZQIDIH(I,K)
     &              + ZQITS(I,K)
            ENDDO
         ENDDO


C        VERTIKALADVEKTION UND -DIFFUSION
         DO I = IAH , IEH
            ZAGCM      = .5*ZETAS(I,2,JM2)*ZEDDPQ(I,1)
            ZTMKVHM    = ZTKVZ*ZTMKVH(I  ,2,JM3) + ZTKVL*
     &           ( ZTMKVH(I+1,2,JM3) + ZTMKVH(I,2,JO3) +
     &             ZTMKVH(I-1,2,JM3) + ZTMKVH(I,2,JU3) )
            ZAGCT      = - ZTMKVHM*ZEDDPQ(I,1)
            AGC(I,1,1) = ZAGCM*ZA1A + ZAGCT*A1T(2)
            AGB(I,1,1) = AGB(I,1,1) - AGC(I,1,1)
            AGD(I,1,1) = AGD(I,1,1) - ( ZA2A*ZAGCM + A2T(2)*ZAGCT )*
     &           ( T(I,J,2,NA) - T(I,J,1,NA) )
     &           -WCPR* ZAGCT*( ZFIHF(I,2,JM5) - ZFIHF(I,1,JM5) )
            AGD(I,1,2) = AGD(I,1,2) - ( ZA2A*ZAGCM + A2T(2)*ZAGCT )*
     &           ( QD(I,J,2,NA) - QD(I,J,1,NA) )
            AGD(I,1,3) = AGD(I,1,3) - ( ZA2A*ZAGCM + A2T(2)*ZAGCT )*
     &           ( QW(I,J,2,NA) - QW(I,J,1,NA) )
            AGD(I,1,4) = AGD(I,1,4) - ( ZA2A*ZAGCM + A2T(2)*ZAGCT )*
     &           ( QI(I,J,2,NA) - QI(I,J,1,NA) )
            AGB(I,1,2) = AGB(I,1,1)
            AGC(I,1,2) = AGC(I,1,1)
            AGB(I,1,3) = AGB(I,1,1)
            AGC(I,1,3) = AGC(I,1,1)
            AGB(I,1,4) = AGB(I,1,1)
            AGC(I,1,4) = AGC(I,1,1)
         ENDDO

         DO K = 2 , KE-1
            DO I  = IAH , IEH
               ZAGAM      = -.5*ZETAS(I,K  ,JM2)*ZEDDPQ(I,K)
               ZAGCM      =  .5*ZETAS(I,K+1,JM2)*ZEDDPQ(I,K)
               ZTMKVHM    = ZTKVZ*ZTMKVH(I  ,K,JM3) + ZTKVL*
     &              ( ZTMKVH(I+1,K,JM3) + ZTMKVH(I  ,K,JO3) +
     &                ZTMKVH(I-1,K,JM3) + ZTMKVH(I  ,K,JU3) )
               ZAGAT      = - ZTMKVHM*ZEDDPQ(I,K)
               ZTMKVHM    = ZTKVZ*ZTMKVH(I  ,K+1,JM3) + ZTKVL*
     &              ( ZTMKVH(I+1,K+1,JM3) + ZTMKVH(I  ,K+1,JO3) +
     &                ZTMKVH(I-1,K+1,JM3) + ZTMKVH(I  ,K+1,JU3) )
               ZAGCT      = - ZTMKVHM*ZEDDPQ(I,K)
               AGA(I,K,1) = ZAGAM*ZA1A + ZAGAT*A1T(K)
               AGC(I,K,1) = ZAGCM*ZA1A + ZAGCT*A1T(K+1)
               AGB(I,K,1) = AGB(I,K,1) - AGA(I,K,1) - AGC(I,K,1)
               AGD(I,K,1) = AGD(I,K,1) - ( ZA2A*ZAGAM + A2T(K)*ZAGAT )*
     &              ( T(I,J,K-1,NA) - T(I,J,K,NA) )
     &              - ( ZA2A*ZAGCM + A2T(K+1)*ZAGCT )*
     &              ( T(I,J,K+1,NA) - T(I,J,K,NA) )
     &              -WCPR* ZAGAT*( ZFIHF(I,K-1,JM5) - ZFIHF(I,K,JM5))
     &              -WCPR* ZAGCT*( ZFIHF(I,K+1,JM5) - ZFIHF(I,K,JM5))
               AGD(I,K,2) = AGD(I,K,2) - ( ZA2A*ZAGAM + A2T(K)*ZAGAT )*
     &              (QD(I,J,K-1,NA) - QD(I,J,K,NA) )
     &              - ( ZA2A*ZAGCM + A2T(K+1)*ZAGCT )*
     &              (QD(I,J,K+1,NA) - QD(I,J,K,NA) )
               AGD(I,K,3) = AGD(I,K,3) - ( ZA2A*ZAGAM + A2T(K)*ZAGAT )*
     &              (QW(I,J,K-1,NA) - QW(I,J,K,NA) )
     &              - ( ZA2A*ZAGCM + A2T(K+1)*ZAGCT )*
     &              (QW(I,J,K+1,NA) - QW(I,J,K,NA) )
               AGD(I,K,4) = AGD(I,K,4) - ( ZA2A*ZAGAM + A2T(K)*ZAGAT )*
     &              (QI(I,J,K-1,NA) - QI(I,J,K,NA) )
     &              - ( ZA2A*ZAGCM + A2T(K+1)*ZAGCT )*
     &              (QI(I,J,K+1,NA) - QI(I,J,K,NA) )
               AGA(I,K,2) = AGA(I,K,1)
               AGB(I,K,2) = AGB(I,K,1)
               AGC(I,K,2) = AGC(I,K,1)
               AGA(I,K,3) = AGA(I,K,1)
               AGB(I,K,3) = AGB(I,K,1)
               AGC(I,K,3) = AGC(I,K,1)
               AGA(I,K,4) = AGA(I,K,1)
               AGB(I,K,4) = AGB(I,K,1)
               AGC(I,K,4) = AGC(I,K,1)
            ENDDO
         ENDDO
C
         DO I = IAH , IEH
            ZAGAM       = -.5*ZETAS(I,KE,JM2)*ZEDDPQ(I,KE)
            ZTMKVHM     = ZTKVZ*ZTMKVH(I  ,KE,JM3) + ZTKVL*
     &           ( ZTMKVH(I+1,KE,JM3) + ZTMKVH(I  ,KE,JO3) +
     &             ZTMKVH(I-1,KE,JM3) + ZTMKVH(I  ,KE,JU3) )
            ZAGAT       = -  ZTMKVHM *ZEDDPQ(I,KE)
            ZAGCT       = - TMCH(I,J)*ZEDDPQ(I,KE)
            AGA(I,KE,1) = ZAGAM*ZA1A  + ZAGAT*A1T(KE)
            AGB(I,KE,1) = AGB(I,KE,1) - AGA(I,KE,1) - ZAGCT*A1T(KE1)
            AGD(I,KE,1) = AGD(I,KE,1) - ( ZA2A*ZAGAM + A2T(KE)*ZAGAT )*
     &           ( T(I,J,KE-1,NA) - T(I,J,KE,NA) )
     &           - A2T(KE1)*ZAGCT*( TG(I,J,NA) - T(I,J,KE,NA) )
     &           - A1T(KE1)*ZAGCT*  TG(I,J,NA)
     &           -WCPR* ZAGAT*( ZFIHF(I,KE-1,JM5) - ZFIHF(I,KE,JM5) )
     &           -WCPR* ZAGCT*( FIB  (I,J)        - ZFIHF(I,KE,JM5) )
            AGD(I,KE,2) = AGD(I,KE,2) - ( ZA2A*ZAGAM + A2T(KE)*ZAGAT )*
     &           ( QD(I,J,KE-1,NA) - QD(I,J,KE,NA) )
     &           - A2T(KE1)*ZAGCT*(QDB(I,J,NA)  - QD(I,J,KE,NA) )
     &           - A1T(KE1)*ZAGCT* QDB(I,J,NA)
            AGD(I,KE,3) = AGD(I,KE,3) - ( ZA2A*ZAGAM + A2T(KE)*ZAGAT )*
     &           ( QW(I,J,KE-1,NA) - QW(I,J,KE,NA) )
     &             + A2T(KE1)*ZAGCT* QW(I,J,KE,NA)
            AGD(I,KE,4) = AGD(I,KE,4) - ( ZA2A*ZAGAM + A2T(KE)*ZAGAT )*
     &           ( QI(I,J,KE-1,NA) - QI(I,J,KE,NA) )
     &             + A2T(KE1)*ZAGCT* QI(I,J,KE,NA)
            AGA(I,KE,2) = AGA(I,KE,1)
            AGB(I,KE,2) = AGB(I,KE,1)
            AGA(I,KE,3) = AGA(I,KE,1)
            AGB(I,KE,3) = AGB(I,KE,1)
            AGA(I,KE,4) = AGA(I,KE,1)
            AGB(I,KE,4) = AGB(I,KE,1)
         ENDDO

C        GAUSS - ELIMINATION UND H-QDW-ZERLEGUNG

         CALL GAUSS ( IAH , IEH , INDEX1, AGA,AGB,AGC,AGD,AGE)
         CALL GAUSS ( IAH , IEH , INDEX2, AGA,AGB,AGC,AGD,AGE)
         CALL GAUSS ( IAH , IEH , INDEX3, AGA,AGB,AGC,AGD,AGE)
         CALL GAUSS ( IAH , IEH , INDEX4, AGA,AGB,AGC,AGD,AGE)

         DO K = 1 , KE
            DO I = IAH , IEH
               HE(I,K)= AGE(I,K,1)*WCP + WLK * AGE(I,K,2)*ZTRC
               QDWE(I,K)= (AGE(I,K,2) + AGE(I,K,3) + AGE(I,K,4))*ZTRC
               QDLE(I,K)= (AGE(I,K,2) + AGE(I,K,3))*ZTRC
               QDIE(I,K)= AGE(I,K,4)*ZTRC
               TSTART(I,K) = T(I,J,K,NJ)
               GQDSTA(I,K) = ZGQD(I,K,JM3)
               PHFSTA(I,K) = ZPHF(I,K,JM5)
               PHFE  (I,K) = AKH(K) + BKH(K)*PS(I,J,NE)
            ENDDO
         ENDDO
C
C        BERECHNUNG DES FLUESSES LATENTER WAERME AM BODEN;
C        DIESER FLUESS WIRD AUFSUMMIERT
C
         DO I  = IAH , IEH
            IF (INFRL(I,J).GT.0) THEN
               BFLQDSL(I,J) = BFLQDSL(I,J) - TMCHL(I,J)*EDG*(WLK*
     &              (A2T(KE1)*(QDBL(I,J,NA) - (QD(I,J,KE,NA)+
     &              (QW(I,J,KE,NA))) ) +
     &              A1T(KE1)*(QDBL(I,J,NE) - QDLE (I,KE  )))-WLS*(
     &              A2T(KE1)*QI(I,J,KE,NA)+A1T(KE1)*QDIE(I,KE)))*DTDEH
            ENDIF
            IF (INFRW(I,J).GT.0) THEN
               BFLQDSW(I,J) = BFLQDSW(I,J) - TMCHW(I,J)*EDG*(WLK*
     &              (A2T(KE1)*(QDBW(I,J,NA) - (QD(I,J,KE,NA)+
     &              (QW(I,J,KE,NA))) ) +
     &              A1T(KE1)*(QDBW(I,J,NE) - QDLE (I,KE  )))-WLS*(
     &              A2T(KE1)*QI(I,J,KE,NA)+A1T(KE1)*QDIE(I,KE)))*DTDEH
            ENDIF
            IF (INFRI(I,J).GT.0) THEN
               BFLQDSI(I,J) = BFLQDSI(I,J) - TMCHI(I,J)*EDG*(WLK*
     &              (A2T(KE1)*(QDBI(I,J,NA) - (QD(I,J,KE,NA)+
     &              (QW(I,J,KE,NA))) ) +
     &              A1T(KE1)*(QDBI(I,J,NE) - QDLE (I,KE  )))-WLS*(
     &              A2T(KE1)*QI(I,J,KE,NA)+A1T(KE1)*QDIE(I,KE)))*DTDEH
            ENDIF
            BFLQDS(I,J) = (FLOAT(INFRL(I,J))*BFLQDSL(I,J)
     &                  +  FLOAT(INFRW(I,J))*BFLQDSW(I,J)
     &                  +  FLOAT(INFRI(I,J))*BFLQDSI(I,J))*EDFAKINF
         ENDDO
C
C        MASSENFLUSS-KORREKTURSCHEMA  (FUER QDW<0)
         IF(LMASSF) THEN
            DO I = IAH , IEH
               ZQKOR(I) = 0.
            ENDDO
            DO K = 1 , KE
               DO I = IAH , IEH
                  ZQKOR(I)   = QDWE(I,K) + ZQKOR(I)*ZEDDPQ(I,K)
                  IF(ZQKOR(I).LT.0.) THEN
                     QDWE(I,K) = 0.
                     ZQKOR(I)  = ZQKOR(I)*ZDP(I,K,JM3)
                  ELSE
                     QDWE(I,K) = ZQKOR(I)
                     ZQKOR(I)  = 0.
                  END IF
               ENDDO
            ENDDO
         END IF

C        BEI ADIABATISCHER INITIALISIERUNG: T MIT UNVERAENDERTEM QD BERECH-
C        NEN; SONST IST DIE SKALIGE KONDENSATIONSRATE UNGLEICH NULL.
         IF (LAISTEP) THEN
            DO K = 1 , KE
               DO I = IAH , IEH
                  QD(I,J,K,NE) = QD(I,J,K,NA)
                  QW(I,J,K,NE) = QW(I,J,K,NA)
                  T (I,J,K,NE) = WCPR*HE(I,K)
                  QI(I,J,K,NE) = QI(I,J,K,NA)
               ENDDO
            ENDDO

         ELSE
            DO K = 1 , KE
               DO I = IAH , IEH
                  QD(I,J,K,NE) = AGE(I,K,2)
                  QW(I,J,K,NE) = AGE(I,K,3)
                  T (I,J,K,NE) = AGE(I,K,1)
                  QI(I,J,K,NE) = AGE(I,K,4)
               ENDDO
            ENDDO
            DO I = IAH , IEH
               ZQKOR(I) = 0.
            ENDDO
            DO K = 1 , KE
               DO I = IAH , IEH
                  ZQW1 = QW(I,J,K,NE)
                  ZQW2 = QI(I,J,K,NE)
                  IF (ZQW1.LT.0.) THEN
                     T (I,J,K,NE) = T (I,J,K,NE) - WLK*QW(I,J,K,NE)*WCPR
                     QD(I,J,K,NE) = QD(I,J,K,NE) + QW(I,J,K,NE)
                     QW(I,J,K,NE) = 0.
                  ENDIF
                  IF (ZQW2.LT.0.) THEN
                     T (I,J,K,NE) = T (I,J,K,NE) - WLS*QI(I,J,K,NE)*WCPR
                     QD(I,J,K,NE) = QD(I,J,K,NE) + QI(I,J,K,NE)
                     QI(I,J,K,NE) = 0.
                  ENDIF
                  ZQKOR(I)   = QD(I,J,K,NE) + ZQKOR(I)*ZEDDPQ(I,K)
                  IF(ZQKOR(I).LT.0.) THEN
                     QD(I,J,K,NE) = 0.
                     QW(I,J,K,NE) = 0.
                     QI(I,J,K,NE) = 0.
                     ZQKOR(I)  = ZQKOR(I)*ZDP(I,K,JM3)
                  ELSE
                     QD(I,J,K,NE) = ZQKOR(I)
                     ZQKOR(I)  = 0.
                  END IF
               ENDDO
            ENDDO
C
C     BERECHNUNG DES FLUESSES SENSIBLER WAERME AM BODEN;
C
            DO I = IAH , IEH
               IF (INFRL(I,J).GT.0) THEN
                  BFLHSL(I,J)= BFLHSL(I,J)- TMCHL(I,J)*EDG*WCP*
     &                 (A2T(KE1)*(TGL(I,J,NA) - T    (I,J,KE,NA)) +
     &                 A1T(KE1)*(TGL(I,J,NE) - T    (I,J,KE,NE)) +
     &                 WCPR*    (FIB(I,J)    - ZFIHF(I,  KE,JM5)))*DTDEH
               ENDIF
               IF (INFRW(I,J).GT.0) THEN
                  BFLHSW(I,J)= BFLHSW(I,J)- TMCHW(I,J)*EDG*WCP*
     &                 (A2T(KE1)*(TGW(I,J,NA) - T    (I,J,KE,NA)) +
     &                 A1T(KE1)*(TGW(I,J,NE) - T    (I,J,KE,NE)) +
     &                 WCPR*    (FIB(I,J)    - ZFIHF(I,  KE,JM5)))*DTDEH
               ENDIF
               IF (INFRI(I,J).GT.0) THEN
                  BFLHSI(I,J)= BFLHSI(I,J)- TMCHI(I,J)*EDG*WCP*
     &                 (A2T(KE1)*(TGI(I,J,NA) - T    (I,J,KE,NA)) +
     &                 A1T(KE1)*(TGI(I,J,NE) - T    (I,J,KE,NE)) +
     &                 WCPR*    (FIB(I,J)    - ZFIHF(I,  KE,JM5)))*DTDEH
               ENDIF
               BFLHS(I,J) = (FLOAT(INFRL(I,J))*BFLHSL(I,J)
     &              +  FLOAT(INFRW(I,J))*BFLHSW(I,J)
     &              +  FLOAT(INFRI(I,J))*BFLHSI(I,J))*EDFAKINF
            ENDDO
C
C         OMEGA-WERTE BESTIMMTER MODELLFLAECHEN SUMMIEREN (PROGCHK)
            DO I = IAH , IEH
               ZOM850M   = ZOM850M +
     &    ABS(ZALPOM(I,KFL850)*ZPHF(I,KFL850,JM5)/(R*ZTV(I,KFL850,JM5)))
               ZOM500M   = ZOM500M +
     &    ABS(ZALPOM(I,KFL500)*ZPHF(I,KFL500,JM5)/(R*ZTV(I,KFL500,JM5)))
               ZOM300M   = ZOM300M +
     &    ABS(ZALPOM(I,KFL300)*ZPHF(I,KFL300,JM5)/(R*ZTV(I,KFL300,JM5)))
            ENDDO
         ENDIF ! LAISTEP


C        3. U - PROGNOSE
C        ---------------

         ZX1   = ACPHIR(J,1)*EDDLAM
         ZFVZO = 0.125*CPHI(J  ,2)/CPHI(J,1)
         ZFVZU = 0.125*CPHI(J-1,2)/CPHI(J,1)

         DO K = 1 , KE
            DO I = IAU , IEU
               ZGRAD(I,K) = -ZPLAM(I,K) - ZX1*
     &             ( (ZFIH(I+1,K,JM3) - ZFIH(I,K,JM3)) *
     &             0.5 * (DWDT(I,J,K,NJ)+DWDT(I+1,J,K,NJ))
     &             + ZEKIN(I+1,K,JM2) - ZEKIN(I,K,JM2) )
               ZVZEQ(I,K) = ( ZZETA(I,K,JO2) + ZZETA(I,K,JM2) )*
     &              (  ZFVZO*( ZGV(I+1,K,JM3) + ZGV(I,K,JM3) )
     &               + ZFVZU*( ZGV(I+1,K,JU3) + ZGV(I,K,JU3) ) )
               ZEDDPQ(I,K) =  2./( ZDP(I+1,K,JM3) + ZDP(I,K,JM3) )
               AGB(I,K,1) = ED2DT
               AGD(I,K,1) =(ED2DT-ZEPSRAY)*U(I,J,K,NA) + ZGRAD(I,K) +
     &              ZVZEQ(I,K) + ZUDIFH(I,K) +
     &              0.5*( ZUTK(I+1,K,JM2) + ZUTK(I,K,JM2) )
            ENDDO
         ENDDO

C        VERTIKALDIFFUSION UND VERTIKALADVEKTION
         DO I = IAU , IEU
            ZAGCM = .25*( ZETAS (I,2,JM2) + ZETAS(I+1,2,JM2) )*
     &           ZEDDPQ(I,1)
            ZAGCT = -0.5*(ZTMKVM(I,2,JM3)+ ZTMKVM(I+1,2,JM3) )*
     &           ZEDDPQ(I,1)
            AGC(I,1,1) = ZAGCM*ZA1A + ZAGCT*A1T(2)
            AGB(I,1,1) = AGB(I,1,1) - AGC(I,1,1)
            AGD(I,1,1) = AGD(I,1,1) - ( ZA2A*ZAGCM + A2T(2)*ZAGCT )*
     &           ( U(I,J,2,NA) - U(I,J,1,NA ) )
         ENDDO

         DO K = 2 , KE-1
            DO I = IAU , IEU
               ZAGAM = - 0.25*(ZETAS(I,K  ,JM2)+ZETAS(I+1,K  ,JM2))*
     &              ZEDDPQ(I,K)
               ZAGCM =   0.25*(ZETAS(I,K+1,JM2)+ZETAS(I+1,K+1,JM2))*
     &              ZEDDPQ(I,K)
               ZAGAT = - 0.5 *(ZTMKVM(I,K  ,JM3)+ZTMKVM(I+1,K  ,JM3))*
     &              ZEDDPQ(I,K)
               ZAGCT = - 0.5 *(ZTMKVM(I,K+1,JM3)+ZTMKVM(I+1,K+1,JM3))*
     &              ZEDDPQ(I,K)
               AGA(I,K,1) = ZAGAM*ZA1A + ZAGAT*A1T(K)
               AGC(I,K,1) = ZAGCM*ZA1A + ZAGCT*A1T(K+1)
               AGB(I,K,1) = AGB(I,K,1) - AGA(I,K,1) - AGC(I,K,1)
               AGD(I,K,1) = AGD(I,K,1) - ( ZA2A*ZAGAM + A2T(K)*ZAGAT )*
     &              ( U(I,J,K-1,NA) - U(I,J,K,NA) )
     &              - ( ZA2A*ZAGCM + A2T(K+1)*ZAGCT )*
     &              ( U(I,J,K+1,NA) - U(I,J,K,NA) )
            ENDDO
         ENDDO

         DO I = IAU , IEU
            ZAGAM = - 0.25*(ZETAS  (I,KE,JM2)+ZETAS (I+1,KE,JM2))*
     &           ZEDDPQ(I,KE)
            ZAGAT = - 0.5 *( ZTMKVM(I,KE,JM3)+ZTMKVM(I+1,KE,JM3))*
     &           ZEDDPQ(I,KE)
            ZAGCT = - 0.5 *( TMCM (I,J   ) + TMCM (I+1,J  ) )    *
     &           ZEDDPQ(I,KE)
            AGA(I,KE,1) = ZAGAM*ZA1A  + ZAGAT*A1T(KE)
            AGB(I,KE,1) = AGB(I,KE,1) - AGA(I,KE,1) - ZAGCT*A1T(KE1)
            AGD(I,KE,1) = AGD(I,KE,1) - ( ZA2A*ZAGAM + A2T(KE)*ZAGAT )*
     &           ( U(I,J,KE-1,NA) - U(I,J,KE,NA) )
     &             + A2T(KE1)*ZAGCT*U(I,J,KE,NA)
         ENDDO

C        GAUSS-ELIMINATION UND U-ENDE-WERTE AUF U(I,J,K,NE) ABLEGEN
         CALL GAUSS ( IAU , IEU , INDEX1, AGA,AGB,AGC,AGD,AGE)

         DO K = 1 , KE
            DO I = IAU , IEU
               U(I,J,K,NE) = AGE(I,K,1)
            ENDDO
         ENDDO

C        BERECHNUNG DES U-IMPULSFLUSSES AM BODEN;
C        DIESER FLUESS WIRD AUFSUMMIERT
         DO I = IAU , IEU
            BFLUS(I,J) = BFLUS(I,J) + 0.5*(TMCM(I,J) + TMCM(I+1,J))*EDG*
     &           (A2T(KE1)*U(I,J,KE,NA) +
     &            A1T(KE1)*U(I,J,KE,NE))*DTDEH
         ENDDO

C        4. V - PROGNOSE
C        ---------------

         IF(J.LE.JEV) THEN

            DO K = 1 , KE
               DO I = IAV , IEV
                  ZGRAD(I,K) =
     &                 - ZPPHI(I,K,JO2) - EDADPHI*
     &                ( (ZFIH(I,K,JO3) - ZFIH(I,K,JM3) ) *
     &                0.5*(DWDT(I,J,K,NJ)+DWDT(I,J+1,K,NJ))
     &                + ZEKIN(I,K,JO2) - ZEKIN(I,K,JM2) )
                  ZUZEQ(I,K) =
     &                 - 0.125*( ZZETA(I-1,K,JO2) + ZZETA(I,K,JO2) )*
     &                          (  ZGU(I-1,K,JO2) + ZGU  (I,K,JO2)
     &                           + ZGU(I-1,K,JM2) + ZGU  (I,K,JM2) )
                  ZEDDPQ(I,K) =  2./( ZDP(I,K,JO3) + ZDP(I,K,JM3) )
                  AGB(I,K,2) = ED2DT
                  AGD(I,K,2) =(ED2DT-ZEPSRAY)*V(I,J,K,NA) +
     &                 ZGRAD(I,K) + ZUZEQ(I,K) + ZVDIFH(I,K) +
     &                 0.5*( ZVTK(I,K,JO2) + ZVTK(I,K,JM2) )
               ENDDO
            ENDDO

C           VERTIKALDIFFUSION UND VERTIKALADVEKTION
            DO I = IAV , IEV
               ZAGCM =  0.25*( ZETAS(I,2,JM2) +  ZETAS(I,2,JO2) )*
     &              ZEDDPQ(I,1)
               ZAGCT = - 0.5*(ZTMKVM(I,2,JM3) + ZTMKVM(I,2,JO3) )*
     &              ZEDDPQ(I,1)
               AGC(I,1,2) = ZAGCM*ZA1A + ZAGCT*A1T(2)
               AGB(I,1,2) = AGB(I,1,2) - AGC(I,1,2)
               AGD(I,1,2) = AGD(I,1,2) - ( ZA2A*ZAGCM + A2T(2)*ZAGCT )*
     &              ( V(I,J,2,NA) - V(I,J,1,NA ) )
            ENDDO

            DO K = 2 , KE-1
               DO I = IAV , IEV
                  ZAGAM = - 0.25*(ZETAS(I,K  ,JM2)+ZETAS(I,K  ,JO2))*
     &                 ZEDDPQ(I,K)
                  ZAGCM =   0.25*(ZETAS(I,K+1,JM2)+ZETAS(I,K+1,JO2))*
     &                 ZEDDPQ(I,K)
                  ZAGAT = 
     &              - 0.5 *( ZTMKVM(I,K  ,JM3) + ZTMKVM(I  ,K  ,JO3) ) *
     &              ZEDDPQ(I,K)
                  ZAGCT =
     &              - 0.5 *( ZTMKVM(I,K+1,JM3) + ZTMKVM(I  ,K+1,JO3) ) *
     &              ZEDDPQ(I,K)
                  AGA(I,K,2) = ZAGAM*ZA1A + ZAGAT*A1T(K)
                  AGC(I,K,2) = ZAGCM*ZA1A + ZAGCT*A1T(K+1)
                  AGB(I,K,2) = AGB(I,K,2) - AGA(I,K,2) - AGC(I,K,2)
                  AGD(I,K,2) = AGD(I,K,2) -
     &                 ( ZA2A*ZAGAM + A2T(K)*ZAGAT )*
     &                 ( V(I,J,K-1,NA) - V(I,J,K,NA) )-
     &                 ( ZA2A*ZAGCM + A2T(K+1)*ZAGCT )*
     &                 ( V(I,J,K+1,NA) - V(I,J,K,NA) )
               ENDDO
            ENDDO

            DO I = IAV , IEV
               ZAGAM = - 0.25*(ZETAS (I,KE,JM2)+ZETAS (I,KE,JO2))*
     &              ZEDDPQ(I,KE)
               ZAGAT = - 0.5 *(ZTMKVM(I,KE,JM3)+ZTMKVM(I,KE,JO3))*
     &              ZEDDPQ(I,KE)
               ZAGCT = - 0.5 *( TMCM (I,J   ) + TMCM (I,J+1 ) )  *
     &              ZEDDPQ(I,KE)
               AGA(I,KE,2) = ZAGAM*ZA1A  + ZAGAT*A1T(KE)
               AGB(I,KE,2) = AGB(I,KE,2) - AGA(I,KE,2) - ZAGCT*A1T(KE1)
               AGD(I,KE,2) = AGD(I,KE,2) -
     &              ( ZA2A*ZAGAM + A2T(KE)*ZAGAT )*
     &              ( V(I,J,KE-1,NA) - V(I,J,KE,NA) )
     &                + A2T(KE1)*ZAGCT*V(I,J,KE,NA)
            ENDDO

C           GAUSS ELIMINATION UND V-ENDE-WERTE AUF V(I,J,K,NE) ABLEGEN
            CALL GAUSS ( IAV , IEV , INDEX2, AGA,AGB,AGC,AGD,AGE)

            DO K = 1 , KE
               DO I = IAV , IEV
                  V(I,J,K,NE) = AGE(I,K,2)
               ENDDO
            ENDDO

C           BERECHNUNG DES V-IMPULSFLUSSES AM BODEN;
C           DIESER FLUESS WIRD AUFSUMMIERT
            DO I = IAV , IEV
               BFLVS(I,J) = BFLVS(I,J) +
     &              0.5*(TMCM(I,J) + TMCM(I,J+1))*EDG*
     &              (A2T(KE1)*V(I,J,KE,NA) +
     &               A1T(KE1)*V(I,J,KE,NE))*DTDEH
            ENDDO

         ENDIF ! IF(J.LE.JEV)

C        ZYKLISCHE VERTAUSCHUNG DER SCHEIBENINDIZES

         JM2 = 3 - JM2
         JO2 = 3 - JO2

         JSP = JU3
         JU3 = JM3
         JM3 = JO3
         JO3 = JSP

         JSP = JS4
         JS4 = JU4
         JU4 = JM4
         JM4 = JO4
         JO4 = JSP

         JSP = JS5
         JS5 = JU5
         JU5 = JM5
         JM5 = JO5
         JO5 = JN5
         JN5 = JSP

      ENDDO ! J = JAH, JEH

C     DIESER FLUESS WIRD AUFSUMMIERT
C     OMEGA-MITTELWERTE BESTIMMTER MODELLFLAECHEN BILDEN
      ZDELTF  = FLOAT( (IEHGG-IAHGG+1)*(JEHGG-JAHGG+1) )
      ZOM850M = ZOM850M/ZDELTF
      ZOM500M = ZOM500M/ZDELTF
      ZOM300M = ZOM300M/ZDELTF

C-----------------------------------------------------------------------
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
      END SUBROUTINE PROGEXP
