C
C     SUBROUTINE ECTIED
C
C**** ECTIED   -   UP: BERECHNUNG VON OMEGA UND FEUCHTE-TENDENZEN
C**   AUFRUF   :   CALL ECTIED  ( JATKON, JETKON,
C**                               GMVERT, GMQTEN, GMQHFL )
C**   ENTRIES  :      ---
C**   ZWECK    :   BERECHNUNG DER FEUCHTEKONVERGENZ UND BEREITSTELLUNG
C**                DER FUER DIE TIEDKE-KONVEKTION BENOETIGTEN FELDER
C**   VERSIONS-
C**   DATUM    :   01.06.1995
C**
C**   EXTERNALS:   CUCALL, GAUSS, COPYRE
C**
C**   EINGABE-
C**   PARAMETER:   JATPROG: ANFANGS-J-INDEX FUER KONVEKTIONSRECHNUNG;
C**                JETPROG: END    -J-INDEX FUER KONVEKTIONSRECHNUNG;
C**                         DIE FESTLEGUNG ERFOLGT IN UP *ECKONT*.
C**   AUSGABE-
C**   PARAMETER:   GMVERT, GMQDQT, GMQHFL
C**
C**   COMMON-
C**   BLOECKE  :   PARAM  , ORG, COMDYN, COMPYH, PHYKON, HIGKON, PARKON,
C**                PROGCHK, COMDIA, COMPCST, COMPMLF, COMPSLF, COMPGP1,
C**                COMPGP3, COMSDS
C**
C**   METHODE  :   ANALOG ZU UP *PROGEXP* WIRD DIE FEUCHTETENDENZ
C**                BERECHNET ALS SUMME AUS DREIDIMENSIONALER LARGE-SCALE
C**                FEUCHTEKONVERGENZ UND DER DIVERGENZ DES VERTIKALEN
C**                TURBULENTEN FLUSSES.
C**   FEHLERBE-
C**   HANDLUNG :      ---
C**   VERFASSER:   D. MAJEWSKI
C
      SUBROUTINE ECTIED(
     &    JATKON, JETKON, GMVERT, GMQTEN,
     &    AK , BK , AKH   ,
     &    BKH   , DAK   , DBK   , ALPHABOUND, A1T, A2T, ACPHIR,
     &    CPHI  , PS    , QDB   , U     , V  , T  , QD    ,
     &    TMKVMH, TMCH)
C
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
      INCLUDE "org.h"
      INCLUDE "comdyn.h"
      INCLUDE "phykon.h"
      INCLUDE "higkon.h"
C
      INCLUDE "comphy.h"
      INCLUDE "parkon.h"
C
C     Local Variables
C
      INTEGER, INTENT(IN) :: JATKON, JETKON
C
C-----------------------------------------------------------------------
C     ERFORDERLICHE EM-FELDER AUS DEM LANGZEITSPEICHER DIMENSIONIEREN
C     ---------------------------------------------------------------
C     VERTIKAL-KOORDINATEN-PARAMETER
      REAL, INTENT(IN) :: AK(KE1), BK(KE1), AKH(KE), 
     &                    BKH(KE), DAK(KE), DBK(KE)
C
C     EXTERNE PARAMETER
      REAL, INTENT(IN) :: ALPHABOUND(IE,JE,3)
C
C     VERTIKAL VARIIERENDER IMPLIZITHEITSGRAD DER V-DIFFUSION
      REAL, INTENT(IN) :: A1T(KE1), A2T(KE1)
C
C     EXTERNE PARAMETER
      REAL, INTENT(IN) :: ACPHIR(JE,2), CPHI(JE,2)
C
C     PROGNOSTISCHE BODENFELDER
      REAL, INTENT(IN) :: PS(IE,JE,3), QDB(IE,JE,3)
C
C     ATMOSPHAEREN-FELDER
      REAL, INTENT(IN) :: U(IE,JE,KE,3), V (IE,JE,KE,3),
     &                    T(IE,JE,KE,3), QD(IE,JE,KE,3)
C
C     FELDER DER KOEFFIZIENTEN UND FLUESSE (PHYSIKALISCHE UP'S)
      REAL, INTENT(IN) :: TMKVMH(IE*(KE-1),JE,2)
C
      REAL, INTENT(IN) :: TMCH(IE,JE)
C
C     LOKALE FELDER DIMENSIONIEREN
C     ----------------------------
      REAL    :: ZTV(IE,KE,3), ZPHF(IE,KE   ),
     &           ZDP(IE,KE,3), ZPNF(IE,KE1,3)
C
      REAL    :: ZTMKVM(IE,2:KE,3), ZTMKVH(IE,2:KE,3)
C
      REAL    :: ZALOPN(IE,KE1,3), ZALOPH(IE,KE ,3),
     &           ZBETAK(IE,KE ,3),
     &           ZPPHI (IE,KE ,2), ZPLAM (IE,KE   ),
     &           ZGU   (IE,KE ,2), ZGV   (IE,KE ,3),
     &           ZETAS (IE,KE1  ), ZSDIV (IE,KE1  )
C
      REAL    :: ZALPOM(IE,KE),  ZQADV(IE,KE),
     &           ZEDDPQ(IE,KE)
C
C     FELDER FUER SUBROUTINE GAUSS
      REAL    :: AGA(IE,KE,4), AGB(IE,KE,4),
     &           AGC(IE,KE,4), AGD(IE,KE,4),
     &           AGE(IE,KE,4)
C
C     UEBERGABEFELDER
      REAL    :: GMQTEN(IE,JE,KE),GMVERT(IE,JE,KE)
C
      REAL    :: ZFDIVO,ZFADVY,ZFADVX,ZALOG2,ZAGCT,ZAGCM,ZAGAT,ZAGAM,
     &           ZA2A,ZQ3,ZQ2,ZQ1,ZMYA,ZGVCU,ZGVCO,ZFDIVX,ZFDIVU,ZX2,
     &           ZTMKVHM,ZTKVZ,ZTKVL,ZQ4,ZA2,ZA1A,ZA1,ZA3,ZX1
      INTEGER :: I,INDEX2,J,JM2,JM3,JO2,JO3,JSP,JU3,K,JZU,JZM
C
!DIR$ NOTASK
C
C-----------------------------------------------------------------------
C     VORBEREITENDE MASSNAHMEN
C     PARAMETER ZUR BERECHNUNG DER DIFFUSION
C     UND LN(2) SETZEN.
      ALCNVA  = 0.5
      ZA1A    = ALCNVA
      IF(ZA1A.LT.0.5) ZA1A = 0.5
      ZA2A    = (1.-ZA1A)
      ZALOG2  = ALOG  (2.)
C
C     ZTKVZ , ZTKVL : GEWICHTSFAKTOREN ZUR HORIZONTALEN MITTELUNG DER
C                     VERTIKALEN TRANSPORTKOEFFIZIENTEN
      ZTKVZ = 0.9
      ZTKVL = (1.-ZTKVZ)*0.25
C
C-----------------------------------------------------------------------
C
C     DIE RECHNUNG ERFOLGT 'SCHEIBENWEISE' VON J = JATKON BIS JETKON.
C     IN EINEM TASK WERDEN JEWEILS ZUSAMMENHAENGENDE BEREICHE VON
C     J = JATKON BIS JETKON BEHANDELT; DIE TASK-EINTEILUNG UND
C     STEUERUNG ERFOLGT IM UP *KONORG*; MAXIMAL KOENNEN VIER TASKS
C     PARALLEL LAUFEN; D.H. DER BEREICH JAH BIS JEH WIRD DANN IN VIER
C     UNTERBEREICHE GETEILT.
C
C***********************************************************************
C*                                                                     *
C*    ZUORDNUNG VON SCHEIBENINDIZES  'SUED'--------->-----------'NORD' *
C*    -----------------------------                                    *
C*                                              J-1      J      J+1    *
C*    U, T, QD, PHY.UP'S                         +       +       +     *
C*    V                                              +       +         *
C*                                                                     *
C*    ZALOPN(K+1) ALOG( P(K+1/2) )              JU3     JM3     JO3    *
C*    ZDP (K)                                   JU3     JM3     JO3    *
C*    ZPNF(K)                                   JU3     JM3     JO3    *
C*    ZTV (K)                                   JU3     JM3     JO3    *
C*    ZALOPH(K)   ALOG( P(K) )                  JU3     JM3     JO3    *
C*    ZBETAK(K)   HILFSGROESSE F. FIH, ALOPN    JU3     JM3     JO3    *
C*    ZGU(K)      .5*(DP(I+1)+DP(I))*U                  JM2     JO2    *
C*    ZGV(K)      .5*(DP(J+1)+DP(J))*V              JU3     JM3     JO3*
C*    ZPPHI(K)    R*TV*GRADY(LN(P))                 JM2     JO2        *
C*    ZSDIV(K+1)  VERT. SUMME DER DIVERGENZEN           JM2            *
C*    ZETAS(K+1)  MODIF. VERTIKALGESCHWINDIGK.          JM2            *
C*    ZTMKVH(K)   TMKVH(K)                      JU3     JM3     JO3    *
C*                                                                     *
C*    ZPHF(K)                                            +             *
C*    ZPLAM(K)    R*TV*GRADX(LN(P))                      +             *
C*    ZALPOM(K)   ALPHA*OMEGA                            +             *
C*    ZPLAM(K)    DRUCKGRADIENT                          +             *
C*    ZUTK(K)     U-TENDENZ (KONVEKTIV)                  +             *
C*    ZVTK(K)     V-TENDENZ (KONVEKTIV)                  +             *
C*    SOWIE WEITERE LOKALE FELDER OHNE SCHEI-            +             *
C*    BENINDEX SIND IN DER J-FLAECHE DEFINIERT           +             *
C*                                                                     *
C*                                              J-1      J      J+1    *
C*                                                                     *
C***********************************************************************
C
C     VORBESETZUNG LOKALER FELDER AM OBER- UND UNTERRAND
C     --------------------------------------------------
      DO I = IAA , IEA
         ZPNF  (I,1,1) = 0.
         ZPNF  (I,1,2) = 0.
         ZPNF  (I,1,3) = 0.
         ZALOPN(I,1,1) = 0.
         ZALOPN(I,1,2) = 0.
         ZALOPN(I,1,3) = 0.
         ZSDIV (I,1  ) = 0.
         ZETAS (I,1  ) = 0.
         ZBETAK(I,1,1) = ZALOG2
         ZBETAK(I,1,2) = ZALOG2
         ZBETAK(I,1,3) = ZALOG2
         ZETAS(I,KE1 ) = 0.
         AGA(I, 1,2)   = 0.
         AGC(I,KE,2)   = 0.
      ENDDO
C
C     ANFANGSINDIZES FUER ZYKLISCHES UMSPEICHERN SETZEN
      JM2 = 1
      JO2 = 2
C
      JU3 = 1
      JM3 = 2
      JO3 = 3
C
C     LOKALE FELDER AM 'SUEDRAND' (J=JATKON, JATKON-1) VORBESETZEN
C     ------------------------------------------------------------
C
C     ACHTUNG:
C     UM BEIM AUTOTASKING (PARALLER DURCHLAUF VON UP *EMTIED*) IDENTI-
C     SCHE ERGEBNISSE ZU ERZIELEN, MUSS DARAUF GEACHTET WERDEN, DASS
C     DIE FORMULIERUNG DER TERME FUER DEN 'SUEDRAND' UND FUER DAS INNE-
C     RE DES MODELLGEBIETES EXAKT IDENTISCH SIND. SONST KANN DURCH DIE
C     UNVERMEIDLICHEN RUNDUNGSFEHLER DAS ERGEBNIS JE NACH ANZAHL DER
C     TASKS VERSCHIEDEN SEIN.
C
      JZU = MAX (1, JATKON-1)
      JZM = JATKON
C
C     GEPACKTE FELDER DER PHYSIKALISCHEN TENDENZEN UND KOEFFIZIENTEN
C     AUSPACKEN
      CALL COPYRE(TMKVMH(1,JZU,1),ZTMKVM(1,2,JU3),IE*(KE-1))
      CALL COPYRE(TMKVMH(1,JZU,2),ZTMKVH(1,2,JU3),IE*(KE-1))
      CALL COPYRE(TMKVMH(1,JZM,1),ZTMKVM(1,2,JM3),IE*(KE-1))
      CALL COPYRE(TMKVMH(1,JZM,2),ZTMKVH(1,2,JM3),IE*(KE-1))
C
      ZX2 = .5*RERD*(CPHI(JZU,1)+CPHI(JZM,1))
C
      DO K = 1 , KE
         DO I = IAA , IEA
            ZDP   (I,K  ,JU3) = DAK(K  ) + DBK(K  )*PS(I,JZU,NJ)
            ZDP   (I,K  ,JM3) = DAK(K  ) + DBK(K  )*PS(I,JZM,NJ)
            ZPNF  (I,K+1,JU3) = AK (K+1) + BK (K+1)*PS(I,JZU,NJ)
            ZPNF  (I,K+1,JM3) = AK (K+1) + BK (K+1)*PS(I,JZM,NJ)
            ZALOPN(I,K+1,JU3) = ALOG  (ZPNF(I,K+1,JU3))
            ZALOPN(I,K+1,JM3) = ALOG  (ZPNF(I,K+1,JM3))
            ZTV   (I,K  ,JU3) =     T(I,JZU,K,NJ) * ( 1.0 + RDDRM1*
     &                             QD(I,JZU,K,NJ))
            ZTV   (I,K  ,JM3) =     T(I,JZM,K,NJ) * ( 1.0 + RDDRM1*
     &                             QD(I,JZM,K,NJ))
         ENDDO
      ENDDO

      DO I = IAA, IEA
         ZALOPH(I,1 ,JU3) = ZALOPN(I,2,JU3) - ZBETAK(I,1 ,JU3)
         ZALOPH(I,1 ,JM3) = ZALOPN(I,2,JM3) - ZBETAK(I,1 ,JM3)
      ENDDO

      DO K = 2, KE-1
         DO I = IAA , IEA
            ZBETAK(I,K,JU3) = 1. - ZPNF(I,K  ,JU3)/ZDP   (I,K,JU3)
     &                        * (ZALOPN(I,K+1,JU3)-ZALOPN(I,K,JU3))
            ZBETAK(I,K,JM3) = 1. - ZPNF(I,K  ,JM3)/ZDP   (I,K,JM3)
     &                        * (ZALOPN(I,K+1,JM3)-ZALOPN(I,K,JM3))
            ZALOPH(I,K,JU3) = ZALOPN(I,K+1,JU3) - ZBETAK(I,K,JU3)
            ZALOPH(I,K,JM3) = ZALOPN(I,K+1,JM3) - ZBETAK(I,K,JM3)
         ENDDO
      ENDDO

      DO I = IAA , IEA
         ZBETAK(I,KE,JU3) = 1. - ZPNF(I,KE ,JU3) / ZDP   (I,KE,JU3)
     &                      * (ZALOPN(I,KE1,JU3) - ZALOPN(I,KE,JU3))
         ZBETAK(I,KE,JM3) = 1. - ZPNF(I,KE ,JM3) / ZDP   (I,KE,JM3)
     &                      * (ZALOPN(I,KE1,JM3) - ZALOPN(I,KE,JM3))
         ZALOPH(I,KE,JU3) =    ZALOPN(I,KE1,JU3) - ZBETAK(I,KE,JU3)
         ZALOPH(I,KE,JM3) =    ZALOPN(I,KE1,JM3) - ZBETAK(I,KE,JM3)
      ENDDO

      DO K = 1, KE
         DO I = IAA , IEH
            ZGU(I,K,JM2) = .5*(ZDP(I,K,JM3)+ZDP(I+1,K,JM3))
     &           *U(I,JZM,K,NJ)
         ENDDO

C     WEITERE LOKALE FELDER AN H/V-GITTERPUNKTEN VORBESETZEN
         ZX1    = 0.5*R*EDADPHI
         ZFDIVX = ACPHIR(JZM,1)*EDDLAM
         ZFDIVO = ACPHIR(JZM,1)*EDDPHI*CPHI(JZM,2)
         ZFDIVU = ACPHIR(JZM,1)*EDDPHI*CPHI(JZU,2)

         DO I = IAH , IEH
            ZGV(I,K,JU3) =  .5*(ZDP(I,K,JU3)+
     &           (DAK(K)+DBK(K)*PS(I,JZM  ,NJ)))*
     &                           V(I,JZU,K,NJ)
            ZGV(I,K,JM3) =  .5*(ZDP(I,K,JM3)+
     &           (DAK(K)+DBK(K)*PS(I,JZM+1,NJ)))*
     &                           V(I,JZM,K,NJ)
            ZPPHI(I,K,JM2) = ZX1*(ZTV(I,K,JM3) + ZTV   (I,K,JU3))*
     &                        (ZALOPH(I,K,JM3) - ZALOPH(I,K,JU3))
         ENDDO
      ENDDO

C     BEGINN DER SCHLEIFE IN J-RICHTUNG
C     ---------------------------------

C     ACHTUNG:
C     UM BEIM AUTOTASKING (PARALLER DURCHLAUF VON UP *EMTIED*) IDENTI-
C     SCHE ERGEBNISSE ZU ERZIELEN, MUSS DARAUF GEACHTET WERDEN, DASS
C     DIE FORMULIERUNG DER TERME FUER DEN 'SUEDRAND' UND FUER DAS INNE-
C     RE DES MODELLGEBIETES EXAKT IDENTISCH SIND. SONST KANN DURCH DIE
C     UNVERMEIDLICHEN RUNDUNGSFEHLER DAS ERGEBNIS JE NACH ANZAHL DER
C     TASKS VERSCHIEDEN SEIN.

      DO J = JATKON, JETKON

C        LOKALE FELDER FUER J+1 BESETZEN
C        GEPACKTE FELDER DER PHYSIKALISCHEN TENDENZEN UND KOEFFIZIENTEN
C        AUSPACKEN
         CALL COPYRE(TMKVMH(1,J+1,1),ZTMKVM(1,2,JO3),IE*(KE-1))
         CALL COPYRE(TMKVMH(1,J+1,2),ZTMKVH(1,2,JO3),IE*(KE-1))

         ZX1 = .5*R*ACPHIR(J,1)*EDDLAM
         ZX2 = .5*RERD*(CPHI(J,1)+CPHI(J+1,1))

         DO K = 1 , KE
            DO I = IAA , IEA
               ZDP   (I,K  ,JO3) = DAK(K  ) + DBK(K  )*PS(I,J+1,NJ)
               ZPHF  (I,K      ) = AKH(K  ) + BKH(K  )*PS(I,J  ,NJ)
               ZPNF  (I,K+1,JO3) = AK (K+1) + BK (K+1)*PS(I,J+1,NJ)
               ZALOPN(I,K+1,JO3) = ALOG  (ZPNF(I,K+1,JO3))
               ZTV   (I,K  ,JO3) = T(I,J+1,K,NJ) * ( 1.0 + RDDRM1*
     &                            QD(I,J+1,K,NJ))
            ENDDO
         ENDDO

         DO I = IAA , IEA
            ZALOPH(I,1,JO3) = ZALOPN(I,2,JO3) -ZBETAK(I,1,JO3)
         ENDDO

         DO K  = 2 , KE-1
            DO I = IAA , IEA
               ZBETAK(I,K,JO3) = 1. - ZPNF(I,K  ,JO3)/ZDP   (I,K,JO3)
     &                           * (ZALOPN(I,K+1,JO3)-ZALOPN(I,K,JO3))
               ZALOPH(I,K,JO3) = ZALOPN(I,K+1,JO3) - ZBETAK(I,K,JO3)
            ENDDO
         ENDDO

         DO I = IAA , IEA
            ZBETAK(I,KE,JO3) = 1. - ZPNF(I,KE ,JO3)  /ZDP   (I,KE,JO3)
     &                         * (ZALOPN(I,KE1,JO3) - ZALOPN(I,KE,JO3))
            ZALOPH(I,KE,JO3) =    ZALOPN(I,KE1,JO3) - ZBETAK(I,KE,JO3)
         ENDDO

         DO K = 1 , KE
            DO I = IAA , IEH
               ZGU(I,K,JO2) = .5*(ZDP(I,K,JO3)+ZDP(I+1,K,JO3))
     &              *U(I,J+1,K,NJ)
               ZPLAM(I,K)   = ZX1*( ZTV(I+1,K,JM3) +    ZTV(I,K,JM3) )
     &                         *(ZALOPH(I+1,K,JM3) - ZALOPH(I,K,JM3))
            ENDDO
            DO I = IAH , IEH
               ZSDIV(I,K+1)   = ZSDIV(I,K) +
     &              ZFDIVX*( ZGU(I,K,JM2) - ZGU(I-1,K,JM2) ) +
     &              ZFDIVO*  ZGV(I,K,JM3) - ZFDIVU*ZGV(I,K,JU3)
            ENDDO
         ENDDO

         DO K = 2 , KE
            DO I = IAH , IEH
               ZETAS(I,K) = BK(K)*ZSDIV(I,KE1) - ZSDIV(I,K)
            ENDDO
         ENDDO

         ZX1 = 0.5*R*EDADPHI
         IF(J.LT.JEH) THEN
            ZFDIVX = ACPHIR(J+1,1)*EDDLAM
            ZFDIVO = ACPHIR(J+1,1)*EDDPHI*CPHI(J+1,2)
            ZFDIVU = ACPHIR(J+1,1)*EDDPHI*CPHI(J  ,2)

            DO K = 1 , KE
               DO I = IAH , IEH
                  ZGV(I,K,JO3) = .5*(ZDP(I,K,JO3)+
     &                 (DAK(K)+DBK(K)*PS(I,J+2,NJ)))*V(I,J+1,K,NJ)
                  ZPPHI(I,K,JO2) =  ZX1*(ZTV(I,K,JO3) +    ZTV(I,K,JM3))
     &                              *(ZALOPH(I,K,JO3) - ZALOPH(I,K,JM3))
               ENDDO
            ENDDO

         ELSE
            DO K = 1 , KE
               DO I = IAH , IEH
                  ZPPHI(I,K,JO2) =  ZX1*(ZTV(I,K,JO3) +    ZTV(I,K,JM3))
     &                              *(ZALOPH(I,K,JO3) - ZALOPH(I,K,JM3))
               ENDDO
            ENDDO
         ENDIF

C        FEUCHTEKONVERGENZ UND OMEGA BERECHNEN FUER KONVEKTION
C        -----------------------------------------------------

         ZX2   =  RERD*ACPHIR(J,1)
         ZFADVX= .5*ACPHIR(J,1)*EDDLAM
         ZFADVY= .5*ACPHIR(J,1)*EDDPHI

         DO K = 1 , KE
            DO I = IAH , IEH
               ZEDDPQ(I,K) = 1./ZDP(I,K,JM3)
               ZGVCO       = ZGV(I,K,JM3)*CPHI(J  ,2)
               ZGVCU       = ZGV(I,K,JU3)*CPHI(J-1,2)
               ZQ1         = ZGU(I  ,K,JM2)
     &              *( QD(I+1,J  ,K,NJ)-QD(I  ,J  ,K,NJ) )
               ZQ2         = ZGU(I-1,K,JM2)
     &              *( QD(I  ,J  ,K,NJ)-QD(I-1,J  ,K,NJ) )
               ZQ3         = ZGVCO*( QD(I,J+1,K,NJ)-QD(I,J  ,K,NJ) )
               ZQ4         = ZGVCU*( QD(I,J  ,K,NJ)-QD(I,J-1,K,NJ) )
               ZQADV(I,K)  = -( ZFADVX*(ZQ1+ZQ2) + ZFADVY*(ZQ3+ZQ4) )
     &              *ZEDDPQ(I,K)

C              ALPHA*OMEGA
               ZA1 =  ZGU(I,K,JM2)*ZPLAM(I,K) + ZGU(I-1,K,JM2)
     &              *ZPLAM(I-1,K)
               ZA2 = ( ZGVCO*ZPPHI(I,K,JO2) + ZGVCU*ZPPHI(I,K,JM2) )*ZX2
               ZA3 = - R*ZTV(I,K,JM3)*ZEDDPQ(I,K)*
     &         (   (ZALOPN(I,K+1,JM3)-ZALOPN(I,K,JM3))*ZSDIV (I,K    )
     &           + (ZSDIV (I,K+1    )-ZSDIV (I,K    ))*ZBETAK(I,K,JM3) )
               ZALPOM(I,K) = .5*( ZA1 + ZA2 )*ZEDDPQ(I,K) + ZA3

               AGB(I,K,2) = ED2DT
               AGD(I,K,2) = ED2DT*QD(I,J,K,NA) + ZQADV(I,K)
            ENDDO
         ENDDO

C     VERTIKALADVEKTION UND -DIFFUSION (IMPLIZIT)
         DO I = IAH , IEH
            ZAGCM      = .5*ZETAS(I,2)*ZEDDPQ(I,1)
            ZTMKVHM    = ZTKVZ*ZTMKVH(I  ,2,JM3) + ZTKVL*
     &                       ( ZTMKVH(I+1,2,JM3) + ZTMKVH(I,2,JO3) +
     &                         ZTMKVH(I-1,2,JM3) + ZTMKVH(I,2,JU3) )
            ZAGCT      = - ZTMKVHM*ZEDDPQ(I,1)
            AGC(I,1,2) = ZAGCM*ZA1A + ZAGCT*A1T(2)
            AGB(I,1,2) = AGB(I,1,2) - AGC(I,1,2)
            AGD(I,1,2) = AGD(I,1,2) - ( ZA2A*ZAGCM + A2T(2)*ZAGCT )*
     &           ( QD(I,J,2,NA) - QD(I,J,1,NA) )
         ENDDO

         DO K = 2 , KE-1
            DO I = IAH , IEH
               ZAGAM      = -.5*ZETAS(I,K  )*ZEDDPQ(I,K)
               ZAGCM      =  .5*ZETAS(I,K+1)*ZEDDPQ(I,K)
               ZTMKVHM = ZTKVZ*ZTMKVH(I  ,K,JM3) + ZTKVL*
     &                       ( ZTMKVH(I+1,K,JM3) + ZTMKVH(I,K,JO3)
     &                      +  ZTMKVH(I-1,K,JM3) + ZTMKVH(I,K,JU3) )
               ZAGAT      = - ZTMKVHM*ZEDDPQ(I,K)
               ZTMKVHM = ZTKVZ*ZTMKVH(I  ,K+1,JM3) + ZTKVL*
     &                       ( ZTMKVH(I+1,K+1,JM3) + ZTMKVH(I,K+1,JO3)
     &                      +  ZTMKVH(I-1,K+1,JM3) + ZTMKVH(I,K+1,JU3) )
               ZAGCT      = - ZTMKVHM*ZEDDPQ(I,K)
               AGA(I,K,2) = ZAGAM*ZA1A + ZAGAT*A1T(K)
               AGC(I,K,2) = ZAGCM*ZA1A + ZAGCT*A1T(K+1)
               AGB(I,K,2) = AGB(I,K,2) - AGA(I,K,2) - AGC(I,K,2)
               AGD(I,K,2) = AGD(I,K,2) - ( ZA2A*ZAGAM + A2T(K  )*ZAGAT )
     &                                  *(QD(I,J,K-1,NA) - QD(I,J,K,NA))
     &                                 - ( ZA2A*ZAGCM + A2T(K+1)*ZAGCT )
     &                                  *(QD(I,J,K+1,NA) - QD(I,J,K,NA))
            ENDDO
         ENDDO

         DO I = IAH , IEH
            ZAGAM     = -.5*ZETAS(I,KE)*ZEDDPQ(I,KE)
            ZTMKVHM   = ZTKVZ*ZTMKVH(I  ,KE,JM3) + ZTKVL*
     &                      ( ZTMKVH(I+1,KE,JM3) + ZTMKVH(I  ,KE,JO3) +
     &                        ZTMKVH(I-1,KE,JM3) + ZTMKVH(I  ,KE,JU3) )
            ZAGAT       = -  ZTMKVHM *ZEDDPQ(I,KE)
            ZAGCT       = - TMCH(I,J)*ZEDDPQ(I,KE)
            AGA(I,KE,2) = ZAGAM*ZA1A  + ZAGAT*A1T(KE)
            AGB(I,KE,2) = AGB(I,KE,2) - AGA(I,KE,2) - ZAGCT*A1T(KE1)
            AGD(I,KE,2) = AGD(I,KE,2) - ( ZA2A*ZAGAM + A2T(KE)*ZAGAT )*
     &                        ( QD(I,J,KE-1,NA) - QD(I,J,KE,NA) )
     &           - A2T(KE1)*ZAGCT*(QDB(I,J,NA)  - QD(I,J,KE,NA) )
     &           - A1T(KE1)*ZAGCT* QDB(I,J,NA)
         ENDDO

C     GAUSS - ELIMINATION UND BERECHNUNG DER FEUCHTETENDENZ
C     DER ZEITPUNKT N+1 (=NE) STEHT IN AGE
         INDEX2 = 2
         CALL GAUSS (IAH, IEH, INDEX2, AGA, AGB, AGC, AGD, AGE)

C     VORBESETZUNG DER TENDENZEN UND NIEDERSCHLAGSRATEN
         DO K = 1, KE
            DO I = 1, IE
               GMQTEN(I,J,K)= 0.0
               GMVERT(I,J,K)= 0.0
            ENDDO
         ENDDO

C     WENN IN DIESER SCHEIBE KEINE KONVEKTION, DANN NAECHSTE SCHEIBE
C     BERECHNEN

         DO K = 1,KE
            DO I = IAH,IEH
               ZMYA         = 1.0 - ALPHABOUND(I,J,1)
               GMQTEN(I,J,K)= ( AGE(I,K,2) - QD(I,J,K,NA) )*ED2DT*ZMYA
               GMVERT(I,J,K)= ZALPOM(I,K)*ZPHF(I,K)/(R*ZTV(I,K,JM3))
            ENDDO
         ENDDO

C     ZYKLISCHE VERTAUSCHUNG DER SCHEIBENINDIZES
         JM2 = 3 - JM2
         JO2 = 3 - JO2

         JSP = JU3
         JU3 = JM3
         JM3 = JO3
         JO3 = JSP

C     ENDE DER SCHLEIFE UEBER DIE REIHEN VON N - S
      ENDDO

      RETURN
      END SUBROUTINE ECTIED
