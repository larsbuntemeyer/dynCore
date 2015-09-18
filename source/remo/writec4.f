      SUBROUTINE WRITEC4(YTYP  , YXVARN , NZVX   , NZTX   , NX    ,
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
C     SUBROUTINE WRITEC4
C
C**** WRITEC4  -   UP:DIE GEWUENSCHTEN FELDER IN DIE ERGEBNISDATEI
C****                 SCHREIBEN; PDB UND GDB BESETZEN.
C**   AUFRUF   :   CALL WRITEC4(YTYP, YXVARN, NZVX, NZTX, NX)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   DIE IN DER LISTE 'YXVARN(NZVX)' ANGEGEBENEN FELDER
C**                AUS DEM LANGZEIT-SPEICHER HOLEN UND IN DIE ERGEBNIS-
C**                DATEI SCHREIBEN. EVT. ZUSAETZLICHE SONDERFELDER
C**                BERECHNEN.
C**   VERSIONS-
C**   DATUM    :   07.03.05
C**                2007
C**
C**   EXTERNALS:   MAKEPDB, MAKEGDB, PRCV, PUTEC4
C**
C**   EINGABE-
C**   PARAMETER:   YTYP :  TYP DER ERGEBNISDATEN:
C**                       'E': 2-D ZEITREIHENDATEN    FUER GESAMTGEBIET
C**                       'D':         ERGEBNISDATEN  FUER TEILGEBIET
C**                       'F': FORTSETZUNGSDATEN      FUER GESAMTGEBIET
C**                       'T': NORMALE ERGEBNISDATEN  FUER GESAMTGEBIET
C**                YXVARN: TABELLE DER FELDNAMEN, DIE AUSGEGEBEN WERDEN
C**                        SOLLEN
C**                NZVX  : LAENGE DER TABELLE YXVARN
C**                NZTX  : ZEITSCHRITT (FUER PDB)
C**                NX    : ZEITEBENE FUER FELDER MIT DREI ZEITEBENEN
C**   AUSGABE-
C**   PARAMETER:   KEINE
C**
C**   COMMON-
C**   BLOECKE  :   ORG, CORG, EMGBCH, EMGBRI, UNITCH, UNITNR,
C**                SUMPAR, UNITZR, GRID, HIGKON
C**
C**   METHODE  :   ERSTELLUNG VON PDB UND GDB; AUSGABE DER GEPACKTEN
C**                FELDER AUF DIE ENTSPRECHENDE EINHEIT
C**   FEHLERBE-
C**   HANDLUNG :   STOP IM FEHLERFALL
C**   VERFASSER:   R.PODZUN
C
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
      INCLUDE "org.h"
      INCLUDE "corg.h"
      INCLUDE "emgbch.h"
      INCLUDE "emgbri.h"
      INCLUDE "unitch.h"
      INCLUDE "unitnr.h"
      INCLUDE "unitzr.h"
      INCLUDE "grid.h"
      INCLUDE "higkon.h"
C
C     Dummy Arguments
C
      INTEGER,   INTENT(IN) :: NZVX   , NZTX   , NX
C
      CHARACTER, INTENT(IN) :: YXVARN(NZVX)*(*), YTYP*(*)
C
      REAL, DIMENSION(KE1), INTENT(IN) :: AK, BK
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
     & TSLECH (IEJE,3), TSWECH(IEJE,3), TSIECH(IEJE,3), QDBL  (IEJE,3),
     & QDBW   (IEJE,3), QDBI  (IEJE,3), SN    (IEJE,3), TD    (IEJE,3),
     & TDCL   (IEJE,3), WL    (IEJE,3), TSN   (IEJE,3), TD3   (IEJE,3),
     & TD4    (IEJE,3), TD5   (IEJE,3)
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
     & FIBOXS (IEJE  )
C
      REAL, INTENT(IN) ::
     & VGRAT  (IEJE  ), VAROR (IEJE  ), VLT   (IEJE  ), T2MAX (IEJE  ),
     & SRAD0U (IEJE  ), SRADSU(IEJE  ), TRADSU(IEJE  ), T2MIN (IEJE  ),
     & SEAICE (IEJE  ), SICED (IEJE  ), FOREST(IEJE  ), TEFF  (IEJE  ),
     & TSMAX  (IEJE  ), TSMIN (IEJE  ), WIMAX (IEJE  ), TOPMAX(IEJE  ),
     & SNMEL  (IEJE  ), TSLIN (IEJE  ), DSNAC (IEJE  ), FAO   (IEJE  ),
     & RGCGN  (IEJE  ), WSMX  (IEJE  ), QVI   (IEJE  ), 
     & ALWCVI (IEJE  ), GLAC  (IEJE  ), DRAIN (IEJE  ), TLAMBDA(IEJE  ),
     & DLAMBDA(IEJE  ), PORVOL(IEJE  ), FCAP  (IEJE  ), WI3   (IEJE,3),
     & WI4    (IEJE,3), WI5   (IEJE,3), WI    (IEJE,3), WICL  (IEJE,3),
     & GHPBL  (IEJE  ), BETA  (IEJE  ), WMINLOK(IEJE ), WMAXLOK(IEJE) ,
     & VBM10M (IEJE  ), CAPE  (IEJE  )
CSP
      REAL, INTENT(IN) ::
     & QI  (IEJE,KE,3), QIVI  (IEJE  ), QIBOXS (IEJE ), 
     & RPRAC (IEJE,KE)
C
      REAL, INTENT(IN) ::
     & WS1    (IEJE  ), WS2   (IEJE  ), WS3   (IEJE  ), WS4   (IEJE  ),
     & WS5    (IEJE  ), DZR   (IEJE  ), DZS   (IEJE  ), FKSAT (IEJE  ),
     & FMPOT  (IEJE  ), BCLAPP(IEJE  ), VPOR  (IEJE  ), ETRANS(IEJE  ),
     & EBSOIL (IEJE  ), ESNOW (IEJE  ), ESKIN (IEJE  ), ERES  (IEJE  )

      REAL ::   PINT(IEJE,KE1,3), DWDT(IEJE,KE ,3),
     &          W   (IEJE,KE1,3)
C
C     Local Variables
C
      REAL :: FTKVM (IEJE,KE)     , FTKVH (IEJE,KE)
      REAL :: ZFTKVM(IE*(KE-1),JE), ZFTKVH(IE*(KE-1),JE)
      REAL :: U10ER(IEJE), V10ER(IEJE), PSRED(IEJE)
C
      REAL :: UFELDFIB(IEJE  ), UFELDP(IEJE  )
      REAL :: GFELDFIB(MOIEJE), GFELDP(MOIEJE)
      REAL :: ZPSHOR  (MOIEJE), VFELD(IEJE)
C
      LOGICAL :: LYXVARN
      REAL    :: ZPOLPHI,ZPHI,ZD2,ZD1,ZCOSPOL,ZBETA,ZARG2,ZARG1,ZARG,
     &           ZSINPOL,ZPOLLAM,ZRLA,ZRLAS,ZZRLA,POLLAMD,FIBTEST
      INTEGER :: NUNIT,NTRI,NTAB2,NTAB1,K,J,IUN,IJ,I,IVLOC
C
C     DATEI-TYP FESTSTELLEN
      IF (YTYP.EQ.'D') THEN
         NUNIT = NUDDAT
      ELSE IF (YTYP.EQ.'F') THEN
         NUNIT = NUFDAT
      ELSE IF (YTYP.EQ.'T') THEN
         NUNIT = NUTDAT
      ELSE IF (YTYP.EQ.'M') THEN
         NUNIT = NUMDAT
      ENDIF
C
      IF (YTYP.EQ.'D'.OR.YTYP.EQ.'F'.OR.YTYP.EQ.'T') THEN
         IF (MYID .EQ. 0) REWIND NUNIT
      ENDIF
C
      CALL SETRA(FTKVM,IEJEKE,0.)
      CALL SETRA(FTKVH,IEJEKE,0.)
      DO J = 1,JE
         CALL COPYRE(TMKVMH(1,J,1),ZFTKVM(1,J),IE*(KE-1))
         CALL COPYRE(TMKVMH(1,J,2),ZFTKVH(1,J),IE*(KE-1))
      ENDDO
      DO K    = 2,KE
         DO J    = 1,JE
            DO I    = 1,IE
               IJ          = I + (J-1)*IE
               FTKVM(IJ,K) = ZFTKVM(I+(K-2)*IE,J)
               FTKVH(IJ,K) = ZFTKVH(I+(K-2)*IE,J)
            ENDDO
         ENDDO
      ENDDO
      IF (YTYP.NE.'F') THEN
C
C        KONSTANTEN BELEGEN
C
C        AUF NN REDUZIERTEN BODENDRUCK BERECHNEN
         CALL PSREDBER(AK, BK,
     &                 T , PS, FIB , PSRED, NX)

         DO IJ=1, IEJE
            UFELDFIB(IJ) = FIB  (IJ)
            UFELDP  (IJ) = PSRED(IJ)
         ENDDO

C        DATEN SAMMELN
         CALL COLLECTDATA(GFELDFIB,UFELDFIB,200,0)
         CALL COLLECTDATA(GFELDP,UFELDP,201,0)

         CALL PSTOP

         IF(MYID.EQ.0) THEN

C           REDUZIERTEN BODENDRUCK UMSPEICHERN
            DO IJ=1, MOIEJE
               ZPSHOR(IJ) = GFELDP(IJ)
            ENDDO

            FIBTEST=9.80665 * 800.
C           REDUZIERTEN BODENDRUCK FILTERN
            CALL HORINT(ZPSHOR, GFELDFIB, FIBTEST)
            CALL FILTER(ZPSHOR, GFELDP)

C           REDUZIERTEN BODENDRUCK VERTEILEN
            CALL SSGB(GFELDP,VFELD)
         ENDIF

C        REDUZIERTEN BODENDRUCK EINSAMMELN
C
C        KONSTANTEN BELEGEN
         TAGCOUNT=1
         TAGTABLE(1)=1
         TYPE=1

         IF(MYID.NE.0) THEN
            CALL PTEST
            CALL PRECVR(VFELD(1))
         ENDIF

         CALL PSTOP

         DO IJ=1, IEJE
            PSRED(IJ)=VFELD(IJ)
         ENDDO
C
C        U10 UND V10 ENTROTIEREN
         ZPOLPHI = POLPHI*DEGRAD
         ZPOLLAM = POLLAM*DEGRAD
         POLLAMD = POLLAM
         ZSINPOL = SIN(ZPOLPHI)
         ZCOSPOL = COS(ZPOLPHI)
         IF (POLLAM.LT.0.0) POLLAMD = 360.0 + POLLAM
C
         DO IJ=1,IEJE
C
            ZRLA    = RLA(IJ)
            ZPHI    = PHI(IJ)
C
C           LAENGE IM ROTIERTEN SYSTEM BERECHNEN
C
            ZZRLA=RLA(IJ)*RADDEG
            IF (ZZRLA.GT.180.0) ZZRLA = ZZRLA - 360.0
            ZZRLA=ZZRLA*DEGRAD
            ZD1=ZZRLA-ZPOLLAM
            ZD2=ZRLA-ZPOLLAM
            ZARG1= - SIN(ZD1)*COS(ZPHI)
            ZARG2= - ZSINPOL*COS(ZPHI)*COS(ZD1)+ZCOSPOL*SIN(ZPHI)
            IF(ABS(ZARG2).LT.1.0 E-20) ZARG2 = 1.0 E-20
            ZRLAS = ATAN2(ZARG1,ZARG2)
C
C           WINKEL ZBETA BERECHEN (SCHNITTWINKEL DER BREITENKREISE)
C
            ZARG = - SIN(ZPOLPHI)*SIN(ZD2)*SIN(ZRLAS) -
     &                            COS(ZD2)*COS(ZRLAS)
            ZARG  = MIN(1.0,  ZARG)
            ZARG  = MAX(-1.0, ZARG)
CHG       AENDERUNG FUER POL
CHG       ZBETA = ACOS(ZARG)
CHG       ZBETA = SIGN(ZBETA, -((RLA(IJ)*RADDEG) - (POLLAMD-180.0)))
            ZBETA = ABS(ACOS(ZARG))
            IF ( -((RLA(IJ)*RADDEG) - (POLLAMD-180.0)) .LT. 0.   .AND.
     &           -((RLA(IJ)*RADDEG) - (POLLAMD-180.0)) .GE. -180.) THEN
               ZBETA = -ZBETA
            ENDIF
C           U10 - WIND TRANSFORMIEREN
            U10ER(IJ)=U10(IJ)*COS(ZBETA)-V10(IJ)*SIN(ZBETA)
C           V10 - WIND TRANSFORMIEREN
            V10ER(IJ)=U10(IJ)*SIN(ZBETA)+V10(IJ)*COS(ZBETA)
         ENDDO
      ENDIF
C
C     HERAUSSUCHEN DER GEWUENSCHTEN FELDER AUS DEM LANGZEITSPEICHER
C     UND AUSSCHREIBEN AUF GRIB-CODE-DATEI
      DO NTAB1 = 1,NZVX
C
C        FELDNAME IN DER TABELLE *YEMNAME* SUCHEN
         LYXVARN=.FALSE.
         DO NTAB2 = 1,NFLDEM
            IF (YXVARN(NTAB1).EQ.YEMNAME(NTAB2)) THEN
               IVLOC = NTAB2
               IF (YTYP.EQ.'E') THEN
                  IUN=MOD(NTAB1,INUMZRM)
                  IF (IUN.EQ.0) IUN=INUMZRM
                  NUNIT=NUEDAT(IUN)
               ENDIF
               LYXVARN=.TRUE.
               EXIT
            ENDIF
         ENDDO
C
C        FELDNAME NICHT IN DER TABELLE DER EM-FELDNAMEN; DIESES FELD
C        KANN NICHT ABGESPEICHERT WERDEN
         IF (.NOT. LYXVARN) THEN
            IF (MYID .EQ. 0) THEN
               CALL PRCV('WRITEC4', 'YXVARN(NTAB1)' , YXVARN(NTAB1))
               CALL PRCV('WRITEC4', 'FEHLERMELDUNG:',
     &              'O.A.FELD UNBEKANNT')
               CALL PRCV('WRITEC4', 'AKTION        ', 
     &              'FELD NICHT GESPEICHERT')
            ENDIF
            CYCLE
         ENDIF
C
C        FELD IN TABELLE DER EM-NAMEN; FELD AUSSCHREIBEN AUF GRIB-CODE-
C        DATEI; TIME RANGE INDICATOR (NTRI) BEACHTEN
C
C        PDB, GDB VORBESETZEN
C        FUER ZEITPUNKT 0 INITIALISIERTE ANALYSE; PDB-SONDERBEHANDLUNG
         IF(NZTX.EQ.0) THEN
            NTRI = 1
         ELSE
C
C        FUER FORTSETZUNGSLAEUFE; PDB-SONDERBEHANDLUNG
            IF(YTYP.EQ.'F') THEN
               NTRI = 10
            ELSE
               NTRI = NGBTRI(IVLOC)
            ENDIF
         ENDIF
         CALL MAKEPDB(YTYP, NZTX, NTRI)
         CALL MAKEGDB(YTYP, YEMNAME(IVLOC), MOIE, MOJE, KE1)
C
         CALL PUTEC4
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

      ENDDO
C
      IF (YTYP.EQ.'D'.OR.YTYP.EQ.'F'.OR.YTYP.EQ.'T') THEN
         IF (MYID .EQ. 0) CLOSE(NUNIT)
      ENDIF
C
      END SUBROUTINE WRITEC4
