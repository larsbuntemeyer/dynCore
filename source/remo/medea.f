      SUBROUTINE MEDEA(MOIE  , MOJE , MOKE, MOIEJE , MOIEJEKE,
     &                 MOIEKE, MOKE1, MYID, NPROCXM, NPROCYM)
      !
      IMPLICIT NONE
      !
      !C
      !C**** MEDEA    -   UP:EINGABEDATEN FUER REMO LESEN UND IN DEN
      !C****                 COMMON-BLOECKEN BEREITSTELLEN
      !C**   AUFRUF   :   CALL MEDEA IN HP: REMORG
      !C**   ZWECK    :   NAMELIST-DATEN FUER DIE COMMON-BLOECKE BEREITSTELLEN
      !C**   VERSIONS-
      !C**   DATUM    :   16.08.05
      !C**                2007
      !C**
      !C**   EXTERNALS:   LENCH
      !C**
      !C**   EINGABE-
      !C**   PARAMETER:   KEINE
      !C**   AUSGABE-
      !C**   PARAMETER:   KEINE
      !C**
      !C**   COMMON-
      !C**   BLOECKE  :   ORG     , CORG  , CORG2, COMNMI, COMPHY, COMDIA,
      !C**                COMECPHY, COMDYN, GRID , UNITNR
      !C**   FEHLERBE-
      !C**   HANDLUNG :   DURCH STOP ERROR
      !C**
      !C**   VERFASSER:   R.PODZUN/ V.WOHLGEMUTH/ D.MAJEWSKI
      !C
      INCLUDE "org.h"
      INCLUDE "corg.h"
      INCLUDE "corg2.h"
      INCLUDE "comnmi.h"
      INCLUDE "comphy.h"
      INCLUDE "comdia.h"
      INCLUDE "comecphy.h"
      INCLUDE "comdyn.h"
      INCLUDE "grid.h"
      INCLUDE "unitnr.h"
      INCLUDE "YOMAER"
      INCLUDE "comconv.h"
      INCLUDE "COMCON"
      INCLUDE "comhyd.h"
      !
      ! Formal Parameters
      !
      INTEGER, INTENT(OUT) :: MOIE,MOJE,MOKE
      INTEGER, INTENT(OUT) :: MOIEJE,MOIEJEKE,MOIEKE,MOKE1
      INTEGER, INTENT(IN)  :: MYID
      INTEGER, INTENT(OUT) :: NPROCXM, NPROCYM
      !
      ! Local Parameters
      !
      INTEGER, PARAMETER   :: NZVDD = 68
      INTEGER, PARAMETER   :: NZVTD = 156
      INTEGER, PARAMETER   :: NZVMD = 123
      INTEGER, PARAMETER   :: NZVND = 12
      !
      !
      ! Local Variables
      !
      CHARACTER(132)       :: YZEILE
      INTEGER              :: I,II,III
      INTEGER              :: ILYDCAT,IZAIA,NHANF,NHDDA,NHDEA
      INTEGER              :: NHDFOR,NHDR,NHDAA,NHDTA,NHEAA,NHENDE
      INTEGER              :: NHFORA,NHORPHD,NHTAA,NZVD2,NZVM2,NZVN2
      INTEGER              :: NZVT2
      REAL                 :: VBMXVD
      !
      INTEGER, EXTERNAL    :: LENCH
C------  NAMELIST INITIALISIERUNG  (DEFAULTS) -------------------------

C    /PARCTL/
      INTEGER :: NPROCXD, NPROCYD

C    /EMGRID/
      INTEGER :: MOIED,MOJED,MOKED

      REAL    :: POLLAMD,POLPHID,DLAMD,DPHID,RLALUD,PHILUD



C    /RUNCTL/
      CHARACTER(10) :: YADATD

      LOGICAL   LQWRD  , LMOMOND, LMOMITD, LSCEND , LGMOND ,
     &          LTAMITD, LNOHALOD
      INTEGER   NANFD  , NENDED , NEAAD  , NHANFD ,
     &          NHENDED, NDEAD  , NHEAAD , NHDEAD , NFORAD ,
     &          NDFORD , NHFORAD, NHDFORD, NDAAD  , NDDAD  ,
     &          NHDAAD , NHDDAD , NDRD   , NHDRD  , NTAAD  ,
     &          NDTAD  , NHTAAD , NHDTAD , NHDMXND,
     &          IADGBD , IEDGBD , JADGBD , JEDGBD
      REAL      DTD    , DLANDD , DNOPRCD

C    /DYNCTL/
      LOGICAL   LSITSD , LVNESTD
      REAL      EPSASSD

C    /HYDCTL/
      LOGICAL :: LHYDROD

C    /PHYCTL/
      LOGICAL   LPHYD  , LECGADD, LAEROZD, L5LAYD , LWDIFD
      LOGICAL   LAKKUD , LNEARD , LVEGD  , LINBOXD, LSICEDD
      LOGICAL   LECDIUD, LECRADD, LECSOLD, LECAERD, LECCFCD
      LOGICAL   LECVDID, LECSURD, LECGWDD, LECCOVD, LECCODD
      LOGICAL   LCOLDCD
      REAL      HDRADD
      INTEGER   IEXCD  , NDFAERD(12)

C    /NMICTL/
      LOGICAL   LANMID , LDNMID
      INTEGER   NVMD   , NITNMID
      REAL      DTNMID

C    /PRICTL/
      LOGICAL   LDIAD

C    /DATEN/
      CHARACTER YADEND * 3, YADCATD*56,
     &          YRDEND * 3, YRDCATD*56,
     &          YEDEND * 3, YEDCATD*56,
     &          YDDEND * 3, YDDCATD*56,
     &          YFDEND * 3, YFDCATD*56,
     &          YTDEND * 3, YTDCATD*56,
     &          YBDCATD*64, YBDNAMD(NZMXVB)*64,
     &          YNDCATD*56, YMDCATD*56,
     &          YUSERAD* 3, YUSERED* 3,
     &          YGDCATD *64, YGDNAMD *64,
     &          YHDCATD *64, YHDNAMD *64,
     &          YO3DCATD*64, YO3DNAMD*64,
     &          YSADCATD*64, YSADNAMD*64,
     &          YSNDCATD*64, YSNDNAMD*64

C---- NAMELIST ---------------------------------------------------------

      NAMELIST  /PARCTL/ NPROCXM, NPROCYM

      NAMELIST  /EMGRID/ MOIE  , MOJE  , MOKE,
     &                   POLLAM, POLPHI, DLAM, DPHI, RLALU, PHILU


      NAMELIST  /RUNCTL/ NANF  , NENDE , NHANF , NHENDE,
     &                   YADAT , NEAA  , NDEA  , NHEAA , NHDEA ,
     &                   NFORA , NDFOR , NHFORA, NHDFOR,
     &                   NDAA  , NDDA  , NHDAA , NHDDA ,
     &                   NTAA  , NDTA  , NHTAA , NHDTA ,
     &                   LMOMON, DT    , LQWR  , LMOMIT, LSCEN ,
     &                   NDR   , NHDR  , NHDMXN, LGMON ,
     &                   IADGB , IEDGB , JADGB , JEDGB ,
     &                   DLAND , DNOPRC, LTAMIT, LNOHALO

      NAMELIST  /DYNCTL/ LSITS , EPSASS, LVNEST, VBMXV

      NAMELIST  /HYDCTL/ LHYDRO

      NAMELIST  /PHYCTL/ LPHY  , LECGAD, LAEROZ, L5LAY , LWDIF ,
     &                   LAKKU , LNEAR , LVEG  , LINBOX, LSICED,
     &                   LECDIU, LECRAD, LECSOL, LECAER, LECCFC,
     &                   LECVDI, LECSUR, LECGWD, LECCOV, LECCOD,
     &                   LCOLDC, HDRAD , IEXC  , NDFAER

      NAMELIST  /NMICTL/ LANMI , LDNMI , NVM   , NITNMI, DTNMI

      NAMELIST  /PRICTL/ LDIA

      NAMELIST  /DATEN/  YADEN , YADCAT,
     &                   YRDEN , YRDCAT,
     &                   YEDEN , YEDCAT,
     &                   YDDEN , YDDCAT, YDVARN,
     &                   YTDEN , YTDCAT, YTVARN,
     &                   YFDEN , YFDCAT,
     &                   YMDCAT, YMVARN,
     &                   YNDCAT, YNVARN,
     &                   YUSERA, YUSERE,
     &                   YBDCAT, YBDNAM,
     &                   YGDCAT, YGDNAM,
     &                   YHDCAT, YHDNAM,
     &                   YO3DCAT, YO3DNAM,
     &                   YSADCAT, YSADNAM,
     &                   YSNDCAT, YSNDNAM
C---- NAMELIST DEFAULTS-------------------------------------------------

C    /PARCTL/
      NPROCXD  =    1
      NPROCYD  =    1

C    /PHYCTL/
      DATA NDFAERD / 12*0 /

C    /EMGRID/    EM3-GITTER
      MOIED    =  129
      MOJED    =  101
      MOKED    =  27
      POLLAMD  = -170.0
      POLPHID  =   32.5
      DLAMD    =    0.5
      DPHID    =    0.5
      RLALUD   =  -57.0
      PHILUD   =  -32.0

C    /RUNCTL/
      NANFD    =       0
      NENDED   =       0
      NHANFD   =       0
      NHENDED  =       0
      NEAAD    = 9999999
      NDEAD    = 9999999
      NHEAAD   =       0
      NHDEAD   =       0
      NFORAD   = 9999999
      NDFORD   = 9999999
      NHFORAD  = 9999999
      NHDFORD  = 9999999
      NDAAD    = 9999999
      NDDAD    = 9999999
      NHDAAD   = 9999999
      NHDDAD   = 9999999
      NTAAD    = 9999999
      NDTAD    = 9999999
      NHTAAD   =       0
      NHDTAD   =       0
      NHDMXND  =       0
      IADGBD   =       1
      IEDGBD   =    MOIE
      JADGBD   =       1
      JEDGBD   =    MOJE
      DTD      =     300.0
      DLANDD   =    7500.0
      DNOPRCD  =    7500.0
      NDRD     =      72
      NHDRD    =       0
      YADATD   = ' '
      LQWRD    =  .TRUE.
      LMOMOND  = .FALSE.
      LMOMITD  = .TRUE.
      LTAMITD  = .TRUE.
      LGMOND   = .TRUE.
      LSCEND   = .FALSE.
      LNOHALOD = .FALSE.

C    /DYNCTL/
      LSITSD   =  .TRUE.
      EPSASSD  =    0.15
      LVNESTD  = .FALSE.
      VBMXVD   =   10.0

C    /HYDCTL/
      LHYDROD  =  .TRUE.

C    /PHYCTL/
      HDRADD   =    1.0
      IEXCD    =    5
      LPHYD    =  .TRUE.
      LAKKUD   = .FALSE.
      LNEARD   = .FALSE.
      LVEGD    =  .TRUE.
      LINBOXD  =  .TRUE.
      LSICEDD  = .FALSE.
      LECDIUD  =  .TRUE.
      LECRADD  =  .TRUE.
      LECSOLD  = .FALSE.
      LECAERD  =  .TRUE.
      LECCFCD  = .FALSE.
      LECVDID  =  .TRUE.
      LECSURD  =  .TRUE.
      LECGWDD  = .FALSE.
      LECCOVD  =  .TRUE.
      LECCODD  =  .TRUE.
      LECGADD  =  .FALSE.
      LAEROZD  =  .FALSE.
      LCOLDCD  =  .TRUE.
      L5LAYD   =  .TRUE.
      LWDIFD   =  .TRUE.

C    /NMICTL/
      NVMD     =    3
      NITNMID  =    2
      DTNMID   =   45.0
      LANMID   =  .TRUE.
      LDNMID   = .FALSE.

C    /PRICTL/
      LDIAD    = .FALSE.

C    /DATEN/
      YADEND   = ' '
      YADCATD  = ' '
      YRDEND   = ' '
      YRDCATD  = ' '
      YEDEND   = ' '
      YEDCATD  = ' '
      YDDEND   = ' '
      YDDCATD  = ' '
      YFDEND   = ' '
      YFDCATD  = ' '
      YTDEND   = ' '
      YTDCATD  = ' '
      YMDCATD  = ' '
      YNDCATD  = ' '
      YUSERAD  = '001'
      YUSERED  = '001'
      YBDCATD  = ' '
      DO III=1,NZMXVB
         YBDNAMD(III) = ' '
      ENDDO
      YGDCATD  = ' '
      YGDNAMD  = ' '
      YHDCATD  = ' '
      YHDNAMD  = ' '
      YO3DCATD  = ' '
      YO3DNAMD  = ' '
      YSADCATD  = ' '
      YSADNAMD  = ' '
      YSNDCATD  = ' '
      YSNDNAMD  = ' '


C---- NAMELIST-VARIABLEN MIT DEN DEFAULTS BESETZEN ---------------------

C    /PARCTL/
      NPROCXM = NPROCXD
      NPROCYM = NPROCYD

C    /EMGRID/
      MOIE   = MOIED
      MOJE   = MOJED
      MOKE   = MOKED
      POLLAM = POLLAMD
      POLPHI = POLPHID
      DLAM   = DLAMD
      DPHI   = DPHID
      RLALU  = RLALUD
      PHILU  = PHILUD

C    /RUNCTL/
      NANF   = NANFD
      NENDE  = NENDED
      NHANF  = NHANFD
      NHENDE = NHENDED
      YADAT  = YADATD

      NEAA   = NEAAD
      NDEA   = NDEAD
      NHEAA  = NHEAAD
      NHDEA  = NHDEAD

      NFORA  = NFORAD
      NDFOR  = NDFORD
      NHFORA = NHFORAD
      NHDFOR = NHDFORD

      NDAA   = NDAAD
      NDDA   = NDDAD
      NHDAA  = NHDAAD
      NHDDA  = NHDDAD

      NTAA   = NTAAD
      NDTA   = NDTAD
      NHTAA  = NHTAAD
      NHDTA  = NHDTAD

      NHDMXN = NHDMXND

      DT     = DTD
      DLAND  = DLANDD
      DNOPRC = DNOPRCD
      LQWR   = LQWRD
      LMOMON = LMOMOND
      LMOMIT = LMOMITD
      LTAMIT = LTAMITD
      LGMON  = LGMOND
      LSCEN  = LSCEND
      LNOHALO= LNOHALOD
      NDR    = NDRD
      NHDR   = NHDRD

      IADGB  = IADGBD
      IEDGB  = IEDGBD
      JADGB  = JADGBD
      JEDGB  = JEDGBD

C    /DYNCTL/
      LSITS  = LSITSD
      EPSASS = EPSASSD
      LVNEST = LVNESTD
      VBMXV  = VBMXVD

C    /HYDCTL/
      LHYDRO  = LHYDROD

C    /PHYCTL/
      LPHY   = LPHYD
      LECDIU = LECDIUD
      LECRAD = LECRADD
      LECSOL = LECSOLD
      LECAER = LECAERD
      LECCFC = LECCFCD
      LECVDI = LECVDID
      LECSUR = LECSURD
      LECGWD = LECGWDD
      LECCOV = LECCOVD
      LECCOD = LECCODD
      LECGAD = LECGADD
      LAEROZ = LAEROZD
      L5LAY  = L5LAYD
      LWDIF  = LWDIFD
      LAKKU  = LAKKUD
      LNEAR  = LNEARD
      LVEG   = LVEGD
      LINBOX = LINBOXD
      LSICED = LSICEDD
      LCOLDC = LCOLDCD
      HDRAD  = HDRADD
      IEXC   = IEXCD

      DO I=1,12
         NDFAER(I) = NDFAERD(I)
      ENDDO

C    /NMICTL/
      LANMI  = LANMID
      LDNMI  = LDNMID
      NVM    = NVMD
      NITNMI = NITNMID
      DTNMI  = DTNMID

C    /PRICTL/
      LDIA    = LDIAD

C    /DATEN/
      YADEN   = YADEND
      YADCAT  = YADCATD

      YRDEN   = YRDEND
      YRDCAT  = YRDCATD

      YEDEN   = YEDEND
      YEDCAT  = YEDCATD

      YDDEN   = YDDEND
      YDDCAT  = YDDCATD

      YFDEN   = YFDEND
      YFDCAT  = YFDCATD

      YTDEN   = YTDEND
      YTDCAT  = YTDCATD

      YMDCAT  = YMDCATD
      YNDCAT  = YNDCATD

      YUSERA  = YUSERAD
      YUSERE  = YUSERED

      YBDCAT  = YBDCATD
      DO III=1,NZMXVB
         YBDNAM(III) = YBDNAMD(III)
      ENDDO

      YGDCAT  = YGDCATD
      YGDNAM  = YGDNAMD

      YHDCAT  = YHDCATD
      YHDNAM  = YHDNAMD

      YO3DCAT  = YO3DCATD
      YO3DNAM  = YO3DNAMD

      YSADCAT  = YSADCATD
      YSADNAM  = YSADNAMD

      YSNDCAT  = YSNDCATD
      YSNDNAM  = YSNDNAMD

C---  EINLESEN DER NAMELIST-KARTEN -------------------------------------
C    /PARCTL/
      READ (NUIN,PARCTL)

C    /EMGRID/
      READ (NUIN,EMGRID)

      NTASKS = 1

      MOIEJE   = MOIE*MOJE
      MOIEJEKE = MOIE*MOJE*MOKE
      MOIEKE   = MOIE*MOKE
      MOKE1    = MOKE+1

C    /RUNCTL/
      READ (NUIN,RUNCTL)

C     CHECK: ANGABE IN STUNDEN ODER IN ZEITSCHRITTEN ?

      IF (NHENDE.NE.NHENDED) THEN
         NANF  = MAX(0, NINT(NHANF *3600./DT) - 1)
         NENDE =        NINT(NHENDE*3600./DT)
      ENDIF
      IF (NHDEA.NE.NHDEAD) THEN
         NEAA  = NINT(NHEAA*3600./DT)
         NDEA  = NINT(NHDEA*3600./DT)
      ENDIF
      IF (NHFORA.NE.NHFORAD) THEN
         NFORA = NINT(NHFORA*3600./DT)
         NDFOR = NINT(NHDFOR*3600./DT)
      ENDIF
      IF (NHDDA.NE.NHDDAD) THEN
         NDAA  = NINT(NHDAA*3600./DT)
         NDDA  = NINT(NHDDA*3600./DT)
      ENDIF
      IF (NHDTA.NE.NHDTAD) THEN
         NTAA  = NINT(NHTAA*3600./DT)
         NDTA  = NINT(NHDTA*3600./DT)
      ENDIF
      IF (NHDR.NE.NHDRD) THEN
         NDR   = NINT(NHDR*3600./DT)
      ENDIF
      IF (NHDMXN.EQ.0) THEN
      IF (MYID .EQ. 0) THEN
         PRINT *,'FEHLER IN NAMELIST RUNCTL:'
         PRINT *,'NHDMXN=',NHDMXN
         PRINT *,'KEIN OUTPUT-INTERVALL GESETZT !'
      ENDIF
         STOP 'ERROR'
      ELSE
         NDMXN   = NINT(NHDMXN*3600./DT)
      ENDIF

C    /DYNCTL/
      READ (NUIN,DYNCTL)

C    /HYDCTL/
      READ (NUIN,HYDCTL)

C    /PHYCTL/
      READ (NUIN,PHYCTL)

C    /NMICTL/
      READ (NUIN,NMICTL)

C    /PRICTL/
      READ (NUIN,PRICTL)

C    /DATEN/
      READ (NUIN,DATEN)

C     LAENGE DER CATALOG-NAMEN FESTSTELLEN; EVT. "/" ERGAENZEN
      ILYDCAT = LENCH(YADCAT)
      IF( ILYDCAT.GT.0 ) THEN
          IF(YADCAT(ILYDCAT:ILYDCAT).NE.'/') THEN
              IF( ILYDCAT.LT.LEN(YADCAT) ) THEN
                  YADCAT = YADCAT(1:ILYDCAT)//'/'
              ELSE
                  CALL REMARK (' YADCAT ZU LANG; STOP ')
                  STOP
              ENDIF
          ENDIF
      ENDIF
      ILYDCAT = LENCH(YRDCAT)
      IF( ILYDCAT.GT.0 ) THEN
          IF(YRDCAT(ILYDCAT:ILYDCAT).NE.'/') THEN
              IF( ILYDCAT.LT.LEN(YRDCAT) ) THEN
                  YRDCAT = YRDCAT(1:ILYDCAT)//'/'
              ELSE
                  CALL REMARK (' YRDCAT ZU LANG; STOP ')
                  STOP
              ENDIF
          ENDIF
      ENDIF
      ILYDCAT = LENCH(YEDCAT)
      IF( ILYDCAT.GT.0 ) THEN
          IF(YEDCAT(ILYDCAT:ILYDCAT).NE.'/') THEN
              IF( ILYDCAT.LT.LEN(YEDCAT) ) THEN
                  YEDCAT = YEDCAT(1:ILYDCAT)//'/'
              ELSE
                  CALL REMARK (' YEDCAT ZU LANG; STOP ')
                  STOP
              ENDIF
          ENDIF
      ENDIF
      ILYDCAT = LENCH(YDDCAT)
      IF( ILYDCAT.GT.0 ) THEN
          IF(YDDCAT(ILYDCAT:ILYDCAT).NE.'/') THEN
              IF( ILYDCAT.LT.LEN(YDDCAT) ) THEN
                  YDDCAT = YDDCAT(1:ILYDCAT)//'/'
              ELSE
                  CALL REMARK (' YDDCAT ZU LANG; STOP ')
                  STOP
              ENDIF
          ENDIF
      ENDIF
      ILYDCAT = LENCH(YFDCAT)
      IF( ILYDCAT.GT.0 ) THEN
          IF(YFDCAT(ILYDCAT:ILYDCAT).NE.'/') THEN
              IF( ILYDCAT.LT.LEN(YADCAT) ) THEN
                  YFDCAT = YFDCAT(1:ILYDCAT)//'/'
              ELSE
                  CALL REMARK (' YFDCAT ZU LANG; STOP ')
                  STOP
              ENDIF
          ENDIF
      ENDIF
      ILYDCAT = LENCH(YTDCAT)
      IF( ILYDCAT.GT.0 ) THEN
          IF(YTDCAT(ILYDCAT:ILYDCAT).NE.'/') THEN
              IF( ILYDCAT.LT.LEN(YTDCAT) ) THEN
                  YTDCAT = YTDCAT(1:ILYDCAT)//'/'
              ELSE
                  CALL REMARK (' YTDCAT ZU LANG; STOP ')
                  STOP
              ENDIF
          ENDIF
      ENDIF

      ILYDCAT = LENCH(YMDCAT)
      IF( ILYDCAT.GT.0 ) THEN
          IF(YMDCAT(ILYDCAT:ILYDCAT).NE.'/') THEN
              IF( ILYDCAT.LT.LEN(YMDCAT) ) THEN
                  YMDCAT = YMDCAT(1:ILYDCAT)//'/'
              ELSE
                  CALL REMARK (' YMDCAT ZU LANG; STOP ')
                  STOP
              ENDIF
          ENDIF
      ENDIF

      ILYDCAT = LENCH(YNDCAT)
      IF( ILYDCAT.GT.0 ) THEN
          IF(YNDCAT(ILYDCAT:ILYDCAT).NE.'/') THEN
              IF( ILYDCAT.LT.LEN(YNDCAT) ) THEN
                  YNDCAT = YNDCAT(1:ILYDCAT)//'/'
              ELSE
                  CALL REMARK (' YNDCAT ZU LANG; STOP ')
                  STOP
              ENDIF
          ENDIF
      ENDIF

      ILYDCAT = LENCH(YHDCAT)
      IF( ILYDCAT.GT.0 ) THEN
          IF(YHDCAT(ILYDCAT:ILYDCAT).NE.'/') THEN
              IF( ILYDCAT.LT.LEN(YHDCAT) ) THEN
                  YHDCAT = YHDCAT(1:ILYDCAT)//'/'
              ELSE
                  CALL REMARK (' YHDCAT ZU LANG; STOP ')
                  STOP
              ENDIF
          ENDIF
      ENDIF
C
      IF (LSCEN) THEN
C
         IF (YGDCAT.EQ.YGDCATD) THEN
            CALL REMARK ('YGDCAT NICHT BELEGT; ERROR')
            STOP
         ELSE
            ILYDCAT = LENCH(YGDCAT)
            IF (ILYDCAT.GT.0) THEN
               IF (YGDCAT(ILYDCAT:ILYDCAT).NE.'/') THEN
                  IF (ILYDCAT.LT.LEN(YGDCAT)) THEN
                     YGDCAT = YGDCAT(1:ILYDCAT)//'/'
                  ELSE
                     CALL REMARK (' YGDCAT ZU LANG; STOP ')
                     STOP
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         IF (YGDNAM.EQ.YGDNAMD) THEN
            CALL REMARK ('YGDNAM NICHT BELEGT; ERROR')
            STOP
         ENDIF
C
      ENDIF
C
      IF (LAEROZ) THEN
C
         IF (YO3DCAT.EQ.YO3DCATD) THEN
            CALL REMARK ('YO3DCAT NICHT BELEGT; ERROR')
            STOP
         ELSE
            ILYDCAT = LENCH(YO3DCAT)
            IF (ILYDCAT.GT.0) THEN
               IF (YO3DCAT(ILYDCAT:ILYDCAT).NE.'/') THEN
                  IF (ILYDCAT.LT.LEN(YO3DCAT)) THEN
                     YO3DCAT = YO3DCAT(1:ILYDCAT)//'/'
                  ELSE
                     CALL REMARK (' YO3DCAT ZU LANG; STOP ')
                     STOP
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         IF (YO3DNAM.EQ.YO3DNAMD) THEN
            CALL REMARK ('YO3DNAM NICHT BELEGT; ERROR')
            STOP
         ENDIF
C
         IF (YSADCAT.EQ.YSADCATD) THEN
            CALL REMARK ('YSADCAT NICHT BELEGT; ERROR')
            STOP
         ELSE
            ILYDCAT = LENCH(YSADCAT)
            IF (ILYDCAT.GT.0) THEN
               IF (YSADCAT(ILYDCAT:ILYDCAT).NE.'/') THEN
                  IF (ILYDCAT.LT.LEN(YSADCAT)) THEN
                     YSADCAT = YSADCAT(1:ILYDCAT)//'/'
                  ELSE
                     CALL REMARK (' YSADCAT ZU LANG; STOP ')
                     STOP
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         IF (YSADNAM.EQ.YSADNAMD) THEN
            CALL REMARK ('YSADNAM NICHT BELEGT; ERROR')
            STOP
         ENDIF
C
         IF (YSNDCAT.EQ.YSNDCATD) THEN
            CALL REMARK ('YSNDCAT NICHT BELEGT; ERROR')
            STOP
         ELSE
            ILYDCAT = LENCH(YSNDCAT)
            IF (ILYDCAT.GT.0) THEN
               IF (YSNDCAT(ILYDCAT:ILYDCAT).NE.'/') THEN
                  IF (ILYDCAT.LT.LEN(YSNDCAT)) THEN
                     YSNDCAT = YSNDCAT(1:ILYDCAT)//'/'
                  ELSE
                     CALL REMARK (' YSNDCAT ZU LANG; STOP ')
                     STOP
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         IF (YSNDNAM.EQ.YSNDNAMD) THEN
            CALL REMARK ('YSNDNAM NICHT BELEGT; ERROR')
            STOP
         ENDIF
C
      ENDIF
C
      IF (LVEG) THEN
C
         IF (YBDCAT.EQ.YBDCATD) THEN
            CALL REMARK ('YBDCAT NICHT BELEGT; ERROR')
            STOP
         ELSE
            ILYDCAT = LENCH(YBDCAT)
            IF (ILYDCAT.GT.0) THEN
               IF (YBDCAT(ILYDCAT:ILYDCAT).NE.'/') THEN
                  IF (ILYDCAT.LT.LEN(YBDCAT)) THEN
                     YBDCAT = YBDCAT(1:ILYDCAT)//'/'
                  ELSE
                     CALL REMARK (' YBDCAT ZU LANG; STOP ')
                     STOP
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         IZAIA=0
         DO III=1,NZMXVB
            IF (YBDNAM(III).EQ.YBDNAMD(III))  THEN
               IZAIA=IZAIA+1
               WRITE(*,*) 'YBDNAM(',III,') NICHT BELEGT'
            ENDIF
         ENDDO
         IF (IZAIA.GT.0) THEN
            CALL REMARK ('YBDNAM NICHT BELEGT; ERROR')
            STOP
         ENDIF
C
      ENDIF

C     CHECK : SIND DIE NAMEN DER FELDER , DIE IN DIE D- UND T-DATEN
C     SOLLEN, BESETZT , SONST DEFAULT !
      NZVD2 = 0
      IZAIA = 0

      DO III = 1,NZMXVD
         IF ( YDVARN(III).EQ.'        ') THEN
            IZAIA = IZAIA+1
         ELSE
            NZVD2 = NZVD2+1
         ENDIF
      ENDDO
      IF (IZAIA.EQ.NZMXVD) THEN
         NZVD  = NZVDD
         DO III  = 1,NZVD
            YDVARN(III) = YDVARND(III)
         ENDDO
      ELSE
         NZVD  = NZVD2
      ENDIF

      NZVT2 = 0
      IZAIA = 0

      DO III  = 1,NZMXVT
         IF ( YTVARN(III).EQ.'        ') THEN
            IZAIA = IZAIA+1
         ELSE
            NZVT2 = NZVT2+1
         ENDIF
      ENDDO
      IF (IZAIA.EQ.NZMXVT) THEN
         NZVT  = NZVTD
         DO III  = 1,NZVT
            YTVARN(III) = YTVARND(III)
         ENDDO
      ELSE
         NZVT  = NZVT2
      ENDIF

      NZVM2 = 0
      IZAIA = 0

      DO III = 1,NZMXVM
         IF ( YMVARN(III).EQ.'        ') THEN
            IZAIA = IZAIA+1
         ELSE
            NZVM2 = NZVM2+1
         ENDIF
      ENDDO
      IF (IZAIA.EQ.NZMXVM) THEN
         NZVM  = NZVMD
         DO III  = 1,NZVM
            YMVARN(III) = YMVARND(III)
         ENDDO
      ELSE
         NZVM  = NZVM2
      ENDIF

      NZVN2 = 0
      IZAIA = 0

      DO III = 1,NZMXVN
         IF ( YNVARN(III).EQ.'        ') THEN
            IZAIA = IZAIA+1
         ELSE
            NZVN2 = NZVN2+1
         ENDIF
      ENDDO
      IF (IZAIA.EQ.NZMXVN) THEN
         NZVN  = NZVND
         DO III  = 1,NZVN
            YNVARN(III) = YNVARND(III)
         ENDDO
      ELSE
         NZVN  = NZVN2
      ENDIF

      IFEL=NZVM
      JFEL=NZVN
      INUMZRM=IFEL
      INUMZRN=JFEL

C      IED=IE
C      JED=JE
C      KED=KE
      NHORPHD=NHORPH
      IF (LDIA) THEN
      IF (MYID .EQ. 0) THEN
C---- AUSGABE AUF UNIT : NUAUFTR VORBEREITEN----------------------------

      YZEILE = ' '
      WRITE(NUAUFTR,'(A)') YZEILE
      WRITE(NUAUFTR,'(A1)') '1'
      WRITE(NUAUFTR,9903) '0',' DIE DATENKARTEN FUER DEN MODELL-LAUF',
     &                        ' ENTHALTEN FOLGENDE WERTE : '

C    /PARCTL/

      WRITE(NUAUFTR,9905) '0',' NAMELIST : PARCTL ',
     & ' ----------------- ','VARIABLENNAME','AKTUELLER WERT',
     & 'DEFAULT','FORMAT'
      WRITE(NUAUFTR,9910) 'NPROCXM',NPROCXM,NPROCXD,' I '
      IF ( NPROCXM.NE.NPROCXD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'NPROCYM',NPROCYM,NPROCYD,' I '
      IF ( NPROCYM.NE.NPROCYD ) WRITE (NUAUFTR,9920) '*'

C    /EMGRID/

      WRITE(NUAUFTR,9905) '0',' NAMELIST : EMGRID ',
     & ' ----------------- ','VARIABLENNAME','AKTUELLER WERT',
     & 'DEFAULT','FORMAT'
      WRITE(NUAUFTR,9910) 'MOIE   ',MOIE   ,MOIED  ,' I '
      IF ( MOIE    .NE.MOIED ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'MOJE   ',MOJE   ,MOJED  ,' I '
      IF ( MOJE    .NE.MOJED ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9912) 'POLLAM ',POLLAM ,POLLAMD,' R '
      IF ( POLLAM.NE.POLLAMD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9912) 'POLPHI ',POLPHI ,POLPHID,' R '
      IF ( POLPHI.NE.POLPHID ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9912) 'DLAM   ',DLAM   ,DLAMD  ,' R '
      IF ( DLAM  .NE.DLAMD   ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9912) 'DPHI   ',DPHI   ,DPHID  ,' R '
      IF ( DPHI  .NE.DPHID   ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9912) 'RLALU  ',RLALU  ,RLALUD ,' R '
      IF ( RLALU .NE.RLALUD  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9912) 'PHILU  ',PHILU  ,PHILUD ,' R '
      IF ( PHILU .NE.PHILUD  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'MOKE   ',MOKE   ,MOKED  ,' I '
      IF ( MOKE    .NE.MOKED ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'NHORPH ',NHORPH ,NHORPHD,' I '
      IF ( NHORPH.NE.NHORPHD ) WRITE (NUAUFTR,9920) '*'


C    /RUNCTL/

      WRITE(NUAUFTR,9905) '0',' NAMELIST : RUNCTL ',
     & ' ----------------- ','VARIABLENNAME','AKTUELLER WERT',
     & 'DEFAULT','FORMAT'
      WRITE(NUAUFTR,9910) 'NANF   ',NANF   ,NANFD  ,' I '
      IF ( NANF  .NE.NANFD   ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'NENDE  ',NENDE  ,NENDED ,' I '
      IF ( NENDE .NE.NENDED  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'NHANF  ',NHANF  ,NHANFD ,' I '
      IF ( NHANF .NE.NHANFD  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'NHENDE ',NHENDE ,NHENDED,' I '
      IF ( NHENDE.NE.NHENDED ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9914) 'YADAT ' ,YADAT  ,YADATD ,'C*8'
      IF ( YADAT .NE.YADATD  ) WRITE (NUAUFTR,9920) '*'

      WRITE(NUAUFTR,9910) 'NEAA   ',NEAA   ,NEAAD  ,' I '
      IF ( NEAA  .NE.NEAAD   ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'NDEA   ',NDEA   ,NDEAD  ,' I '
      IF ( NDEA  .NE.NDEAD   ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'NHEAA  ',NHEAA  ,NHEAAD ,' I '
      IF ( NHEAA .NE.NHEAAD  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'NHDEA  ',NHDEA  ,NHDEAD ,' I '
      IF ( NHDEA .NE.NHDEAD  ) WRITE (NUAUFTR,9920) '*'

      WRITE(NUAUFTR,9910) 'NFORA  ',NFORA  ,NFORAD ,' I '
      IF ( NFORA .NE.NFORAD  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'NDFOR  ',NDFOR  ,NDFORD ,' I '
      IF ( NDFOR .NE.NDFORD  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'NHFORA ',NHFORA ,NHFORAD,' I '
      IF ( NHFORA.NE.NHFORAD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'NHDFOR ',NHDFOR ,NHDFORD,' I '
      IF ( NHDFOR.NE.NHDFORD ) WRITE (NUAUFTR,9920) '*'

      WRITE(NUAUFTR,9910) 'NDAA   ',NDAA   ,NDAAD  ,' I '
      IF ( NDAA  .NE.NDAAD   ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'NDDA   ',NDDA   ,NDDAD  ,' I '
      IF ( NDDA  .NE.NDDAD   ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'NHDAA  ',NHDAA  ,NHDAAD ,' I '
      IF ( NHDAA .NE.NHDAAD  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'NHDDA  ',NHDDA  ,NHDDAD ,' I '
      IF ( NHDDA .NE.NHDDAD  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'IADGB  ',IADGB  ,IADGBD ,' I '
      IF ( IADGB .NE.IADGBD  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'IEDGB  ',IEDGB  ,IEDGBD ,' I '
      IF ( IEDGB .NE.IEDGBD  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'JADGB  ',JADGB  ,JADGBD ,' I '
      IF ( JADGB .NE.JADGBD  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'JEDGB  ',JEDGB  ,JEDGBD ,' I '
      IF ( JEDGB .NE.JEDGBD  ) WRITE (NUAUFTR,9920) '*'


      WRITE(NUAUFTR,9910) 'NTAA   ',NTAA   ,NTAAD  ,' I '
      IF ( NTAA  .NE.NTAAD   ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'NDTA   ',NDTA   ,NDTAD  ,' I '
      IF ( NDTA  .NE.NDTAD   ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'NHTAA  ',NHTAA  ,NHTAAD ,' I '
      IF ( NHTAA .NE.NHTAAD  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'NHDTA  ',NHDTA  ,NHDTAD ,' I '
      IF ( NHDTA .NE.NHDTAD  ) WRITE (NUAUFTR,9920) '*'

      WRITE(NUAUFTR,9910) 'NHDMXN ',NHDMXN ,NHDMXND,' I '
      IF ( NHDMXN.NE.NHDMXND ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9939) 'NDMXN  ',NDMXN  ,'        ',' I '

      WRITE(NUAUFTR,9910) 'NDR    ',NDR    ,NDRD   ,' I '
      IF ( NDR   .NE.NDRD    ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'NHDR   ',NHDR   ,NHDRD  ,' I '
      IF ( NHDR  .NE.NHDRD   ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9939) 'NTASKS ',NTASKS,'        ',' I '
      WRITE(NUAUFTR,9912) 'DT     ',DT     ,DTD    ,' R '
      IF ( DT    .NE.DTD     ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9912) 'DLAND  ',DLAND  ,DLANDD ,' R '
      IF ( DLAND .NE.DLANDD  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9912) 'DNOPRC ',DNOPRC ,DNOPRCD,' R '
      IF ( DNOPRC.NE.DNOPRCD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LQWR   ',LQWR   ,LQWRD  ,' L '
      IF ( LQWR  .NEQV.LQWRD   ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LMOMON ',LMOMON ,LMOMOND,' L '
      IF ( LMOMON.NEQV.LMOMOND ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LMOMIT ',LMOMIT ,LMOMITD,' L '
      IF ( LMOMIT.NEQV.LMOMITD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LTAMIT ',LTAMIT ,LTAMITD,' L '
      IF ( LTAMIT.NEQV.LTAMITD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LGMON  ',LGMON  ,LGMOND,' L '
      IF ( LGMON .NEQV.LGMOND  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LSCEN  ',LSCEN  ,LSCEND,' L '
      IF ( LSCEN .NEQV.LSCEND  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LNOHALO',LSCEN  ,LNOHALOD,' L '
      IF ( LNOHALO .NEQV.LNOHALOD  ) WRITE (NUAUFTR,9920) '*'

C    /DYNCTL/

      WRITE(NUAUFTR,9905) '0',' NAMELIST : DYNCTL ',
     & ' ----------------- ','VARIABLENNAME','AKTUELLER WERT',
     & 'DEFAULT','FORMAT'
      WRITE(NUAUFTR,9916) 'LSITS  ',LSITS  ,LSITSD ,' L '
      IF ( LSITS .NEQV.LSITSD  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9912) 'EPSASS ',EPSASS ,EPSASSD,' R '
      IF ( EPSASS.NE.EPSASSD   ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LVNEST ',LVNEST ,LVNESTD,' L '
      IF ( LVNEST.NEQV.LVNESTD ) WRITE (NUAUFTR,9920) '*'

C    /HYDCTL/

      WRITE(NUAUFTR,9905) '0',' NAMELIST : HYDCTL ',
     & ' ----------------- ','VARIABLENNAME','AKTUELLER WERT',
     & 'DEFAULT','FORMAT'
      WRITE(NUAUFTR,9916) 'LHYDRO  ',    LHYDRO ,LHYDROD,' L '
      IF ( LVNEST.NEQV.LHYDROD )      WRITE (NUAUFTR,9920) '*'

C    /PHYCTL/

      WRITE(NUAUFTR,9905) '1',' NAMELIST : PHYCTL ',
     & ' ----------------- ','VARIABLENNAME','AKTUELLER WERT',
     & 'DEFAULT','FORMAT'
      WRITE(NUAUFTR,9916) 'LPHY   ',LPHY   ,LPHYD  ,' L '
      IF ( LPHY  .NEQV.LPHYD   ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LECDIU ',LECDIU ,LECDIUD,' L '
      IF ( LECDIU.NEQV.LECDIUD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LECRAD ',LECRAD ,LECRADD,' L '
      IF ( LECRAD.NEQV.LECRADD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LECSOL ',LECSOL ,LECSOLD,' L '
      IF ( LECSOL.NEQV.LECSOLD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LECAER ',LECAER ,LECAERD,' L '
      IF ( LECAER.NEQV.LECAERD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LECCFC ',LECCFC ,LECCFCD,' L '
      IF ( LECCFC.NEQV.LECCFCD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LECVDI ',LECVDI ,LECVDID,' L '
      IF ( LECVDI.NEQV.LECVDID ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LECSUR ',LECSUR ,LECSURD,' L '
      IF ( LECSUR.NEQV.LECSURD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LECGWD ',LECGWD ,LECGWDD,' L '
      IF ( LECGWD.NEQV.LECGWDD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LECCOV ',LECCOV ,LECCOVD,' L '
      IF ( LECCOV.NEQV.LECCOVD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LECCOD ',LECCOD ,LECCODD,' L '
      IF ( LECCOD.NEQV.LECCODD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LECGAD ',LECGAD ,LECGADD,' L '
      IF ( LECGAD.NEQV.LECGADD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LAEROZ ',LAEROZ ,LAEROZD,' L '
      IF ( LAEROZ.NEQV.LAEROZD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'L5LAY  ',L5LAY  ,L5LAYD ,' L '
      IF ( L5LAY.NEQV.L5LAYD   ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LWDIF  ',LWDIF  ,LWDIFD ,' L '
      IF ( LWDIF.NEQV.LWDIFD   ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LAKKU  ',LAKKU  ,LAKKUD ,' L '
      IF ( LAKKU.NEQV.LAKKUD   ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LNEAR  ',LNEAR  ,LNEARD ,' L '
      IF ( LNEAR.NEQV.LNEARD   ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LVEG   ',LVEG   ,LVEGD  ,' L '
      IF ( LVEG .NEQV.LVEGD    ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LINBOX ',LINBOX ,LINBOXD,' L '
      IF ( LINBOX.NEQV.LINBOXD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LSICED ',LSICED ,LSICEDD,' L '
      IF ( LSICED.NEQV.LSICEDD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LCOLDC ',LCOLDC ,LCOLDCD,' L '
      IF ( LCOLDC.NEQV.LCOLDCD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9912) 'HDRAD  ',HDRAD  ,HDRADD ,' R '
      IF ( HDRAD .NE.HDRADD  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'IEXC   ',IEXC   ,IEXCD  ,' I '
      IF ( IEXC  .NE.IEXCD   ) WRITE (NUAUFTR,9920) '*'
      DO I=1,12
         WRITE(NUAUFTR,9910) 'NDFAER ',NDFAER(I) ,NDFAERD(I),' I '
         IF ( NDFAER(I).NE.NDFAERD(I) ) WRITE (NUAUFTR,9920) '*'
      ENDDO

C    /NMICTL/

      WRITE(NUAUFTR,9905) '0',' NAMELIST : NMICTL ',
     & ' ----------------- ','VARIABLENNAME','AKTUELLER WERT',
     & 'DEFAULT','FORMAT'
      WRITE(NUAUFTR,9916) 'LANMI  ',LANMI  ,LANMID ,' L '
      IF ( LANMI .NEQV.LANMID  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9916) 'LDNMI  ',LDNMI  ,LDNMID ,' L '
      IF ( LDNMI .NEQV.LDNMID  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9912) 'DTNMI  ',DTNMI  ,DTNMID ,' R '
      IF ( DTNMI .NE.DTNMID    ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'NVM    ',NVM    ,NVMD   ,' I '
      IF ( NVM   .NE.NVMD      ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9910) 'NITNMI ',NITNMI ,NITNMID,' I '
      IF ( NITNMI.NE.NITNMID   ) WRITE (NUAUFTR,9920) '*'

C    /PRICTL/

      WRITE(NUAUFTR,9905) '1',' NAMELIST : PRICTL ',
     & ' ----------------- ','VARIABLENNAME','AKTUELLER WERT',
     & 'DEFAULT','FORMAT'
      WRITE(NUAUFTR,9916) 'LDIA   ',LDIA   ,LDIAD  ,' L '
      IF ( LDIA  .NEQV.LDIAD   ) WRITE (NUAUFTR,9920) '*'

C    /DATEN/

      WRITE(NUAUFTR,9905) '0',' NAMELIST : DATEN  ',
     & ' ----------------- ','VARIABLENNAME','AKTUELLER WERT',
     & 'DEFAULT','FORMAT'
      WRITE(NUAUFTR,9914) 'YADEN  ',YADEN  ,YADEND ,'C* 3'
      IF ( YADEN .NE.YADEND  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9914) 'YADCAT ',YADCAT ,YADCATD,'C*56'
      IF ( YADCAT.NE.YADCATD ) WRITE (NUAUFTR,9920) '*'

      WRITE(NUAUFTR,9914) 'YRDEN  ',YRDEN  ,YRDEND ,'C* 3'
      IF ( YRDEN .NE.YRDEND  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9914) 'YRDCAT ',YRDCAT ,YRDCATD,'C*56'
      IF ( YRDCAT.NE.YRDCATD ) WRITE (NUAUFTR,9920) '*'

      WRITE(NUAUFTR,9914) 'YEDEN  ',YEDEN  ,YEDEND ,'C* 3'
      IF ( YEDEN .NE.YEDEND  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9914) 'YEDCAT ',YEDCAT ,YEDCATD,'C*56'
      IF ( YEDCAT.NE.YEDCATD ) WRITE (NUAUFTR,9920) '*'

      WRITE(NUAUFTR,9914) 'YDDEN  ',YDDEN  ,YDDEND ,'C* 3'
      IF ( YDDEN .NE.YDDEND  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9914) 'YDDCAT ',YDDCAT ,YDDCATD,'C*56'
      IF ( YDDCAT.NE.YDDCATD ) WRITE (NUAUFTR,9920) '*'

      WRITE(NUAUFTR,9914) 'YFDEN  ',YFDEN  ,YFDEND ,'C* 3'
      IF ( YFDEN .NE.YFDEND  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9914) 'YFDCAT ',YFDCAT ,YFDCATD,'C*56'
      IF ( YFDCAT.NE.YFDCATD ) WRITE (NUAUFTR,9920) '*'

      WRITE(NUAUFTR,9914) 'YTDEN  ',YTDEN  ,YTDEND ,'C* 3'
      IF ( YTDEN .NE.YTDEND  ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9914) 'YTDCAT ',YTDCAT ,YTDCATD,'C*56'
      IF ( YTDCAT.NE.YTDCATD ) WRITE (NUAUFTR,9920) '*'

      WRITE(NUAUFTR,9914) 'YMDCAT ',YMDCAT ,YMDCATD,'C*56'
      IF ( YMDCAT.NE.YMDCATD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9914) 'YNDCAT ',YNDCAT ,YNDCATD,'C*56'
      IF ( YNDCAT.NE.YNDCATD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9914) 'YHDCAT ',YHDCAT ,YHDCATD,'C*56'
      IF ( YHDCAT.NE.YHDCATD ) WRITE (NUAUFTR,9920) '*'

      WRITE(NUAUFTR,9914) 'YUSERA ',YUSERA ,YUSERAD,'C* 3'
      IF ( YUSERA.NE.YUSERAD ) WRITE (NUAUFTR,9920) '*'
      WRITE(NUAUFTR,9914) 'YUSERE ',YUSERE ,YUSERED,'C* 3'
      IF ( YUSERE.NE.YUSERED ) WRITE (NUAUFTR,9920) '*'

C     AUSGABE DER VARIABLEN FUER DIE D- UND T-DATEIEN
      WRITE(NUAUFTR,9903) '0','  VARIABLEN IN DER D-',
     &                        'DATEI : '

      ENDIF

      YZEILE = ' '
      DO III = 1,NZVD,10
         IF((NZVD-III).GE.10) THEN
            WRITE(YZEILE,'(10X,10A10)') YDVARN(III),YDVARN(III+1),
     &           YDVARN(III+2),YDVARN(III+3),YDVARN(III+4),
     &           YDVARN(III+5),YDVARN(III+6),YDVARN(III+7),
     &           YDVARN(III+8),YDVARN(III+9)
            IF (MYID .EQ. 0)WRITE(NUAUFTR,'(A)') YZEILE
         ELSE
            YZEILE = ' '
            IF((NZVD-III).NE.0) THEN
               DO II = 1,NZVD-III+1
                  WRITE(YZEILE(13+(II-1)*10:20+(II-1)*10),'(A8)')
     &                 YDVARN(III-1+II)
               ENDDO
               IF (MYID .EQ. 0) WRITE(NUAUFTR,'(A)') YZEILE
            ELSE
               IF (MYID .EQ. 0) WRITE(NUAUFTR,'(12X,A8)') YDVARN(NZVD)
            ENDIF
            EXIT
         ENDIF
      ENDDO

      IF (MYID .EQ. 0) WRITE(NUAUFTR,9903) '0','  VARIABLEN IN DER T-',
     &                        'DATEI : '

      YZEILE = ' '
      DO III = 1,NZVT,10
         IF((NZVT-III).GE.10) THEN
            WRITE(YZEILE,'(10X,10A10)') YTVARN(III),YTVARN(III+1),
     &           YTVARN(III+2),YTVARN(III+3),YTVARN(III+4),
     &           YTVARN(III+5),YTVARN(III+6),YTVARN(III+7),
     &           YTVARN(III+8),YTVARN(III+9)
            IF (MYID .EQ. 0) WRITE(NUAUFTR,'(A)') YZEILE
         ELSE
            YZEILE = ' '
            IF((NZVT-III).NE.0) THEN
               DO II = 1,NZVT-III+1
                  WRITE(YZEILE(13+(II-1)*10:20+(II-1)*10),'(A8)')
     &                 YTVARN(III-1+II)
               ENDDO
               IF (MYID .EQ. 0) WRITE(NUAUFTR,'(A)') YZEILE
            ELSE
               IF (MYID .EQ. 0) WRITE(NUAUFTR,'(12X,A8)') YTVARN(NZVT)
            ENDIF
            EXIT
         ENDIF
      ENDDO

      ENDIF

 9903 FORMAT(A1,T4,A,A)

 9905 FORMAT(A1,T4 ,A/T4 ,A//T6,A,T21,A,T82,A,T114,A)

 9910 FORMAT(T8,A,T21,I8,T82,I8,T115,A3)

 9912 FORMAT(T8,A,T21,E12.4,T82,E12.4,T115,A3)

 9914 FORMAT(T8,A,T21,A,T82,A,T115,A4)

 9916 FORMAT(T8,A,T21,L1,T82,L1,T115,A3)

 9920 FORMAT('+',T79,A1)

 9939 FORMAT(T8,A7,T20,I8,T82,A7,T115,A3)

      CLOSE(NUIN)

      RETURN
      END SUBROUTINE MEDEA
