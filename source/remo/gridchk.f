C
C     SUBROUTINE GRIDCHK 
C
C**** GRIDCHK  -   UP:GITTER AUS GDB MIT CB *PARAM*, *GRID* VERGLEICHEN
C**   AUFRUF   :   CALL GRIDCHK(IE,JE,KE,ISTAT)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   UBERPRUEFEN, OB DIE GEGEBENEN WERTE FUER IE, JE, KE,
C**                POLPHI, POLLAM, DPHI, DLAM, PHILU, RLALU MIT DEN
C**                ANGABEN IM GDB UEBEREINSTIMMEN
C**   VERSIONS-
C**   DATUM    :   08.02.89
C**
C**   EXTERNALS:   PRIV, PRRV
C**
C**   EINGABE-
C**   PARAMETER:   KEINE
C**   AUSGABE-
C**   PARAMETER:   ISTAT: STATUSPARAMETER; ISTAT=0 : ALLES O.K.
C**
C**   COMMON-
C**   BLOECKE  :   PARAM, GRID UND GRIB
C**
C**   METHODE  :   AUS GDB DIE WERTE FUER IE, JE, KE HOLEN UND MIT DEN
C**                WERTEN IM CB *PARAM* VERGLEICHEN
C**                AUS GDB DIE WERTE FUER POLPHI, POLLAM, DLAM, DPHI,
C**                PHILU, RLALU HOLEN UND MIT DEN WERTEN IM CB *GRID*
C**                VERGLEICHEN
C**   FEHLERBE-
C**   HANDLUNG :   FEHLERSCHLUESSEL SETZEN
C**   VERFASSER:   D.MAJEWSKI
C
      SUBROUTINE GRIDCHK (MOIE,MOJE,KE,ISTAT)
C
      IMPLICIT NONE
C
      INCLUDE "grid.h"
      INCLUDE "grib.h"
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN)  :: MOIE,MOJE,KE
      INTEGER, INTENT(OUT) :: ISTAT
C
C     Local Variables
C
      INTEGER :: IEGB,JEGB,KEGB,NDGB,MMGB,JJGB
      REAL    :: RLAROGB,RLALUGB,POLPHIG,POLLAMG,
     &           PHIROGB,PHILUGB,DPHIGB,DLAMGB
      REAL    :: FACRES,RESFAC
C
      IEGB  =  IGDB(5)
      JEGB  =  IGDB(6)
      KEGB  = (IGDB(1) - 42)/8 - 1

      ISTAT = 0

C     VERGLEICHEN MIT DEN WERTEN IN CB *PARAM*
      IF(IEGB.NE.MOIE) THEN
         CALL PRIV('GRIDCHK','MOIE  ', MOIE)
         CALL PRIV('GRIDCHK','IEGB', IEGB)
         ISTAT = 1
      ENDIF

      IF(JEGB.NE.MOJE) THEN
         CALL PRIV('GRIDCHK','MOJE  ', MOJE)
         CALL PRIV('GRIDCHK','JEGB', JEGB)
         ISTAT = 2
      ENDIF

      IF(KEGB.NE.KE) THEN

C     WENN DER FLAECHENTYP = 100 IST, SO HANDELT ES SICH UM AUF P-FLAE-
C     CHEN INTERPOLIERTE FELDER, DIE OHNE VERTIKALKOORDINATEN ABGELEGT
C     SIND; AUF NN REDUZIERTE FELDER HABEN DEN FLAECHENTYP 102
         IF(IPDB(8).EQ.100 .OR. IPDB(8).EQ.102) RETURN

C     WENN DAS DATUM D=000000, SO HANDELT ES SICH UM KONSTANTE BODEN-
C     FELDER, DIE OHNE VERTIKALKOORDINATEN ABGELEGT SIND
         JJGB = IPDB(11)
         MMGB = IPDB(12)
         NDGB = IPDB(13)
         IF(JJGB+MMGB+NDGB.GT.0) THEN
            CALL PRIV('GRIDCHK','KE  ', KE)
            CALL PRIV('GRIDCHK','KEGB', KEGB)
            ISTAT = 3
         ENDIF
      ENDIF
CKS
CKS   CHECK FOR RESOLUTION FACTOR
CKS
      IF ( IGDB(20).EQ.0) THEN
         IRESFAC = 1000
      ELSE
         IRESFAC = IGDB(20)
      ENDIF
      RESFAC = FLOAT( IRESFAC )
      FACRES = 1./ FLOAT( IRESFAC )
CKS
C     POLPHI, POLLAM, DLAM, DPHI, PHILU, RLALU AUS GDB HOLEN
      PHILUGB  =  IGDB( 7)*FACRES
      RLALUGB  =  IGDB( 8)*FACRES
      PHIROGB  =  IGDB(10)*FACRES
      RLAROGB  =  IGDB(11)*FACRES
      POLPHIG  = -IGDB(17)*FACRES
      POLLAMG  =  IGDB(18)*FACRES - 180.0
      DPHIGB   =  (PHIROGB-PHILUGB)/FLOAT(JEGB-1)
      DLAMGB   =  (RLAROGB-RLALUGB)/FLOAT(IEGB-1)

C     VERGLEICHEN MIT DEN WERTEN IN CB *GRID*
      IF(ABS(PHILUGB - PHILU).GT.1.0 E-5) THEN
         IF(IPDB(7).NE.132 .AND. IPDB(7).NE.181) THEN
            CALL PRRV('GRIDCHK','PHILU  ', PHILU  )
            CALL PRRV('GRIDCHK','PHILUGB', PHILUGB)
            ISTAT = 4
         ENDIF
      ENDIF

      IF(ABS(RLALUGB - RLALU).GT.1.0 E-5) THEN
         IF(IPDB(7).NE.131 .AND. IPDB(7).NE.180) THEN
            CALL PRRV('GRIDCHK','RLALU  ', RLALU  )
            CALL PRRV('GRIDCHK','RLALUGB', RLALUGB)
            ISTAT = 5
         ENDIF
      ENDIF

      IF(ABS(POLPHIG - POLPHI).GT.1.0 E-5) THEN
         CALL PRRV('GRIDCHK','POLPHI ', POLPHI )
         CALL PRRV('GRIDCHK','POLPHIG', POLPHIG)
         ISTAT = 6
      ENDIF

      IF(ABS(POLLAMG - POLLAM).GT.1.0 E-5) THEN
         CALL PRRV('GRIDCHK','POLLAM ', POLLAM )
         CALL PRRV('GRIDCHK','POLLAMG', POLLAMG)
         ISTAT = 7
      ENDIF

      IF(IABS(NINT(RESFAC*DPHIGB) - NINT(RESFAC*DPHI)).GT.0) THEN
         CALL PRRV('GRIDCHK','DPHI   ', DPHI   )
         CALL PRRV('GRIDCHK','DPHIGB ', DPHIGB )
         ISTAT = 8
      ENDIF

      IF(IABS(NINT(RESFAC*DLAMGB) - NINT(RESFAC*DLAM)).GT.0) THEN
         CALL PRRV('GRIDCHK','DLAM   ', DLAM   )
         CALL PRRV('GRIDCHK','DLAMGB ', DLAMGB )
         ISTAT = 9
      ENDIF

      RETURN
      END SUBROUTINE GRIDCHK