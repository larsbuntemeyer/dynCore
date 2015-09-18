      SUBROUTINE PLATTE(MYID)
      !
      IMPLICIT NONE
      !
      !@(#) BIBLIOTHEK EM: EUROPA-MODELL / DEUTSCHLAND-MODELL
      !@(#) MODUL PLATTE.F, V2.12 VOM 4/27/94, EXTRAHIERT AM 4/28/94
      !
      !**** PLATTE   -   UP: ZUORDNUNG UND EROEFFNUNG DER EM-DATEIEN
      !**   AUFRUF   :   CALL PLATTE IN HP *EMORG*
      !**   ENTRIES  :   KEINE
      !**   ZWECK    :   ALLEN IM EM VERWENDETEN FORMATTIERTEN DATEIEN WERDEN
      !**                DATEINAME UND UNIT-NUMMER ZUGEORDNET. DIE DATEIEN
      !**                WERDEN EROEFFNET.
      !**   VERSIONS-
      !**   DATUM    :   30.01.89
      !**                2007
      !**
      !**   EXTERNALS:   KEINE
      !**
      !**   EINGABE-
      !**   PARAMETER:   KEINE
      !**   AUSGABE-
      !**   PARAMETER:   KEINE
      !**
      !**   COMMON-
      !**   BLOECKE  :   UNITCH, UNITNR, UNITDT, ORG
      !**
      !**   METHODE  :   OPEN (STANDARD FORTRAN)
      !**   FEHLERBE-
      !**   HANDLUNG :   ABBRUCH BEI FEHLER
      !**   VERFASSER:   D.MAJEWSKI
      !
      ! not used INCLUDE "org.h"
      INCLUDE "unitch.h"
      INCLUDE "unitnr.h"
      !
      ! Formal Parameters
      !
      INTEGER, INTENT(IN) :: MYID
      !
      ! Local Variables
      !
      INTEGER             :: NIOSTAT
      !
      !
      OPEN(NUIN   , FILE=YINPUT, FORM='FORMATTED', STATUS='UNKNOWN',
     &     IOSTAT=NIOSTAT)
      IF (NIOSTAT.NE.0) THEN
         CALL REMARK( 'PLATTE: OPEN INPUT' )
         STOP
      ENDIF
      IF (MYID .GT. 0) RETURN
      !
      !
      OPEN(NUAUFTR, FILE=YUAUFTR, FORM='FORMATTED', STATUS='UNKNOWN',
     &     IOSTAT=NIOSTAT)
      !
      IF (NIOSTAT.NE.0) THEN
         CALL REMARK( 'PLATTE: OPEN YUAUFTR' )
         STOP
      ENDIF
      !
      REWIND NUAUFTR
      !
      RETURN
      !
      END SUBROUTINE PLATTE
