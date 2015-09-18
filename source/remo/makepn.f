      SUBROUTINE MAKEPN (YTYP, NZEIT)
      !
      IMPLICIT NONE
      !
      !
      !**** MAKEPN   -   UP:PERMANENTEN DATEINAMEN ERZEUGEN FUER REMO-DATEIEN
      !**   AUFRUF   :   CALL MAKEPN(YTYP, NZEIT)
      !**   ENTRIES  :   KEINE
      !**   ZWECK    :   ERZEUGUNG DER PERMANENTEN DATEINAMEN FUER REMO-DATEIEN
      !**                ABSPEICHERUNG DES NAMENS IN DEN COMMON *CORG*
      !**   VERSIONS-
      !**   DATUM    :   26.04.07
      !**
      !**   EXTERNALS:   DATUTC
      !**
      !**   EINGABE-
      !**   PARAMETER:   YTYP:  TYP DER REMO-DATEI: 'A', 'R', 'E', 'D', 'F' ,'T'
      !**                                           'H', 'M', 'S'
      !**                NZEIT: ZEITSCHRITT
      !**   AUSGABE-
      !**   PARAMETER:   DATEINAME IM COMMON *CORG*
      !**
      !**   COMMON-
      !**   BLOECKE  :   ORG, CORG, COMDYN, UNITCH
      !**
      !**   METHODE  :   AUS USER-, EXPERIMENTNUMMER, DATUM UND ZEITSCHRITT
      !**                DEN DATEI-NAMEN BILDEN
      !**   FEHLERBE-
      !**   HANDLUNG :   STOP IM FEHLERFALLE
      !**
      !**   VERFASSER:   R.PODZUN
      !
      ! Dummy arguments
      !
      CHARACTER, INTENT(IN) :: YTYP*(*)
      INTEGER  , INTENT(IN) :: NZEIT
      !
      ! Local variables
      !
      CHARACTER :: YAKDRD1*10, YAKDRD2*21,
     &             YCODE*1   , YNZT*4
      INTEGER   :: NAKTARD
      REAL      :: AKHHRD
      !
      INCLUDE "org.h"
      INCLUDE "corg.h"
      INCLUDE "comdyn.h"
      INCLUDE "unitch.h"
      !
CCM
C      WRITE (*,*) "YTYP:",YTYP
      YCODE='0'
      !
      WRITE(YNZT(1:4),'(I4.4)') NZT
      ! BELEGUNG DER CODENUMMERN
      IF (YTYP.EQ.'A'.OR.YTYP.EQ.'R') YCODE='a'
      IF (YTYP.EQ.'D') YCODE='d'
      IF (YTYP.EQ.'F') YCODE='f'
      IF (YTYP.EQ.'H') YCODE='h'
      IF (YTYP.EQ.'T') YCODE='t'
      IF (YTYP.EQ.'M') YCODE='m'
      IF (YTYP.EQ.'S') YCODE='s'
C      WRITE (*,*) "YCODE=",YCODE
      CALL DATUTC(NZEIT, YADAT, DT, YAKDRD1, YAKDRD2, NAKTARD,
     &            AKHHRD)
      !
      YADNAM = 'a'//YUSERA//YADEN//YCODE//YAKDRD1
      YRDNAM = 'a'//YUSERA//YRDEN//YCODE//YAKDRD1
      YDDNAM = 'e'//YUSERE//YDDEN//YCODE//YAKDRD1
      YFDNAM = 'e'//YUSERE//YFDEN//YCODE//YAKDRD1
      YHDNAM = 'e'//YUSERE//YEDEN//YCODE//YAKDRD1
      YTDNAM = 'e'//YUSERE//YTDEN//YCODE//YAKDRD1
      YMDNAM = 'e'//YUSERE//YEDEN//YCODE//YAKDRD1(1:6)
      YSDNAM = 'e'//YUSERE//YEDEN//YCODE//YAKDRD1(1:6)
C      WRITE (*,*) "MAKEPN: mark1"
C     FALSCHER TYP VON DATEN
      IF (YCODE.EQ.'0') THEN
         CALL PRCV('MAKEPN','YTYP',YTYP)
         CALL REMARK( 'MAKEPN:  FALSCHER DATEI-TYP (YTYP)' )
         STOP
      ENDIF
C      WRITE (*,*) "MAKEPN: end"
      RETURN
      END SUBROUTINE MAKEPN 
