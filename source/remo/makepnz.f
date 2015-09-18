      SUBROUTINE MAKEPNZ(NZEIT, NUN, YYTYP, ICODE, ISPLIT)
      !
      IMPLICIT NONE
      !
      !C
      !C**** MAKEPNZ  -   UP:DATEINAMEN FUER ZEITREIHEN-DATEIEN ERZEUGEN
      !C**
      !C**   AUFRUF   :   CALL MAKEPNZ(NZEIT, NCODE, NUN)
      !C**
      !C**   ENTRIES  :   KEINE
      !C**
      !C**   ZWECK    :   ERZEUGUNG DER DATEINAMEN FUER ZEITREIHEN-DATEIEN
      !C**                ABSPEICHERUNG DES NAMENS IN DEN COMMON *CORG*
      !C**   VERSIONS-
      !C**   DATUM    :   19.12.03
      !C**
      !C**   EXTERNALS:   DATUTC
      !C**
      !C**   EINGABE-
      !C**   PARAMETER:   NZEIT: ZEITSCHRITT
      !C**                NCODE: CODE-NUMMER DER ZEITREIHEN-DATEI
      !C**                NUN  : INDEX DES AKTUELLEN DATEINAMEN-FELDES
      !C**   AUSGABE-
      !C**   PARAMETER:   DATEINAME IM COMMON *CORG*
      !C**
      !C**   COMMON-
      !C**   BLOECKE  :   ORG, COMDYN
      !C**
      !C**   METHODE  :   -
      !C**
      !C**   FEHLERBE-
      !C**   HANDLUNG :   KEINE
      !C**
      !C**   VERFASSER:   R. PODZUN
      !C
      INCLUDE "corg.h"
      INCLUDE "comdyn.h"
      !
      ! Formal Parameters
      !
      INTEGER     , INTENT(IN) :: NZEIT, NUN, ICODE, ISPLIT
      CHARACTER(1), INTENT(IN) :: YYTYP
      !
      ! Local Variables
      !
      CHARACTER :: YAKDRD1*10, YAKDRD2*21, YMONDAT*6
      CHARACTER :: YCODE*3, YFILE*2, YTYP*1
      REAL      :: AKHHRD
      INTEGER   :: NAKTARD
      !
      !
      ! BERECHNUNG DES AKTUELLEN TERMINS
      !
      CALL DATUTC(NZEIT, YADAT, DT, YAKDRD1, YAKDRD2, NAKTARD,
     &            AKHHRD)
      !
      YMONDAT=YAKDRD1(1:6)
      !
      IF (YYTYP.EQ.'E') THEN
         !
         YTYP='e'
         !
         IF (ISPLIT.EQ.1) THEN
            !
            WRITE(YFILE,'(I2.2)') NUN
            YEDNAM(NUN)=
     &          'e'//YUSERE//YEDEN//YTYP//'_zrfile'//YFILE//'_'//YMONDAT
            !
         ELSE
            !
            WRITE(YCODE,'(I3.3)') ICODE
            YEDNAM(NUN)=
     &          'e'//YUSERE//YEDEN//YTYP//'_c'//YCODE//'_'//YMONDAT
            !
         ENDIF
         !
      ENDIF
      !
      IF (YYTYP.EQ.'N') THEN
         !
         YTYP='n'
         ! 
         IF (ISPLIT.EQ.1) THEN
            !
            WRITE(YFILE,'(I2.2)') NUN
            YNDNAM(NUN)=
     &          'e'//YUSERE//YEDEN//YTYP//'_zrfile'//YFILE//'_'//YMONDAT
            !
         ELSE
            !
            WRITE(YCODE,'(I3.3)') ICODE
            YNDNAM(NUN)=
     &          'e'//YUSERE//YEDEN//YTYP//'_c'//YCODE//'_'//YMONDAT
            ! 
         ENDIF
         !
      ENDIF
      !
      RETURN
      !
      END SUBROUTINE MAKEPNZ
