      SUBROUTINE INIZR(NZTANF)
      !
      IMPLICIT NONE
      !
C
C**** INIZR    -   UP: VORBEREITUNG DES ABSPEICHERNS IN
C****                  ZEITREIHEN-DATEIEN
C**
C**   AUFRUF   :   CALL INIZR
C**
C**   ENTRIES  :   KEINE
C**
C**   ZWECK    :   VORBEREITUNG DES ABSPEICHERNS IN ZEITREIHEN-DATEIEN
C**
C**   VERSIONS-
C**   DATUM    :   19.12.03
C**
C**   EXTERNALS:   MAKEPNZ, SEND
C**
C**   EINGABE-
C**   PARAMETER:   KEINE
C**
C**   AUSGABE-
C**   PARAMETER:   UNITNUMERN IN COMMON *UNITZR*
C**
C**   COMMON-
C**   BLOECKE  :   ORG, CORG, EMGBCH, EMGBRI, SUMPAR, UNITZR
C**
C**   METHODE  :   -
C**
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**
C**   VERFASSER:   R. PODZUN
C
      INCLUDE "org.h"
      INCLUDE "corg.h"
      INCLUDE "emgbch.h"
      INCLUDE "emgbri.h"
      INCLUDE "unitzr.h"
      !
      ! Formal Parameters
      !
      INTEGER   ::  NZTANF
      !
      ! Local Variables
      !
      INTEGER      ::  ICODE
      INTEGER      ::  MSPLIT,NSPLIT
      INTEGER      ::  IUN,NTAB2
      CHARACTER(1) ::  YTYP
C
      ICODE=0
      MSPLIT=1
      NSPLIT=1
      IF (INUMZRM.EQ.IFEL) MSPLIT=0
      IF (INUMZRN.EQ.JFEL) NSPLIT=0
C
      DO IUN=1,INUMZRM
C
C     DATEINAME ERZEUGEN UND UNIT OEFFNEN
C
         IF (MSPLIT.EQ.0) THEN
C
            DO NTAB2 = 1,NFLDEM
               IF (YMVARN(IUN).EQ.YEMNAME(NTAB2)) THEN
                  ICODE=NEMGBNR(NTAB2)
                  EXIT
               ENDIF
            ENDDO
C
         ENDIF
C
         YTYP='E'
         NUEDAT(IUN)=IUN+100
         CALL MAKEPNZ(NZTANF, IUN, YTYP, ICODE, MSPLIT)
         CALL SEND(NUEDAT(IUN), YEDNAM(IUN), YEDCAT)
C
      ENDDO
C
      IF (.NOT. LTAMIT) RETURN
      ICODE=0
      DO IUN=1,INUMZRN
C
C     DATEINAME ERZEUGEN UND UNIT OEFFNEN
C
         IF (NSPLIT.EQ.0) THEN
C
            DO NTAB2 = 1,NFLDEM
               IF (YNVARN(IUN).EQ.YEMNAME(NTAB2)) THEN
                  ICODE=NEMGBNR(NTAB2)
                  EXIT
               ENDIF
            ENDDO
C
         ENDIF
C
         YTYP='N'
         NUNDAT(IUN)=IUN+300
         CALL MAKEPNZ(NZTANF, IUN, YTYP, ICODE, NSPLIT)
         CALL SEND(NUNDAT(IUN), YNDNAM(IUN), YNDCAT)
C
      ENDDO
C
      RETURN
      !
      END SUBROUTINE INIZR
