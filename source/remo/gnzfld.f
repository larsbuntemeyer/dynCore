C
C     SUBROUTINE GNZFLD
C
C**** GNZFLD   -   UP:ANZAHL DER ZU LESENDEN A- ODER R-FELDER FESTLEGEN
C**   AUFRUF   :   CALL GNZFLD(YTYP, NZGFELD)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   ANZAHL DER ZU LESENDEN A- ODER R-FELDER FESTLEGEN
C**   VERSIONS-
C**   DATUM    :   09.02.89
C**
C**   EXTERNALS:   KEINE
C**
C**   EINGABE-
C**   PARAMETER:   YTYP: TYP DER DATEN ('A' ODER 'R')
C**   AUSGABE-
C**   PARAMETER:   NZGFELD: ANZAHL DER ZU LESENDEN FELDER
C**
C**   COMMON-
C**   BLOECKE  :   ORG, CORG, EMGBCH, EMGBRI
C**
C**   METHODE  :   JEDES ZU LESENDE BODENFELD ERHOEHT NZGFELD UM 1,
C**                JEDES ATMOSPHAERENFELD UM KAKE(NTAB,2)-KAKE(NTAB,1)+1
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   D.MAJEWSKI
C
C
      SUBROUTINE GNZFLD ( YTYP, NZGFELD )
C
      IMPLICIT NONE
C
      INCLUDE "org.h"
      INCLUDE "corg.h"
      INCLUDE "emgbch.h"
      INCLUDE "emgbri.h"
C
      CHARACTER, INTENT(IN)  ::    YTYP*(*)
      INTEGER,   INTENT(OUT) ::    NZGFELD
C      
      INTEGER :: NTAB1,NTAB2
      
      NZGFELD = 0

      IF(YTYP.EQ.'A') THEN
         OUTERA: DO NTAB1 = 1,NZVA
            DO NTAB2 = 1,NFLDEM
               IF(YAVARN(NTAB1).EQ.YEMNAME(NTAB2)) THEN
                  NZGFELD = NZGFELD + KAKE(NTAB2,2) - KAKE(NTAB2,1) + 1
                  CYCLE OUTERA
               ENDIF
            ENDDO
         ENDDO OUTERA

      ELSE IF(YTYP.EQ.'R') THEN
         OUTERR: DO NTAB1 = 1,NZVR
            DO NTAB2 = 1,NFLDEM
               IF(YRVARN(NTAB1).EQ.YEMNAME(NTAB2)) THEN
                  NZGFELD = NZGFELD + KAKE(NTAB2,2) - KAKE(NTAB2,1) + 1
                  CYCLE OUTERR
               ENDIF
            ENDDO
         ENDDO OUTERR

      ELSE
         CALL REMARK( 'GNZFLD: FALSCHER TYP (YTYP)' )
         STOP
      ENDIF

      RETURN
      END SUBROUTINE GNZFLD