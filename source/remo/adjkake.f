      SUBROUTINE ADJKAKE(KE)
      !
      IMPLICIT NONE
      !
C
C**** ADJKAKE  -   UP:TABELLE KAKE AUS CB *EMGRIB* ANPASSEN, WENN KE
C****                 UNGLEICH 20 IST.
C**   AUFRUF   :   CALL ADJKAKE
C**   ENTRIES  :   KEINE
C**   ZWECK    :   TABELLE KAKE AUS CB *EMGRIB* ANPASSEN, WENN KE UN-
C**                GLEICH 20 IST.
C**   VERSIONS-
C**   DATUM    :   09.02.89
C**
C**   EXTERNALS:   KEINE
C**
C**   EINGABE-
C**   PARAMETER:   KEINE
C**   AUSGABE-
C**   PARAMETER:   KEINE
C**
C**   COMMON-
C**   BLOECKE  :   PARAM, EMGBCH, EMGBRI
C**
C**   METHODE  :   DIFFERENZ KD = KE - KAKE(1,2) BERECHNEN UND DIE
C**                TABELLE KAKE MIT KD KORRIGIEREN. KAKE WIRD IN
C**                CB *EMGBDT* BESETZT.
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   D.MAJEWSKI
      !
      ! Dummy argumtens
      !
      INTEGER, INTENT(IN) :: KE
      !
      ! Local Variables
      !
      INTEGER :: KD,NTAB
      !
      INCLUDE "emgbch.h"
      INCLUDE "emgbri.h"
      !
      KD = KE - 20
      !
      DO NTAB = 1,NFLDEM
         IF(KAKE(NTAB,1).EQ.KAKE(NTAB,2) .AND. KAKE(NTAB,1).GT.20) THEN
            KAKE(NTAB,1) = KAKE(NTAB,1) + KD
         ENDIF
         IF(KAKE(NTAB,2).GE.20) THEN
            KAKE(NTAB,2) = KAKE(NTAB,2) + KD
         ENDIF
      ENDDO
      !
      RETURN
      !
      END SUBROUTINE ADJKAKE
