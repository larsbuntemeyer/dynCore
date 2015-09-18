      SUBROUTINE SCALGB ( FELDIN, IEJE, FAK, BIAS, FELDOUT )
C
C@(#) BIBLIOTHEK EM: EUROPA-MODELL / DEUTSCHLAND-MODELL
C@(#) MODUL SCALGB.F, V2.12 VOM 4/27/94, EXTRAHIERT AM 4/28/94
C
C**** SCALGB   -   UP:SKALIEREN DES FELDES 'FELDIN' MIT FAKTOR 'FAK' UND
C****                 BIAS 'BIAS' UND ABSPEICHERN AUF 'FELDOUT'
C**   AUFRUF   :   CALL SCALGB(FELDIN, IEJE, FAK, BIAS, FELDOUT)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   SKALIEREN DES FELDES 'FELDIN' MIT FAKTOR 'FAK' UND
C**                BIAS 'BIAS' UND ABSPEICHERN AUF 'FELDOUT'
C**   VERSIONS-
C**   DATUM    :   14.02.89
C**
C**   EXTERNALS:   KEINE
C**
C**   EINGABE-
C**   PARAMETER:   FELDIN:  ZU SKALIERENDES FELD
C**                IEJE:    DIMENSION VON FELDIN: FELDIN(IEJE)
C**                FAK:     SKALIERUNGSFAKTOR
C**                BIAS:    ADDITIVER BIAS
C**   AUSGABE-
C**   PARAMETER:   FELDOUT: SKALIERTES FELD
C**
C**   COMMON-
C**   BLOECKE  :   KEINE
C**
C**   METHODE  :   SKALIEREN MIT FAKTOR UND BIAS
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   D.MAJEWSKI
C
      IMPLICIT NONE
C      
      INTEGER, INTENT(IN)    :: IEJE
      REAL,    INTENT(IN)    :: FAK, BIAS
      REAL,    INTENT(IN)    :: FELDIN(IEJE)
      REAL,    INTENT(INOUT) :: FELDOUT(IEJE)
C
      INTEGER :: IJ
C
      DO IJ  = 1,IEJE
         FELDOUT(IJ) = (FELDIN(IJ) + BIAS)*FAK
      ENDDO
      DO IJ  = 1,IEJE
         IF (ABS(FELDOUT(IJ)).LT.1.E-36) FELDOUT(IJ)=0.0
      ENDDO

      RETURN
      END SUBROUTINE SCALGB