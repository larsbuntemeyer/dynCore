      SUBROUTINE PRTLOG ( YTEXT )
C
C@(#) BIBLIOTHEK EM: EUROPA-MODELL / DEUTSCHLAND-MODELL
C@(#) MODUL PRTLOG.F, V2.16 VOM 9/6/94, EXTRAHIERT AM 9/21/94
C
C**** PRTLOG   -   UP:SCHREIBEN VON TEXTEN IN DEN LOGFILE MIT REMARK
C**   AUFRUF   :   SUBROUTINE PRTLOG ( YTEXT )
C**   ENTRIES  :   KEINE
C**   ZWECK    :   SCHREIBEN VON TEXTEN IN DEN LOGFILE MIT REMARK,
C**                JEWEILS 70 ZEICHEN IN EINE ZEILE.
C**   VERSIONS-
C**   DATUM    :   20.06.94
C**
C**   EXTERNALS:   LENCH, REMARK
C**
C**   EINGABE-
C**   PARAMETER:   YTEXT:   TEXT (CHARACTER-STRING), DER GESCHRIEBEN
C**                         WERDEN SOLL
C**   AUSGABE-
C**   PARAMETER:   KEINE
C**
C**   COMMON-
C**   BLOECKE  :   KEINE
C**
C**   METHODE  :   AUFRUF VON REMARK
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   D.MAJEWSKI
C
      IMPLICIT NONE
C
      CHARACTER*(*), INTENT(IN) ::  YTEXT
C
      INTEGER, EXTERNAL :: LENCH
      INTEGER           :: IYA, IYE, IY, ITER, ILYTEXT
C
C     LAENGE VON YTEXT FESTSTELLEN, D.H. ERSTES BLANK VON RECHTS
C     SUCHEN
      ILYTEXT = LENCH ( YTEXT )

C     TEXT AUFBRECHEN IN ZEILEN VON 70 ZEICHEN LAENGE
      IF ( ILYTEXT.LE.70 ) THEN
         CALL REMARK ( YTEXT(1:ILYTEXT) )
      ELSE
         ITER = ILYTEXT/70 + 1
         DO IY = 1,ITER
            IYA  = (IY - 1)*70 + 1
            IYE  = MIN ( IYA + 69, ILYTEXT )
            CALL REMARK ( YTEXT(IYA:IYE) )
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE PRTLOG