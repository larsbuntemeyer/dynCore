      INTEGER FUNCTION LENCH (YTEXT)
      !
      IMPLICIT NONE
      !
      !C
      !C@(#) BIBLIOTHEK EM: EUROPA-MODELL / DEUTSCHLAND-MODELL
      !C@(#) MODUL LENCH.F, V2.12 VOM 4/27/94, EXTRAHIERT AM 4/28/94
      !C
      !C**** LENCH    -   FC:TEXTLAENGE: POSITION DES ERSTEN BLANK VON RECHTS
      !C**   AUFRUF   :   ILEN = LENCH(YTEXT)
      !C**   ENTRIES  :   KEINE
      !C**   ZWECK    :   DIE LAENGE EINES TEXTSTRINGS BESTIMMEN. DIE LAENGE
      !C**                IST DEFINIERT ALS DIE ZAHL DER CHARACTER VON LINKS
      !C**                VOR DEM ERSTEN BLANK BZW. NULLBYTE
      !C**   VERSIONS-
      !C**   DATUM    :   07.02.89
      !C**
      !C**   EXTERNALS:   LEN
      !C**
      !C**   EINGABE-
      !C**   PARAMETER:   YTEXT: TEXTSTRING, DESSEN LAENGE BESTIMMT WERDEN SOLL
      !C**   AUSGABE-
      !C**   PARAMETER:   LAENGE DES TEXTES ALS WERT DER FUNCTION
      !C**
      !C**   COMMON-
      !C**   BLOECKE  :   KEINE
      !C**
      !C**   METHODE  :   ERSTES BLANK BZW. NULLBYTE VON RECHTS BESTIMMEN
      !C**   FEHLERBE-
      !C**   HANDLUNG :   KEINE
      !C**   VERFASSER:   D.MAJEWSKI
      !
      ! Formal Parameters
      !
      CHARACTER, INTENT(IN) :: YTEXT*(*)
      !
      LENCH = LEN(YTEXT)
      !
      DO
CKS      WHITESPACE FOUND!         
         IF(YTEXT(LENCH:LENCH).EQ.' ') THEN
            LENCH = LENCH - 1
            IF(LENCH.EQ.0) THEN
               EXIT
            ENDIF
CKS      ZERO BYTE FOUND!
         ELSE IF(YTEXT(LENCH:LENCH).EQ.CHAR(0)) THEN
            LENCH = LENCH - 1
            IF(LENCH.EQ.0) THEN
               EXIT
            ENDIF
CKS      CHARACTER LENGTH FOUND!
         ELSE
            EXIT
         ENDIF
      ENDDO
      !
      RETURN
      !
      END FUNCTION LENCH
