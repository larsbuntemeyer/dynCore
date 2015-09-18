C
C@(#) BIBLIOTHEK EM: EUROPA-MODELL / DEUTSCHLAND-MODELL
C@(#) MODUL PROGCHK.H, V2.12 VOM 4/27/94, EXTRAHIERT AM 4/28/94
C
C     These fields are initiated by progec4
C
      REAL    ::  OM850M, OM500M, OM300M
C
      INTEGER ::  KFL850, KFL500, KFL300
C
      COMMON  / PROGCHK / OM850M, OM500M, OM300M, KFL850, KFL500, KFL300

C**** PROGCHK  -   CB:MITTLERE VERTIKALBEWEGUNG IN 3 NIVEAUS
C**
C**   BESCHREIBUNG DER VARIABLEN:
C**   KFL850   :   K-INDEX DER EM-HAUPTFLAECHE, DEREN DRUCK ETWA 850 HPA
C**                BETRAEGT
C**   KFL500   :   K-INDEX DER EM-HAUPTFLAECHE, DEREN DRUCK ETWA 500 HPA
C**                BETRAEGT
C**   KFL300   :   K-INDEX DER EM-HAUPTFLAECHE, DEREN DRUCK ETWA 300 HPA
C**                BETRAEGT
C**   OM850M   :   FLAECHENMITTEL DES ABSOLUTBETRAGES DER VERTIKALBEWE-
C**                GUNG 'OMEGA' IN DER FLAECHE K=KFL850
C**   OM500M   :   FLAECHENMITTEL DES ABSOLUTBETRAGES DER VERTIKALBEWE-
C**                GUNG 'OMEGA' IN DER FLAECHE K=KFL500
C**   OM300M   :   FLAECHENMITTEL DES ABSOLUTBETRAGES DER VERTIKALBEWE-
C**                GUNG 'OMEGA' IN DER FLAECHE K=KFL300
C**
C**   VERFASSER:   D.MAJEWSKI