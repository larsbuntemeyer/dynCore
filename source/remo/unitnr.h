C
      INTEGER ::         NUIN   , NUAUFTR,
     &                   NUADAT , NURDAT , NUDDAT , NUFDAT ,
     &                   NUTDAT , NUMDAT , NUSDAT , NUHDAT
C
      COMMON / UNITNR /  NUIN   , NUAUFTR,
     &                   NUADAT , NURDAT , NUDDAT , NUFDAT ,
     &                   NUTDAT , NUMDAT , NUSDAT , NUHDAT
C
C**** UNITNR   -   CB:NAMEN UND LOGICAL UNITS DER EM-DATEIEN
C****                 *UNITCH* UND *UNITNR* GEHOEREN ZUSAMMEN
C**   BESCHREIBUNG DER VARIABLEN:
C**
C**                LOGISCHE EINHEIT (UNIT-NUMMER) VON:
C**                FORMATTIERTEN DATEIEN:
C**   NUIN     :   INPUT-DATEI FUER NAMELIST-VARIABLEN
C**   NUAUFTR  :   DATEI FUER PROTOKOLLIERUNG DES AUFTRAGS
C**                DATEIEN
C**                LOGISCHE EINHEIT (UNIT-NUMMER) VON:
C**                UNFORMATTIERTEN DATEIEN:
C**   NUADAT   :   DATEI MIT DEN ANFANGSDATEN
C**   NURDAT   :   DATEI MIT DEN RANDDATEN
C**   NUDDAT   :   DATEI MIT DEN ERGEBNISDATEN (TEIL  -GEBIET)
C**   NUFDAT   :   DATEI MIT DEN FORTSETZUNGSDATEN
C**   NUTDAT   :   DATEI MIT DEN DATEN FUER TRAJEKTORIENBERECHNUNG
C**   NUMDAT   :   DATEI MIT DEN MONATSMITTELWERTEN
C**   NUMDAT   :   DATEI MIT STANDARDABWEICHUNGEN BZGL. EINES MONATS
C**   NUHDAT   :   DATEI MIT DEN FORTSETZUNGSDATEN (MITTELBILDUNG)
C**
C**   VERFASSER:   D.MAJEWSKI / R. PODZUN
