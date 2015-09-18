C
C@(#) BIBLIOTHEK EM: EUROPA-MODELL / DEUTSCHLAND-MODELL
C@(#) MODUL PHYKON.H, V2.12 VOM 4/27/94, EXTRAHIERT AM 4/28/94
C
      REAL    ::          T0    , R     , RD   , WCP   , WLK  , WLF   ,
     &                    WLS   , G     , RERD , STAG  , RHF  , SIGMA ,
     &                    SOKO
C
      COMMON  / PHYKON /  T0    , R     , RD   , WCP   , WLK  , WLF   ,
     &                    WLS   , G     , RERD , STAG  , RHF  , SIGMA ,
     &                    SOKO
C**** PHYKON   -   CB:PHYSIKALISCHE KONSTANTEN
C**
C**   BESCHREIBUNG DER VARIABLEN:
C**   T0       :   NULLPUNKT DER TEMPERATUR               (K)
C**   R        :   GASKONSTANTE FUER TROCKENE LUFT        (J/(KG*K))
C**   RD       :   GASKONSTANTE FUER WASSERDAMPF          (J/(KG*K))
C**   WCP      :   SPEZ.WAERME  FUER TROCKENE LUFT        (J/(KG*K))
C**   WLK      :   VERDAMPFUNGSWAERME                     (J/KG)
C**   WLF      :   GEFRIERWAERME                          (J/KG)
C**   WLS      :   SUBLIMATIONSWAERME                     (J/KG)
C**   G        :   ERDBESCHLEUNIGUNG                      (M/S**2)
C**   RERD     :   MITTLERER ERDRADIUS                    (M)
C**   STAG     :   MITTLERER STERNTAG                     (S)
C**   RHF      :   DICHTE DES WASSERS                     (KG/M**3)
C**   SIGMA    :   BOLTZMANN-KONSTANTE                    (W/(M**2*K**4)
C**   SOKO     :   SOLARKONSTANTE                         (W/M**2)
C**
C**   VERFASSER:   D.MAJEWSKI
