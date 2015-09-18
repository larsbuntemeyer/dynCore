C
      LOGICAL :: LECDIU, LECRAD, LECSOL, LECAER, LECCFC, LAEROZ, L5LAY
      LOGICAL :: LECVDI, LECSUR, LECGWD, LECCOV, LECCOD, LECGAD, LWDIF
C
      REAL    :: HDRAD
C
      INTEGER :: IEXC
C
      COMMON  / COMECPHY /LECDIU, LECRAD, LECSOL, LECAER, LECCFC,
     &                    LECVDI, LECSUR, LECGWD, LECCOV, LECCOD,
     &                    LECGAD, LAEROZ, L5LAY , LWDIF , HDRAD ,
     &                    IEXC
C
C**** COMECPHY -   CB:STEUERPARAMETER DER ECHAMPHYSIK
C**
C**   BESCHREIBUNG DER VARIABLEN:
C**   LECDIU   :   TRUE:  MIT TAGESGANG                       (--)
C**                FALSE: OHNE
C**   LECRAD   :   TRUE:  MIT STRAHLUNG                       (--)
C**                FALSE: OHNE
C**   LECSOL   :   TRUE:  MIT CLEAR SKY DIAGNOSTIK            (--)
C**                FALSE: OHNE
C**   LECAER   :   TRUE:  MIT ??                              (--)
C**                FALSE: OHNE
C**   LECCFC   :   TRUE:  MIT ??                              (--)
C**                FALSE: OHNE
C**   LECVDI   :   TRUE:  MIT VERTIKALDIFFUSION               (--)
C**                FALSE: OHNE
C**   LECSUR   :   TRUE:  MIT BODENMODELL                     (--)
C**                FALSE: OHNE
C**   LECGWD   :   TRUE:  MIT GRAVITY WAVE DRAG               (--)
C**                FALSE: OHNE
C**   LECCOV   :   TRUE:  MIT KONVEKTION                      (--)
C**                FALSE: OHNE
C**   LECCOD   :   TRUE:  MIT LARGE SCALE PRECIPITATION       (--)
C**                FALSE: OHNE
C**   LECGAD   :   TRUE:  GADS AEROSOLS DEPEND ON REL. HUMIDITY
C**                FALSE: NO DEPENDENCY ON REL. HUM. (80% FIX)
C**   LAEROZ   :   TRUE:  MIT EXTERNEN AEROSOLEN UND OZON     (--)
C**                FALSE: OHNE
C**   L5LAY    :   TRUE:  MIT 5 WASSERSCHICHTEN IM BODEN      (--)
C**                FALSE: OHNE (ORIGINAL ECHAM4)
C**   LWDIF    :   TRUE:  KEINE FEUCHTIGKEITSABHAENGIGKEIT    (--)
C**                       DER DIFFUSIVITAET UND WAERMEKAPA-
C**                       ZITAET IM BODEN
C**                       (REMO5.0=ORIGINAL ECHAM4)
C**                FALSE: MIT (NACH TIDO SEMMLER)
C**
C**   HDRAD    :   ALLE HDRAD: VOLLE STRAHLUNGSBERECHNUNG    (STD)
C**
C**   IEXC     :   SCHALTER FUER BODENFEUCHTEMODELL           (--)
C**                IEXC=1 AUFRUF STANDARD ARNOSCHEME
C**                IEXC=5 AUFRUF IMPROVED ARNOSCHEME
C**                IEXC=7 AUFRUF OLD REMO SIMULACRUM ARNOSCHEME
C**
C**   VERFASSER:   R.PODZUN
