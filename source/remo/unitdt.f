      BLOCK DATA UNITDT
      IMPLICIT NONE
      INCLUDE "unitch.h"
      INCLUDE "unitnr.h"
C
C**** UNITDT   -   DATA: BESETZEN DER DATEI-NAMEN UND LOGICAL UNITS
C****                    FUER DAS EM
C**   UNITDT   -   DATA: BESETZEN DER DATEI-NAMEN UND LOGICAL UNITS
C**                      FUER DAS EM
C**
C**   VERFASSER:   D.MAJEWSKI /  R.PODZUN

      DATA      YINPUT           / 'INPUT  '             /
      DATA      YUAUFTR          / 'YUAUFTR'             /
C
      DATA      NUIN             /     7                 /
      DATA      NUAUFTR          /     9                 /
      DATA      NUMDAT , NUSDAT  /    36     ,    38     /
      DATA      NUADAT , NURDAT  /    40     ,    42     /
      DATA      NUDDAT , NUHDAT  /    46     ,    44     /
      DATA      NUFDAT , NUTDAT  /    48     ,    50     /
      END BLOCK DATA UNITDT
