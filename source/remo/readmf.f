      SUBROUTINE READMF(NUHDAT)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: NUHDAT
C
      INCLUDE "corg.h"
      INCLUDE "org.h"
      INCLUDE "sumfel.h"
C
      READ(NUHDAT) SUM, SUM2, DSUMX, DSUMN, SUMT, DSUMXT, DSUMNT
      READ(NUHDAT) JPDB, JGDB, KPDB, KGDB
      READ(NUHDAT) ISJAHR, IJAHR  , ISMON  , IMON   ,
     &             ITAG  , ISTD   , IMONDAY, IAKTTER,
     &             IT    , ID     , IC     , IM     ,
     &             IN    , IFIRST ,
     &             JSJAHR, JSMON  , JSTAG  , JTAG   ,
     &             JSTD  , JAKTTER, JT     , JC     ,
     &             JM    , JN     , JFIRST , KFIRST ,
     &             JSTADAT
C
      CLOSE(NUHDAT)
C
      RETURN
      END SUBROUTINE READMF
