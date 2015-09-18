      SUBROUTINE WRITEMF(NUHDAT)
      !
      IMPLICIT NONE
      !
      INCLUDE "corg.h"
      INCLUDE "org.h"
      INCLUDE "sumfel.h"
      !
      INTEGER, INTENT(IN) :: NUHDAT
      !
      WRITE(NUHDAT) SUM, SUM2, DSUMX, DSUMN, SUMT, DSUMXT, DSUMNT
      WRITE(NUHDAT) JPDB, JGDB, KPDB, KGDB
      WRITE(NUHDAT) ISJAHR, IJAHR  , ISMON  , IMON   ,
     &              ITAG  , ISTD   , IMONDAY, IAKTTER,
     &              IT    , ID     , IC     , IM     ,
     &              IN    , IFIRST ,
     &              JSJAHR, JSMON  , JSTAG  , JTAG   ,
     &              JSTD  , JAKTTER, JT     , JC     ,
     &              JM    , JN     , JFIRST , KFIRST ,
     &              JSTADAT
      !
      CLOSE(NUHDAT)
      !
      RETURN
      !
      END SUBROUTINE WRITEMF
