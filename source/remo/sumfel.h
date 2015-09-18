C
C     CB *SUMFEL* FUER DIE ZWISCHENZUSPEICHERENDEN
C     FELDER UND VARIABLEN BEI DER MONATS-TAGESMITTELBILDUNG
C     *SUMFEL* IMMER MIT *CORG* VERWENDEN
C
      INCLUDE "param.h"
C
      DOUBLE PRECISION SUM(IDIMSUM,NZMXVM), SUM2(IDIMSUM,NZMXVM),
     &                 SUMT(IDIMSUM,NZMXVN)
      REAL    ::        DSUMX (IDIMSUM,4), DSUMN (IDIMSUM,2)
      REAL    ::        DSUMXT(IDIMSUM,4), DSUMNT(IDIMSUM,2)
      INTEGER ::        JPDB(37,NZMXVM), JGDB(22,NZMXVN)
      INTEGER ::        KPDB(37,NZMXVM), KGDB(22,NZMXVN)
!      
      INTEGER ::        ISJAHR, IJAHR  , ISMON  , IMON   ,
     &                  ITAG  , ISTD   , IMONDAY, IAKTTER,
     &                  IT    , ID     , IC     , IM     ,
     &                  IN    , IFIRST ,
     &                  JSJAHR, JSMON  , JSTAG  , JTAG   ,
     &                  JSTD  , JAKTTER, JT     , JC     ,
     &                  JM    , JN     , JFIRST , KFIRST ,
     &                  JSTADAT
!
      COMMON / SUMFEL / SUM, SUM2, DSUMX, DSUMN, SUMT, DSUMXT, DSUMNT
!
      COMMON / INTFEL / JPDB, JGDB, KPDB, KGDB
!
      COMMON / INTVAR / ISJAHR, IJAHR  , ISMON  , IMON   ,
     &                  ITAG  , ISTD   , IMONDAY, IAKTTER,
     &                  IT    , ID     , IC     , IM     ,
     &                  IN    , IFIRST ,
     &                  JSJAHR, JSMON  , JSTAG  , JTAG   ,
     &                  JSTD  , JAKTTER, JT     , JC     ,
     &                  JM    , JN     , JFIRST , KFIRST ,
     &                  JSTADAT
!  
! DOUBLE PRECISION SUM(IDIMSUM,NZMXVM), SUM2(IDIMSUM,NZMXVM),
!&                 SUMT(IDIMSUM,NZMXVN)
! DIMENSION DSUMX (IDIMSUM,4), DSUMN (IDIMSUM,2)
! DIMENSION DSUMXT(IDIMSUM,4), DSUMNT(IDIMSUM,2)
! DIMENSION JPDB(37,NZMXVM), JGDB(22,NZMXVN)
! DIMENSION KPDB(37,NZMXVM), KGDB(22,NZMXVN)
! COMMON / SUMFEL / SUM, SUM2, DSUMX, DSUMN, SUMT, DSUMXT, DSUMNT
! COMMON / INTFEL / JPDB, JGDB, KPDB, KGDB
! COMMON / INTVAR / ISJAHR, IJAHR  , ISMON  , IMON   ,
!&                  ITAG  , ISTD   , IMONDAY, IAKTTER,
!&                  IT    , ID     , IC     , IM     ,
!&                  IN    , IFIRST ,
!&                  JSJAHR, JSMON  , JSTAG  , JTAG   ,
!&                  JSTD  , JAKTTER, JT     , JC     ,
!&                  JM    , JN     , JFIRST , KFIRST ,
!&                  JSTADAT
C
C     VERFASSER: R. PODZUN
C
