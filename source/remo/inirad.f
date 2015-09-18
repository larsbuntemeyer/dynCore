C      SUBROUTINE INIRAD(CETAH,CO2FAC,CVDAES,CVDAEL,
C     &                  CVDAEU,CVDAED,NLEVP1)
C
C*** *INIRAD*  PRESET AND MODIFY CONSTANTS IN RADIATION COMMON BLOCKS.
C
C      M.JARRAUD      E.C.M.W.F.     13/12/1982.
C      MODIFIED BY
C      R. PODZUN      DKRZ           30/01/1995
C
C     PURPOSE.
C     --------
C
C             THIS SUBROUTINE PRESET AND MODIFY CONSTANTS IN COMMON
C     BLOCKS USED IN THE PARAMETERISATION OF RADIATIVE PROCESSES.
C
C
C**   INTERFACE.
C     ----------
C
C             *INIRAD* IS CALLED FROM *INIT*.
C
C     EXTERNALS.
C     ----------
C            *AERDIS*    COMPUTE AEROSOL DISTRIBUTIONS.
C
      SUBROUTINE INIRAD(CETAH,CO2FAC,CVDAES,CVDAEL,
     &                  CVDAEU,CVDAED,NLEVP1)
C
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
      INCLUDE "org.h"
      INCLUDE "comecphy.h"
      INCLUDE "COMRAD2"
C
      INTEGER, INTENT(IN) :: NLEVP1
      REAL,    INTENT(IN) :: CO2FAC
      REAL,    INTENT(IN) :: CETAH(NLEVP1)
      REAL,    INTENT(INOUT) :: CVDAES(NLEVP1), 
     &                       CVDAEL(NLEVP1), CVDAEU(NLEVP1), 
     &                       CVDAED(NLEVP1)
C
C     ------------------------------------------------------------
C
C*        1.       PRESET CONSTANTS IN RADIATION
C                  ------ --------- -- ---------
C
C     CALL SURAD
      IF (LSCEN) THEN
         CALL SURADIS(CO2FAC)
      ELSE
         CALL SURADI(CO2FAC)
      ENDIF
      IF (LAEROZ) THEN
         CALL SUAERX
      ENDIF
C     CALL SUAER
C     CALL SUCFC
C     CALL SULW
      CALL SULWX
C     CALL SUSW
      CALL SUSWX
C
C     ------------------------------------------------------------
C
C*        2.   COMPUTE AEROSOL DISTRIBUTION.
C              ------- ------- ------------
C
      CALL ECAERDI(CETAH,CVDAES,CVDAEL,CVDAEU,CVDAED,NLEVP1,
     &             CTRBGA,CVOBGA,CSTBGA,CAEOPS,CAEOPL,CAEOPU,CAEOPD,
     &             CTRPT,CAEADK,CAEADM)
C
      RETURN
      END SUBROUTINE INIRAD
