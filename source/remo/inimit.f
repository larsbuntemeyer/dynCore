      SUBROUTINE INIMIT
      !
      IMPLICIT NONE
      !
      INCLUDE "parorg.h"
      INCLUDE "corg.h"
      INCLUDE "org.h"
      INCLUDE "sumfel.h"
      !
      ! Local Variables
      !
      INTEGER :: I,J,IJ,IW
      !
      DO J=1,IFEL
         DO I=1,MOIEJE
            SUM(I,J)=0.
            SUM2(I,J)=0.
         ENDDO
      ENDDO
      !
      DO J=1,JFEL
         DO I=1,MOIEJE
            SUMT(I,J)=0.
         ENDDO
      ENDDO
      !
      DO IW=1,4
         DO IJ=1,MOIEJE
            DSUMX(IJ,IW)=-1.E-34
            DSUMXT(IJ,IW)=-1.E-34
         ENDDO
      ENDDO
      !
      DO IW=1,2
         DO IJ=1,MOIEJE
            DSUMN(IJ,IW)=1.E34
            DSUMNT(IJ,IW)=1.E34
         ENDDO
      ENDDO
      !
      IFIRST=1
      JFIRST=1
      KFIRST=1
      !
      RETURN
      !
      END SUBROUTINE INIMIT
