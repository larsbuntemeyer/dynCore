      SUBROUTINE FILTER4(FIN, SP, FOUT, IDIM)
C
C**** FILTER4  -   AUFRUF EINES SYMETRISCHEN FILTERS DER LAENGE 4
C**   AUFRUF   :   CALL FILTER4(FIN, IE, SP, FOUT) IN FILTER
C**   ZWECK    :   GLAETTUNG EINES FELDES,DIE FILTERUNG WIRD ME-
C**                MAL DURCHGEFUEHRT.
C**   PARAMETER:   AUFRUF DER PARMETER ERFOLGT UEBER FILTER
C**   VERSIONS-
C**   DATUM    :   03.MAI 1989
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   V:WOHLGEMUTH/D.MAJEWSKI
C     FILTERN VON FIN MIT EINEM SYMMETRISCHEN FILTER DER LAENGE N=4
C     DIE FILTERUNG WIRD MEND-MAL DURCHGEFUEHRT
C
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: IDIM
      REAL,    INTENT(IN)    :: FIN(IDIM)
      REAL,    INTENT(INOUT) :: SP(IDIM),FOUT(IDIM)
C
C     Local Variables
C
      REAL    :: W(5)
      INTEGER :: MEND,M,KDUM,K,I
C
      W(1)=-0.00390625
      W(2)=0.03125
      W(3)=-0.109375
      W(4)=0.21875
      W(5)=0.7265625
C
      DO I = 1,IDIM
         FOUT(I) = FIN(I)
      ENDDO
C
      MEND = 50
      DO M = 1,MEND
C
         DO I = 6,IDIM - 5
            SP(I) = W(5)*FOUT(I) + W(4)*(FOUT(I-1)+FOUT(I+1)) +
     &                             W(3)*(FOUT(I-2)+FOUT(I+2)) +
     &                             W(2)*(FOUT(I-3)+FOUT(I+3)) +
     &                             W(1)*(FOUT(I-4)+FOUT(I+4))
         ENDDO
         DO I = 1,5
            SP(I) = SP(6)
            SP(IDIM-I+1) = SP(IDIM-5)
         ENDDO
         DO I = 6,IDIM-5
            FOUT(I) = W(5)*SP(I) +
     &                W(4)*(SP(I-1)+SP(I+1)) + W(3)*(SP(I-2)+SP(I+2)) +
     &                W(2)*(SP(I-3)+SP(I+3)) + W(1)*(SP(I-4)+SP(I+4))
         ENDDO
         DO I = 1,5
            FOUT(I) = FOUT(6)
            FOUT(IDIM-I+1) = FOUT(IDIM-5)
         ENDDO
C
      ENDDO
C
      DO I=1,5
         FOUT(I) = FIN(I)
         FOUT(IDIM-I+1) = FIN(IDIM-I+1)
      ENDDO
      KDUM=8
      DO K=1,3
         KDUM=KDUM-1
         DO I = 2,KDUM
            FOUT(I)       =0.25*FOUT(I-1)     +0.5*FOUT(I)       +
     &                     0.25*FOUT(I+1)
            FOUT(IDIM-I+1)=0.25*FOUT(IDIM-I+2)+0.5*FOUT(IDIM-I+1)+
     &                     0.25*FOUT(IDIM-I)
         ENDDO
      ENDDO
C
      RETURN
      END SUBROUTINE FILTER4
