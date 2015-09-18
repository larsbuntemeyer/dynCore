      SUBROUTINE EXRANDASS2(ZZMY,
     &                      LSCHALT1,LSCHALT2,LSCHALT3,LSCHALT4,
     &                      LSCHALT5,LSCHALT6,LSCHALT7,LSCHALT8)
      IMPLICIT NONE
      INCLUDE "parorg.h"
C
C     Dummy Arguments
C
      REAL,    INTENT(INOUT) :: ZZMY(IE,JE,KE)

      LOGICAL, INTENT(IN)    :: LSCHALT1(JE*KE),LSCHALT2(JE*KE),
     &                          LSCHALT3(IE*KE),LSCHALT4(IE*KE)

      LOGICAL, INTENT(INOUT) :: LSCHALT5(KE),LSCHALT6(KE),
     &                          LSCHALT7(KE),LSCHALT8(KE)

      LOGICAL :: LSCHALT5TMP(KE),LSCHALT6TMP(KE),
     &           LSCHALT7TMP(KE),LSCHALT8TMP(KE)
C
C     Local Variables
C
      INTEGER :: NPP,K,JJ,J,IMESLEN,II,I
C
      TAGCOUNT    = 1

C****************WESTRAND*******************

      TAGTABLE(1) = 501

      IF (NPPMAX_SEND1 .GT. 0) THEN
         IMESLEN = 0
         DO K = 1,KE
            DO J = JAMIN1,JEMAX1
               IMESLEN = IMESLEN + 1
               IF (LSCHALT1(IMESLEN)) THEN
                  DO I = 1,IEMAX1
                     ZZMY(I,J,K) = 0.
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (NPPMAX_RECV1 .EQ. 1) THEN
         CALL PTEST
         CALL PRECVL(LSCHALT1)
         IMESLEN = 0
         DO K = 1,KE
            DO J = JAMIN1,JEMAX1
               IMESLEN = IMESLEN + 1
               IF (LSCHALT1(IMESLEN)) THEN
                  DO I = 1, IEMAX1
                     ZZMY(I,J,K) = 0.
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      CALL PSTOP

C*************OSTRAND******************

      TAGTABLE(1) = 502

      IF (NPPMAX_SEND2 .GT. 0) THEN
         IMESLEN = 0
         DO K = 1,KE
            DO J = JAMIN2,JEMAX2
               IMESLEN = IMESLEN + 1
               IF (LSCHALT2(IMESLEN)) THEN
                  DO I = IAMIN2,IE
                     ZZMY(I,J,K) = 0.
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (NPPMAX_RECV2 .EQ. 1) THEN
         CALL PTEST
         CALL PRECVL(LSCHALT2)
         IMESLEN = 0
         DO K = 1,KE
            DO J = JAMIN2,JEMAX2
               IMESLEN = IMESLEN + 1
               IF (LSCHALT2(IMESLEN)) THEN
                  DO I = IAMIN2,IE
                     ZZMY(I,J,K) = 0.
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      CALL PSTOP

C*************NORDRAND******************

      TAGTABLE(1) = 503

      IF (NPPMAX_SEND3 .GT. 0) THEN
         IMESLEN = 0
         DO K = 1,KE
            DO I = IAMIN3,IEMAX3
               IMESLEN = IMESLEN + 1
               IF (LSCHALT3(IMESLEN)) THEN
                  DO J = JAMIN3,JE
                     ZZMY(I,J,K) = 0.
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (NPPMAX_RECV3 .EQ. 1) THEN
         CALL PTEST
         CALL PRECVL(LSCHALT3)
         IMESLEN = 0
         DO K = 1,KE
            DO I = IAMIN3,IEMAX3
               IMESLEN = IMESLEN + 1
               IF (LSCHALT3(IMESLEN)) THEN
                  DO J = JAMIN3,JE
                     ZZMY(I,J,K) = 0.
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF


      CALL PSTOP
C*************SUEDRAND******************

      TAGTABLE(1) = 504

      IF (NPPMAX_SEND4 .GT. 0) THEN
         IMESLEN = 0
         DO K = 1,KE
            DO I = IAMIN4,IEMAX4
               IMESLEN = IMESLEN + 1
               IF (LSCHALT4(IMESLEN)) THEN
                  DO J = 1,JEMAX4
                     ZZMY(I,J,K) = 0.
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (NPPMAX_RECV4 .EQ. 1) THEN
         CALL PTEST
         CALL PRECVL(LSCHALT4)
         IMESLEN = 0
         DO K = 1,KE
            DO I = IAMIN4,IEMAX4
               IMESLEN = IMESLEN + 1
               IF (LSCHALT4(IMESLEN)) THEN
                  DO J = 1,JEMAX4
                     ZZMY(I,J,K) = 0.
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      CALL PSTOP

C***********SUED-KANTE UND WEST-KANTE DER SUEDWESTECKE *************

      TAGTABLE(1) = 505

      IF (NPPMAX_RECV5 .GT. 0) THEN
         DO NPP=1,NPPMAX_RECV5

            CALL PTEST
            CALL PRECVL(LSCHALT5TMP)
            DO K = 1,KE
               IF (LSCHALT5TMP(K)) THEN
                  LSCHALT5(K)=LSCHALT5TMP(K)
               ENDIF
            ENDDO

         ENDDO
      ENDIF

      DO K=1,KE
         IF (LSCHALT5(K)) THEN
            DO JJ=JAMIN5,JEMAX5
               DO II=IAMIN5,IEMAX5
                  ZZMY(II,JJ,K) = 0.
               ENDDO
            ENDDO
         ENDIF
      ENDDO

      CALL PSTOP

C***********SUED-KANTE UND OST-KANTE DER SUEDOSTECKE *************

      TAGTABLE(1) = 506

      IF (NPPMAX_RECV6 .GT. 0) THEN
         DO NPP=1,NPPMAX_RECV6

            CALL PTEST
            CALL PRECVL(LSCHALT6TMP)
            DO K = 1,KE
               IF (LSCHALT6TMP(K)) THEN
                  LSCHALT6(K)=LSCHALT6TMP(K)
               ENDIF
            ENDDO

         ENDDO
      ENDIF

      DO K=1,KE
         IF (LSCHALT6(K)) THEN
            DO JJ=JAMIN6,JEMAX6
               DO II=IAMIN6,IEMAX6
                  ZZMY(II,JJ,K) = 0.
               ENDDO
            ENDDO
         ENDIF
      ENDDO

      CALL PSTOP
C***********OST-KANTE UND NORD-KANTE DER NORDOSTECKE *************

      TAGTABLE(1) = 507

      IF (NPPMAX_RECV7 .GT. 0) THEN
         DO NPP=1,NPPMAX_RECV7

            CALL PTEST
            CALL PRECVL(LSCHALT7TMP)
            DO K = 1,KE
               IF (LSCHALT7TMP(K)) THEN
                  LSCHALT7(K)=LSCHALT7TMP(K)
               ENDIF
            ENDDO

         ENDDO
      ENDIF

      DO K=1,KE
         IF (LSCHALT7(K)) THEN
            DO JJ=JAMIN7,JEMAX7
               DO II=IAMIN7,IEMAX7
                  ZZMY(II,JJ,K) = 0.
               ENDDO
            ENDDO
         ENDIF
      ENDDO

      CALL PSTOP

C***********NORD-KANTE UND WEST-KANTE DER NORDWESTECKE *************

      TAGTABLE(1) = 508

      IF (NPPMAX_RECV8 .GT. 0) THEN
         DO NPP=1,NPPMAX_RECV8

            CALL PTEST
            CALL PRECVL(LSCHALT8TMP)
            DO K = 1,KE
               IF (LSCHALT8TMP(K)) THEN
                  LSCHALT8(K)=LSCHALT8TMP(K)
               ENDIF
            ENDDO

         ENDDO
      ENDIF

      DO K=1,KE
         IF (LSCHALT8(K)) THEN
            DO JJ=JAMIN8,JEMAX8
               DO II=IAMIN8,IEMAX8
                  ZZMY(II,JJ,K) = 0.
               ENDDO
            ENDDO
         ENDIF
      ENDDO

      CALL PWAIT

      RETURN
      END SUBROUTINE EXRANDASS2
