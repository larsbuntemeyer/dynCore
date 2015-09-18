      SUBROUTINE EXRANDASS1(U,V,NE,
     &                      LSCHALT1,LSCHALT2,LSCHALT3,LSCHALT4,
     &                      LSCHALT5,LSCHALT6,LSCHALT7,LSCHALT8)
      IMPLICIT NONE
      INCLUDE "parorg.h"
C
C     Dummy Arguments
C
      REAL,    INTENT(IN)    :: U(IEJE,KE,3),V(IEJE,KE,3)

      LOGICAL, INTENT(INOUT) :: LSCHALT1(JE*KE),LSCHALT2(JE*KE),
     &                          LSCHALT3(IE*KE),LSCHALT4(IE*KE)

      LOGICAL, INTENT(INOUT) :: LSCHALT5(KE),LSCHALT6(KE),
     &                          LSCHALT7(KE),LSCHALT8(KE)
C
C     Local Variables
C
      INTEGER :: NPP,M,K,J,IMESLEN,IJ1,I,NE
C
      TAGCOUNT    = 1

C****************WESTRAND*******************

      TAGTABLE(1) = 501

      DO M=1,JE*KE
         LSCHALT1(M)=.FALSE.
      ENDDO

      IF (NPPMAX_SEND1 .GT. 0) THEN
         IMESLEN = 0
         DO K = 1,KE
            DO J = JAMIN1,JEMAX1
               IMESLEN = IMESLEN + 1
               IJ1 = 2 + (J - 1)*IE
               IF (U(IJ1,K,NE) .LT. 0.0) THEN
                  LSCHALT1(IMESLEN) = .TRUE.
               ENDIF
            ENDDO
         ENDDO

         DO NPP = 2,NPPMAX_SEND1
            TYPE  = TAGTABLE(1)
            COUNT = IMESLEN*1
            DEST  = NPPINDEX1(NPP)
            CALL PSENDL(LSCHALT1)
         ENDDO

      ENDIF


C*************OSTRAND******************

      TAGTABLE(1) = 502

      DO M=1,JE*KE
         LSCHALT2(M)=.FALSE.
      ENDDO

      IF (NPPMAX_SEND2 .GT. 0) THEN
         IMESLEN = 0
         DO K = 1,KE
            DO J = JAMIN2,JEMAX2
               IMESLEN = IMESLEN + 1
               IJ1 = IE - 1 + (J - 1)*IE
               IF (U(IJ1,K,NE) .GT. 0.0) THEN
                  LSCHALT2(IMESLEN) = .TRUE.
               ENDIF
            ENDDO
         ENDDO

         DO NPP = 2,NPPMAX_SEND2
            TYPE  = TAGTABLE(1)
            COUNT = IMESLEN*1
            DEST  = NPPINDEX2(NPP)
            CALL PSENDL(LSCHALT2)
         ENDDO

      ENDIF


C*************NORDRAND******************

      TAGTABLE(1) = 503

      DO M=1,IE*KE
         LSCHALT3(M)=.FALSE.
      ENDDO

      IF (NPPMAX_SEND3 .GT. 0) THEN
         IMESLEN = 0
         DO K = 1,KE
            DO I = IAMIN3,IEMAX3
               IMESLEN = IMESLEN + 1
               IJ1 = I + (JE-1 - 1)*IE
               IF (V(IJ1,K,NE) .GT. 0.0) THEN
                  LSCHALT3(IMESLEN) = .TRUE.
               ENDIF
            ENDDO
         ENDDO

         DO NPP = 2,NPPMAX_SEND3
            TYPE  = TAGTABLE(1)
            COUNT = IMESLEN*1
            DEST  = NPPINDEX3(NPP)
            CALL PSENDL(LSCHALT3)
         ENDDO

      ENDIF


C*************SUEDRAND******************

      TAGTABLE(1) = 504

      DO M=1,IE*KE
         LSCHALT4(M)=.FALSE.
      ENDDO

      IF (NPPMAX_SEND4 .GT. 0) THEN
         IMESLEN = 0
         DO K = 1,KE
            DO I = IAMIN4,IEMAX4
               IMESLEN = IMESLEN + 1
               IJ1 = I + (2 - 1)*IE
               IF (V(IJ1,K,NE) .LT. 0.0) THEN
                  LSCHALT4(IMESLEN) = .TRUE.
               ENDIF
            ENDDO
         ENDDO

         DO NPP = 2,NPPMAX_SEND4
            TYPE  = TAGTABLE(1)
            COUNT = IMESLEN*1
            DEST  = NPPINDEX4(NPP)
            CALL PSENDL(LSCHALT4)
         ENDDO

      ENDIF


C***********SUED-KANTE UND WEST-KANTE DER SUEDWESTECKE *************

      TAGTABLE(1) = 505

      DO K=1,KE
         LSCHALT5(K)=.FALSE.
      ENDDO

      IF (NPPMAX_SEND5 .GT. 0) THEN
         DO K = 1,KE

            IF (IS5 .EQ. 1) THEN
               DO I = IAMIN5,IEMAX5
                  IJ1 = I + (2 - 1)*IE
                  IF (V(IJ1,K,NE) .LT. 0.0) THEN
                     LSCHALT5(K) = .TRUE.
                     EXIT
                  ENDIF
               ENDDO
            ENDIF

            IF ((IW5 .EQ. 1) .AND. (.NOT. LSCHALT5(K))) THEN
               DO J = JAMIN5,JEMAX5
                  IJ1 = 2 + (J - 1)*IE
                  IF (U(IJ1,K,NE) .LT. 0.0) THEN
                     LSCHALT5(K) = .TRUE.
                     EXIT
                  ENDIF
               ENDDO
            ENDIF

         ENDDO

C     NPPMAX_SEND5 IST U.U. GLEICH 1, D.H. SCHLEIFE WIRD DANN 
C     NICHT AUSGEFUEHRT.

         DO NPP = 2,NPPMAX_SEND5
            TYPE  = TAGTABLE(1)
            COUNT = KE*1
            DEST  = NPPINDEX5(NPP)
            CALL PSENDL(LSCHALT5)
         ENDDO

      ENDIF


C***********SUED-KANTE UND OST-KANTE DER SUEDOSTECKE *************

      TAGTABLE(1) = 506

      DO K=1,KE
         LSCHALT6(K)=.FALSE.
      ENDDO

      IF (NPPMAX_SEND6 .GT. 0) THEN
         DO K = 1,KE

            IF (IS6 .EQ. 1) THEN
               DO I = IAMIN6,IEMAX6
                  IJ1 = I + (2 - 1)*IE
                  IF (V(IJ1,K,NE) .LT. 0.0) THEN
                     LSCHALT6(K) = .TRUE.
                     EXIT
                  ENDIF
               ENDDO
            ENDIF

            IF ((IO6 .EQ. 1) .AND. (.NOT. LSCHALT6(K))) THEN
               DO J = JAMIN6,JEMAX6
                  IJ1 = IE - 1 + (J - 1)*IE
                  IF (U(IJ1,K,NE) .GT. 0.0) THEN
                     LSCHALT6(K) = .TRUE.
                     EXIT
                  ENDIF
               ENDDO
            ENDIF

         ENDDO


         DO NPP = 2,NPPMAX_SEND6
            TYPE  = TAGTABLE(1)
            COUNT = KE*1
            DEST  = NPPINDEX6(NPP)
            CALL PSENDL(LSCHALT6)
         ENDDO

      ENDIF


C***********OST-KANTE UND NORD-KANTE DER NORDOSTECKE *************

      TAGTABLE(1) = 507

      DO K=1,KE
         LSCHALT7(K)=.FALSE.
      ENDDO

      IF (NPPMAX_SEND7 .GT. 0) THEN
         DO K = 1,KE

            IF (IN7 .EQ. 1) THEN
               DO I = IAMIN7,IEMAX7
                  IJ1 = I + (JE - 1 - 1)*IE
                  IF (V(IJ1,K,NE) .GT. 0.0) THEN
                     LSCHALT7(K) = .TRUE.
                     EXIT
                  ENDIF
               ENDDO
            ENDIF

            IF ((IO7 .EQ. 1) .AND. (.NOT. LSCHALT7(K))) THEN
               DO J = JAMIN7,JEMAX7
                  IJ1 = IE - 1 + (J - 1)*IE
                  IF (U(IJ1,K,NE) .GT. 0.0) THEN
                     LSCHALT7(K) = .TRUE.
                     EXIT
                  ENDIF
               ENDDO
            ENDIF

         ENDDO


         DO NPP = 2,NPPMAX_SEND7
            TYPE  = TAGTABLE(1)
            COUNT = KE*1
            DEST  = NPPINDEX7(NPP)
            CALL PSENDL(LSCHALT7)
         ENDDO

      ENDIF


C***********NORD-KANTE UND WEST-KANTE DER NORDWESTECKE *************

      TAGTABLE(1) = 508

      DO K=1,KE
         LSCHALT8(K)=.FALSE.
      ENDDO

      IF (NPPMAX_SEND8 .GT. 0) THEN
         DO K = 1,KE

            IF (IN8 .EQ. 1) THEN
               DO I = IAMIN8,IEMAX8
                  IJ1 = I + (JE - 1 - 1)*IE
                  IF (V(IJ1,K,NE) .GT. 0.0) THEN
                     LSCHALT8(K) = .TRUE.
                     EXIT
                  ENDIF
               ENDDO
            ENDIF

            IF ((IW8 .EQ. 1) .AND. (.NOT. LSCHALT8(K))) THEN
               DO J = JAMIN8,JEMAX8
                  IJ1 = 2 + (J - 1)*IE
                  IF (U(IJ1,K,NE) .LT. 0.0) THEN
                     LSCHALT8(K) = .TRUE.
                     EXIT
                  ENDIF
               ENDDO
            ENDIF

         ENDDO


         DO NPP = 2,NPPMAX_SEND8
            TYPE  = TAGTABLE(1)
            COUNT = KE*1
            DEST  = NPPINDEX8(NPP)
            CALL PSENDL(LSCHALT8)
         ENDDO

      ENDIF

      CALL PWAIT

      RETURN
      END SUBROUTINE EXRANDASS1
