      SUBROUTINE EXRANDASSSETUP
      !
      IMPLICIT NONE
      !
      INCLUDE "parorg.h"

C     ES WIRD VORAUSGESETZT, DASS ES MINDESTENS (2*8+1)=17 PUNKTE IN I-RICHTUNG
C     UND MINDESTENS (2*8+1)=17 PUNKTE IN J-RICHTUNG GIBT.
C     (INSBESONDERE MUSS ALSO Z.B. GESICHERT SEIN, DASS I(OSTRAND) > 8) 

C     WEITERE VORAUSSETZUNG: KEIN PROZESSOR HAT NUR RANDPUNKTE

C     NPPINDEX(1) = MYID;
C     ES WIRD AN (NPPMAX_SEND - 1) PROZESSOREN GESENDET (NICHT AN SICH SELBST).

C     WESTRAND: NPPINDEX, NPPMAX_SEND UND NPPMAX_RECV MUESSTEN EIGENTLICH NUR
C     FUER EIN J BERECHNET WERDEN; DIE SCHLEIFE UEBER ALLE J IST ABER WEGEN DER
C     BESTIMMUNG VON JAMIN UND JEMAX NOTWENDIG. ANALOG ANDERE RAENDER.

C     ZU BEACHTEN: ISUBPOS ENTHAELT KEINE RANDPUNKTE

C     ES SIND NUR REGULAERE GEBIETSZERLEGUNGEN ZUGELASSEN.
      !
      ! Local Variables
      !
      INTEGER :: I,J,M,IG,JG,NP,NPP
      INTEGER :: IANZP5,IANZP6,IANZP7,IANZP8
      INTEGER :: IMESLEN,IEXIT
      !
      !
      !
      NP = MYID + 1

C ********************* WESTRAND******************

      DO M=1,MAXPRC
         NPPINDEX1(M) = -1
      ENDDO

      NPPMAX_SEND1 = 0
      NPPMAX_RECV1 = 0
      JAMIN1 = JE + 1
      JEMAX1 = 0


      DO JG = ISUBPOS(NP,2) , ISUBPOS(NP,4)
         IMESLEN = 1
         IF ((JG .GE. 9) .AND. (JG .LE. (MOJE-8))) THEN
            J = JG - MYPOSGRD(2) + 3

            IF (NEIGHBOR(1) .EQ. -1) THEN

C     EIGENE DATEN
               IF (J .LT. JAMIN1) JAMIN1 = J
               IF (J .GT. JEMAX1) JEMAX1 = J
               NPPMAX_SEND1 = 1
               NPPINDEX1(IMESLEN) = NP - 1

C     FUER ANDERE PROZESSOREN
               DO NPP = 1,NPROC
                  IF (NPP .NE. NP) THEN
                     IF ((JG .GE. ISUBPOS(NPP,2)) .AND.
     &                   (JG .LE. ISUBPOS(NPP,4))) THEN

C     EIGENTLICH IG=2,8 AUSREICHEND
                        DO IG = 1,8
                           IF((IG .GE. ISUBPOS(NPP,1)) .AND.
     &                        (IG .LE. ISUBPOS(NPP,3))) THEN
                              IMESLEN = IMESLEN + 1
                              NPPINDEX1(IMESLEN) = NPP - 1
                              NPPMAX_SEND1   = NPPMAX_SEND1 + 1
                              EXIT
                           ENDIF
                        ENDDO

                     ENDIF
                  ENDIF
               ENDDO

            ELSE 
               
C     EIGENTLICH IG=2,8 AUSREICHEND
               DO IG = 1,8
                  IF((IG .GE. ISUBPOS(NP,1)) .AND.
     &               (IG .LE. ISUBPOS(NP,3))) THEN

C     JAMIN1 UND JEMAX1 DUERFEN NUR GESETZT WERDEN, WENN DER JEWEILIGE
C     PROZESSOR EINEN DER SPEZIFIZIERTEN JG-WERTE (JG-SCHLEIFE) UND EINEN
C     DER SPEZIFIZIERTEN IG-WERTE (IG-SCHLEIFE) BESITZT, D.H. DIE BEIDEN
C     WERTE DUERFEN NUR INNERHALB DER IG-SCHLEIFE GESETZT WERDEN.

                     IF (J .LT. JAMIN1) JAMIN1 = J
                     IF (J .GT. JEMAX1) JEMAX1 = J
                     NPPMAX_RECV1 = 1
                     EXIT
                  ENDIF
               ENDDO
                     
            ENDIF

         ENDIF
      ENDDO



      IEMAX1 = 0
C     IAMIN1 WIRD EIGENTLICH NICHT BENOETIGT
      IAMIN1 = 1

      IF (NEIGHBOR(1) .EQ. -1) THEN
         IEMAX1 = 2
      ENDIF
      
C     RANDWERTE WERDEN VON ISUBPOS NICHT ERFASST, ALSO IG=2,8
      DO IG = 2,8
         IF((IG .GE. ISUBPOS(NP,1)) .AND. (IG .LE. ISUBPOS(NP,3))) THEN
            I = IG - MYPOSGRD(1) + 3
            IEMAX1 = I
         ENDIF
      ENDDO


C ********************* OSTRAND******************

      DO M=1,MAXPRC
         NPPINDEX2(M) = -1
      ENDDO

      NPPMAX_SEND2 = 0
      NPPMAX_RECV2 = 0
      JAMIN2 = JE + 1
      JEMAX2 = 0


      DO JG = ISUBPOS(NP,2) , ISUBPOS(NP,4)
         IMESLEN = 1
         IF ((JG .GE. 9) .AND. (JG .LE. (MOJE-8))) THEN
            J = JG - MYPOSGRD(2) + 3

            IF (NEIGHBOR(3) .EQ. -1) THEN

C     EIGENE DATEN
               IF (J .LT. JAMIN2) JAMIN2 = J
               IF (J .GT. JEMAX2) JEMAX2 = J
               NPPMAX_SEND2 = 1
               NPPINDEX2(IMESLEN) = NP - 1

C     FUER ANDERE PROZESSOREN
               DO NPP = 1,NPROC
                  IF (NPP .NE. NP) THEN
                     IF ((JG .GE. ISUBPOS(NPP,2)) .AND.
     &                   (JG .LE. ISUBPOS(NPP,4))) THEN

C     EIGENTLICH IG=MOIE-7,MOIE-1 AUSREICHEND
                        DO IG = MOIE-7,MOIE
                           IF((IG .GE. ISUBPOS(NPP,1)) .AND.
     &                        (IG .LE. ISUBPOS(NPP,3))) THEN
                              IMESLEN = IMESLEN + 1
                              NPPINDEX2(IMESLEN) = NPP - 1
                              NPPMAX_SEND2   = NPPMAX_SEND2 + 1
                              EXIT
                           ENDIF
                        ENDDO

                     ENDIF
                  ENDIF
               ENDDO

            ELSE 
               
C     EIGENTLICH IG=MOIE-7,MOIE-1 AUSREICHEND
               DO IG = MOIE-7,MOIE
                  IF((IG .GE. ISUBPOS(NP,1)) .AND.
     &               (IG .LE. ISUBPOS(NP,3))) THEN
                     IF (J .LT. JAMIN2) JAMIN2 = J
                     IF (J .GT. JEMAX2) JEMAX2 = J
                     NPPMAX_RECV2 = 1
                     EXIT
                  ENDIF
               ENDDO
                     
            ENDIF

         ENDIF
      ENDDO



      IAMIN2 = IE + 1
      IEMAX2 = IE

      DO IG = MOIE-7,MOIE-1
         IF((IG .GE. ISUBPOS(NP,1)) .AND. (IG .LE. ISUBPOS(NP,3))) THEN
            I = IG - MYPOSGRD(1) + 3
            IF (I .LT. IAMIN2) IAMIN2 = I
         ENDIF
      ENDDO

      IF (NEIGHBOR(3) .EQ. -1) THEN
         IF (IAMIN2 .EQ. (IE+1)) IAMIN2 = IE - 1
      ENDIF


C ********************* NORDRAND******************

      DO M=1,MAXPRC
         NPPINDEX3(M) = -1
      ENDDO

      NPPMAX_SEND3 = 0
      NPPMAX_RECV3 = 0
      IAMIN3 = IE + 1
      IEMAX3 = 0


      DO IG = ISUBPOS(NP,1) , ISUBPOS(NP,3)
         IMESLEN = 1
         IF ((IG .GE. 9) .AND. (IG .LE. (MOIE-8))) THEN
            I = IG - MYPOSGRD(1) + 3

            IF (NEIGHBOR(2) .EQ. -1) THEN

C     EIGENE DATEN
               IF (I .LT. IAMIN3) IAMIN3 = I
               IF (I .GT. IEMAX3) IEMAX3 = I
               NPPMAX_SEND3 = 1
               NPPINDEX3(IMESLEN) = NP - 1

C     FUER ANDERE PROZESSOREN
               DO NPP = 1,NPROC
                  IF (NPP .NE. NP) THEN
                     IF ((IG .GE. ISUBPOS(NPP,1)) .AND.
     &                   (IG .LE. ISUBPOS(NPP,3))) THEN

C     EIGENTLICH JG=MOJE-7,MOJE-1 AUSREICHEND
                        DO JG = MOJE-7,MOJE
                           IF((JG .GE. ISUBPOS(NPP,2)) .AND.
     &                        (JG .LE. ISUBPOS(NPP,4))) THEN
                              IMESLEN = IMESLEN + 1
                              NPPINDEX3(IMESLEN) = NPP - 1
                              NPPMAX_SEND3   = NPPMAX_SEND3 + 1
                              EXIT
                           ENDIF
                        ENDDO

                     ENDIF
                  ENDIF
               ENDDO

            ELSE 
               
C     EIGENTLICH JG=MOJE-7,MOJE-1 AUSREICHEND
               DO JG = MOJE-7,MOJE
                  IF((JG .GE. ISUBPOS(NP,2)) .AND.
     &               (JG .LE. ISUBPOS(NP,4))) THEN
                     IF (I .LT. IAMIN3) IAMIN3 = I
                     IF (I .GT. IEMAX3) IEMAX3 = I
                     NPPMAX_RECV3 = 1
                     EXIT
                  ENDIF
               ENDDO
                     
            ENDIF

         ENDIF
      ENDDO



      JAMIN3 = JE + 1
      JEMAX3 = JE

      DO JG = MOJE-7,MOJE-1
         IF((JG .GE. ISUBPOS(NP,2)) .AND. (JG .LE. ISUBPOS(NP,4))) THEN
            J = JG - MYPOSGRD(2) + 3
            IF (J .LT. JAMIN3) JAMIN3 = J
         ENDIF
      ENDDO

      IF (NEIGHBOR(2) .EQ. -1) THEN
         IF (JAMIN3 .EQ. (JE+1)) JAMIN3 = JE - 1
      ENDIF

C ********************* SUEDRAND******************

      DO M=1,MAXPRC
         NPPINDEX4(M) = -1
      ENDDO

      NPPMAX_SEND4 = 0
      NPPMAX_RECV4 = 0
      IAMIN4 = IE + 1
      IEMAX4 = 0


      DO IG = ISUBPOS(NP,1) , ISUBPOS(NP,3)
         IMESLEN = 1
         IF ((IG .GE. 9) .AND. (IG .LE. (MOIE-8))) THEN
            I = IG - MYPOSGRD(1) + 3

            IF (NEIGHBOR(4) .EQ. -1) THEN

C     EIGENE DATEN
               IF (I .LT. IAMIN4) IAMIN4 = I
               IF (I .GT. IEMAX4) IEMAX4 = I
               NPPMAX_SEND4 = 1
               NPPINDEX4(IMESLEN) = NP - 1

C     FUER ANDERE PROZESSOREN
               DO NPP = 1,NPROC
                  IF (NPP .NE. NP) THEN
                     IF ((IG .GE. ISUBPOS(NPP,1)) .AND.
     &                   (IG .LE. ISUBPOS(NPP,3))) THEN

C     EIGENTLICH JG=2,8 AUSREICHEND
                        DO JG = 1,8
                           IF((JG .GE. ISUBPOS(NPP,2)) .AND.
     &                        (JG .LE. ISUBPOS(NPP,4))) THEN
                              IMESLEN = IMESLEN + 1
                              NPPINDEX4(IMESLEN) = NPP - 1
                              NPPMAX_SEND4   = NPPMAX_SEND4 + 1
                              EXIT
                           ENDIF
                        ENDDO

                     ENDIF
                  ENDIF
               ENDDO

            ELSE 
               
C     EIGENTLICH JG=2,8 AUSREICHEND
               DO JG = 1,8
                  IF((JG .GE. ISUBPOS(NP,2)) .AND.
     &               (JG .LE. ISUBPOS(NP,4))) THEN
                     IF (I .LT. IAMIN4) IAMIN4 = I
                     IF (I .GT. IEMAX4) IEMAX4 = I
                     NPPMAX_RECV4 = 1
                     EXIT
                  ENDIF
               ENDDO
                     
            ENDIF

         ENDIF
      ENDDO



      JEMAX4 = 0
      JAMIN4 = 1

      IF (NEIGHBOR(4) .EQ. -1) THEN
         JEMAX4 = 2
      ENDIF
      
      DO JG = 2,8
         IF((JG .GE. ISUBPOS(NP,2)) .AND. (JG .LE. ISUBPOS(NP,4))) THEN
            J = JG - MYPOSGRD(2) + 3
            JEMAX4 = J
         ENDIF
      ENDDO

C****************SUED-KANTE UND WEST-KANTE DER SUEDWESTECKE*****************

C     NPPMAX_SEND5: HABE ENTWEDER AUF SUEDKANTE MINDESTENS EINEN PUNKT 
C                   ODER AUF WESTKANTE ODER AUF BEIDEN.

      DO M=1,MAXPRC
         NPPINDEX5(M) = -1
      ENDDO

      NPPMAX_SEND5 = 0
      NPPMAX_RECV5 = 0
      IAMIN5 = IE + 1
      IEMAX5 = 0
      JAMIN5 = JE + 1
      JEMAX5 = 0
      IANZP5 = 0
      IW5 = 0
      IS5 = 0

C     EIGENE DATEN

      DO JG = ISUBPOS(NP,2) , ISUBPOS(NP,4)
         DO IG = ISUBPOS(NP,1) , ISUBPOS(NP,3)
            IF ((IG .GE. 2) .AND. (IG .LE. 8) .AND. 
     &          (JG .GE. 2) .AND. (JG .LE. 8)) THEN

               I = IG - MYPOSGRD(1) + 3
C              WENN EINEM PROZESSOR IG=2 GEHOERT, GEHOERT IHM AUCH IG=1,
C              ALSO I=2 (UND I=1).
               IF (IG .EQ. 2) I = 1

               J = JG - MYPOSGRD(2) + 3
               IF (JG .EQ. 2) J = 1

               IF (I .LT. IAMIN5) IAMIN5 = I
               IF (I .GT. IEMAX5) IEMAX5 = I

               IF (J .LT. JAMIN5) JAMIN5 = J
               IF (J .GT. JEMAX5) JEMAX5 = J

               IF (NEIGHBOR(1) .EQ. -1) THEN
                  IMESLEN = 1
                  NPPMAX_SEND5 = 1
                  NPPINDEX5(IMESLEN) = NP - 1
                  IW5 = 1
               ENDIF

               IF (NEIGHBOR(4) .EQ. -1) THEN
                  IMESLEN = 1
                  NPPMAX_SEND5 = 1
                  NPPINDEX5(IMESLEN) = NP - 1
                  IS5 = 1
               ENDIF

               IANZP5 = 1

            ENDIF
         ENDDO
      ENDDO

      IF (IANZP5 .NE. 0) THEN

C        SENDEN AN ANDERE PROZESSOREN

         IF (NPPMAX_SEND5 .GT. 0) THEN
            DO NPP = 1,NPROC
               IF (NPP .NE. NP) THEN
               
                  IEXIT = 0
                  DO JG = ISUBPOS(NPP,2) , ISUBPOS(NPP,4)
                     DO IG = ISUBPOS(NPP,1) , ISUBPOS(NPP,3)
                        IF ((IG .GE. 2) .AND. (IG .LE. 8) .AND. 
     &                      (JG .GE. 2) .AND. (JG .LE. 8)) THEN
                        
                           IMESLEN = IMESLEN + 1
                           NPPINDEX5(IMESLEN) = NPP - 1
                           NPPMAX_SEND5   = NPPMAX_SEND5 + 1
                           IEXIT = 1
                           EXIT
                        
                        ENDIF
                     ENDDO
                     IF (IEXIT .EQ. 1) EXIT
                  ENDDO
               
               ENDIF
            ENDDO
         ENDIF


C        EMPFANGEN

         DO NPP = 1,NPROC
            IF (NPP .NE. NP) THEN

               IEXIT = 0
               DO JG = ISUBPOS(NPP,2) , ISUBPOS(NPP,4)
                  DO IG = ISUBPOS(NPP,1) , ISUBPOS(NPP,3)
                     IF ((IG .GE. 2) .AND. (IG .LE. 8) .AND. 
     &                   (JG .GE. 2) .AND. (JG .LE. 8) .AND. 
     &                   ((ISUBNEIGH(NPP,1) .EQ. -1) .OR.
     &                   (ISUBNEIGH(NPP,4) .EQ. -1))) THEN

                        NPPMAX_RECV5 = NPPMAX_RECV5 + 1
                        IEXIT = 1
                        EXIT

                     ENDIF
                  ENDDO
                  IF (IEXIT .EQ. 1) EXIT
               ENDDO

            ENDIF
         ENDDO

      ENDIF                     !IF (IANZP5 .NE. 0)

C****************SUED-KANTE UND OST-KANTE DER SUEDOSTECKE*****************


      DO M=1,MAXPRC
         NPPINDEX6(M) = -1
      ENDDO

      NPPMAX_SEND6 = 0
      NPPMAX_RECV6 = 0
      IAMIN6 = IE + 1
      IEMAX6 = 0
      JAMIN6 = JE + 1
      JEMAX6 = 0
      IANZP6  = 0
      IS6 = 0
      IO6 = 0

C     EIGENE DATEN

      DO JG = ISUBPOS(NP,2) , ISUBPOS(NP,4)
         DO IG = ISUBPOS(NP,1) , ISUBPOS(NP,3)
            IF ((IG .GE. MOIE-7) .AND. (IG .LE. MOIE-1) .AND. 
     &          (JG .GE. 2) .AND. (JG .LE. 8)) THEN

               I = IG - MYPOSGRD(1) + 3
               IF (IG .EQ. MOIE-1) I = IE

               J = JG - MYPOSGRD(2) + 3
               IF (JG .EQ. 2) J = 1

               IF (I .LT. IAMIN6) IAMIN6 = I
               IF (I .GT. IEMAX6) IEMAX6 = I

               IF (J .LT. JAMIN6) JAMIN6 = J
               IF (J .GT. JEMAX6) JEMAX6 = J

               IF (NEIGHBOR(4) .EQ. -1) THEN
                  IMESLEN = 1
                  NPPMAX_SEND6 = 1
                  NPPINDEX6(IMESLEN) = NP - 1
                  IS6 = 1
               ENDIF

               IF (NEIGHBOR(3) .EQ. -1) THEN
                  IMESLEN = 1
                  NPPMAX_SEND6 = 1
                  NPPINDEX6(IMESLEN) = NP - 1
                  IO6 = 1
               ENDIF

               IANZP6 = 1

            ENDIF
         ENDDO
      ENDDO

      IF (IANZP6 .NE. 0) THEN

C        SENDEN AN ANDERE PROZESSOREN

         IF (NPPMAX_SEND6 .GT. 0) THEN
            DO NPP = 1,NPROC
               IF (NPP .NE. NP) THEN
               
                  IEXIT = 0
                  DO JG = ISUBPOS(NPP,2) , ISUBPOS(NPP,4)
                     DO IG = ISUBPOS(NPP,1) , ISUBPOS(NPP,3)
                        IF ((IG .GE. MOIE-7) .AND. (IG .LE. MOIE-1)
     &                       .AND. (JG .GE. 2) .AND. (JG .LE. 8)) THEN
                        
                           IMESLEN = IMESLEN + 1
                           NPPINDEX6(IMESLEN) = NPP - 1
                           NPPMAX_SEND6   = NPPMAX_SEND6 + 1
                           IEXIT = 1
                           EXIT
                        
                        ENDIF
                     ENDDO
                     IF (IEXIT .EQ. 1) EXIT
                  ENDDO
               
               ENDIF
            ENDDO
         ENDIF


C        EMPFANGEN

         DO NPP = 1,NPROC
            IF (NPP .NE. NP) THEN

               IEXIT = 0
               DO JG = ISUBPOS(NPP,2) , ISUBPOS(NPP,4)
                  DO IG = ISUBPOS(NPP,1) , ISUBPOS(NPP,3)
                     IF ((IG .GE. MOIE-7) .AND. (IG .LE. MOIE-1) .AND. 
     &                    (JG .GE. 2) .AND. (JG .LE. 8) .AND. 
     &                    ((ISUBNEIGH(NPP,4) .EQ. -1) .OR.
     &                    (ISUBNEIGH(NPP,3) .EQ. -1))) THEN

                        NPPMAX_RECV6 = NPPMAX_RECV6 + 1
                        IEXIT = 1
                        EXIT

                     ENDIF
                  ENDDO
                  IF (IEXIT .EQ. 1) EXIT
               ENDDO

            ENDIF
         ENDDO


      ENDIF                     !IF (IANZP6 .NE. 0)

C****************OST-KANTE UND NORD-KANTE DER NORDOSTECKE*****************


      DO M=1,MAXPRC
         NPPINDEX7(M) = -1
      ENDDO

      NPPMAX_SEND7 = 0
      NPPMAX_RECV7 = 0
      IAMIN7 = IE + 1
      IEMAX7 = 0
      JAMIN7 = JE + 1
      JEMAX7 = 0
      IANZP7  = 0
      IO7 = 0
      IN7 = 0

C     EIGENE DATEN

      DO JG = ISUBPOS(NP,2) , ISUBPOS(NP,4)
         DO IG = ISUBPOS(NP,1) , ISUBPOS(NP,3)
            IF ((IG .GE. MOIE-7) .AND. (IG .LE. MOIE-1) .AND. 
     &          (JG .GE. MOJE-7) .AND. (JG .LE. MOJE-1)) THEN

               I = IG - MYPOSGRD(1) + 3
               IF (IG .EQ. MOIE-1) I = IE

               J = JG - MYPOSGRD(2) + 3
               IF (JG .EQ. MOJE-1) J = JE 

               IF (I .LT. IAMIN7) IAMIN7 = I
               IF (I .GT. IEMAX7) IEMAX7 = I

               IF (J .LT. JAMIN7) JAMIN7 = J
               IF (J .GT. JEMAX7) JEMAX7 = J

               IF (NEIGHBOR(3) .EQ. -1) THEN
                  IMESLEN = 1
                  NPPMAX_SEND7 = 1
                  NPPINDEX7(IMESLEN) = NP - 1
                  IO7 = 1
               ENDIF

               IF (NEIGHBOR(2) .EQ. -1) THEN
                  IMESLEN = 1
                  NPPMAX_SEND7 = 1
                  NPPINDEX7(IMESLEN) = NP - 1
                  IN7 = 1
               ENDIF

               IANZP7 = 1

            ENDIF
         ENDDO
      ENDDO

      IF (IANZP7 .NE. 0) THEN

C        SENDEN AN ANDERE PROZESSOREN

         IF (NPPMAX_SEND7 .GT. 0) THEN
            DO NPP = 1,NPROC
               IF (NPP .NE. NP) THEN
               
                  IEXIT = 0
                  DO JG = ISUBPOS(NPP,2) , ISUBPOS(NPP,4)
                     DO IG = ISUBPOS(NPP,1) , ISUBPOS(NPP,3)
                        IF ((IG .GE. MOIE-7) .AND. (IG .LE. MOIE-1)
     &                       .AND. (JG .GE. MOJE-7) 
     &                       .AND. (JG .LE. MOJE-1)) THEN
                        
                           IMESLEN = IMESLEN + 1
                           NPPINDEX7(IMESLEN) = NPP - 1
                           NPPMAX_SEND7   = NPPMAX_SEND7 + 1
                           IEXIT = 1
                           EXIT
                        
                        ENDIF
                     ENDDO
                     IF (IEXIT .EQ. 1) EXIT
                  ENDDO
               
               ENDIF
            ENDDO
         ENDIF


C        EMPFANGEN

         DO NPP = 1,NPROC
            IF (NPP .NE. NP) THEN

               IEXIT = 0
               DO JG = ISUBPOS(NPP,2) , ISUBPOS(NPP,4)
                  DO IG = ISUBPOS(NPP,1) , ISUBPOS(NPP,3)
                     IF ((IG .GE. MOIE-7) .AND. (IG .LE. MOIE-1) .AND. 
     &                    (JG .GE. MOJE-7) .AND. (JG .LE. MOJE-1) .AND. 
     &                    ((ISUBNEIGH(NPP,3) .EQ. -1) .OR.
     &                    (ISUBNEIGH(NPP,2) .EQ. -1))) THEN

                        NPPMAX_RECV7 = NPPMAX_RECV7 + 1
                        IEXIT = 1
                        EXIT

                     ENDIF
                  ENDDO
                  IF (IEXIT .EQ. 1) EXIT
               ENDDO

            ENDIF
         ENDDO

      ENDIF                     !IF (IANZP7 .NE. 0)

C****************NORD-KANTE UND WEST-KANTE DER NORDWESTECKE*****************


      DO M=1,MAXPRC
         NPPINDEX8(M) = -1
      ENDDO

      NPPMAX_SEND8 = 0
      NPPMAX_RECV8 = 0
      IAMIN8 = IE + 1
      IEMAX8 = 0
      JAMIN8 = JE + 1
      JEMAX8 = 0
      IANZP8  = 0
      IN8 = 0
      IW8 = 0

C     EIGENE DATEN

      DO JG = ISUBPOS(NP,2) , ISUBPOS(NP,4)
         DO IG = ISUBPOS(NP,1) , ISUBPOS(NP,3)
            IF ((IG .GE. 2) .AND. (IG .LE. 8) .AND. 
     &          (JG .GE. MOJE-7) .AND. (JG .LE. MOJE-1)) THEN

               I = IG - MYPOSGRD(1) + 3
               IF (IG .EQ. 2) I = 1

               J = JG - MYPOSGRD(2) + 3
               IF (JG .EQ. MOJE-1) J = JE 

               IF (I .LT. IAMIN8) IAMIN8 = I
               IF (I .GT. IEMAX8) IEMAX8 = I

               IF (J .LT. JAMIN8) JAMIN8 = J
               IF (J .GT. JEMAX8) JEMAX8 = J

               IF (NEIGHBOR(2) .EQ. -1) THEN
                  IMESLEN = 1
                  NPPMAX_SEND8 = 1
                  NPPINDEX8(IMESLEN) = NP - 1
                  IN8 = 1
               ENDIF

               IF (NEIGHBOR(1) .EQ. -1) THEN
                  IMESLEN = 1
                  NPPMAX_SEND8 = 1
                  NPPINDEX8(IMESLEN) = NP - 1
                  IW8 = 1
               ENDIF

               IANZP8 = 1

            ENDIF
         ENDDO
      ENDDO

      IF (IANZP8 .NE. 0) THEN

C        SENDEN AN ANDERE PROZESSOREN

         IF (NPPMAX_SEND8 .GT. 0) THEN
            DO NPP = 1,NPROC
               IF (NPP .NE. NP) THEN
               
                  IEXIT = 0
                  DO JG = ISUBPOS(NPP,2) , ISUBPOS(NPP,4)
                     DO IG = ISUBPOS(NPP,1) , ISUBPOS(NPP,3)
                        IF ((IG .GE. 2) .AND. (IG .LE. 8) .AND. 
     &                       (JG .GE. MOJE-7)
     &                       .AND. (JG .LE. MOJE-1)) THEN
                        
                           IMESLEN = IMESLEN + 1
                           NPPINDEX8(IMESLEN) = NPP - 1
                           NPPMAX_SEND8   = NPPMAX_SEND8 + 1
                           IEXIT = 1
                           EXIT
                        
                        ENDIF
                     ENDDO
                     IF (IEXIT .EQ. 1) EXIT
                  ENDDO
               
               ENDIF
            ENDDO
         ENDIF


C        EMPFANGEN

         DO NPP = 1,NPROC
            IF (NPP .NE. NP) THEN

               IEXIT = 0
               DO JG = ISUBPOS(NPP,2) , ISUBPOS(NPP,4)
                  DO IG = ISUBPOS(NPP,1) , ISUBPOS(NPP,3)
                     IF ((IG .GE. 2) .AND. (IG .LE. 8) .AND. 
     &                    (JG .GE. MOJE-7) .AND. (JG .LE. MOJE-1) .AND. 
     &                    ((ISUBNEIGH(NPP,2) .EQ. -1) .OR.
     &                    (ISUBNEIGH(NPP,1) .EQ. -1))) THEN

                        NPPMAX_RECV8 = NPPMAX_RECV8 + 1
                        IEXIT = 1
                        EXIT

                     ENDIF
                  ENDDO
                  IF (IEXIT .EQ. 1) EXIT
               ENDDO

            ENDIF
         ENDDO


      ENDIF                     !IF (IANZP8 .NE. 0)

c      call pstop
c
c      if (myid .eq. 0) then
c      open (unit=300, file='/pf/k/k202057/zzmy0',form='formatted')
c      write(300,*) nppindex1(1),nppindex1(2),nppindex1(3),nppindex1(4)
c      write(300,*) nppmax_send1,nppmax_recv1
c      write(300,*) iamin1,iemax1
c      write(300,*) jamin1,jemax1
c      write(300,*) next
c
c      write(300,*) nppindex2(1),nppindex2(2),nppindex2(3),nppindex2(4)
c      write(300,*) nppmax_send2,nppmax_recv2
c      write(300,*) iamin2,iemax2
c      write(300,*) jamin2,jemax2
c      write(300,*) next
c
c      write(300,*) nppindex3(1),nppindex3(2),nppindex3(3),nppindex3(4)
c      write(300,*) nppmax_send3,nppmax_recv3
c      write(300,*) iamin3,iemax3
c      write(300,*) jamin3,jemax3
c      write(300,*) next
c
c      write(300,*) nppindex4(1),nppindex4(2),nppindex4(3),nppindex4(4)
c      write(300,*) nppmax_send4,nppmax_recv4
c      write(300,*) iamin4,iemax4
c      write(300,*) jamin4,jemax4
c      write(300,*) next
c
c      write(300,*) nppindex5(1),nppindex5(2),nppindex5(3),nppindex5(4)
c      write(300,*) nppmax_send5,nppmax_recv5
c      write(300,*) iamin5,iemax5
c      write(300,*) jamin5,jemax5
c      write(300,*) is5,iw5
c      write(300,*) next
c
c      write(300,*) nppindex6(1),nppindex6(2),nppindex6(3),nppindex6(4)
c      write(300,*) nppmax_send6,nppmax_recv6
c      write(300,*) iamin6,iemax6
c      write(300,*) jamin6,jemax6
c      write(300,*) is6,io6
c      write(300,*) next
c
c      write(300,*) nppindex7(1),nppindex7(2),nppindex7(3),nppindex7(4)
c      write(300,*) nppmax_send7,nppmax_recv7
c      write(300,*) iamin7,iemax7
c      write(300,*) jamin7,jemax7
c      write(300,*) in7,io7
c      write(300,*) next
c
c      write(300,*) nppindex8(1),nppindex8(2),nppindex8(3),nppindex8(4)
c      write(300,*) nppmax_send8,nppmax_recv8
c      write(300,*) iamin8,iemax8
c      write(300,*) jamin8,jemax8
c      write(300,*) in8,iw8
c      write(300,*) next
c
c      close(300)
c      endif
c
c      call pstop
c
c      if (myid .eq. 1) then
c      open (unit=300, file='/pf/k/k202057/zzmy1',form='formatted')
c      write(300,*) nppindex1(1),nppindex1(2),nppindex1(3),nppindex1(4)
c      write(300,*) nppmax_send1,nppmax_recv1
c      write(300,*) iamin1,iemax1
c      write(300,*) jamin1,jemax1
c      write(300,*) next
c
c      write(300,*) nppindex2(1),nppindex2(2),nppindex2(3),nppindex2(4)
c      write(300,*) nppmax_send2,nppmax_recv2
c      write(300,*) iamin2,iemax2
c      write(300,*) jamin2,jemax2
c      write(300,*) next
c
c      write(300,*) nppindex3(1),nppindex3(2),nppindex3(3),nppindex3(4)
c      write(300,*) nppmax_send3,nppmax_recv3
c      write(300,*) iamin3,iemax3
c      write(300,*) jamin3,jemax3
c      write(300,*) next
c
c      write(300,*) nppindex4(1),nppindex4(2),nppindex4(3),nppindex4(4)
c      write(300,*) nppmax_send4,nppmax_recv4
c      write(300,*) iamin4,iemax4
c      write(300,*) jamin4,jemax4
c      write(300,*) next
c
c      write(300,*) nppindex5(1),nppindex5(2),nppindex5(3),nppindex5(4)
c      write(300,*) nppmax_send5,nppmax_recv5
c      write(300,*) iamin5,iemax5
c      write(300,*) jamin5,jemax5
c      write(300,*) is5,iw5
c      write(300,*) next
c
c      write(300,*) nppindex6(1),nppindex6(2),nppindex6(3),nppindex6(4)
c      write(300,*) nppmax_send6,nppmax_recv6
c      write(300,*) iamin6,iemax6
c      write(300,*) jamin6,jemax6
c      write(300,*) is6,io6
c      write(300,*) next
c
c      write(300,*) nppindex7(1),nppindex7(2),nppindex7(3),nppindex7(4)
c      write(300,*) nppmax_send7,nppmax_recv7
c      write(300,*) iamin7,iemax7
c      write(300,*) jamin7,jemax7
c      write(300,*) in7,io7
c      write(300,*) next
c
c      write(300,*) nppindex8(1),nppindex8(2),nppindex8(3),nppindex8(4)
c      write(300,*) nppmax_send8,nppmax_recv8
c      write(300,*) iamin8,iemax8
c      write(300,*) jamin8,jemax8
c      write(300,*) in8,iw8
c      write(300,*) next
c
c      close(300)
c      endif
c
c      call pstop
c
c      if (myid .eq. 2) then
c      open (unit=300, file='/pf/k/k202057/zzmy2',form='formatted')
c      write(300,*) nppindex1(1),nppindex1(2),nppindex1(3),nppindex1(4)
c      write(300,*) nppmax_send1,nppmax_recv1
c      write(300,*) iamin1,iemax1
c      write(300,*) jamin1,jemax1
c      write(300,*) next
c
c      write(300,*) nppindex2(1),nppindex2(2),nppindex2(3),nppindex2(4)
c      write(300,*) nppmax_send2,nppmax_recv2
c      write(300,*) iamin2,iemax2
c      write(300,*) jamin2,jemax2
c      write(300,*) next
c
c      write(300,*) nppindex3(1),nppindex3(2),nppindex3(3),nppindex3(4)
c      write(300,*) nppmax_send3,nppmax_recv3
c      write(300,*) iamin3,iemax3
c      write(300,*) jamin3,jemax3
c      write(300,*) next
c
c      write(300,*) nppindex4(1),nppindex4(2),nppindex4(3),nppindex4(4)
c      write(300,*) nppmax_send4,nppmax_recv4
c      write(300,*) iamin4,iemax4
c      write(300,*) jamin4,jemax4
c      write(300,*) next
c
c      write(300,*) nppindex5(1),nppindex5(2),nppindex5(3),nppindex5(4)
c      write(300,*) nppmax_send5,nppmax_recv5
c      write(300,*) iamin5,iemax5
c      write(300,*) jamin5,jemax5
c      write(300,*) is5,iw5
c      write(300,*) next
c
c      write(300,*) nppindex6(1),nppindex6(2),nppindex6(3),nppindex6(4)
c      write(300,*) nppmax_send6,nppmax_recv6
c      write(300,*) iamin6,iemax6
c      write(300,*) jamin6,jemax6
c      write(300,*) is6,io6
c      write(300,*) next
c
c      write(300,*) nppindex7(1),nppindex7(2),nppindex7(3),nppindex7(4)
c      write(300,*) nppmax_send7,nppmax_recv7
c      write(300,*) iamin7,iemax7
c      write(300,*) jamin7,jemax7
c      write(300,*) in7,io7
c      write(300,*) next
c
c      write(300,*) nppindex8(1),nppindex8(2),nppindex8(3),nppindex8(4)
c      write(300,*) nppmax_send8,nppmax_recv8
c      write(300,*) iamin8,iemax8
c      write(300,*) jamin8,jemax8
c      write(300,*) in8,iw8
c      write(300,*) next
c
c      close(300)
c      endif
c
c      call pstop
c
c      if (myid .eq. 3) then
c      open (unit=300, file='/pf/k/k202057/zzmy3',form='formatted')
c      write(300,*) nppindex1(1),nppindex1(2),nppindex1(3),nppindex1(4)
c      write(300,*) nppmax_send1,nppmax_recv1
c      write(300,*) iamin1,iemax1
c      write(300,*) jamin1,jemax1
c      write(300,*) next
c
c      write(300,*) nppindex2(1),nppindex2(2),nppindex2(3),nppindex2(4)
c      write(300,*) nppmax_send2,nppmax_recv2
c      write(300,*) iamin2,iemax2
c      write(300,*) jamin2,jemax2
c      write(300,*) next
c
c      write(300,*) nppindex3(1),nppindex3(2),nppindex3(3),nppindex3(4)
c      write(300,*) nppmax_send3,nppmax_recv3
c      write(300,*) iamin3,iemax3
c      write(300,*) jamin3,jemax3
c      write(300,*) next
c
c      write(300,*) nppindex4(1),nppindex4(2),nppindex4(3),nppindex4(4)
c      write(300,*) nppmax_send4,nppmax_recv4
c      write(300,*) iamin4,iemax4
c      write(300,*) jamin4,jemax4
c      write(300,*) next
c
c      write(300,*) nppindex5(1),nppindex5(2),nppindex5(3),nppindex5(4)
c      write(300,*) nppmax_send5,nppmax_recv5
c      write(300,*) iamin5,iemax5
c      write(300,*) jamin5,jemax5
c      write(300,*) is5,iw5
c      write(300,*) next
c
c      write(300,*) nppindex6(1),nppindex6(2),nppindex6(3),nppindex6(4)
c      write(300,*) nppmax_send6,nppmax_recv6
c      write(300,*) iamin6,iemax6
c      write(300,*) jamin6,jemax6
c      write(300,*) is6,io6
c      write(300,*) next
c
c      write(300,*) nppindex7(1),nppindex7(2),nppindex7(3),nppindex7(4)
c      write(300,*) nppmax_send7,nppmax_recv7
c      write(300,*) iamin7,iemax7
c      write(300,*) jamin7,jemax7
c      write(300,*) in7,io7
c      write(300,*) next
c
c      write(300,*) nppindex8(1),nppindex8(2),nppindex8(3),nppindex8(4)
c      write(300,*) nppmax_send8,nppmax_recv8
c      write(300,*) iamin8,iemax8
c      write(300,*) jamin8,jemax8
c      write(300,*) in8,iw8
c      write(300,*) next
c
c      close(300)
c      endif
c
c
c
c      stop


      RETURN
      END

