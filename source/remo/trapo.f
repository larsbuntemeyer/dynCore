C
C     SUBROUTINE TRAPO
C
C**********************************************************************
C
C     NAME      :  NTRAPOFAST2
C     CALL      :  CALL NTRAPOFAST2 FROM NPARSOLVE3D
C     PURPOSE   :  TRANSPOSES A 3D-FIELD BETWEEN DIFFERENT DECOMPOSITIONS
C     AUTHOR    :  ULRICH SCHAETTLER
C     DATE      :  23.12.94
C
C     EXTERNALS :                 -----
C
C     PARAMATERLIST INPUT:      FIELDIN, ALL INTEGER INDICES
C
C     PARAMETERLIST OUTPUT:     FIELDOUT
C
C     HEADER FILES:  NPARAM.H, NORG.H, NSENDBUF.H, NPARORG.H
C
C     VARIABLES OF NAMED COMMON-BLOCKS THAT ARE CHANGED:
C                                 -----
C
C     METHOD:       PARMACS-MACROS
C                   EVERY PROCESSOR SENDS THE DATA NECESSARY TO
C                   EACH OTHER PROCESSOR. 4 PARTS ARE IMPLEMENTED FOR
C                   THE 4 DIFFERENT TRANSPOSITIONS.
C
C                   THE DIFFERENCE TO NTRAPOFAST IS THAT EVERY
C                   PROCESSOR SENDS TO ONE OTHER PROCESSOR AND
C                   THEN IMMEDIATELY RECEIVES THE DATA SENT TO IT.
C
C**********************************************************************
C
      SUBROUTINE TRAPO (FIELDIN,IFIN,JFIN,KFIN,
     &                  FIELDOUT,IFOUT,JFOUT,KFOUT,INDEX)
C
      IMPLICIT NONE
C     (* INCLUDING HEADER-FILES *)
C
      INCLUDE "org.h"
      INCLUDE "parorg.h"
C
C**********************************************************************
C
C     TYPES OF THE PARAMETERS
      INTEGER, INTENT(IN)    :: IFIN, JFIN, KFIN, IFOUT, 
     &                          JFOUT, KFOUT, INDEX

      REAL,    INTENT(IN)    :: FIELDIN(IFIN,JFIN,KFIN)
      REAL,    INTENT(INOUT) :: FIELDOUT(IFOUT,JFOUT,KFOUT)                        
C
C     LOCAL VARIABLES
      REAL          ::   REBUF(IE*JE*KE,NPROC)
      INTEGER       ::   I, J, K, L, IMESLEN,
     &                   ILOC, JLOC, JLOC1,
     &                   JLO,JUP,IPARTNER,DL
      DATA DL/0/

C**********************************************************************


C      DL=0
C      CALL PSTOP

C     1. TRANSPOSITION:  2D ---> FFT
      IF (INDEX .EQ. 1) THEN
         TAGCOUNT=1
         TAGTABLE(1)=40000+DL

!        LOOP OVER NUMBER OF PROCESSORS IN X-DIRECTION
         DO L = 1 , NPROCX-1

C        SEND THE DATA TO THE NEXT PROCESSOR

C        DETERMINATION OF THE RECEIVER
            IF (MYPOS(1)+L .LE. NPROCX) THEN
               IPARTNER = MYID + L+1
            ELSE
               IPARTNER = MYID + L - NPROCX+1
            ENDIF

C        PUT THE DATA INTO THE BUFFER
C        J-INDICES FOR INPUT-FIELD
            JLO = JAH + NFFTJLO(IPARTNER) - MYPOSGRD(2)
            JUP = JAH + NFFTJUP(IPARTNER) - MYPOSGRD(2)
            IMESLEN = 0
            DO K = 1 , KFIN
               DO J = JLO , JUP
                  DO I = IAH , IEH
                     IMESLEN = IMESLEN + 1
                     REBUF(IMESLEN,L) = FIELDIN(I,J,K)
                  ENDDO
               ENDDO
            ENDDO

C        SEND THE DATA
            TYPE=TAGTABLE(1)
            COUNT=IMESLEN*1
            DEST=IPARTNER-1
            CALL PSENDR(REBUF(1,L))

C        RECEIVE DATA FROM OTHER PROCESSORS
            CALL PTEST
            CALL PRECVR(REBUF(1,NPROC))
            IPARTNER=SOURCE+1

C        FETCH DATA FROM BUFFER
            IMESLEN = 0
            DO K = 1 , KFIN
               DO J = 2 , JFOUT-1
                  DO I = ISUBPOS(IPARTNER,1) , ISUBPOS(IPARTNER,3)
                     IMESLEN = IMESLEN + 1
                     FIELDOUT(I,J,K) = REBUF(IMESLEN,NPROC)
                  ENDDO
               ENDDO
            ENDDO

         ENDDO !DO L = 1 , NPROCX-1

C      PUT OWN VALUES TO OUTPUT-FIELD
         JLOC1 = JAH + NFFTJLO(MYID+1) - MYPOSGRD(2) - 1
         DO K = 1 , KFIN
            JLOC = JLOC1
            DO J = 2 , JFOUT-1
               JLOC = JLOC + 1
               ILOC = 2
               DO I = ISUBPOS(MYID+1,1) , ISUBPOS(MYID+1,3)
                  ILOC = ILOC + 1
                  FIELDOUT(I,J,K) = FIELDIN(ILOC,JLOC,K)
               ENDDO
            ENDDO
         ENDDO

      ENDIF
C     END OF 1. TRANSPOSITION 2D ---> FFT

C----------------------------------------------------------------------

C     2. TRANSPOSITION:  FFT ---> GAUSS
      IF (INDEX .EQ. 2) THEN
         TAGCOUNT=1
         TAGTABLE(1)=40000+DL+1

         DO L = 1 , NPROC-1

C     SEND OWN DATA TO OTHER PROCESSORS

C     DETERMINATION OF THE RECEIVER
            IF (MYID+L+1 .LE. NPROC) THEN
               IPARTNER = MYID + L+1
            ELSE
               IPARTNER = MYID + L - NPROC+1
            ENDIF

C     PUT OWN DATA INTO THE BUFFER
            IMESLEN = 0
            DO K = 1 , KFIN
               DO J = 2 , JFIN-1
                  DO I = NGAUSSILO(IPARTNER) , NGAUSSIUP(IPARTNER)
                     IMESLEN = IMESLEN + 1
                     REBUF(IMESLEN,L) = FIELDIN(I,J,K)
                  ENDDO
               ENDDO
            ENDDO

C     SEND THE DATA
            TYPE=TAGTABLE(1)
            COUNT=IMESLEN*1
            DEST=IPARTNER-1
            CALL PSENDR(REBUF(1,L))

C     RECEIVE THE DATA FROM OTHER PROCESSORS
            CALL PTEST
            CALL PRECVR(REBUF(1,NPROC))

C     FETCH DATA FROM BUFFER
            IPARTNER=SOURCE+1
            IMESLEN = 0
            DO K = 1 , KFIN
               DO J = NFFTJLO(IPARTNER) , NFFTJUP(IPARTNER)
                  DO I = 2 , IFOUT-1
                     IMESLEN = IMESLEN + 1
                     FIELDOUT(I,J,K) = REBUF(IMESLEN,NPROC)
                  ENDDO
               ENDDO
            ENDDO

         ENDDO

C      PUT THE OWN VALUES INTO THE OUTPUT-FIELD
         DO K = 1 , KFIN
            JLOC = 1
            DO J = NFFTJLO(MYID+1) , NFFTJUP(MYID+1)
               JLOC = JLOC + 1
               ILOC = NGAUSSILO(MYID+1) - 1
               DO I = 2 , IFOUT-1
                  ILOC = ILOC + 1
                  FIELDOUT(I,J,K) = FIELDIN(ILOC,JLOC,K)
               ENDDO
            ENDDO
         ENDDO

      ENDIF
C     END OF 2. TRANSPOSITION FFT ---> GAUSS

C----------------------------------------------------------------------

C     3. TRANSPOSITION:  GAUSS ---> FFT
      IF (INDEX .EQ. 3) THEN
         TAGCOUNT=1
         TAGTABLE(1)=40000+DL+2

         DO L = 1 , NPROC-1

C        SEND OWN DATA TO OTHER PROCESSORS

C        DETERMINATION OF THE RECEIVER
            IF (MYID+L+1 .LE. NPROC) THEN
               IPARTNER = MYID + L +1
            ELSE
               IPARTNER = MYID + L - NPROC +1
            ENDIF

C        PUT OWN DATA INTO THE BUFFER
            IMESLEN = 0
            DO K = 1 , KFIN
               DO J = NFFTJLO(IPARTNER) , NFFTJUP(IPARTNER)
                  DO I = 2 , IFIN-1
                     IMESLEN = IMESLEN + 1
                     REBUF(IMESLEN,L) = FIELDIN(I,J,K)
                  ENDDO
               ENDDO
            ENDDO

C        SEND THE DATA
            TYPE=TAGTABLE(1)
            COUNT=IMESLEN*1
            DEST=IPARTNER-1
            CALL PSENDR(REBUF(1,L))

C        RECEIVE THE DATA FROM OTHER PROCESSORS
            CALL PTEST
            CALL PRECVR(REBUF(1,NPROC))
            IPARTNER=SOURCE+1

C        FETCH DATA FROM BUFFER
            IMESLEN = 0
            DO K = 1 , KFIN
               DO J = 2 , JFOUT-1
                  DO I = NGAUSSILO(IPARTNER) , NGAUSSIUP(IPARTNER)
                     IMESLEN = IMESLEN + 1
                     FIELDOUT(I,J,K) = REBUF(IMESLEN,NPROC)
                  ENDDO
               ENDDO
            ENDDO

         ENDDO

C      PUT OWN VALUES INTO THE OUTPUT-FIELD
         DO K = 1 , KFIN
            JLOC = NFFTJLO(MYID+1) - 1
            DO J = 2 , JFOUT-1
               JLOC = JLOC + 1
               ILOC = 1
               DO I = NGAUSSILO(MYID+1) , NGAUSSIUP(MYID+1)
                  ILOC = ILOC + 1
                  FIELDOUT(I,J,K) = FIELDIN(ILOC,JLOC,K)
               ENDDO
            ENDDO
         ENDDO

      ENDIF
C     END OF 3. TRANSPOSITION GAUSS ---> FFT

C----------------------------------------------------------------------

C     4. TRANSPOSITION:  FFT ---> 2D
      IF (INDEX .EQ. 4) THEN
         TAGCOUNT=1
         TAGTABLE(1)=40000+DL+3

         DO L = 1 , NPROCX-1

C        SEND OWN DATA TO OTHER PROCESSORS
C        DETERMINATION OF THE RECEIVER
            IF (MYPOS(1)+L .LE. NPROCX) THEN
               IPARTNER = MYID + L+1
            ELSE
               IPARTNER = MYID + L - NPROCX +1
            ENDIF

C        PUT OWN DATA INTO THE BUFFER
            IMESLEN = 0
            DO K = 1 , KFIN
               DO J = 2 , JFIN-1
                  DO I = ISUBPOS(IPARTNER,1) , ISUBPOS(IPARTNER,3)
                     IMESLEN = IMESLEN + 1
                     REBUF(IMESLEN,L) = FIELDIN(I,J,K)
                  ENDDO
               ENDDO
            ENDDO

C        SEND THE DATA
            TYPE=TAGTABLE(1)
            COUNT=IMESLEN*1
            DEST=IPARTNER-1

            CALL PSENDR(REBUF(1,L))

C        RECEIVE DATA FROM OTHER PROCESSORS
            CALL PTEST
            CALL PRECVR(REBUF(1,NPROC))
            IPARTNER=SOURCE+1

C        FETCH DATA FROM BUFFER
C        J-INDICES FOR INPUT-FIELD
            JLO = JAH + NFFTJLO(IPARTNER) - MYPOSGRD(2)
            JUP = JAH + NFFTJUP(IPARTNER) - MYPOSGRD(2)
            IMESLEN = 0
            DO K = 1 , KFIN
               DO J = JLO , JUP
                  DO I = IAH , IEH
                     IMESLEN = IMESLEN + 1
                     FIELDOUT(I,J,K) = REBUF(IMESLEN,NPROC)
                  ENDDO
               ENDDO
            ENDDO

       ENDDO

C      PUT OWN VALUES INTO THE OUTPUT-FIELD
       JLOC1 = JAH + NFFTJLO(MYID+1) - MYPOSGRD(2) - 1
       DO K = 1 , KFIN
         JLOC = JLOC1
         DO J = 2 , JFIN-1
            JLOC = JLOC + 1
            ILOC = ISUBPOS(MYID+1,1) - 1
            DO I = IAH , IEH
               ILOC = ILOC + 1
               FIELDOUT(I,JLOC,K) = FIELDIN(ILOC,J,K)
            ENDDO
         ENDDO
      ENDDO

      ENDIF
C     END OF 4. TRANSPOSITION FFT ---> 2D

C----------------------------------------------------------------------

      CALL PWAIT
      DL=DL+4

C     MPI ON T3D: MAXTAG=65534
      IF (40000+DL+3 .GT. 65534) THEN
         DL=0
      ENDIF

      CALL PSTOP
      RETURN
      END SUBROUTINE TRAPO
