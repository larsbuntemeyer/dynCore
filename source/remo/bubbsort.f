      SUBROUTINE BUBBSORT(VECTOR,FELD,M,N)
C
      IMPLICIT NONE
C
C     ZWECK:
C     BUBBSORT SORTIERT DEN VECTOR NACH DEM BUBBLESORT-ALGORITHMUS
C     FALLS ZWEI ELEMENTE DES VECTORS GETAUSCH WURDEN, WIRD ANALOG
C     AUCH ZWEI VECTOREN DES FELDES GETAUSCH.
C     VECTOR    (INPUT/OUTPUT) : ZU SORTIERENDER VECTOR
C     FELD      (INPUT/OUTPUT) : ZU SORTIERENDES FELD
C     N         (INPUT)        : ANZAHL DER ELEMENTE VON VECTOR
C     M         (INPUT)        : 2. DIMENSION VON FELD
C     FEHLERBEHANDLUNG: FALLS N GROESSER ALS DIE DIMENSIONIERUNG
C                       VON FDUM IST, WIRD DAS PROGRAMM ABGEBROCHEN.
C     AUTHOR: AXEL DIEHL
C     DATUM : FRI MAR 10 09:29:37 UTC 1995
C
      INTEGER, INTENT(IN)    :: M,N
      REAL,    INTENT(INOUT) :: VECTOR(N), FELD(M,N)
      REAL    :: VDUM, FDUM(50)
      LOGICAL :: TAUSCH
      INTEGER :: I, J
C
      IF (N.GT.50) THEN
         PRINT *, 'FEHLER IN BUBBSORT'
         STOP
      END IF
C
C     BUBBLESORT-ALGORITHMUS:
C     EINE SCHLEIFE LAEUFT UEBER DEN GESAMTEN VECTOR. ES WERDEN JEWEILS
C     ZWEI BENACHBARTE ELEMENTE VERGLICHEN. FALLS ELEMENT VECTOR(I)
C     KLEINER ALS VECTOR(I+1) IST, WERDEN DIE BEIDEN ELEMENTE VERTAUSCHT.
C     DIESER VORGANG WIRD SOOFT WIEDERHOLT, WIE ES VERTAUSCHUNGEN GIBT.
C     ZUM SCHLUSS SIND ALLE ELEMENTE IN ABSTEIGENDER REIHENFOLGE SORTIERT.
C
      TAUSCH = .TRUE.
      DO WHILE (TAUSCH)
         DO J=1,M-1
            IF (VECTOR(J).LT.VECTOR(J+1)) THEN
               VDUM = VECTOR(J)
               DO I=1,N
                  FDUM(I) = FELD(I,J)
               ENDDO
               VECTOR(J) = VECTOR(J+1)
               DO I=1,N
                  FELD(I,J) = FELD(I,J+1)
               ENDDO
               VECTOR(J+1) = VDUM
               DO I=1,N
                  FELD(I,J+1) = FDUM(I)
               ENDDO
               TAUSCH = .TRUE.
            ELSE
               TAUSCH = .FALSE.
            ENDIF
         ENDDO
      ENDDO
C
      RETURN
      END SUBROUTINE BUBBSORT
