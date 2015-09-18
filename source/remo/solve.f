C
C     SUBROUTINE SOLVE
C
C**********************************************************************
C
C**** SOLVE    -   UP:LOESUNG EINER ZWEIDIMENSIONALEN HELMHOLTZ-
C****                 GLEICHUNG MIT FFT IN W-E- UND GAUSS-ELIMINATION
C****                 IN S-N-RICHTUNG.
C**   AUFRUF       CALL SOLVE(X,F,COEF,BPT,W1,W2,M,N,AP,TRIGS,MT,
C**                           M1,N1,N2,Z1,Z2,IFAX)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   LOESUNG EINER ZWEIDIMENSIONALEN HELMHOLTZ-
C**                GLEICHUNG MIT FFT IN W-E- UND GAUSS-ELIMINATION
C**                IN S-N-RICHTUNG; FORTRAN77-VERSION.
C**   VERSIONS-
C**   DATUM    :   05.04.89
C**
C**   EXTERNALS:   FFT
C**
C**   EINGABE-
C**   PARAMETER:   F(M,N) : RECHTE SEITE DER HELMHOLTZ-GLEICHUNG
C**   AUSGABE-
C**   PARAMETER:   X(M,N) : LOESUNGSFELD
C**
C**   COMMON-
C**   BLOECKE  :   KEINE
C**
C**   METHODE  :   FFT IN W-E-RICHTUNG; GAUSS-ELIMINATION IN S-N
C**     SCHNELLE FOURIERTRANSFORMATION IN X-RICHTUNG
C**     TRIDIAGONALE LOESUNGEN MIT HILFE DER GAUSS SCHEN ELIMINATION
C**     FFT / TRIDIAGONALE LOESUNG / INVERSE FFT
C**     M1 = M - 1 BESCHRAENKT AUF PRODUKT VON POTENZEN VON 2, 3, 5, 6
C**     N1 = N - 1 BELIEBIG
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   P.PROHL
C**
C**********  PARALLELE VERSION  ***************************************
C**
C**   DATUM:   09.06.94
C**   AUTOR:   ULRICH SCHAETTLER
C**
C**********************************************************************
      SUBROUTINE SOLVE (X2D,F2D,F_FFT,W1FFT,W2FFT,X_GAUSS,COEF,
     &                  BPT,AP,TRIGS,Z1,Z2,IFAX,MYNFFT,MYMGAUSS,MYMT)
C
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
C-----------------------------------------------------------------------
C     Dummy Arguments
C
C     TYPES OF THE PARAMETERS
      INTEGER, INTENT(IN)    :: MYNFFT, MYMGAUSS, MYMT, IFAX(11)
      REAL,    INTENT(IN)    :: F2D(IE,JE), 
     &                          COEF(MYMGAUSS,MOJE),BPT(MYMGAUSS,MOJE),
     &                          AP(MOJE), TRIGS(MYMT),  Z1(6),  Z2(6)
      REAL,    INTENT(INOUT) :: X2D(IE,JE),  X_GAUSS(MYMGAUSS,MOJE),
     &                          F_FFT(MOIE,MYNFFT),
     &                          W1FFT(MOIE,MYNFFT),
     &                          W2FFT(MOIE,MYNFFT)
C-----------------------------------------------------------------------
C     Local Declarations
C
C     LOCAL VARIABLES
      INTEGER    ::       I, J

C**********************************************************************
 
C     INITIALISIERUNGEN
      DO J = 1 , MYNFFT
        DO I = 1 , MOIE
          F_FFT(I,J) = 0.0
        ENDDO
      ENDDO

      DO J = 1,MOJE
        DO I = 1,MYMGAUSS
          X_GAUSS(I,J) = 0.0
        ENDDO
      ENDDO 

      DO J = 1,JE
        DO I = 1,IE
          X2D(I,J) = 0.0
        ENDDO
      ENDDO 

C---------------------------------------------------------------------- 
C     1. TRANSPOSITION:   F2D ---> F_FFT
C----------------------------------------------------------------------


      CALL TRAPO (F2D,IE,JE,1,
     &            F_FFT,MOIE,MYNFFT,1,1)


C----------------------------------------------------------------------
C     SINUSTRANSFORMATION FUER DIE N-2 REIHEN DES INNEREN DES GEBIETES
C     BERECHNET WIRD G(DACH) = S * F
C     SOLUTION:   W1FFT
C----------------------------------------------------------------------

      CALL FFT (F_FFT,W1FFT,W2FFT,MOIE,MYNFFT,TRIGS,MYMT,
     &          +1,Z1,Z2,IFAX)



C----------------------------------------------------------------------
C     2. TRANSPOSITION:   W1FFT ---> X_GAUSS
C----------------------------------------------------------------------
 
      CALL TRAPO (W1FFT,MOIE,MYNFFT,1,
     &            X_GAUSS,MYMGAUSS,MOJE,1,2)


C----------------------------------------------------------------------
C     1. UND 2. TEIL DER GAUSS SCHEN ELIMINATION
C     BERECHNET WIRD Y(DACH)
C     SOLUTION:  X_GAUSS
C----------------------------------------------------------------------
 
      DO I = 2 , MYMGAUSS-1
        X_GAUSS(I,MYGAUSS_JLO) = X_GAUSS(I,MYGAUSS_JLO)
     &             * BPT(I,MYGAUSS_JLO)
      END DO
 
      DO J = MYGAUSS_JLO+1 , MYGAUSS_JUP
        DO I = 2 , MYMGAUSS-1
          X_GAUSS(I,J) = BPT(I,J) * 
     &                        (X_GAUSS(I,J) - AP(J) * X_GAUSS(I,J-1))
        END DO
      END DO
 
      DO J = MYGAUSS_JUP-1 , MYGAUSS_JLO , -1
         DO I = 2 , MYMGAUSS-1
            X_GAUSS(I,J) = X_GAUSS(I,J) - COEF(I,J) * X_GAUSS(I,J+1)
          END DO
       END DO

C----------------------------------------------------------------------
C     3. TRANSPOSITION:   X_GAUSS --->  F_FFT
C----------------------------------------------------------------------

      CALL TRAPO (X_GAUSS,MYMGAUSS,MOJE,1,
     &            F_FFT,MOIE,MYNFFT,1,3)

 
C----------------------------------------------------------------------
C     INVERSE SINUSTRANSFORMATION
C     BERECHNET WIRD Y = S**(-1) * Y(DACH)
C     SOLUTION:   W1FFT
C----------------------------------------------------------------------
 
      CALL FFT (F_FFT,W1FFT,W2FFT,MOIE,MYNFFT,TRIGS,MYMT,
     &          -1,Z1,Z2,IFAX)
 
C----------------------------------------------------------------------
C     4. TRANSPOSITION:   W1FFT  --->  X2D
C----------------------------------------------------------------------


      CALL TRAPO (W1FFT,MOIE,MYNFFT,1,
     &            X2D,IE,JE,1,4)


C----------------------------------------------------------------------

      RETURN
      END SUBROUTINE SOLVE
