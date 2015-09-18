C
C     SUBROUTINE SOLVE3D
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
C**********  PARALLELE VERSION MIT 3D TRANSPOSITION *******************
C**
C**   DATUM:   01.11.94
C**   AUTOR:   ULRICH SCHAETTLER
C**
C**********************************************************************
      SUBROUTINE SOLVE3D (X3D,F3D,
     &           AP,BP,CP,DI,PI,TRIGS,Z1,Z2,IFAX,
     &           NFFT,MGAUSS,MT)
C
      IMPLICIT NONE
C     INCLUDING OF HEADER-FILES
      INCLUDE "parorg.h"
C-----------------------------------------------------------------------
C     Dummy Arguments
C
C     TYPES OF THE PARAMETERS
      INTEGER, INTENT(IN)    ::  NFFT, MGAUSS, MT, IFAX(11)
      REAL,    INTENT(INOUT) ::  X3D(IE,JE,KE)
      REAL,    INTENT(IN)    ::  F3D(IE,JE,KE), AP(MOJE), BP(MOJE,KE), 
     &                           CP(MOJE), TRIGS(MT),
     &                           Z1(6),Z2(6),PI
      REAL,    INTENT(INOUT) ::  DI(MOIE)
C-----------------------------------------------------------------------
C     Local Declarations
C
C     LOCAL VARIABLES
      INTEGER :: I,J,K
      REAL    :: F_FFT(MOIE,NFFT,KE), W1FFT(MOIE,NFFT),
     &           W2FFT(MOIE,NFFT), X_GAUSS(MGAUSS,MOJE,KE),
     &           COEF(MGAUSS,MOJE),BPT(MGAUSS,MOJE)
C
C**********************************************************************


C     INITIALISIERUNGEN
      DO K = 1 , KE
         DO J = 1 , NFFT
            DO I = 1 , MOIE
               F_FFT(I,J,K) = 0.0
            ENDDO
         ENDDO

         DO J = 1,MOJE
            DO I = 1,MGAUSS
               X_GAUSS(I,J,K) = 0.0
            ENDDO
         ENDDO

         DO J = 1,JE
            DO I = 1,IE
               X3D(I,J,K) = 0.0
            ENDDO
         ENDDO

      ENDDO

C----------------------------------------------------------------------
C     1. TRANSPOSITION:   F3D ---> F_FFT   (PROZESSOREN-ZEILE)
C----------------------------------------------------------------------


      CALL TRAPO (F3D,IE,JE,KE,
     &     F_FFT,MOIE,NFFT,KE,1)



C----------------------------------------------------------------------
C     SINUSTRANSFORMATION FUER DIE N-2 REIHEN DES INNEREN DES GEBIETES
C     BERECHNET WIRD G(DACH) = S * F
C     SOLUTION:   W1FFT
C----------------------------------------------------------------------

      DO K = 1 , KE
C        BERECHNEN DER TRANSFORMATIONEN
         CALL FFT (F_FFT(1,1,K),W1FFT,W2FFT,MOIE,NFFT,TRIGS,MT,
     &        +1,Z1,Z2,IFAX)


C        UMSPEICHERN DER LOESUNG AUF F_FFT
         DO J = 2 , NFFT-1
            DO I = 2 , MOIE-1
               F_FFT(I,J,K) = W1FFT(I,J)
            ENDDO
         ENDDO
      ENDDO


C----------------------------------------------------------------------
C     2. TRANSPOSITION:   F_FFT ---> F3D ---> X_GAUSS
C----------------------------------------------------------------------

      CALL TRAPO (F_FFT,MOIE,NFFT,KE,
     &     X_GAUSS,MGAUSS,MOJE,KE,2)


C----------------------------------------------------------------------
C     1. UND 2. TEIL DER GAUSS SCHEN ELIMINATION
C     BERECHNET WIRD Y(DACH)
C     SOLUTION:  X_GAUSS
C----------------------------------------------------------------------

      DO K = 1 , KE
C       BERECHNEN DER KOEFFIZIENTEN
         CALL GAELKO (COEF,BPT,MYGAUSS_ILO,MYGAUSS_IUP,MYGAUSS_JLO,
     &        MYGAUSS_JUP,MOIE,MOJE,AP,BP(1,K),CP,DI,PI)

         DO I = 2 , MGAUSS-1
            X_GAUSS(I,MYGAUSS_JLO,K) = X_GAUSS(I,MYGAUSS_JLO,K)
     &           * BPT(I,MYGAUSS_JLO)
         ENDDO

         DO J = MYGAUSS_JLO+1 , MYGAUSS_JUP
            DO I = 2 , MGAUSS-1
               X_GAUSS(I,J,K) = BPT(I,J) *
     &              (X_GAUSS(I,J,K) - AP(J) * X_GAUSS(I,J-1,K))
            ENDDO
         ENDDO

         DO J = MYGAUSS_JUP-1 , MYGAUSS_JLO , -1
            DO I = 2 , MGAUSS-1
               X_GAUSS(I,J,K) = X_GAUSS(I,J,K) - COEF(I,J)
     &              * X_GAUSS(I,J+1,K)
            ENDDO
         ENDDO
      ENDDO

C----------------------------------------------------------------------
C     3. TRANSPOSITION:   X_GAUSS --->  F3D ---> F_FFT
C----------------------------------------------------------------------


      CALL TRAPO (X_GAUSS,MGAUSS,MOJE,KE,
     &     F_FFT,MOIE,NFFT,KE,3)


C----------------------------------------------------------------------
C     INVERSE SINUSTRANSFORMATION
C     BERECHNET WIRD Y = S**(-1) * Y(DACH)
C     SOLUTION:   W1FFT
C----------------------------------------------------------------------

      DO K = 1 , KE
         CALL FFT (F_FFT(1,1,K),W1FFT,W2FFT,MOIE,NFFT,TRIGS,MT,
     &        -1,Z1,Z2,IFAX)

         DO J = 2 , NFFT-1
            DO I = 2 , MOIE-1
               F_FFT(I,J,K) = W1FFT(I,J)
            ENDDO
         ENDDO
      ENDDO

C----------------------------------------------------------------------
C     4. TRANSPOSITION:   F_FFT  --->  X3D
C----------------------------------------------------------------------


      CALL TRAPO (F_FFT,MOIE,NFFT,KE,
     &     X3D,IE,JE,KE,4)

C----------------------------------------------------------------------

      RETURN
      END SUBROUTINE SOLVE3D
