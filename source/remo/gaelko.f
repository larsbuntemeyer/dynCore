C**********************************************************************
C
C                  SUBROUTINE GAELKO 
C
C**** GAELKO   -   UP:BERECHNUNG DER KOEFFIZIENTEN ZUR LOESUNG DER
C****                 HELMHOLTZ-GLEICHUNG IN N-S-RICHTUNG MIT
C****                 GAUSSSCHER ELIMINATION; FORTRAN77-VERSION
C**   AUFRUF   :   CALL GAELKO (COEF,BPT,M,N,AP,BP,CP,DI,M1,N1,PI)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   BERECHNUNG DER KOEFFIZIENTEN ZUR LOESUNG DER
C**                HELMHOLTZ-GLEICHUNG IN N-S-RICHTUNG MIT
C**                GAUSSSCHER ELIMINATION
C**   VERSIONS-
C**   DATUM    :   05.04.89
C**
C**   EXTERNALS:   KEINE
C**
C**   EINGABE-
C**   PARAMETER:   M    : ANZAHL DER GITTERPUNKTE IN I-RICHTUNG
C**                N    : ANZAHL DER GITTERPUNKTE IN J-RICHTUNG
C**                M1   : M - 1
C**                N1   : N - 1
C**                PI   : ZAHL PI (3.141592654)
C**   AUSGABE-
C**   PARAMETER:   BP,DI: DIAGONALEN
C**                AP,CP: NEBENDIAGONALEN
C**                COEF : KOEFFIZIENTENMATRIX
C**
C**   COMMON-
C**   BLOECKE  :   KEINE
C**
C**   METHODE  :   LAM = S * A * S**(-1)
C**                LAM = BP(J) + 2 * COS ((I-1)*PI/(M-1)), I = 1,...,M-1
C**
C**                GAUSS SCHE ELIMINATION EINER TRIDIAGONALEN MATRIX
C**                DIAGONALE BP(J) + DI(I), NEBENDIAGONALEN AP BZW. CP
C**               (SIEHE VARGA, MATRIX ITERATIVE ANALYSIS, S. 195 FF.)
C**
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
      SUBROUTINE GAELKO (COEF,BPT,ILO,IUP,JLO,JUP,M,N,AP,BP,CP,DI,PI)

      IMPLICIT NONE

C**********************************************************************

      INCLUDE "parorg.h"

C     VARIABLEN DER PARAMETERLISTE
      INTEGER, INTENT(IN)               ::  ILO, IUP, JLO, JUP, M, N
      REAL,    INTENT(IN), DIMENSION(N) ::  AP, BP, CP
      REAL,    INTENT(INOUT)            ::
     &                        COEF((ILO-1):(IUP+1),(JLO-1):(JUP+1)), 
     &                        BPT ((ILO-1):(IUP+1),(JLO-1):(JUP+1)), 
     &                        DI  ((ILO-1):(IUP+1))

      REAL  PI

C     LOKALE VARIABLE 
      REAL  PIR
      INTEGER   I, J

C**********************************************************************

      PIR = PI / FLOAT (M-1)
 
      DO I = ILO , IUP
         DI(I) =  2. * COS (PIR * FLOAT (I-1))
      ENDDO
 
      DO J = JLO , JUP
         DO I = ILO , IUP
            BPT(I,J) = BP(J) + DI(I)
         ENDDO
      ENDDO
 
      DO I = ILO , IUP
         BPT(I,JLO)  = 1. / BPT(I,JLO)
         COEF(I,JLO) = CP(2) * BPT(I,JLO)
      ENDDO
 
      IF (JLO .EQ. JUP)  RETURN
      DO J = JLO+1 , JUP
         DO I = ILO , IUP
            BPT(I,J)  = 1. / (BPT(I,J) - AP(J) * COEF(I,J-1))
            COEF(I,J) = CP(J) * BPT(I,J)
         ENDDO
      ENDDO


      RETURN
      END SUBROUTINE GAELKO
