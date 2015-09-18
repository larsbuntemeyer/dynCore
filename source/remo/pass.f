C
C
C@(#) BIBLIOTHEK EM: EUROPA-MODELL / DEUTSCHLAND-MODELL
C@(#) MODUL PASS.F, V2.12 VOM 4/27/94, EXTRAHIERT AM 4/28/94
C
C**** PASS     -   UP:DURCHFUEHRUNG DER FAST FOURIER TRANSFORM
C**   AUFRUF       CALL PASS(U,X,M,N,TRIGS,MT,INC,NH,IFAC,LA,
C**                          F1,F2)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   DURCHFUEHRUNG DER FAST FOURIER TRANSFORM OHNE PRE-
C**                UND POSTPROCESSING, VERSION IN FORTRAN77
C**   VERSIONS-
C**   DATUM    :   05.04.89
C**
C**   EXTERNALS:   KEINE
C**
C**   EINGABE-
C**   PARAMETER:   U(M,N) : INPUTFELD
C**                INC    : ADDRESSENINCREMENT FUER U(1) UND X(1)
C**                NH     : LAENGE DER KOMPLEXEN VEKTOREN
C**                IFAC   : LAUFENDER FAKTOR
C**                MC     : ZUSCHLAG FUER COSINUSTABELLE
C**                LA     : PRODUKT DER VORHERGEHENDEN FAKTOREN
C**   AUSGABE-
C**   PARAMETER:   X(M,N) : OUTPUTFELD
C**
C**   COMMON-
C**   BLOECKE  :   KEINE
C**
C**   METHODE  :   FFT NACH TEMPERTON
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   P.PROHL
C----------------------------------------------------------------------
      SUBROUTINE PASS (U,X,M,N,TRIGS,MT,INC,NH,IFAC,LA,
     &                 F1,F2)
C
      IMPLICIT NONE
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: M,N,MT,INC,NH,IFAC,LA
      REAL,    INTENT(IN)    :: F1(6),F2(6)
      REAL,    INTENT(IN)    :: U(M,N)
      REAL,    INTENT(IN)    :: TRIGS(MT)
      REAL,    INTENT(INOUT) :: X(M,N)
C
C     U(1), U(2) INPUT -VEKTOR (REAL- UND IMAGINAERTEIL)
C     X(1), X(2) OUTPUT-VEKTOR (REAL- UND IMAGINAERTEIL)
C     INC  ADDRESSENINCREMENT FUER U(1) UND X(1)
C     NH   LAENGE DER KOMPLEXEN VEKTOREN
C     IFAC LAUFENDER FAKTOR
C     MC   ZUSCHLAG FUER COSINUSTABELLE
C     LA   PRODUKT DER VORHERGEHENDEN FAKTOREN
C
C     Local Variables
C
      REAL :: A0, A00, A01, A1, A10, A11, A2, A20, A21, A3, A30, A31, 
     &        A40, A41, A5, A50, A51, A6, A7, A4
      REAL :: A8, A9, B0, B00, B01, B1, B10, B11, B2, B20, B21, B3, B30, 
     &        B31, B4, B40, B41, B5, B50, B51
      REAL :: B6, B7, B8, B9, CO1, CO2, CO3, CO4, CO5, CO6, CO7, SI1, 
     &        SI3, SI4, SI5, SI6, SI7, Z11, Z12, SI2
      REAL :: Z13, Z14, Z21, Z22, Z23, ZWSP1, ZWSP2, ZWSPI1, ZWSPI2, 
     &        ZWSPI3, ZWSPI4, ZWSPI5, ZWSPI6, ZWSPI7, ZWSPR1, ZWSPR2, 
     &        ZWSPR3, ZWSPR4, ZWSPR5, ZWSPR6
      REAL :: ZWSPR7
      INTEGER :: ID3, IU1, IU2, IU3, IU4, IU5, IU6, IU7, IU8, IX1, IX2, 
     &           IX3, IX4, IX5, IX6, IX7, IX8, J, JD3, JUMP
      INTEGER :: K, KB, KC, KD, KE, KF, KG, KH, L, MC, MI, N1
C----------------------------------------------------------------------
      N1   = N -1
      MI   = NH / IFAC
      ID3  = MI * INC
      JD3  = LA * INC
      JUMP = (IFAC  - 1) * JD3
      MC   = 2 * NH

      IU1 = 1
      IX1 = 1
      IU2 = IU1 + ID3
      IX2 = IX1 + JD3

CKS      IGO = IFAC - 1
CKS      IF (IFAC .EQ. 8)  IGO = IGO - 1
CKS      GO TO (10,110,210,310,410,510), IGO
      SELECT CASE (IFAC)

C-----------------------------------------------------------------------

C     CODING FUER FAKTOR 2
      CASE (2)
CKS 10      DO K = 1, MI, LA
         DO K = 1, MI, LA
            IF (K .NE. 1) THEN
               KB = 2 * K - 1
               SI1 = TRIGS (KB)
               CO1 = TRIGS (KB+MC)
            ENDIF

            DO L = 1, LA
               IF (K .EQ. 1) THEN
                  DO J = 2 , N1
                     X(IX1,J)   = U(IU1,J)   + U(IU2,J)
                     X(IX1+1,J) = U(IU1+1,J) + U(IU2+1,J)
                     X(IX2,J)   = U(IU1,J)   - U(IU2,J)
                     X(IX2+1,J) = U(IU1+1,J) - U(IU2+1,J)
                  ENDDO
               ELSE
                  DO J = 2 , N1
                     X(IX1,J)   = U(IU1,J)   + U(IU2,J)
                     X(IX1+1,J) = U(IU1+1,J) + U(IU2+1,J)
                     ZWSP1 = U(IU1,J)   - U(IU2,J)
                     ZWSP2 = U(IU1+1,J) - U(IU2+1,J)
                     X(IX2,J)   = CO1 * ZWSP1 - SI1 * ZWSP2
                     X(IX2+1,J) = SI1 * ZWSP1 + CO1 * ZWSP2
                  ENDDO
               ENDIF
               IU1 = IU1 + INC
               IU2 = IU2 + INC
               IX1 = IX1 + INC
               IX2 = IX2 + INC
            ENDDO
            IX1 = IX1 + JUMP
            IX2 = IX2 + JUMP
         ENDDO

C-----------------------------------------------------------------------

C     CODING FUER FAKTOR 3
      CASE (3)
CKS  110 IU3 = IU2 + ID3
         IU3 = IU2 + ID3
         IX3 = IX2 + JD3
         Z11 = F1(1)
         Z21 = F2(1)
         DO K = 1, MI, LA
            IF (K .NE. 1) THEN
               KB = 2 * K - 1
               KC = 2 * KB - 1
               SI1 = TRIGS(KB)
               CO1 = TRIGS(KB+MC)
               SI2 = TRIGS(KC)
               CO2 = TRIGS(KC+MC)
            ENDIF

            DO L = 1, LA
               DO J = 2 , N1
                  A1 = U(IU2,J) + U(IU3,J)
                  A2 = U(IU1,J) - Z11 * A1
                  A3 = Z21 * (U(IU2,J) - U(IU3,J))
                  B1 = U(IU2+1,J) + U(IU3+1,J)
                  B2 = U(IU1+1,J) - Z11 * B1
                  B3 = Z21 * (U(IU2+1,J) - U(IU3+1,J))
                  X(IX1,J)   = U(IU1,J)   + A1
                  X(IX1+1,J) = U(IU1+1,J) + B1
                  X(IX3,J)   = A2 + B3
                  X(IX2,J)   = A2 - B3
                  X(IX2+1,J) = B2 + A3
                  X(IX3+1,J) = B2 - A3
               ENDDO

               IF (K .NE. 1) THEN
                  DO J = 2 , N1
                     ZWSPR1 = CO1 * X(IX2,J) - SI1 * X(IX2+1,J)
                     ZWSPI1 = SI1 * X(IX2,J) + CO1 * X(IX2+1,J)
                     X(IX2,J) = ZWSPR1
                     X(IX2+1,J) = ZWSPI1
                     ZWSPR2 = CO2 * X(IX3,J) - SI2 * X(IX3+1,J)
                     ZWSPI2 = SI2 * X(IX3,J) + CO2 * X(IX3+1,J)
                     X(IX3,J)   = ZWSPR2
                     X(IX3+1,J) = ZWSPI2
                  ENDDO
               ENDIF
               IU1 = IU1 + INC
               IU2 = IU2 + INC
               IU3 = IU3 + INC
               IX1 = IX1 + INC
               IX2 = IX2 + INC
               IX3 = IX3 + INC
            ENDDO
            IX1 = IX1 + JUMP
            IX2 = IX2 + JUMP
            IX3 = IX3 + JUMP
         ENDDO

C-----------------------------------------------------------------------

C     CODING FUER FAKTOR 4
      CASE (4)
CKS  210 IU3 = IU2 + ID3
         IU3 = IU2 + ID3
         IU4 = IU3 + ID3
         IX3 = IX2 + JD3
         IX4 = IX3 + JD3
         DO K = 1, MI, LA
            IF (K .NE. 1) THEN
               KB = 2 * K - 1
               KC = 2 * KB - 1
               KD = KB + KC - 1
               SI1 = TRIGS(KB)
               CO1 = TRIGS(KB+MC)
               SI2 = TRIGS(KC)
               CO2 = TRIGS(KC+MC)
               SI3 = TRIGS(KD)
               CO3 = TRIGS(KD+MC)
            ENDIF

            DO L = 1, LA
               DO J = 2 , N1
                  A0 = U(IU1,J) + U(IU3,J)
                  A1 = U(IU2,J) + U(IU4,J)
                  A2 = U(IU1,J) - U(IU3,J)
                  A3 = U(IU2,J) - U(IU4,J)
                  B0 = U(IU1+1,J) + U(IU3+1,J)
                  B1 = U(IU2+1,J) + U(IU4+1,J)
                  B2 = U(IU1+1,J) - U(IU3+1,J)
                  B3 = U(IU2+1,J) - U(IU4+1,J)
                  X(IX1,J)   = A0 + A1
                  X(IX3,J)   = A0 - A1
                  X(IX1+1,J) = B0 + B1
                  X(IX3+1,J) = B0 - B1
                  X(IX2,J)   = A2 - B3
                  X(IX4,J)   = A2 + B3
                  X(IX2+1,J) = B2 + A3
                  X(IX4+1,J) = B2 - A3
               ENDDO

               IF (K .NE. 1) THEN
                  DO J = 2 , N1
                     ZWSPR1 = CO1 * X(IX2,J) - SI1 * X(IX2+1,J)
                     ZWSPI1 = SI1 * X(IX2,J) + CO1 * X(IX2+1,J)
                     X(IX2,J)   = ZWSPR1
                     X(IX2+1,J) = ZWSPI1
                     ZWSPR2 = CO2 * X(IX3,J) - SI2 * X(IX3+1,J)
                     ZWSPI2 = SI2 * X(IX3,J) + CO2 * X(IX3+1,J)
                     X(IX3,J)   = ZWSPR2
                     X(IX3+1,J) = ZWSPI2
                     ZWSPR3 = CO3 * X(IX4,J) - SI3 * X(IX4+1,J)
                     ZWSPI3 = SI3 * X(IX4,J) + CO3 * X(IX4+1,J)
                     X(IX4,J)   = ZWSPR3
                     X(IX4+1,J) = ZWSPI3
                  ENDDO
               ENDIF
               IU1 = IU1 + INC
               IU2 = IU2 + INC
               IU3 = IU3 + INC
               IU4 = IU4 + INC
               IX1 = IX1 + INC
               IX2 = IX2 + INC
               IX3 = IX3 + INC
               IX4 = IX4 + INC
            ENDDO
            IX1 = IX1 + JUMP
            IX2 = IX2 + JUMP
            IX3 = IX3 + JUMP
            IX4 = IX4 + JUMP
         ENDDO

C-----------------------------------------------------------------------

C     CODING FUER FAKTOR 5
      CASE (5)
CKS  310 IU3 = IU2 + ID3
         IU3 = IU2 + ID3
         IU4 = IU3 + ID3
         IU5 = IU4 + ID3
         IX3 = IX2 + JD3
         IX4 = IX3 + JD3
         IX5 = IX4 + JD3
         Z12 = F1(2)
         Z13 = F1(3)
         Z22 = F2(2)
         Z23 = F2(3)
         DO K = 1, MI, LA
            IF (K .NE. 1) THEN
               KB = 2*K - 1
               KC = 2*KB - 1
               KD = KB + KC - 1
               KE = KB + KD - 1
               SI1 = TRIGS(KB)
               CO1 = TRIGS(KB+MC)
               SI2 = TRIGS(KC)
               CO2 = TRIGS(KC+MC)
               SI3 = TRIGS(KD)
               CO3 = TRIGS(KD+MC)
               SI4 = TRIGS(KE)
               CO4 = TRIGS(KE+MC)
            ENDIF

            DO L = 1, LA
               DO J = 2 , N1
                  A1 = U(IU2,J) + U(IU5,J)
                  A2 = U(IU3,J) + U(IU4,J)
                  A3 = U(IU2,J) - U(IU5,J)
                  A4 = U(IU3,J) - U(IU4,J)
                  A5 = A1 + A2
                  A6 = Z12 * (A1 - A2)
                  A7 = U(IU1,J) - Z13 * A5
                  B1 = U(IU2+1,J) + U(IU5+1,J)
                  B2 = U(IU3+1,J) + U(IU4+1,J)
                  B3 = U(IU2+1,J) - U(IU5+1,J)
                  B4 = U(IU3+1,J) - U(IU4+1,J)
                  B5 = B1 + B2
                  B6 = Z12 * (B1 - B2)
                  B7 = U(IU1+1,J) - Z13 * B5
                  A10 = A7 + A6
                  A11 = B3 * Z23 + B4 * Z22
                  A20 = A7 - A6
                  A21 = B3 * Z22 - B4 * Z23
                  B10 = B7 + B6
                  B11 = A3 * Z23 + A4 * Z22
                  B20 = B7 - B6
                  B21 = A3 * Z22 - A4 * Z23
                  X(IX1,J)   = U(IU1,J) + A5
                  X(IX1+1,J) = U(IU1+1,J) + B5
                  X(IX2,J)   = A10 - A11
                  X(IX2+1,J) = B10 + B11
                  X(IX3,J)   = A20 - A21
                  X(IX3+1,J) = B20 + B21
                  X(IX4,J)   = A20 + A21
                  X(IX4+1,J) = B20 - B21
                  X(IX5,J)   = A10 + A11
                  X(IX5+1,J) = B10 - B11
               ENDDO

               IF (K .NE. 1) THEN
                  DO J = 2 , N1
                     ZWSPR1 = CO1 * X(IX2,J) - SI1 * X(IX2+1,J)
                     ZWSPI1 = SI1 * X(IX2,J) + CO1 * X(IX2+1,J)
                     X(IX2,J)   = ZWSPR1
                     X(IX2+1,J) = ZWSPI1
                     ZWSPR2 = CO2 * X(IX3,J) - SI2 * X(IX3+1,J)
                     ZWSPI2 = SI2 * X(IX3,J) + CO2 * X(IX3+1,J)
                     X(IX3,J)   = ZWSPR2
                     X(IX3+1,J) = ZWSPI2
                     ZWSPR3 = CO3 * X(IX4,J) - SI3 * X(IX4+1,J)
                     ZWSPI3 = SI3 * X(IX4,J) + CO3 * X(IX4+1,J)
                     X(IX4,J)   = ZWSPR3
                     X(IX4+1,J) = ZWSPI3
                     ZWSPR4 = CO4 * X(IX5,J) - SI4 * X(IX5+1,J)
                     ZWSPI4 = SI4 * X(IX5,J) + CO4 * X(IX5+1,J)
                     X(IX5,J)   = ZWSPR4
                     X(IX5+1,J) = ZWSPI4
                  ENDDO
               ENDIF
               IU1 = IU1 + INC
               IU2 = IU2 + INC
               IU3 = IU3 + INC
               IU4 = IU4 + INC
               IU5 = IU5 + INC
               IX1 = IX1 + INC
               IX2 = IX2 + INC
               IX3 = IX3 + INC
               IX4 = IX4 + INC
               IX5 = IX5 + INC
            ENDDO
            IX1 = IX1 + JUMP
            IX2 = IX2 + JUMP
            IX3 = IX3 + JUMP
            IX4 = IX4 + JUMP
            IX5 = IX5 + JUMP
         ENDDO

C-----------------------------------------------------------------------

C     CODING FUER FAKTOR 6
      CASE (6)
CKS  410 IU3 = IU2 + ID3
         IU3 = IU2 + ID3
         IU4 = IU3 + ID3
         IU5 = IU4 + ID3
         IU6 = IU5 + ID3
         IX3 = IX2 + JD3
         IX4 = IX3 + JD3
         IX5 = IX4 + JD3
         IX6 = IX5 + JD3
         Z11 = F1(1)
         Z21 = F2(1)
         DO K = 1, MI, LA
            IF (K .NE. 1) THEN
               KB = 2*K - 1
               KC = 2*KB - 1
               KD = KB + KC - 1
               KE = KB + KD - 1
               KF = KB + KE - 1
               SI1 = TRIGS(KB)
               CO1 = TRIGS(KB+MC)
               SI2 = TRIGS(KC)
               CO2 = TRIGS(KC+MC)
               SI3 = TRIGS(KD)
               CO3 = TRIGS(KD+MC)
               SI4 = TRIGS(KE)
               CO4 = TRIGS(KE+MC)
               SI5 = TRIGS(KF)
               CO5 = TRIGS(KF+MC)
            ENDIF

            DO L = 1, LA
               DO J = 2 , N1
                  A10 = U(IU1,J) + U(IU4,J)
                  A11 = U(IU1,J) - U(IU4,J)
                  A20 = U(IU3,J) + U(IU6,J)
                  A21 = U(IU3,J) - U(IU6,J)
                  A30 = U(IU5,J) + U(IU2,J)
                  A31 = U(IU5,J) - U(IU2,J)
                  A40 = A20 + A30
                  A50 = A21 + A31
                  A1 = A10 - Z11 * A40
                  A2 = A11 - Z11 * A50
                  A3 = Z21 * (A20 - A30)
                  A4 = Z21 * (A21 - A31)
                  B10 = U(IU1+1,J) + U(IU4+1,J)
                  B11 = U(IU1+1,J) - U(IU4+1,J)
                  B20 = U(IU3+1,J) + U(IU6+1,J)
                  B21 = U(IU3+1,J) - U(IU6+1,J)
                  B30 = U(IU5+1,J) + U(IU2+1,J)
                  B31 = U(IU5+1,J) - U(IU2+1,J)
                  B40 = B20 + B30
                  B50 = B21 + B31
                  B1 = B10 - Z11 * B40
                  B2 = B11 - Z11 * B50
                  B3 = Z21 * (B20 - B30)
                  B4 = Z21 * (B21 - B31)
                  X(IX1,J)   = A10 + A40
                  X(IX1+1,J) = B10 + B40
                  X(IX4,J)   = A11 + A50
                  X(IX4+1,J) = B11 + B50
                  X(IX3,J)   = A1 + B3
                  X(IX3+1,J) = B1 - A3
                  X(IX5,J)   = A1 - B3
                  X(IX5+1,J) = B1 + A3
                  X(IX2,J)   = A2 - B4
                  X(IX2+1,J) = B2 + A4
                  X(IX6,J)   = A2 + B4
                  X(IX6+1,J) = B2 - A4
               ENDDO

               IF (K .NE. 1) THEN
                  DO J = 2 , N1
                     ZWSPR1 = CO1 * X(IX2,J) - SI1 * X(IX2+1,J)
                     ZWSPI1 = SI1 * X(IX2,J) + CO1 * X(IX2+1,J)
                     X(IX2,J)   = ZWSPR1
                     X(IX2+1,J) = ZWSPI1
                     ZWSPR2 = CO2 * X(IX3,J) - SI2 * X(IX3+1,J)
                     ZWSPI2 = SI2 * X(IX3,J) + CO2 * X(IX3+1,J)
                     X(IX3,J)   = ZWSPR2
                     X(IX3+1,J) = ZWSPI2
                     ZWSPR3 = CO3 * X(IX4,J) - SI3 * X(IX4+1,J)
                     ZWSPI3 = SI3 * X(IX4,J) + CO3 * X(IX4+1,J)
                     X(IX4,J)   = ZWSPR3
                     X(IX4+1,J) = ZWSPI3
                     ZWSPR4 = CO4 * X(IX5,J) - SI4 * X(IX5+1,J)
                     ZWSPI4 = SI4 * X(IX5,J) + CO4 * X(IX5+1,J)
                     X(IX5,J)   = ZWSPR4
                     X(IX5+1,J) = ZWSPI4
                     ZWSPR5 = CO5 * X(IX6,J) - SI5 * X(IX6+1,J)
                     ZWSPI5 = SI5 * X(IX6,J) + CO5 * X(IX6+1,J)
                     X(IX6,J)   = ZWSPR5
                     X(IX6+1,J) = ZWSPI5
                  ENDDO
               ENDIF
               IU1 = IU1 + INC
               IU2 = IU2 + INC
               IU3 = IU3 + INC
               IU4 = IU4 + INC
               IU5 = IU5 + INC
               IU6 = IU6 + INC
               IX1 = IX1 + INC
               IX2 = IX2 + INC
               IX3 = IX3 + INC
               IX4 = IX4 + INC
               IX5 = IX5 + INC
               IX6 = IX6 + INC
            ENDDO
            IX1 = IX1 + JUMP
            IX2 = IX2 + JUMP
            IX3 = IX3 + JUMP
            IX4 = IX4 + JUMP
            IX5 = IX5 + JUMP
            IX6 = IX6 + JUMP
         ENDDO

C-----------------------------------------------------------------------

C     CODING FUER FAKTOR 8
      CASE (8)
CKS  510 IU3 = IU2 + ID3
         IU3 = IU2 + ID3
         IU4 = IU3 + ID3
         IU5 = IU4 + ID3
         IU6 = IU5 + ID3
         IU7 = IU6 + ID3
         IU8 = IU7 + ID3
         IX3 = IX2 + JD3
         IX4 = IX3 + JD3
         IX5 = IX4 + JD3
         IX6 = IX5 + JD3
         IX7 = IX6 + JD3
         IX8 = IX7 + JD3
         Z14 = F1(4)
         DO K = 1, MI, LA
            IF (K .NE. 1) THEN
               KB = 2*K - 1
               KC = 2*KB - 1
               KD = KB + KC - 1
               KE = KB + KD - 1
               KF = KB + KE - 1
               KG = KB + KF - 1
               KH = KB + KG - 1
               SI1 = TRIGS(KB)
               CO1 = TRIGS(KB+MC)
               SI2 = TRIGS(KC)
               CO2 = TRIGS(KC+MC)
               SI3 = TRIGS(KD)
               CO3 = TRIGS(KD+MC)
               SI4 = TRIGS(KE)
               CO4 = TRIGS(KE+MC)
               SI5 = TRIGS(KF)
               CO5 = TRIGS(KF+MC)
               SI6 = TRIGS(KG)
               CO6 = TRIGS(KG+MC)
               SI7 = TRIGS(KH)
               CO7 = TRIGS(KH+MC)
            ENDIF

            DO L = 1, LA
               DO J = 2 , N1
                  A00 = U(IU1,J) + U(IU5,J)
                  A01 = U(IU3,J) + U(IU7,J)
                  A10 = U(IU2,J) + U(IU6,J)
                  A11 = U(IU4,J) + U(IU8,J)
                  A0 = A00 + A01
                  A2 = A00 - A01
                  A1 = A10 + A11
                  A3 = A10 - A11
                  B00 = U(IU1+1,J) + U(IU5+1,J)
                  B01 = U(IU3+1,J) + U(IU7+1,J)
                  B10 = U(IU2+1,J) + U(IU6+1,J)
                  B11 = U(IU4+1,J) + U(IU8+1,J)
                  B0 = B00 + B01
                  B2 = B00 - B01
                  B1 = B10 + B11
                  B3 = B10 - B11
                  X(IX1,J)   = A0 + A1
                  X(IX1+1,J) = B0 + B1
                  X(IX5,J)   = A0 - A1
                  X(IX5+1,J) = B0 - B1
                  X(IX3,J)   = A2 - B3
                  X(IX3+1,J) = B2 + A3
                  X(IX7,J)   = A2 + B3
                  X(IX7+1,J) = B2 - A3
                  A40 = U(IU1,J) - U(IU5,J)
                  A41 = U(IU3,J) - U(IU7,J)
                  B40 = U(IU1+1,J) - U(IU5+1,J)
                  B41 = U(IU3+1,J) - U(IU7+1,J)
                  A4 = A40 - B41
                  A6 = A40 + B41
                  B4 = B40 + A41
                  B6 = B40 - A41
                  A50 = U(IU2,J) - U(IU6,J)
                  A51 = U(IU4,J) - U(IU8,J)
                  A8 = A50 + A51
                  A9 = A50 - A51
                  B50 = U(IU2+1,J) - U(IU6+1,J)
                  B51 = U(IU4+1,J) - U(IU8+1,J)
                  B8 = B50 + B51
                  B9 = B50 - B51
                  A5 = Z14 * (A9 - B8)
                  A7 =-Z14 * (A9 + B8)
                  B5 = Z14 * (A8 + B9)
                  B7 = Z14 * (A8 - B9)
                  X(IX2,J)   = A4 + A5
                  X(IX2+1,J) = B4 + B5
                  X(IX6,J)   = A4 - A5
                  X(IX6+1,J) = B4 - B5
                  X(IX4,J)   = A6 + A7
                  X(IX4+1,J) = B6 + B7
                  X(IX8,J)   = A6 - A7
                  X(IX8+1,J) = B6 - B7
               ENDDO

               IF (K .NE. 1) THEN
                  DO J = 2 , N1
                     ZWSPR1 = CO1 * X(IX2,J) - SI1 * X(IX2+1,J)
                     ZWSPI1 = SI1 * X(IX2,J) + CO1 * X(IX2+1,J)
                     X(IX2,J)   = ZWSPR1
                     X(IX2+1,J) = ZWSPI1
                     ZWSPR2 = CO2 * X(IX3,J) - SI2 * X(IX3+1,J)
                     ZWSPI2 = SI2 * X(IX3,J) + CO2 * X(IX3+1,J)
                     X(IX3,J)   = ZWSPR2
                     X(IX3+1,J) = ZWSPI2
                     ZWSPR3 = CO3 * X(IX4,J) - SI3 * X(IX4+1,J)
                     ZWSPI3 = SI3 * X(IX4,J) + CO3 * X(IX4+1,J)
                     X(IX4,J)   = ZWSPR3
                     X(IX4+1,J) = ZWSPI3
                     ZWSPR4 = CO4 * X(IX5,J) - SI4 * X(IX5+1,J)
                     ZWSPI4 = SI4 * X(IX5,J) + CO4 * X(IX5+1,J)
                     X(IX5,J)   = ZWSPR4
                     X(IX5+1,J) = ZWSPI4
                     ZWSPR5 = CO5 * X(IX6,J) - SI5 * X(IX6+1,J)
                     ZWSPI5 = SI5 * X(IX6,J) + CO5 * X(IX6+1,J)
                     X(IX6,J)   = ZWSPR5
                     X(IX6+1,J) = ZWSPI5
                     ZWSPR6 = CO6 * X(IX7,J) - SI6 * X(IX7+1,J)
                     ZWSPI6 = SI6 * X(IX7,J) + CO6 * X(IX7+1,J)
                     X(IX7,J)   = ZWSPR6
                     X(IX7+1,J) = ZWSPI6
                     ZWSPR7 = CO7 * X(IX8,J) - SI7 * X(IX8+1,J)
                     ZWSPI7 = SI7 * X(IX8,J) + CO7 * X(IX8+1,J)
                     X(IX8,J)   = ZWSPR7
                     X(IX8+1,J) = ZWSPI7
                  ENDDO
               ENDIF
               IU1 = IU1 + INC
               IU2 = IU2 + INC
               IU3 = IU3 + INC
               IU4 = IU4 + INC
               IU5 = IU5 + INC
               IU6 = IU6 + INC
               IU7 = IU7 + INC
               IU8 = IU8 + INC
               IX1 = IX1 + INC
               IX2 = IX2 + INC
               IX3 = IX3 + INC
               IX4 = IX4 + INC
               IX5 = IX5 + INC
               IX6 = IX6 + INC
               IX7 = IX7 + INC
               IX8 = IX8 + INC
            ENDDO
            IX1 = IX1 + JUMP
            IX2 = IX2 + JUMP
            IX3 = IX3 + JUMP
            IX4 = IX4 + JUMP
            IX5 = IX5 + JUMP
            IX6 = IX6 + JUMP
            IX7 = IX7 + JUMP
            IX8 = IX8 + JUMP
         ENDDO

      END SELECT

      RETURN
      END SUBROUTINE PASS
