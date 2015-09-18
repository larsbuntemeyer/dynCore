C
C     SUBROUTINE FFT
C
C@(#) BIBLIOTHEK EM: EUROPA-MODELL / DEUTSCHLAND-MODELL
C@(#) MODUL FFT.F, V2.12 VOM 4/27/94, EXTRAHIERT AM 4/28/94
C
C**** FFT      -   UP:BERECHNUNG DER FAST FOURIER TRANSFORM (HIN- UND
C****                 RUECKTRANSFORMATION), VERSION IN FORTRAN77
C**   AUFRUF   :   CALL FFT(A,B,W,M,N,TRIGS,MT,MODE,Z1,Z2,IFAX)
C**
C**   ENTRIES  :   KEINE
C**   ZWECK    :   BERECHNUNG DER FAST FOURIER TRANSFORM (HIN- UND
C**                RUECKTRANSFORMATION)
C**   VERSIONS-
C**   DATUM    :   13.03.89
C**
C**   EXTERNALS:   PASS
C**
C**   EINGABE-
C**   PARAMETER:   A(M,N) : EINGABE-FELD (DATEN / FOURIERKOEFFIZIENTEN)
C**                MODE = +1 INPUT=DATEN, OUTPUT = FOURIERKOEFFIZIENTEN
C**                MODE = -1 INPUT=FOURIERKOEFFIZIENTEN, OUTPUT = DATEN
C**                IFAX(1):  ANZAHL DER FAKTOREN
C**                IFAX(2) - IFAX(IFAX(1)+1) : FAKTOREN
C**   AUSGABE-
C**   PARAMETER:   B(N,M) : AUSGABE-FELD (FOURIERKOEFFIZIENTEN / DATEN)
C**
C**   COMMON-
C**   BLOECKE  :   KEINE
C**
C**   METHODE  :   FFT NACH TEMPERTON
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   P.PROHL
C
      SUBROUTINE FFT (A,B,W,M,N,TRIGS,MT,MODE,Z1,Z2,IFAX)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN)    :: M,N,MT,MODE
      REAL,    INTENT(IN)    :: A(M,N),TRIGS(MT)
      REAL,    INTENT(INOUT) :: B(M,N),W(M,N)
      REAL,    INTENT(IN)    :: Z1(6),Z2(6)
      INTEGER, INTENT(IN)    :: IFAX(11)
C
      INTEGER :: N1,MP2,MP1,MH2,MH1,MH,MD,M4,M2,M1,NFAX,
     &           K,J,IGO,LA,IFAC,ID2,I
      REAL    :: ZWSP2,SI1,SCALE,SCA4,ZWSP1,CO1,A3,A2,A1,A0
C
      M1  = M - 1
      M2  = M - 2
      M4  = M - 4
      MP1 = M + 1
      MP2 = M + 2
      MH  = M1 / 2
      MH1 = MH + 1
      MH2 = MH + 2
      MD  = M1 * 2

      N1  = N - 1
      ID2 = 2

C*****ERKLAERUNG PRE- UND POSTPROCESSING SIEHE COOLEY, LEWIS UND WELCH
C     PLUS ERGAENZUNGEN
C     PREPROCESSING

      SCALE = 0.25
      IF (MODE .GT. 0)  SCALE = 0.25 / FLOAT (MH)
      SCA4 = 4. * SCALE
      DO J = 2 , N1
         W(1,J)   = 0.
         W(M,J)   = 0.
      ENDDO

      DO J = 2 , N1
         W(MH1,J) = SCA4 * A(MH1,J)
      ENDDO

      DO I = 2, MH
         SI1 = TRIGS(MD+I)
         DO J = 2 , N1
            A0 = (A(I,J) + A(MP1-I,J)) * SI1
            A1 =  A(I,J) - A(MP1-I,J)
            W(I,J)     = SCALE * (A0 + A1)
            W(MP1-I,J) = SCALE * (A0 - A1)
         ENDDO
      ENDDO
C-----------------------------------------------------------------------

C     COMPLEX FFT
C     MH IN IFAX(1) FAKTOREN IFAX(2), ... ZERLEGT
C     LA = MH / (IFAX(K)*...*IFAX(IFAX(1)+1)), K = 2,..,IFAX(1)+1

      NFAX = IFAX(1) + 1
      LA = 1
      IGO = +1
      DO K = 2, NFAX
         IFAC = IFAX(K)
         IF (IGO .EQ. -1) THEN
            CALL PASS (B,W,M,N,TRIGS,MT,ID2,MH,IFAC,LA,Z1,Z2)
         ELSE
            CALL PASS (W,B,M,N,TRIGS,MT,ID2,MH,IFAC,LA,Z1,Z2)
         ENDIF
         LA = IFAC * LA
         IGO = -IGO
      ENDDO

C     EVT. UMSPEICHERN

      IF (MOD(IFAX(1),2) .NE. 1) THEN
         DO I = 1 , M
            DO J = 1 , N
               B(I,J) = W(I,J)
            ENDDO
         ENDDO
      ENDIF

C-----------------------------------------------------------------------

C     POSTPROCESSING, ERSTE STUFE

      DO J = 2 , N1
         W(1,J) = B(1,J) + B(2,J)
         W(2,J) = B(1,J) - B(2,J)
      ENDDO

CKS      IF (MOD(MH,2) .EQ. 1)  GO TO 125
      IF (MOD(MH,2) .NE. 1) THEN
         DO J = 2 , N1
            W(MH1,J) = 2. * B(MH1,J)
            W(MH2,J) = 2. * B(MH2,J)
         ENDDO
      ENDIF

CKS      IF (M1 .EQ. 4)  GO TO 140
      IF (MOD(MH,2) .EQ. 1 .OR. M1 .NE. 4) THEN

CKS  125 K = 2
         K = 2
         DO I = 3, MH, 2
            SI1  = TRIGS(K)
            CO1  = TRIGS(K+M1)
            K   = K + 1
            DO J = 2 , N1
               A0 = B(I,J)   + B(MP1-I,J)
               A1 = B(I,J)   - B(MP1-I,J)
               A2 = B(I+1,J) + B(MP2-I,J)
               A3 = B(I+1,J) - B(MP2-I,J)
               ZWSP1 = A1 * SI1 + A2 * CO1
               W(I,J)     = A0 + ZWSP1
               W(MP1-I,J) = A0 - ZWSP1
               ZWSP2 = A2 * SI1 - A1 * CO1
               W(I+1,J)   = ZWSP2 + A3
               W(MP2-I,J) = ZWSP2 - A3
            ENDDO
         ENDDO
      ENDIF

C     POSTPROCESSING, ZWEITE STUFE

CKS  140 DO J = 2 , N1
      DO J = 2 , N1
         B(1,J)  =  0.
         B(2,J)  =  W(1,J)
         B(M2,J) =  W(M1,J)
         B(M1,J) = -W(2,J)
      ENDDO

      IF (M1 .EQ. 4)  RETURN
      DO I = 3 , M4, 2
         DO J = 2 , N1
            B(I,J)   = W(I+1,J)
            B(I+1,J) = B(I-1,J) + W(I,J)
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE FFT
