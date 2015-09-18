C
C
C@(#) BIBLIOTHEK EM: EUROPA-MODELL / DEUTSCHLAND-MODELL
C@(#) MODUL PASSKO.F, V2.12 VOM 4/27/94, EXTRAHIERT AM 4/28/94
C
C**** PASSKO   -   UP:SETZEN VON TRIGONOMETRISCHEN FUNKTIONSWERTEN FUER
C****                 EINE KOMPLEXE FFT DER LAENGE (M-1)/2
C**                   VERSION IN FORTRAN77
C**   AUFRUF   :   CALL PASSKO (TRIGS,MT,M1,PI,Z1,Z2,IFAX)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   SETZEN VON TRIGONOMETRISCHEN FUNKTIONSWERTEN FUER
C**                EINE KOMPLEXE FFT DER LAENGE (M-1)/2
C**   VERSIONS-
C**   DATUM    :   05.04.89
C**
C**   EXTERNALS:   KEINE
C**
C**   EINGABE-
C**   PARAMETER:   MT : 5*(M-1)/2
C**
C**   AUSGABE-
C**   PARAMETER:   TRIGS(MT) MIT:
C**                TRIGS(1,...,M1)   = SIN (2*J*PI/M1)  , J=0,.,M1-1
C**                TRIGS(M1+1,.,2M1) = COS (2*J*PI/M1)  , J=0,.,M1-1
C**                TRIGS(2M1+1,.,MT) = 2 * SIN (J*PI/M1), J=0,.,(M1/2)-1
C**
C**   COMMON-
C**   BLOECKE  :   KEINE
C**
C**   METHODE  :   BERECHNUNG DER SINUS- UND COSINUS-TERME
C**   FEHLERBE-
C**   HANDLUNG :   STOP BEI FALSCHER ANZAHL VON GITTERPUNKTEN FUER FFT
C**   VERFASSER:   P.PROHL
C
      SUBROUTINE PASSKO (TRIGS,MT,M1,PI,Z1,Z2,IFAX)
C
      IMPLICIT NONE
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: MT,M1
      REAL,    INTENT(INOUT) :: TRIGS(MT)
      REAL,    INTENT(INOUT) :: Z1(6),Z2(6)
      INTEGER, INTENT(INOUT) :: IFAX(11)
C
C     Local Variables
C
      INTEGER :: JFAX(10),LFAX(8)
      INTEGER :: I, IFAC, J, MD, MH, MX, NFAX
      REAL    :: PI, PI2R, PI36, PI45, PI60, PIR
C
      DATA  LFAX   /8,6,5,4,3,2,1,0/
C
      MH   = M1 / 2
      MD   = M1 * 2
      PIR  = PI / FLOAT (M1)
      PI2R = 2. * PIR

      DO I = 1 , M1
         TRIGS(I)    = SIN (PI2R * FLOAT (I-1))
         TRIGS(I+M1) = COS (PI2R * FLOAT (I-1))
      ENDDO

      DO I = 1, MH
         TRIGS(I+MD) = 2. * SIN (PIR * FLOAT (I-1))
      ENDDO

C----------------------------------------------------------------------
C     FUER DIE FFT ZERLEGUNG VON (M-1)/2 IN FAKTOREN
C     (M-1)/2 = LFAX(I)*LFAX(J)*LFAX(K)..., FAKTOREN GROESSER 1
C     GLEICHE FAKTOREN KOENNEN MEHRFACH VERWENDUNG FINDEN (Z.B. 4*4)
C     IFAX(1) ANZAHL DER FAKTOREN
C     IFAX(2) BIS IFAX(IFAX(1)+1) FAKTOREN (WERTE ABNEHMEND)

      IF (MOD(M1,2) .NE. 0) THEN
         WRITE (2,100)
 100     FORMAT (' ', 'M-1 UNGERADE, STOP')
         CALL REMARK( 'PASSKO: M-1 UNGERADE' )
         STOP
      ENDIF

      MX = M1 / 2
      IFAC = LFAX(1)
      I = 0
      J = 1
      DO
         IF (MOD(MX,IFAC) .NE. 0) THEN
            J = J +1
            IFAC = LFAX(J)
            IF (IFAC .GT. 1) THEN
               CYCLE
            ELSE
               WRITE (2,140)
 140           FORMAT (' ', '(M-1)/2 ENTHAELT ILLEGALE FAKTOREN, STOP')
               CALL REMARK(
     &              'PASSKO: (M-1)/2 ENTHAELT ILLEGALE FAKTOREN')
               STOP
            ENDIF
         ELSE
            I = I + 1
            JFAX(I) = IFAC
            MX = MX / IFAC
            IF (MX .EQ. 1) EXIT
            CYCLE
         ENDIF
      ENDDO

      NFAX =I
      IFAX(1) = NFAX
      DO I = 1 , NFAX
         IFAX(I+1) = JFAX(I)
      ENDDO

C----------------------------------------------------------------------

C     KONSTANTEN FUER SUBROUTINE PASS

      PI60  = PI / 3.
      PI45  = PI / 4.
      PI36  = PI / 5.
      Z1(1) = COS (PI60)
      Z2(1) = SIN (PI60)
      Z1(2) = COS (PI36) - 0.25
      Z2(2) = SIN (PI36)
      Z1(3) = 0.25
      Z2(3) = SIN (2. * PI36)
      Z1(4) = COS (PI45)

      RETURN
      END SUBROUTINE PASSKO
