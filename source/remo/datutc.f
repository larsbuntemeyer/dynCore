      SUBROUTINE DATUTC(NZTNR,YADAT,DT, YAKDAT1, YAKDAT2, NAKJATA, AKHH)
      !
      IMPLICIT NONE
      !
      !C
      !C**** DATUTCM  -   UP:BERECHNUNG DES AKTUELLEN TERMINS
      !C**   AUFRUF   :   CALL DATUTC(NZTNR,YADAT,DT,
      !C**                            YAKDAT1,YAKDAT2,NAKJATA,AKHH)
      !C**                IN REMORG
      !C**   ENTRIES  :   KEINE
      !C**   ZWECK    :   BERECHNUNG DES AKT. TERMINS AUS YADAT, DT, NZTNR; DER
      !C**                AKT. TERMIN WIRD IN DREI VERSCHIEDENEN FORMEN DARGE-
      !C**                STELLT. (SOWOHL FUER 30 TAGE MONATE WIE IN ECHAM,
      !C**                ALS AUCH FUER KALENDERMONATE WIE BEI ANALYSEN)
      !C**   VERSIONS-
      !C**   DATUM    :   30.01.89
      !C**                08.03.96
      !C**
      !C**   EXTERNALS:   KEINE
      !C**
      !C**   EINGABE-
      !C**   PARAMETER:   NZTNR:   ZEITSCHRITTNUMMER
      !C**                YADAT: ANFANGSDATUM DER VORHERSAGE IN DER FORM:
      !C**                      'DDMMYYHH' (TAG, MONAT, JAHR, TERMIN)
      !C**                DT:    ZEITSCHRITT (IN S)
      !C**   AUSGABE-
      !C**   PARAMETER:   YAKDAT1: AKTUELLER TERMIN IN DER FORM:
      !C**                      'DDMMYYHH' (TAG, MONAT, JAHR, TERMIN)
      !C**                YAKDAT2: AKTUELLER TERMIN IN DER FORM:
      !C**                      'WD  DD.MM.YY  HH UTC (WD: WOCHENTAG)
      !C**                NAKJATA: JAHRESTAG DES AKTUELLEN TERMINS
      !C**                AKHH:    STUNDE DES AKTUELLEN TERMINS
      !C**
      !C**   COMMON-
      !C**   BLOECKE  :   ORG
      !C**
      !C**   METHODE  :   BERECHNUNG UNTER BERUECKSICHTIGUNG VON SCHALTJAHREN
      !C**   FEHLERBE-
      !C**   HANDLUNG :   KEINE
      !C**   VERFASSER:   I.JACOBSEN/ D.MAJEWSKI/ R.PODZUN
      !C
      INCLUDE "comdia.h"
      !
      ! Formal Parameters
      !
      INTEGER,   INTENT(IN)    :: NZTNR 
      CHARACTER, INTENT(IN)    :: YADAT*10
      REAL,      INTENT(IN)    :: DT
      CHARACTER, INTENT(OUT)   :: YAKDAT1*10,YAKDAT2*21
      INTEGER,   INTENT(OUT)   :: NAKJATA
      REAL,      INTENT(OUT)   :: AKHH
      !
      ! Local Variables
      !
      INTEGER         :: JJJ,MMM,ITA,IHANF,JJ,M,MM
      INTEGER         :: I400,I100,I4,ITAG,IH
      INTEGER         :: ISCHLT,MOSTU
      INTEGER         :: JASTU,IWO
      INTEGER(KIND=8) :: NADDH, IDT, INZT
C
      INTEGER   :: MONATA(12), MONATM(12), MONASU(13)
      CHARACTER :: YWO(7)*2
      LOGICAL   :: LJASTU
C
      DATA         MONATM / 30 ,  30 ,  30 ,  30 ,  30 ,  30 ,
     &                      30 ,  30 ,  30 ,  30 ,  30 ,  30 /
      DATA         MONATA / 31 ,  28 ,  31 ,  30 ,  31 ,  30 ,
     &                      31 ,  31 ,  30 ,  31 ,  30 ,  31 /
      DATA         YWO    /'MO', 'DI', 'MI', 'DO', 'FR', 'SA', 'SO' /
C
      IDT=NINT(DT)
      INZT=NZTNR
C
C     YADAT ZERLEGEN IN TAG, MONAT, JAHR UND STUNDE; MONATSSUMMEN
C     MONATSSUMMEN BERECHNEN
      READ ( YADAT , '(I4,3I2)' ) JJJ, MMM, ITA , IHANF
C
      IF (LMOMON) THEN
C
         MONASU(1) =  0
         DO M =  2 , 13
            MONASU(M) =  MONASU(M-1) + MONATM(M-1)
         ENDDO
C
C        AKHH IM JAHR BERECHNEN
         JASTU = (ITA*24)+MONASU(MMM)*24+(IHANF-24)
         NADDH = (INZT*IDT)/3600
         JASTU = JASTU + INT(NADDH)
         JJ=JJJ
C
         LJASTU = .TRUE.
         DO WHILE (LJASTU)
            LJASTU = .FALSE.
C
C        JAHRESUEBERSCHREITUNG BERUECKSICHTIGEN
            IF (JASTU.LT.0) THEN
               JJ     = JJ-1
               JASTU  = 8640+JASTU
            ELSE IF (JASTU.GE.8640) THEN
               JASTU  = JASTU-8640
               JJ     = JJ+1
               IF (JASTU.GE.8640) THEN
                  LJASTU = .TRUE.
               ENDIF
            ELSE
               JJ     = JJ
            ENDIF
         ENDDO
C
      ELSE
C
C        SCHALTJAHR BESTIMMEN
         ISCHLT=0
         I400=MOD(JJJ,400)
         I100=MOD(JJJ,100)
         I4=MOD(JJJ,4)
         IF (I4.EQ.0) ISCHLT=1
         IF (I100.EQ.0) ISCHLT=0
         IF (I400.EQ.0) ISCHLT=1
C
         MONATA(2) =  28+ISCHLT
         MONASU(1) =  0
         DO M =  2 , 13
            MONASU(M) =  MONASU(M-1) + MONATA(M-1)
         ENDDO
C
C        AKHH IM JAHR BERECHNEN
         JASTU = (ITA*24)+MONASU(MMM)*24+(IHANF-24)
         NADDH = (INZT*IDT)/3600
         JASTU = JASTU + INT(NADDH)
         JJ=JJJ
C
         LJASTU = .TRUE.
         DO WHILE (LJASTU)
            LJASTU = .FALSE.
C
C           JAHRESUEBERSCHREITUNG BERUECKSICHTIGEN
            IF (JASTU.LT.0) THEN
               JJ     = JJ-1
C
C              SCHALTJAHR BESTIMMEN
               ISCHLT=0
               I400=MOD(JJ,400)
               I100=MOD(JJ,100)
               I4=MOD(JJ,4)
               IF (I4.EQ.0) ISCHLT=1
               IF (I100.EQ.0) ISCHLT=0
               IF (I400.EQ.0) ISCHLT=1
               JASTU  = 8760+ISCHLT*24+JASTU
            ELSE IF (JASTU.GE.(8760+ISCHLT*24)) THEN
               JASTU  = JASTU-(8760+ISCHLT*24)
               JJ     = JJ+1
C
C              SCHALTJAHR BESTIMMEN
               ISCHLT=0
               I400=MOD(JJ,400)
               I100=MOD(JJ,100)
               I4=MOD(JJ,4)
               IF (I4.EQ.0) ISCHLT=1
               IF (I100.EQ.0) ISCHLT=0
               IF (I400.EQ.0) ISCHLT=1
               IF (JASTU.GE.(8760+ISCHLT*24)) THEN
                  LJASTU = .TRUE.
               ENDIF
            ELSE
               JJ     = JJ
            ENDIF
         ENDDO !DO WHILE (LJASTU)
C
C        SCHALTJAHR BESTIMMEN
         ISCHLT=0
         I400=MOD(JJ,400)
         I100=MOD(JJ,100)
         I4=MOD(JJ,4)
         IF (I4.EQ.0) ISCHLT=1
         IF (I100.EQ.0) ISCHLT=0
         IF (I400.EQ.0) ISCHLT=1
C
         MONATA(2) =  28+ISCHLT
         MONASU(1) =  0
         DO M =  2 , 13
            MONASU(M) =  MONASU(M-1) + MONATA(M-1)
         ENDDO
C
      ENDIF
C
C     AUS JAHRESSTUNDE WIEDER DATUM BERECHNEN
      DO M = 2,13
         MOSTU   = JASTU-MONASU(M)*24
         IF (MOSTU.LT.0) EXIT
      ENDDO

      MM      = M-1
      MOSTU   = JASTU-MONASU(MM)*24
      ITAG    = MOSTU/24+1
      IH      = MOD(MOSTU,24)
      AKHH    = FLOAT(IH)+DT/3600.*MOD(NZTNR,IFIX(3600./DT+0.01))+0.0001
      IH      = IFIX(AKHH)
      NAKJATA = MONASU(MM)+ITAG+INT(AKHH/24.+0.0001)
      IWO     = 1
C
      WRITE ( YAKDAT1(1:4) , '(I4.4)' ) JJ
      WRITE ( YAKDAT1(5:6) , '(I2.2)' ) MM
      WRITE ( YAKDAT1(7:8) , '(I2.2)' ) ITAG
      WRITE ( YAKDAT1(9:10) , '(I2.2)' ) IH
      YAKDAT2 = YWO(IWO)//' '//YAKDAT1(7:8)//'.'// YAKDAT1(5:6)//'.'//
     &                         YAKDAT1(1:4)//'  '//YAKDAT1(9:10)//' UTC'
C
      RETURN
      !
      END SUBROUTINE DATUTC
