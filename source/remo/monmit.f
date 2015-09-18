      SUBROUTINE MONMIT(IPDB,IGDB,REF,ZAK,ZBK,SFELD,
     &                  KEMAX,YYYAKT,NHDMXN,NZT,IFEL,IEJE)
C
      IMPLICIT NONE
C
      INCLUDE "corg.h"
      INCLUDE "sumfel.h"
      INCLUDE "comdia.h"
      INCLUDE "unitch.h"
      INCLUDE "unitnr.h"
C
C     EINGABEFELDER UND -VARAIBLEN
C
      INTEGER,   INTENT(IN) :: IPDB(37), IGDB(22)
      INTEGER,   INTENT(IN) :: KEMAX,NHDMXN,NZT,IFEL,IEJE
      DOUBLE PRECISION, INTENT(IN) :: REF,ZAK(KEMAX),ZBK(KEMAX),
     &                                SFELD(IEJE)
      CHARACTER, INTENT(IN) :: YYYAKT*10
C
C     AUSGABEFELDER UND -VARAIBLEN
C
      REAL :: SUMS(IEJE,IFEL), SUM2S(IEJE,IFEL)
      REAL :: AK(KEMAX), BK(KEMAX)
C
C     HILFSFELDER UND VARIABLEN
C
      DOUBLE PRECISION :: VAR, F1, F2
      REAL             :: FIN(IEJE)
      INTEGER          :: MONAT(12)
      LOGICAL          :: LMONANF, LNEUTER
      INTEGER :: I, I100, I4, I400, IAKTDAT, ICODE, IHAKT, IHANF, II,  
     &           IJ, ISCHLT, ISTADAT, ITA1, ITA2, IW, JJJ1, JJJ2, KC,
     &           IIV, KL
      INTEGER :: MMM1, MMM2
      REAL    :: RREF
C
C     HILFSGROSSEN BELEGEN
C
      MONAT( 1)=31
      MONAT( 2)=28
      MONAT( 3)=31
      MONAT( 4)=30
      MONAT( 5)=31
      MONAT( 6)=30
      MONAT( 7)=31
      MONAT( 8)=31
      MONAT( 9)=30
      MONAT(10)=31
      MONAT(11)=30
      MONAT(12)=31
C
      IIV=NHDMXN
      RREF=REAL(REF)
      LMONANF=.FALSE.
      LNEUTER=.FALSE.
C
C     AM MONATSANFANG FESTSTELLEN
C
      READ(YANDAT,'(I4,3I2)') JJJ1, MMM1, ITA1 , IHANF
      READ(YYYAKT,'(I4,3I2)') JJJ2, MMM2, ITA2 , IHAKT
      ISTADAT=JJJ1*1000000+MMM1*10000+ITA1*100+IHANF
      IAKTDAT=JJJ2*1000000+MMM2*10000+ITA2*100+IHAKT
      LMONANF=((ISTADAT+IIV).EQ.IAKTDAT).AND.(IFIRST.EQ.1)
C
      IF (LMONANF) THEN
C
         IF (LDIA) PRINT *,'SUBROUTINE MONMIT:'
         IF (LDIA) PRINT *,'YANDAT=',YANDAT
         IF (LDIA) PRINT *,'ISTADAT=',ISTADAT,' IAKTDAT=',IAKTDAT
C
         IFIRST=0
C
C        ZU SICHERNDE VARIABLEN INITIALISIEREN
C
         IAKTTER=IAKTDAT
         IJAHR=IPDB(11)
         ISJAHR=IPDB(11)
         IMON=IPDB(12)
         ISMON=IPDB(12)
         ITAG=IPDB(13)
         ISTD=IPDB(14)
         IF (LMOMON) THEN
            IMONDAY=30
         ELSE
            IMONDAY=MONAT(IPDB(12))
            IF (IPDB(12).EQ.2) THEN
               ISCHLT=0
               I4=MOD(IPDB(11),4)
               I100=MOD(IPDB(11),100)
               I400=MOD(IPDB(11),400)
               IF (I4.EQ.0) ISCHLT=1
               IF (I100.EQ.0) ISCHLT=0
               IF (I400.EQ.0) ISCHLT=1
               IMONDAY=IMONDAY+ISCHLT
            ENDIF
         ENDIF
C
         IT=1
         ID=0
         IC=0
         IM=0
         IN=0
C
      ENDIF
C
C     FESTSTELLEN OB NEUER TERMIN
C
      LNEUTER=IAKTDAT.NE.IAKTTER
C
      IF (LNEUTER) THEN
         IF (LDIA) PRINT *,'IAKTTER=',IAKTTER,' IAKTDAT=',IAKTDAT
         IAKTTER=IAKTDAT
         IF (IIV.EQ.ISTD) THEN
            DO IW=1,4
               DO IJ=1,IEJE
                  DSUMX(IJ,IW)=-1.E-34
               ENDDO
            ENDDO
            DO IW=1,2
               DO IJ=1,IEJE
                  DSUMN(IJ,IW)=1.E34
               ENDDO
            ENDDO
         ENDIF
C
C         ZAEHLER FUER DIE TERMINE EINES MONATS HOCHZAEHLEN
C
         IT=IT+1
         IM=0
         IN=0
         IC=0
      ENDIF
C
C
C     ZAEHLER FUER DIE FELDER EINES TERMINS HOCHZAEHLEN
C
      IC=IC+1
C
C     UMSPEICHERN DER EINGABEFELDER
C
      DO II=1,37
         JPDB(II,IC)=IPDB(II)
      ENDDO
      DO II=1,22
         JGDB(II,IC)=IGDB(II)
      ENDDO
      DO KL=1,KEMAX
         AK(KL)=REAL(ZAK(KL))
         BK(KL)=REAL(ZBK(KL))
      ENDDO
      DO IJ=1,IEJE
         FIN(IJ)=REAL(SFELD(IJ))
      ENDDO
C
      ICODE=JPDB(7,IC)
C
      IF (ICODE.EQ.201.OR.ICODE.EQ.214.OR.
     &    ICODE.EQ.216.OR.ICODE.EQ.217) THEN
         IM=IM+1
         DO IJ=1,IEJE
            DSUMX(IJ,IM)=AMAX1(FIN(IJ),DSUMX(IJ,IM))
         ENDDO
         IF (ISTD.EQ.0) THEN
            DO IJ=1,IEJE
               SUM(IJ,IC)=SUM(IJ,IC)+DSUMX(IJ,IM)
               SUM2(IJ,IC)=SUM2(IJ,IC)+DSUMX(IJ,IM)*DSUMX(IJ,IM)
            ENDDO
         ENDIF
C
      ELSE IF (ICODE.EQ.202.OR.ICODE.EQ.215) THEN
         IN=IN+1
         DO IJ=1,IEJE
            DSUMN(IJ,IN)=AMIN1(FIN(IJ),DSUMN(IJ,IN))
         ENDDO
         IF (ISTD.EQ.0) THEN
            DO IJ=1,IEJE
               SUM(IJ,IC)=SUM(IJ,IC)+DSUMN(IJ,IN)
               SUM2(IJ,IC)=SUM2(IJ,IC)+DSUMN(IJ,IN)*DSUMN(IJ,IN)
            ENDDO
         ENDIF
C
      ELSE
         DO IJ=1,IEJE
            SUM(IJ,IC)=SUM(IJ,IC)+FIN(IJ)
            SUM2(IJ,IC)=SUM2(IJ,IC)+FIN(IJ)*FIN(IJ)
         ENDDO
C
      ENDIF
C
C     SOLANGE NICHT ALLE FELDER EINES TERMINS VERARBEITET
C     SIND: 'RETURN', SONST WEITER UND 'RETURN' ODER
C     MONATSENDE
C
      IF (IC.NE.IFEL) THEN
         RETURN
      ELSE
         IAKTTER=IAKTDAT
         ISTD=ISTD+IIV
         IF (ISTD.EQ.24) THEN
            ISTD=0
            ITAG=ITAG+1
            ID=ID+1
            IF (ITAG.EQ.IMONDAY+1) THEN
               ITAG=1
               IMON=IMON+1
               IF (IMON.GT.12) THEN
                  IMON=1
                  IJAHR=IJAHR+1
               ENDIF
            ENDIF
         ENDIF
C
C        WENN NICHT MONATSENDE, BEENDE ROUTINE
C
         IF (ISTD.GT.0.AND.IMON.NE.ISMON) THEN
            ID=ID-1
         ELSE
            RETURN
         ENDIF
      ENDIF
C
      IF (LDIA) PRINT *,'IT=',IT,' IMONDAY=',IMONDAY,' ID=',ID
C
C     MITTELWERTE UND STANDARDABWEICHUNGEN BERECHNEN
C
      F1=0.
      F2=1.
      DO KC=1,IFEL
C
         ICODE=JPDB(7,KC)
C
         IF (ICODE.EQ.201.OR.ICODE.EQ.202.OR.
     &       ICODE.EQ.214.OR.ICODE.EQ.215.OR.
     &       ICODE.EQ.216.OR.ICODE.EQ.217) THEN
            DO IJ=1,IEJE
               VAR=DMAX1(F1,(F2/DBLE(IMONDAY-1))*(SUM2(IJ,KC)-
     &              ((F2/DBLE(IMONDAY))*SUM(IJ,KC)*SUM(IJ,KC))))
               SUM(IJ,KC)=SUM(IJ,KC)/DBLE(IMONDAY)
               SUM2(IJ,KC)=DSQRT(VAR)
               IF (DABS(SUM(IJ,KC)).LT.1.E-34) SUM(IJ,KC)=0.0
               IF (DABS(SUM2(IJ,KC)).LT.1.E-34) SUM2(IJ,KC)=0.0
               SUMS(IJ,KC)=REAL(SUM(IJ,KC))
               SUM2S(IJ,KC)=REAL(SUM2(IJ,KC))
            ENDDO
         ELSE
C
            DO IJ=1,IEJE
               VAR=DMAX1(F1,(F2/DBLE(IT-1))*(SUM2(IJ,KC)-
     &              ((F2/DBLE(IT))*SUM(IJ,KC)*SUM(IJ,KC))))
               SUM(IJ,KC)=SUM(IJ,KC)/DBLE(IT)
               SUM2(IJ,KC)=DSQRT(VAR)
               IF (DABS(SUM(IJ,KC)).LT.1.E-34) SUM(IJ,KC)=0.0
               IF (DABS(SUM2(IJ,KC)).LT.1.E-34) SUM2(IJ,KC)=0.0
               SUMS(IJ,KC)=REAL(SUM(IJ,KC))
               SUM2S(IJ,KC)=REAL(SUM2(IJ,KC))
            ENDDO
         ENDIF
C
      ENDDO ! KC=1,IFEL
C
C     DATEINAMEN ERZEUGEN UND AUSGABEDATEIEN OEFFNEN
C
      CALL MAKEPN('M', NZT + 1)
      WRITE(YMDNAM(9:12),'(I4.4)') ISJAHR
      WRITE(YMDNAM(13:14),'(I2.2)') ISMON
      CALL SEND2(NUMDAT, YMDNAM, YMDCAT)
      CALL MAKEPN('S', NZT + 1)
      WRITE(YSDNAM(9:12),'(I4.4)') ISJAHR
      WRITE(YSDNAM(13:14),'(I2.2)') ISMON
      CALL SEND2(NUSDAT, YSDNAM, YMDCAT)
C
C     FELDER AUSGEBEN
C
      DO KC=1,IFEL
         JPDB(11,KC)=ISJAHR
         JPDB(12,KC)=ISMON
         JPDB(13,KC)=0
         JPDB(14,KC)=0
         JPDB(20,KC)=1
         WRITE(NUMDAT) (JPDB(I,KC),I=1,37),(JGDB(I,KC),I=1,18),RREF,
     &        (JGDB(I,KC),I=20,22),AK,BK
         WRITE(NUMDAT) (SUMS(IJ,KC),IJ=1,IEJE)
         JPDB(20,KC)=2
         WRITE(NUSDAT) (JPDB(I,KC),I=1,37),(JGDB(I,KC),I=1,18),RREF,
     &        (JGDB(I,KC),I=20,22),AK,BK
         WRITE(NUSDAT) (SUM2S(IJ,KC),IJ=1,IEJE)
      ENDDO
C
C     NEUEN MONAT SETZEN
C
      IFIRST=1
C
C     AUSGABEDATEIEN SCHLIESSEN
C
      CLOSE(NUMDAT)
      CLOSE(NUSDAT)
      RETURN
      END SUBROUTINE MONMIT
