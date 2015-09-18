      SUBROUTINE INIT
     &  (KE     , IEJE  , KE1    ,
     &   NZT    , NANF  , NRADIA , NVDIFF, CEVAPCU, VCT    , FAO    ,
     &   RGCGN  , CVDAES, CVDAEL , CVDAEU, CVDAED , NRADFR , LGADSRH,
        !
        !   WARNING: LRAD differs from actual parameter name when
        !   called in phyec: it's called with                LECRAD
        !                                                     |  
     &   NRADPFR, LSOLC , LAER   , LCFC  , EPS    , LDIUR  , LRAD   ,
     &   LVDIFF , LSURF , LGWDRAG, LCONV , LCOND  , NCBASE , NCDATA ,
     &   NTBASE , NTDATA, NMONTH , NTIMST, TLAMBDA, DLAMBDA, PORVOL ,
     &   FCAP)
        !   All other dummy names correspond to their actual parameters
C
      IMPLICIT NONE
C
      INCLUDE "comdyn.h"
      INCLUDE "comdia.h"
      INCLUDE "corg.h"
      INCLUDE "phykon.h"
      INCLUDE "comecphy.h"
C
      INTEGER, INTENT(IN)  :: KE     , IEJE  , KE1    ,
     &                        NZT    , NANF
      INTEGER, INTENT(OUT) :: NRADIA , NVDIFF, NRADFR, NRADPFR,
     &                        NCBASE , NCDATA ,
     &                        NTBASE , NTDATA, NMONTH , NTIMST
C
      REAL, INTENT(IN)     :: VCT(2*KE1)
      REAL, INTENT(INOUT)  :: CEVAPCU(KE)
      REAL, INTENT(INOUT)  :: CVDAES(KE1),CVDAEL(KE1),
     &                        CVDAEU(KE1),CVDAED(KE1)
      REAL, INTENT(IN)     :: FAO(IEJE)
      REAL, INTENT(INOUT)  :: RGCGN(IEJE)
      REAL, INTENT(INOUT)  :: TLAMBDA(IEJE),DLAMBDA(IEJE),PORVOL(IEJE)
      REAL, INTENT(INOUT)  :: FCAP(IEJE)
      REAL,    INTENT(OUT) :: EPS
c
      LOGICAL, INTENT(OUT) :: LDIUR,LRAD,LSOLC,LAER,LCFC
      LOGICAL, INTENT(OUT) :: LVDIFF,LSURF,LGWDRAG,LCONV,LCOND
      LOGICAL, INTENT(OUT) :: LGADSRH
C
      REAL    :: CETA(KE),CETAH(KE1)
      REAL    :: APZERO, CO2FAC, ZDT, ZETAM, ZETAP
      INTEGER :: I, IDA, IDAY, IFDAT, IFTIM, IHA, IHOUR, IMA, IMONTH, 
     &           IVDAT, IVTIM, IYA, IYEAR, JLEV, KYMD2C, MYMD2C
C
      READ(YAKDAT1,'(I4,3I2)') IYEAR,IMONTH,IDAY,IHOUR
      READ(YADAT,'(I4,3I2)') IYA,IMA,IDA,IHA
      NRADIA=0
      NVDIFF=0
      IF (NZT.EQ.0) THEN
         ZDT=2.*DT
      ELSE
         ZDT=DT
      ENDIF
      NRADFR=NINT(HDRAD*3600./ZDT)
      NRADPFR=10
      CO2FAC=1.
C     SCHALTER DER ECHAMPHYSIK MIT DEN ENTSPRECHENDEN
C     NAMELIST-VARIABLEN BELEGEN
C     LSOLC=.TRUE. FUER CLEAR SKY DIAGNOSTICS
      LSOLC=LECSOL
      LAER=LECAER
      LCFC=LECCFC
      LDIUR=LECDIU
      LRAD=LECRAD
      LVDIFF=LECVDI
      LSURF=LECSUR
      LGWDRAG=LECGWD
      LCONV=LECCOV
      LCOND=LECCOD
      LGADSRH=LECGAD
      EPS=0.05
C
      IFDAT=IYA*10000+IMA*100+IDA
      IFTIM=IHA*10000
      IVDAT=IYEAR*10000+IMONTH*100+IDAY
      IVTIM=IHOUR*10000
      IF (LMOMON) THEN
         NCBASE=MYMD2C(IFDAT)
         NCDATA=MYMD2C(IVDAT)
      ELSE
         NCBASE=KYMD2C(IFDAT)
         NCDATA=KYMD2C(IVDAT)
      ENDIF
CRP   WEGEN ANSCHEINEND FEHLENDEM TAG IM VERGLEICH ZUM DWD-TEIL WIRD
C     EIN TAG AUF NCBASE AUFADDIERT.
      NCBASE=NCBASE+1
      NTBASE=((IFTIM/10000)*60+(MOD(IFTIM,10000)/100))*60+MOD(IFTIM,100)
      NTDATA=((IVTIM/10000)*60+(MOD(IVTIM,10000)/100))*60+MOD(IVTIM,100)
      NTIMST=1
C     DEFAULT: NMONTH=0, PERPETUAL-RUNS NMONTH=IMONTH (LT. U.SCHLESE)
      NMONTH=0
C
      IF ((NZT.EQ.0).OR.((NZT.EQ.NANF+1).AND.(NANF.NE.0))) THEN
C
C     INITIALISE COMMON-BLOCK COMCON
C
         CALL INICON(R,RD,WCP,WLK,WLF,WLS,G,RERD,RHF,SIGMA,SOKO)
C
         CALL INILUC
C
C
         APZERO=101325.
         ZETAM=VCT(1)/APZERO+VCT(KE1+1)
         CETAH(1)=ZETAM
C
         DO JLEV=1,KE
            ZETAP=VCT(JLEV+1)/APZERO+VCT(KE1+1+JLEV)
            CETA(JLEV)=(ZETAM+ZETAP)*.5
            CETAH(JLEV+1)=ZETAP
            ZETAM=ZETAP
         ENDDO
C
         CALL INIPHY(CETA,CEVAPCU,KE)
C
         CALL INIRAD(CETAH,CO2FAC,CVDAES,CVDAEL,CVDAEU,CVDAED,KE1)
C
C     INITIALISE RGCGN,TLAMBDA,DLAMBDA,PORVOL,FCAP
C
CSH   *** CHANGES MADE TO USE SOIL MOISTURE INDEPENDENT DIFFUSIVITY/CAPACITY
C     *** SODIF VALUES PUT TO TLAMBDA
C
         IF (.NOT.LWDIF) THEN
C
            DO I=1,IEJE
               IF ((FAO(I).GT.0.5).AND.(FAO(I).LE.1.5)) THEN
                  RGCGN(I)=1.93E+06
                  TLAMBDA(I)=8.7E-7
                  DLAMBDA(I)=2.40
                  PORVOL(I)=0.364
                  FCAP(I)=0.196
               ELSE IF ((FAO(I).GT.1.5).AND.(FAO(I).LE.2.5)) THEN
                  RGCGN(I)=2.10E+06
                  TLAMBDA(I)=8.0E-7
                  DLAMBDA(I)=2.40
                  PORVOL(I)=0.445
                  FCAP(I)=0.260
               ELSE IF ((FAO(I).GT.2.5).AND.(FAO(I).LE.3.5)) THEN
                  RGCGN(I)=2.25E+06
                  TLAMBDA(I)=7.4E-7
                  DLAMBDA(I)=1.58
                  PORVOL(I)=0.455
                  FCAP(I)=0.340
               ELSE IF ((FAO(I).GT.3.5).AND.(FAO(I).LE.4.5)) THEN
                  RGCGN(I)=2.36E+06
                  TLAMBDA(I)=7.1E-7
                  DLAMBDA(I)=1.55
                  PORVOL(I)=0.475
                  FCAP(I)=0.370
               ELSE IF ((FAO(I).GT.4.5).AND.(FAO(I).LE.5.5)) THEN
                  RGCGN(I)=2.48E+06
                  TLAMBDA(I)=6.7E-7
                  DLAMBDA(I)=1.50
                  PORVOL(I)=0.507
                  FCAP(I)=0.463
               ELSE IF ((FAO(I).GT.5.5).AND.(FAO(I).LE.6.5)) THEN
                  RGCGN(I)=2.59E+06
                  TLAMBDA(I)=6.5E-7
                  DLAMBDA(I)=0.50
                  PORVOL(I)=0.863
                  FCAP(I)=0.763
               ELSE
                  RGCGN(I)=2.25E+06
                  TLAMBDA(I)=7.4E-7
                  DLAMBDA(I)=1.58
                  PORVOL(I)=0.455
                  FCAP(I)=0.340
               ENDIF
            ENDDO
C
         ELSE
C
            DO I=1,IEJE
               IF ((FAO(I).GT.0.5).AND.(FAO(I).LE.1.5)) THEN
                  RGCGN(I)=1.28E+06
                  TLAMBDA(I)=0.30
                  DLAMBDA(I)=2.40
                  PORVOL(I)=0.364
                  FCAP(I)=0.196
               ELSE IF ((FAO(I).GT.1.5).AND.(FAO(I).LE.2.5)) THEN
                  RGCGN(I)=1.35E+06
                  TLAMBDA(I)=0.28
                  DLAMBDA(I)=2.40
                  PORVOL(I)=0.445
                  FCAP(I)=0.260
               ELSE IF ((FAO(I).GT.2.5).AND.(FAO(I).LE.3.5)) THEN
                  RGCGN(I)=1.42E+06
                  TLAMBDA(I)=0.25
                  DLAMBDA(I)=1.58
                  PORVOL(I)=0.455
                  FCAP(I)=0.340
               ELSE IF ((FAO(I).GT.3.5).AND.(FAO(I).LE.4.5)) THEN
                  RGCGN(I)=1.50E+06
                  TLAMBDA(I)=0.21
                  DLAMBDA(I)=1.55
                  PORVOL(I)=0.475
                  FCAP(I)=0.370
               ELSE IF ((FAO(I).GT.4.5).AND.(FAO(I).LE.5.5)) THEN
                  RGCGN(I)=1.63E+06
                  TLAMBDA(I)=0.18
                  DLAMBDA(I)=1.50
                  PORVOL(I)=0.507
                  FCAP(I)=0.463
               ELSE IF ((FAO(I).GT.5.5).AND.(FAO(I).LE.6.5)) THEN
                  RGCGN(I)=5.8E+05
                  TLAMBDA(I)=0.06
                  DLAMBDA(I)=0.50
                  PORVOL(I)=0.863
                  FCAP(I)=0.763
               ELSE
                  RGCGN(I)=1.42E+06
                  TLAMBDA(I)=0.25
                  DLAMBDA(I)=1.58
                  PORVOL(I)=0.455
                  FCAP(I)=0.340
               ENDIF
            ENDDO
C
         ENDIF
C
      ENDIF
C
      RETURN
      END SUBROUTINE INIT
