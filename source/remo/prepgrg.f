C
C
C**** *PREPGRG* - PREPARES GREENHOUSE GASES FOR MODEL USE
C
C      U. SCHLESE   DKRZ - HAMBURG   JUNE-95
C      R. PODZUN    MPI  - HAMBURG   APRIL-08
C
C      METHOD.
C      ------
C
C      1. INTERPOLATION IN TIME.
C      2. MODIFY VALUES TO GIVE THE RIGHT FORCING AND ACCOUNT FOR
C         USING THE 1990 VALUES IN THE CONTROL RUN INSTEAD OF 1860
C      3. CONVERT FROM VOLUME MIXING RATIO TO MASS MIXING RATIO.
C
C      THE RESULTS ARE STORED IN COMMON *YOMRDI* OVERWRITING
C      THE VALUES SET IN *SURADI* USED BY THE CONTROL RUN
C
C
C     INTERFACE.
C     ----------
C
C     *PREPGRG* IS CALLED FROM *STEPON* AT FULL RADIATION TIMESTEPS.
C     NOTE: IT MUST NOT BE CALLED IN A PARALLEL REGION!
C
      SUBROUTINE PREPGRG
C
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
      INCLUDE "org.h"
      INCLUDE "corg.h"
      INCLUDE "COMGTS"
      INCLUDE "YOMRDI"
      INCLUDE "COMCON"
      INCLUDE "comdyn.h"
      INCLUDE "comdia.h"
C
      INTEGER :: MONATA(11), MONASU(12)
      INTEGER :: I, I100, I4, I400, IDAY, IHOUR, IM, IMONTH, ISCHLT,  
     &           ISPD, ISPH, ISPM, IYEAR, IYR, IYRM1, IYRP1, JC, JJJ, 
     &           ISMD, M
      REAL    :: ZAIRMWG, ZAKSEC, ZCH4INT, ZCH4MWG, ZCO2INT, ZCO2MWG, 
     &           ZN2OMWG, ZSNZT, ZW1, ZW2, ZYEARL, ZYRSEC, ZN2OINT
      LOGICAL :: LOPRINT 
C
      DATA         MONATA / 31 ,  28 ,  31 ,  30 ,  31 ,  30 ,
     &                      31 ,  31 ,  30 ,  31 ,  30 /
C
C    CONVERSION FACTORS VOLUME TO MASS MIXING RATIO (FROM *SURADI*)
C
      ZAIRMWG=28.970
      ZCO2MWG=44.011
      ZCH4MWG=16.043
      ZN2OMWG=44.013
C
C     -----------------------------------------------------------------
C
C*      1.       INTERPOLATION IN TIME.
C                ------------- -- ----
C
      CALL DATUTC(NZT, YADAT, DT, YAKDAT1, YAKDAT2, NAKJATA,
     &            AKHH)
C
C     ACTUAL DATE OF THE INTEGRATION
C
      READ(YAKDAT1,'(I4,3I2)') IYEAR,IMONTH,IDAY,IHOUR
C
C     CHECK WHETHER THE TIMESERIES OF GRG-GASES IS LONG ENOUGH
C
      IF (MYID.EQ.0) THEN
         IF (JPGTS1.LT.IYEAR.AND.JPGTS2.LT.IYEAR) THEN
            PRINT *,'ERROR PREPGRG: WRONG GRG-TIME SERIES !'
            PRINT *,'CHECK GRG-GAS FILE !'
            STOP 1
         ENDIF
         IF (IYEAR.GE.JPGTS2.AND.IMONTH.GE.7) THEN
            PRINT *,'ERROR PREPGRG: GRG-TIME SERIES NOT LONG ENOUGH !'
            PRINT *,'CHECK JPGTS2 AND GRG-GAS FILE !'
            STOP 1
         ENDIF
      ENDIF
C
      IF (LMOMON) THEN
C
         MONASU(1) =  0
         DO M = 2 , 12
            MONASU(M) =  MONASU(M-1) + 30
         ENDDO
C
         ZYEARL=360.
C
      ELSE
C
C        ESTIMATION OF LEAPYEAR
C
         JJJ=IYEAR
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
         DO M = 2 , 12
            MONASU(M) =  MONASU(M-1) + MONATA(M-1)
         ENDDO
C
         ZYEARL=365.+FLOAT(ISCHLT)
C
      ENDIF
C
C     SET SOME NECESSARY CONSTANTS AND VARIABLES
C
      ISMD=MONASU(IMONTH)
      ISPH=3600
      ISPD=86400
      ISPM=ISMD*ISPD
C
C     YEARLENGTH IN SECONDS
C
      ZYRSEC=ZYEARL*DAYL
C
C     ACTUAL DATE IN SECONDS
C
      ZSNZT=MOD(FLOAT(NZT)*DT,3600.)
      ZAKSEC=FLOAT(ISPM+(IDAY-1)*ISPD+IHOUR*ISPH)+ZSNZT
C
      IYR=IYEAR
      IYRP1=IYR+1
      IYRM1=IYR-1
      IM=IMONTH
C
C       FIRST HALF OF YEAR
C
      IF (IM.LE.6) THEN
C
         ZW1=(ZAKSEC/ZYRSEC)+0.5
         ZW2=1.-ZW1
         ZCO2INT=ZW1*CO2TS(IYR)+ZW2*CO2TS(IYRM1)
         ZCH4INT=ZW1*CH4TS(IYR)+ZW2*CH4TS(IYRM1)
         ZN2OINT=ZW1*AN2OTS(IYR)+ZW2*AN2OTS(IYRM1)
C
         DO JC=1,2
            ZCFC(JC)=ZW1*CFCTS(JC,IYR)+ZW2*CFCTS(JC,IYRM1)
         ENDDO
C
C      SECOND HALF OF YEAR
C
      ELSE
C
         ZW2=(ZAKSEC/ZYRSEC)-0.5
         ZW1=1.-ZW2
         IF (ZW1.GT.1.0) THEN
            ZW2=0.0
            ZW1=1.0
         ENDIF
         ZCO2INT=ZW1*CO2TS(IYR)+ZW2*CO2TS(IYRP1)
         ZCH4INT=ZW1*CH4TS(IYR)+ZW2*CH4TS(IYRP1)
         ZN2OINT=ZW1*AN2OTS(IYR)+ZW2*AN2OTS(IYRP1)
C
         DO JC=1,2
            ZCFC(JC)=ZW1*CFCTS(JC,IYR)+ZW2*CFCTS(JC,IYRP1)
         ENDDO
C
      ENDIF
C
      If (MYID.EQ.0) THEN
         IF (IDAY.EQ.1.AND.IHOUR.EQ.0.AND.NINT(ZSNZT).EQ.0) THEN
            PRINT *,'YAKDAT1=',YAKDAT1
            PRINT *,'ZAKSEC=',ZAKSEC,' ZYRSEC=',ZYRSEC
            PRINT *,'ZW1=',ZW1,' ZW2=',ZW2
            WRITE(*,9001) NZT,ZCO2INT,ZCH4INT,ZN2OINT
         ENDIF
C
C     ------------------------------------------------------------------
C
C        PRINT DIAGNOSTICS
C
         LOPRINT=(MOD(NZT,8640).EQ.0).OR.(NZT.EQ.0)
         IF (LOPRINT) THEN
            WRITE(*,9001) NZT,ZCO2INT,ZCH4INT,ZN2OINT
            WRITE(*,9004) (ZCFC(I),I=1,16)
         ENDIF
      ENDIF
C
C     ------------------------------------------------------------------
C
C       3.      CONVERT TO MASS MIXING RATIO.
C               ------- -- ---- ------ ------
C
      ZCARDI=ZCO2INT*1.E-06*ZCO2MWG/ZAIRMWG
      ZMETHA=ZCH4INT*1.E-09*ZCH4MWG/ZAIRMWG
      ZNITOX=ZN2OINT*1.E-09*ZN2OMWG/ZAIRMWG
C
      DO JC=1,16
         ZCFC(JC)=ZCFC(JC)*1.E-12
      ENDDO
C
C     ------------------------------------------------------------------
C
 9001 FORMAT(' NSTEP=',I8,  ' CO2=',F9.4,' CH4=',F8.2,' N2O=',F8.3)
 9004 FORMAT(' CFC:',16F7.2)
C
      RETURN
      END SUBROUTINE PREPGRG