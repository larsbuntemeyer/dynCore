C
C
C**** *RADINT* - ORGANISES THE RADIATION FULL COMPUTATIONS.
C
C     U. SCHLESE     DKRZ-HAMBURG   JUL-93
C
C     E. VAN MEIJGAARD KNMI-DE BILT NOV-93  UPDATED FOR HIRHAM4/RACMO
C
C     PURPOSE.
C
C
C          THIS ROUTINE ORGANISES THE INPUT/OUTPUT FOR THE BLACK-BOX
C     RADIATION COMPUTATIONS PERFORMED IN *RADLSW* EVERY TIME THERE IS A
C     FULL RADIATION TIME STEP. INPUT ARE PROGNOSTIC MODEL VARIABLES AT
C     TIME STEP T-1, SURFACE VALUES OF SHORT-WAVE ALBEDO AND LONG-WAVE
C     EMISSIVITY AND CLIMATOLOGICAL VALUES FOR AEROSOLS AND OZONE (TIME
C     OF THE YEAR DEPENDENT). OUTPUT ARE FLUX TRANSMISSIVITIES AND
C     EMISSIVITIES AT ALL THE HALF LEVELS OF THE GRID (RESPECTIVELY
C     RATIO SOLAR FLUX/SOLAR INPUT AND RATIO THERMAL FLUX/LOCAL
C     BLACK-BODY FLUX). THIS OUTPUT WILL BE USED IN *RADHEAT* AT ALL
C     TIME STEPS UNTIL THE NEXT FULL RADIATION TIME STEP.
C
C**   INTERFACE.
C
C
C          *RADINT* IS CALLED FROM *PHYSC*.
C          THE COMMUNICATIONS WITH *RADHEAT* HAVE BEEN DESCRIBED ABOVE
C     AND THOSE WITH *RADLSW* ARE VIA DUMMY LIST.
C
C     METHOD.
C
C
C          A CALL TO SUBROUTINE *SOLANG* GIVES FIELDS OF SOLAR ZENITH
C     ANGLES AND RELATIVE DAY LENGTH (RESULTS DEPENDING ON THE SWITCH ON
C     OR OFF OF THE DIURNAL CYCLE. THE CONSISTENCY OF THESE VALUES WITH
C     THE USE THEY WILL LATER HAVE IN *RADHEAT* WAS ENSURED BY GIVING TO
C     *SOLANG* AN INPUT CORRESPONDING TO THE MIDDLE OF THE PERIOD DURING
C     WHICH THE RESULTS OF *RADINT* WILL BE VALID.
C          *LEGENDRE AND *FOURIER TRANSFORMS ALLOW TO GO FROM THE T5 AND
C     T10 SPECTRAL DEFINITIONS OF OZONE AND AEROSOLS' CLIMATOLOGIES TO
C     THE MODELS' GRID. THE COMPUTATION OF THE *LEGENDRE POLYNOMIALS IS
C     DONE IN SUBROUTINE *LEGTRI* THAT GETS AS INPUT THE SINE OF
C     LATITUDE.
C          TEMPERATURES AND RELATIVE HUMIDITIES ARE INTER-/EXTRAPOLATED
C     FROM THE PROGNOSTIC LEVELS TO THE LEVELS' BOUNDARIES WITH A METHOD
C     (NON LINEAR IN P) CONSISTENT WITH THE ONE USED IN *RADHEAT* FOR
C     THE SAME PURPOSE.
C         THE ACTUAL HANDLING OF THE INPUT/OUTPUT FOR *RADLSW* IS RATHER
C     STRAIGHTFORWARD IF ONE KNOWS THE CHOICES MADE FOR THE VERTICAL
C     DISTRIBUTIONS OF THE RADIATIVE AGENTS: CLOUDS, AEROSOLS, GASES
C     AND TEMPERATURES.
C
C     EXTERNALS.
C
C
C          *SOLANG*, *LEGTRI* AND *RADLSW*
C

C     REFERENCE.
C
C
C          SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION FOR DETAILS
C     ABOUT THE MATHEMATICS OF THIS ROUTINE.
C
C-----------------------------------------------------------------------
      SUBROUTINE RADINT
     &    (KLON, NLEV , NLEVP1, KSTART, KSTOP,
     &     LRAD, LSOLC, LAER  , LCFC  ,
C          -- PROGNOSTIC FIELDS AND TENDENCIES --
     &     TM1, XM1, QM1,
C          -- PRESSURE FIELDS AND MASKS --
     &     APM1, APHM1, LOLAND, LOGLAC,
C          -- GRID COORDINATES --
     &     SINLAT, COSLAT, SINLON, COSLON,
C          -- SURFACE FIELDS  --
     &     TSM1M , SNM1M, ACLC  ,
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
     &     FOREST, ALB   , ALBEDO,
     &     TSLM1M, TSIM1M, INFRL, INFRW, INFRI,
     &     ALSOL , ALSOW , ALSOI ,
CTS 250100
C          --  EMISSIVITIES --
     &     EMTER, TRSOL, EMTEF, TRSOF,
C          -- REMAINING ELEMENTS IN *RADINT* --
     &     ACLCV, CVDAES, CVDAEL, CVDAEU, CVDAED,
     &     XIM1)
C
      IMPLICIT NONE
C
      INCLUDE "COMCON"
      INCLUDE "COMRAD1"
      INCLUDE "COMPH1"
C not used      INCLUDE "COMPH2"
      INCLUDE "COMRAD2"
C
      INCLUDE "YOMRDI"
      INCLUDE "YOTLUC"
      INCLUDE "faktinf.h"
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN) :: KLON, NLEV , NLEVP1, KSTART, KSTOP
      LOGICAL, INTENT(IN) :: LRAD, LSOLC,LAER,LCFC
C
C          -- PROGNOSTIC FIELDS AND TENDENCIES --
      REAL, INTENT(INOUT) ::
     &     TM1(KLON,NLEV), XM1(KLON,NLEV), QM1(KLON,NLEV),
C          -- PRESSURE FIELDS --
     &     APM1(KLON,NLEV), APHM1(KLON,NLEVP1),
C          -- GRID COORDINATES --
     &     SINLAT(KLON), COSLAT(KLON), SINLON(KLON), COSLON(KLON),
C          -- SURFACE FIELDS  --
     &     TSM1M(KLON), SNM1M(KLON) , ACLC(KLON,NLEV), FOREST(KLON),
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
     &     ALB(KLON)   , ALBEDO(KLON),
     &     TSLM1M(KLON), TSIM1M(KLON),
     &     ALSOL(KLON) , ALSOW(KLON) , ALSOI(KLON) ,
CTS 250100
C          -- REMAINING ELEMENTS IN *RADINT* --
     &     ACLCV(KLON)   , CVDAES(NLEVP1), CVDAEL(NLEVP1),
     &     CVDAEU(NLEVP1), CVDAED(NLEVP1),
C          --  EMISSIVITIES AND TRANSMISSIVITIES --
     &     EMTER(KLON,NLEVP1), TRSOL(KLON,NLEVP1),
     &     EMTEF(KLON,2)     , TRSOF(KLON,2)
C          -- LAND MASK + ICE MASK --
      LOGICAL, INTENT(IN) :: LOLAND(KLON), LOGLAC(KLON)
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
      INTEGER, INTENT(IN) :: INFRL(KLON),INFRW(KLON),INFRI(KLON)
CTS 250100
C     ZCLIC: CLOUD ICE
      REAL, INTENT(IN) :: XIM1(KLON,NLEV)
C
C     -------------------------------------------
C     Local Variables and Arrays
C     DECLARATION WORK ARRAYS IN RADIATION MODULE
C
      REAL AMU0(KLON)     , RDAYL(KLON)     ,
     &     ZCLC(KLON,NLEV), ZCLWA(KLON,NLEV), ZMU0(KLON), ZALSO(KLON),
     &     ZQ(KLON,NLEV)  , ZQS(KLON,NLEV)  , ZOZQ(KLON), ZOZH(KLON) ,
     &     ZAES(KLON)     , ZAEL(KLON)      , ZAEU(KLON), ZAED(KLON) ,
C
     &     ZHTI(KLON,NLEVP1) , ZDP(KLON,NLEV)   , ZQOF(KLON,NLEV),
     &     ZSAER(KLON,NLEV,5), ZAETRO(KLON)     , ZAETRN(KLON)   ,
     &     ZFLT(KLON,NLEVP1) , ZFLS(KLON,NLEVP1),
     &     ZFLTC(KLON,2)     , ZFLSC(KLON,2), ZCLIC(KLON,NLEV)
C
C     DECLARE AUXILIARY ARRAYS RELATED TO LEGENDRE EXPANSION
C     OF TRACE GASES CLIMATOLOGIES
C
      INTEGER, PARAMETER :: JPLGOZ= 6, JPLGOZD=2*JPLGOZ-1
      INTEGER, PARAMETER :: JPLGAE=11, JPLGAED=2*JPLGAE-1
      REAL :: ZFOZQ(KLON,JPLGOZD), ZFOZH(KLON,JPLGOZD),
     &        ZFAES(KLON,JPLGAED), ZFAEL(KLON,JPLGAED),
     &        ZFAEU(KLON,JPLGAED), ZFAED(KLON,JPLGAED)
C
      INTEGER, PARAMETER :: JPALP=JPLGAE*(JPLGAE+1)/2
      REAL :: ZALP(KLON,JPALP)
C
C     WORK SPACE FOR ROUTINE *LEGTRI*
C
      REAL    :: WLEGTR(KLON)
      INTEGER :: ILEGTR(KLON)
CC
      REAL :: CAEROS, ZAEQDN, ZAEQDO, ZAEQLN, ZAEQLO, ZAEQSN, ZAEQSO, 
     &        ZAEQUN, ZAEQUO, ZAETR, ZALBMAX, ZALBMIN, ZALBMN0, ZALBMN1, 
     &        ZALBMX0, ZALBMX1, ZALBSN, ZCOS1, ZCOS10, ZCOS2
      REAL :: ZCOS3, ZCOS4, ZCOS5, ZCOS6, ZCOS7, ZCOS8, ZCOS9, ZCPHN3, 
     &        ZCPHO3, ZCRAE, ZDALB, ZQOFN, ZQOFO, ZSCT, ZSDPN3, ZSDPO3, 
     &        ZSIN1, ZSIN10, ZSIN2, ZSIN3
      REAL :: ZSIN4, ZSIN5, ZSIN6, ZSIN7, ZSIN8, ZSIN9, ZTALB, ZTIM1, 
     &        ZTIM2, ZTIM3, ZTSNMELT, ZZOZH
      INTEGER :: IMM, IMNC, IMNS, IT, JAER, JK, JKL, JL, JMM, JNN,
     &           KAER, KCFC, KMODE, NAER, NCP
C
C     -------------------------------------------
C  
C*    DATA STATEMENTS.
C
C
C          *ZALBICE*, *ZALBSEA* AND *ZALBSNO* ARE ALBEDO VALUES FOR
C     FLOATING ICE, OPEN SEA AND THICK SNOW RESPECTIVELY. *ZSNOWAL*
C     IS A VALUE OF SNOW DEPTH (IN EQUIVALENT WATER) FOR WHICH THE
C     SNOW STARTS BEEING CONSIDERED AS THICK. *ZCARDI* IS THE SPECIFIC
C     ATMOSPHERIC CONTENT IN CO2.
C
      ZTALB=TMELT-10.
      DATA ZALBMN0, ZALBMX0, ZALBMN1, ZALBMX1 /.4,.8,.3,.4/
C
C*    CO2 MASS MIXING RATIO *ZCARDI* [KG/KG] (!!!) INITIALIZED IN SURADI
C
C
C*    SECURITY PARAMETER.
C
C
C     *ZEPSEC* AVOIDS 0/0 IN THE DIAGNOSTIC OF TOTAL CLOUD COVER.
C     *ZEPCLC* IS A SECURITY TO AVOID ZERO OR ONE CLOUD COVERS AND
C     *ZEPH2O* IS A SECURITY TO AVOID WATER VAPOUR CONTENT IN A LAYER
C    *         TO BE MORE THEN THE RESPECTIVE VALUE AT SATURATION.
C     *ZEPALB* IS A SECURITY TO AVOID ZERO ALBEDOS.
C
C
C*    COMPUTATIONAL CONSTANTS.
C
C
CMB
      ZSCT=CDISSEM*SOLC
      ZCRAE=CRAE*(CRAE+2.)
C***
      IF(LRAD) THEN
C
C
C*         2.     SOLAR ANGLE AND OZONE/AEROSOL PARAMETERS COMPUTATIONS.
C
C
C     USE MEAN VALUE BETWEEN TWO FULL RADIATION TIME-STEPS FOR ORBITAL
C     PARAMETERS.
C
         ZTIM1=CZEN1M
         ZTIM2=-CZEN2M
         ZTIM3=CZEN3M
C
C
C*         2.2     CALL TO *SOLANG* FOR ZENITH ANGLE AND DAYLENGTH.
C
C***
         CALL SOLANG (KLON  , KSTART, KSTOP , ZTIM1 , ZTIM2, ZTIM3,
     &        COSLAT, SINLAT, COSLON, SINLON, AMU0 , RDAYL)
C***
         DO JL=KSTART,KSTOP
            ZMU0(JL)=CRAE/(SQRT(AMU0(JL)**2+ZCRAE)-AMU0(JL))
         ENDDO
C

C*         2.3     CALL TO LEGTRI.
C
C***
         NCP=6
         CALL LEGTRI(KLON  , KSTART, KSTOP, NCP   ,
     &        SINLAT, COSLAT, ZALP , WLEGTR, ILEGTR)
C***
C
C*         2.4     *LEGENDRE TRANSFORM FOR OZONE.
C
         CALL SETRA(ZFOZQ,KLON*(2*NCP-1),0.0)
         CALL SETRA(ZFOZH,KLON*(2*NCP-1),0.0)
C
         IMM=0
         IMNC=0
         IMNS=0
         DO JMM=1,NCP
            IMM=IMM+1
            DO JNN=JMM,NCP
               IMNC=IMNC+1
               DO JL=KSTART,KSTOP
                  ZFOZQ(JL,IMM)=ZFOZQ(JL,IMM)+ZALP(JL,IMNC)*COZQC(IMNC)
                  ZFOZH(JL,IMM)=ZFOZH(JL,IMM)+ZALP(JL,IMNC)*COZHC(IMNC)
               ENDDO
            ENDDO
            IF(JMM.NE.1) THEN
               IMM=IMM+1
               DO JNN=JMM,NCP
                  IMNS=IMNS+1
                  DO JL=KSTART,KSTOP
                     ZFOZQ(JL,IMM)=ZFOZQ(JL,IMM)+
     &                    ZALP(JL,IMNS+NCP)*COZQS(IMNS)
                     ZFOZH(JL,IMM)=ZFOZH(JL,IMM)+
     &                    ZALP(JL,IMNS+NCP)*COZHS(IMNS)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
C
C***
C         2.5     CALL TO LEGTRI
C
         NCP=11
         CALL LEGTRI(KLON  , KSTART, KSTOP, NCP   ,
     &        SINLAT, COSLAT, ZALP , WLEGTR, ILEGTR)
C
C         2.6     *LEGENDRE TRANSFORM FOR AEROSOLS
C
         CAEROS=1.E-15
         DO JMM=1,2*NCP-1
            DO JL=KSTART,KSTOP
               ZFAES(JL,JMM) = 0.
               ZFAEL(JL,JMM) = 0.
               ZFAEU(JL,JMM) = 0.
               ZFAED(JL,JMM) = 0.
            ENDDO
         ENDDO
         IMM=0
         IMNC=0
         IMNS=0
         DO JMM=1,NCP
            IMM=IMM+1
            DO JNN=JMM,NCP
               IMNC=IMNC+1
               DO JL=KSTART,KSTOP
                  ZFAES(JL,IMM) = ZFAES(JL,IMM)+
     &                 ZALP(JL,IMNC)*CAESC(IMNC)
                  ZFAEL(JL,IMM) = ZFAEL(JL,IMM)+
     &                 ZALP(JL,IMNC)*CAELC(IMNC)
                  ZFAEU(JL,IMM) = ZFAEU(JL,IMM)+
     &                 ZALP(JL,IMNC)*CAEUC(IMNC)
                  ZFAED(JL,IMM) = ZFAED(JL,IMM)+
     &                 ZALP(JL,IMNC)*CAEDC(IMNC)
               ENDDO
            ENDDO
            IF(JMM.NE.1) THEN
               IMM=IMM+1
               DO JNN=JMM,NCP
                  IMNS=IMNS+1
                  DO JL=KSTART,KSTOP
                     ZFAES(JL,IMM) = ZFAES(JL,IMM)+
     &                    ZALP(JL,IMNS+NCP)*CAESS(IMNS)
                     ZFAEL(JL,IMM) = ZFAEL(JL,IMM)+
     &                    ZALP(JL,IMNS+NCP)*CAELS(IMNS)
                     ZFAEU(JL,IMM) = ZFAEU(JL,IMM)+
     &                    ZALP(JL,IMNS+NCP)*CAEUS(IMNS)
                     ZFAED(JL,IMM) = ZFAED(JL,IMM)+
     &                    ZALP(JL,IMNS+NCP)*CAEDS(IMNS)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
C
C
C*         3.     PREPARE INPUT FOR RADIATION.
C
C
C*         3.2    HUMIDITY
C
         DO JK=1,NLEV
            DO JL=KSTART,KSTOP
               IT=INT(TM1(JL,JK)*1000.)
               ZQS(JL,JK)=TLUCUA(IT)/APM1(JL,JK)
               ZQS(JL,JK)=ZQS(JL,JK)/(1.-VTMPC1*ZQS(JL,JK))
               ZQS(JL,JK)=MAX(2.*ZEPH2O,ZQS(JL,JK))
               ZQ(JL,JK)=AMAX1(QM1(JL,JK),ZEPH2O)
            ENDDO
         ENDDO
C
C*         3.4   HALF-LEVEL TEMPERATURES
C
         DO JK=2,NLEV
            DO JL=KSTART,KSTOP
               ZHTI(JL,JK)=
     &            (TM1(JL,JK-1)*APM1(JL,JK-1)*(APM1(JL,JK)-APHM1(JL,JK))
     &            +TM1(JL,JK)*APM1(JL,JK)*(APHM1(JL,JK)-APM1(JL,JK-1)))
     &            *(1./(APHM1(JL,JK)*(APM1(JL,JK)-APM1(JL,JK-1))))
            ENDDO
         ENDDO
C
         DO JL=KSTART,KSTOP
            ZHTI(JL,NLEVP1)=TSM1M(JL)
            ZHTI(JL,1)=TM1(JL,1)-APM1(JL,1)*(TM1(JL,1)-ZHTI(JL,2))
     &           /(APM1(JL,1)-APHM1(JL,2))
         ENDDO
C
C*         3.6     CLOUDS.
C
         DO JK=1,NLEV
            DO JL=KSTART,KSTOP
               ZCLC(JL,JK)=ACLC(JL,JK)
               ZCLC(JL,JK)=MIN(MAX(ZCLC(JL,JK),ZEPCLC),1.-ZEPCLC)
               ZCLWA(JL,JK)=AMAX1(XM1(JL,JK),0.)
               ZCLIC(JL,JK)=AMAX1(XIM1(JL,JK),0.)
            ENDDO
         ENDDO
C
         DO JL=KSTART,KSTOP
            ACLCV(JL)=1.-ACLC(JL,1)
         ENDDO
         DO JK=2,NLEV
            DO JL=KSTART,KSTOP
               ACLCV(JL)=ACLCV(JL)*(1.-AMAX1(ACLC(JL,JK),ACLC(JL,JK-1)))
     &              /(1.-AMIN1(ACLC(JL,JK-1),1.-ZEPSEC))
            ENDDO
         ENDDO
         DO JL=KSTART,KSTOP
            ACLCV(JL)=1.-ACLCV(JL)
         ENDDO
C
C*         3.8     ALBEDO.
C
         DO JL=KSTART,KSTOP
            ZALSO(JL)=ALB(JL)
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
            IF (INFRL(JL).GT.0) THEN
CTS 250100
               ZTSNMELT=TMELT
               IF(LOGLAC(JL)) THEN
                  ZALBMIN=0.6
                  ZALBMAX=0.8
               ELSE
                  ZALBMIN=(1.-FOREST(JL))*ZALBMN0+FOREST(JL)*ZALBMN1
                  ZALBMAX=(1.-FOREST(JL))*ZALBMX0+FOREST(JL)*ZALBMX1
               ENDIF
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
               IF (TSLM1M(JL).GE.ZTSNMELT) THEN
CTS 250100
                  ZALBSN=ZALBMIN
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
               ELSEIF (TSLM1M(JL).LT.ZTALB) THEN
CTS 250100
                  ZALBSN=ZALBMAX
               ELSE
                  ZDALB=(ZALBMAX-ZALBMIN)/(ZTSNMELT-ZTALB)
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
                  ZALBSN=ZALBMIN+ZDALB*(ZTSNMELT-TSLM1M(JL))
CTS 250100
               ENDIF
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
               ALSOL(JL)=ZALSO(JL)+(ZALBSN-ZALSO(JL))
     &              *(SNM1M(JL)/(SNM1M(JL)+ZSNOWAL))
            ENDIF
            IF (INFRI(JL).GT.0) THEN
               ZTSNMELT=TMELT
               ZALBMIN=0.5
               ZALBMAX=0.75
               IF (TSIM1M(JL).GE.ZTSNMELT) THEN
                  ALSOI(JL)=ZALBMIN
               ELSEIF (TSIM1M(JL).LT.ZTALB) THEN
                  ALSOI(JL)=ZALBMAX
               ELSE
                  ZDALB=(ZALBMAX-ZALBMIN)/(ZTSNMELT-ZTALB)
                  ALSOI(JL)=ZALBMIN+ZDALB*(ZTSNMELT-TSIM1M(JL))
               ENDIF
            ENDIF
            IF (INFRW(JL).GT.0) THEN
C           ZWI=ACOS(AMU0(JL))*180./API
C           ZW=COS(API/180.*AMIN1(ZWI,80.))
C           ALSOW(JL)=0.05/(ZW+0.15)
               ALSOW(JL)=ZALBSEA
            ENDIF
            ZALSO(JL)=(FLOAT(INFRL(JL))*ALSOL(JL)
     &               + FLOAT(INFRW(JL))*ALSOW(JL)
     &               + FLOAT(INFRI(JL))*ALSOI(JL))*EDFAKINF
CTS 250100
            ZALSO(JL)=MIN(ZALBSNO,MAX(ZEPALB,ZALSO(JL)))
            ALBEDO(JL)=ZALSO(JL)
         ENDDO
C
C
C*         4.2     FOURIER TRANSFORM FOR OZONE AND AEROSOLS
C
         DO JL=KSTART,KSTOP
            ZCOS1=COSLON(JL)
            ZSIN1=SINLON(JL)
C
            ZCOS2=ZCOS1*ZCOS1-ZSIN1*ZSIN1
            ZSIN2=ZSIN1*ZCOS1+ZCOS1*ZSIN1
            ZCOS3=ZCOS2*ZCOS1-ZSIN2*ZSIN1
            ZSIN3=ZSIN2*ZCOS1+ZCOS2*ZSIN1
            ZCOS4=ZCOS3*ZCOS1-ZSIN3*ZSIN1
            ZSIN4=ZSIN3*ZCOS1+ZCOS3*ZSIN1
            ZCOS5=ZCOS4*ZCOS1-ZSIN4*ZSIN1
            ZSIN5=ZSIN4*ZCOS1+ZCOS4*ZSIN1
C
            ZOZQ(JL)=  ZFOZQ(JL, 1)
     &           +2.*( ZFOZQ(JL, 2)*ZCOS1
     &               + ZFOZQ(JL, 3)*ZSIN1
     &               + ZFOZQ(JL, 4)*ZCOS2
     &               + ZFOZQ(JL, 5)*ZSIN2
     &               + ZFOZQ(JL, 6)*ZCOS3
     &               + ZFOZQ(JL, 7)*ZSIN3
     &               + ZFOZQ(JL, 8)*ZCOS4
     &               + ZFOZQ(JL, 9)*ZSIN4
     &               + ZFOZQ(JL,10)*ZCOS5
     &               + ZFOZQ(JL,11)*ZSIN5  )
C
            ZZOZH   =  ZFOZH(JL, 1)
     &           +2.*( ZFOZH(JL, 2)*ZCOS1
     &               + ZFOZH(JL, 3)*ZSIN1
     &               + ZFOZH(JL, 4)*ZCOS2
     &               + ZFOZH(JL, 5)*ZSIN2
     &               + ZFOZH(JL, 6)*ZCOS3
     &               + ZFOZH(JL, 7)*ZSIN3
     &               + ZFOZH(JL, 8)*ZCOS4
     &               + ZFOZH(JL, 9)*ZSIN4
     &               + ZFOZH(JL,10)*ZCOS5
     &               + ZFOZH(JL,11)*ZSIN5  )
C
            ZOZH(JL) =  SQRT(ZZOZH)**3
C
C
            ZCOS6    = ZCOS5*ZCOS1-ZSIN5*ZSIN1
            ZSIN6    = ZSIN5*ZCOS1+ZCOS5*ZSIN1
            ZCOS7    = ZCOS6*ZCOS1-ZSIN6*ZSIN1
            ZSIN7    = ZSIN6*ZCOS1+ZCOS6*ZSIN1
            ZCOS8    = ZCOS7*ZCOS1-ZSIN7*ZSIN1
            ZSIN8    = ZSIN7*ZCOS1+ZCOS7*ZSIN1
            ZCOS9    = ZCOS8*ZCOS1-ZSIN8*ZSIN1
            ZSIN9    = ZSIN8*ZCOS1+ZCOS8*ZSIN1
            ZCOS10   = ZCOS9*ZCOS1-ZSIN9*ZSIN1
            ZSIN10   = ZSIN9*ZCOS1+ZCOS9*ZSIN1
C
            ZAES(JL) = ZFAES(JL, 1)
     &           +2.*( ZFAES(JL, 2)*ZCOS1
     &               + ZFAES(JL, 3)*ZSIN1
     &               + ZFAES(JL, 4)*ZCOS2
     &               + ZFAES(JL, 5)*ZSIN2
     &               + ZFAES(JL, 6)*ZCOS3
     &               + ZFAES(JL, 7)*ZSIN3
     &               + ZFAES(JL, 8)*ZCOS4
     &               + ZFAES(JL, 9)*ZSIN4
     &               + ZFAES(JL,10)*ZCOS5
     &               + ZFAES(JL,11)*ZSIN5
     &               + ZFAES(JL,12)*ZCOS6
     &               + ZFAES(JL,13)*ZSIN6
     &               + ZFAES(JL,14)*ZCOS7
     &               + ZFAES(JL,15)*ZSIN7
     &               + ZFAES(JL,16)*ZCOS8
     &               + ZFAES(JL,17)*ZSIN8
     &               + ZFAES(JL,18)*ZCOS9
     &               + ZFAES(JL,19)*ZSIN9
     &               + ZFAES(JL,20)*ZCOS10
     &               + ZFAES(JL,21)*ZSIN10 )
C
            ZAEL(JL) = ZFAEL(JL, 1)
     &           +2.*( ZFAEL(JL, 2)*ZCOS1
     &               + ZFAEL(JL, 3)*ZSIN1
     &               + ZFAEL(JL, 4)*ZCOS2
     &               + ZFAEL(JL, 5)*ZSIN2
     &               + ZFAEL(JL, 6)*ZCOS3
     &               + ZFAEL(JL, 7)*ZSIN3
     &               + ZFAEL(JL, 8)*ZCOS4
     &               + ZFAEL(JL, 9)*ZSIN4
     &               + ZFAEL(JL,10)*ZCOS5
     &               + ZFAEL(JL,11)*ZSIN5
     &               + ZFAEL(JL,12)*ZCOS6
     &               + ZFAEL(JL,13)*ZSIN6
     &               + ZFAEL(JL,14)*ZCOS7
     &               + ZFAEL(JL,15)*ZSIN7
     &               + ZFAEL(JL,16)*ZCOS8
     &               + ZFAEL(JL,17)*ZSIN8
     &               + ZFAEL(JL,18)*ZCOS9
     &               + ZFAEL(JL,19)*ZSIN9
     &               + ZFAEL(JL,20)*ZCOS10
     &               + ZFAEL(JL,21)*ZSIN10 )
C
            ZAEU(JL) = ZFAEU(JL, 1)
     &           +2.*( ZFAEU(JL, 2)*ZCOS1
     &               + ZFAEU(JL, 3)*ZSIN1
     &               + ZFAEU(JL, 4)*ZCOS2
     &               + ZFAEU(JL, 5)*ZSIN2
     &               + ZFAEU(JL, 6)*ZCOS3
     &               + ZFAEU(JL, 7)*ZSIN3
     &               + ZFAEU(JL, 8)*ZCOS4
     &               + ZFAEU(JL, 9)*ZSIN4
     &               + ZFAEU(JL,10)*ZCOS5
     &               + ZFAEU(JL,11)*ZSIN5
     &               + ZFAEU(JL,12)*ZCOS6
     &               + ZFAEU(JL,13)*ZSIN6
     &               + ZFAEU(JL,14)*ZCOS7
     &               + ZFAEU(JL,15)*ZSIN7
     &               + ZFAEU(JL,16)*ZCOS8
     &               + ZFAEU(JL,17)*ZSIN8
     &               + ZFAEU(JL,18)*ZCOS9
     &               + ZFAEU(JL,19)*ZSIN9
     &               + ZFAEU(JL,20)*ZCOS10
     &               + ZFAEU(JL,21)*ZSIN10 )
C
            ZAED(JL) = ZFAED(JL, 1)
     &           +2.*( ZFAED(JL, 2)*ZCOS1
     &               + ZFAED(JL, 3)*ZSIN1
     &               + ZFAED(JL, 4)*ZCOS2
     &               + ZFAED(JL, 5)*ZSIN2
     &               + ZFAED(JL, 6)*ZCOS3
     &               + ZFAED(JL, 7)*ZSIN3
     &               + ZFAED(JL, 8)*ZCOS4
     &               + ZFAED(JL, 9)*ZSIN4
     &               + ZFAED(JL,10)*ZCOS5
     &               + ZFAED(JL,11)*ZSIN5
     &               + ZFAED(JL,12)*ZCOS6
     &               + ZFAED(JL,13)*ZSIN6
     &               + ZFAED(JL,14)*ZCOS7
     &               + ZFAED(JL,15)*ZSIN7
     &               + ZFAED(JL,16)*ZCOS8
     &               + ZFAED(JL,17)*ZSIN8
     &               + ZFAED(JL,18)*ZCOS9
     &               + ZFAED(JL,19)*ZSIN9
     &               + ZFAED(JL,20)*ZCOS10
     &               + ZFAED(JL,21)*ZSIN10 )
         ENDDO
C
C
C*         5.1     INPUT: CO2, O3 AND AEROSOLS.
C                  (AEROSOL IS STORED UPSIDE DOWN)
C
         DO JL=KSTART,KSTOP
            ZAETRO(JL)=1.
         ENDDO
         DO JK=1,NLEV
            JKL=NLEV+1-JK
            DO JL=KSTART,KSTOP
               ZCPHO3=APHM1(JL,JK)**3
               ZCPHN3=APHM1(JL,JK+1)**3
               ZSDPO3=SQRT(ZCPHO3)
               ZSDPN3=SQRT(ZCPHN3)
               ZAEQSO=CAEOPS*ZAES(JL)*CVDAES(JK)
               ZAEQSN=CAEOPS*ZAES(JL)*CVDAES(JK+1)
               ZAEQLO=CAEOPL*ZAEL(JL)*CVDAEL(JK)
               ZAEQLN=CAEOPL*ZAEL(JL)*CVDAEL(JK+1)
               ZAEQUO=CAEOPU*ZAEU(JL)*CVDAEU(JK)
               ZAEQUN=CAEOPU*ZAEU(JL)*CVDAEU(JK+1)
               ZAEQDO=CAEOPD*ZAED(JL)*CVDAED(JK)
               ZAEQDN=CAEOPD*ZAED(JL)*CVDAED(JK+1)
               ZAETRN(JL)=ZAETRO(JL)*(MIN(1.,ZHTI(JL,JK)/ZHTI(JL,JK+1)))
     &              **CTRPT
               ZAETR=SQRT(ZAETRN(JL)*ZAETRO(JL))
               ZQOFO=ZOZQ(JL)*ZSDPO3/(ZSDPO3+ZOZH(JL))
               ZQOFN=ZOZQ(JL)*ZSDPN3/(ZSDPN3+ZOZH(JL))
               ZDP(JL,JK)=APHM1(JL,JK+1)-APHM1(JL,JK)
               ZSAER(JL,JKL,1)=(1.-ZAETR)*(CTRBGA*ZDP(JL,JK)
     &              + ZAEQLN-ZAEQLO + ZAEQDN-ZAEQDO)
               ZSAER(JL,JKL,2)=(1.-ZAETR)*(ZAEQSN-ZAEQSO)
               ZSAER(JL,JKL,3)=(1.-ZAETR)*(ZAEQUN-ZAEQUO)
               ZSAER(JL,JKL,4)= ZAETR * CVOBGA*ZDP(JL,JK)
               ZSAER(JL,JKL,5)= ZAETR * CSTBGA*ZDP(JL,JK)
               ZQOF(JL,JK)=ZQOFN-ZQOFO
               ZAETRO(JL)=ZAETRN(JL)
            ENDDO
         ENDDO
C
C    AEROSOL SWITCH
C
CRVD
C*I RADINT.626
C     GUIDE FOR USING EITHER THE 5 TANRE AEROSOL TYPES OR THE 11 (12)

C     COMPONENTS OF THE MUNICH AEROSOL CLIMATOLOGY (MAK1):
C*    USE DECK *SUAER* AND COMDECK *YOMAER* TYPICAL FOR THE CLIMATOLOGY
C*    SET NAER=5 FOR TANRE AEROSOLS (5 TYPES)
C     - ZSAER IS OPTICAL THICKNESS OF AEROSOLS
C*    SET NAER=11 (12) FOR MAK1 11 (12) AEROSOL COMPONENTS
C     - ZSAER IS AEROSOL AMOUNT IN PARTICLE/CM2
C     - ZFAER IS AEROSOL DENSITY IN PARTICLE/CM3
C
         IF(LAER) THEN
            KAER=1
            NAER=5
         ELSE
            KAER=0
            NAER=0
         ENDIF
C
C*   KAER IS THE MULTIPLICATION FACTOR FOR AEROSOL OPTICAL DEPTHS
C
CRVD
C      NAER=0: NO AEROSOLS; ARRAY ZSAER IS NOT USED IN RADIATION CODE
C      NAER=5: USE TANRE ZSAER VALUES EXCEPT FOR VOLCANIC TYPE (INDEX 4
C              BACKGROUND VALUE CAEROS IS USED)
C      NAER=11(12): USE 11 (12) COMPONENTS MAK1
C
C
         IF(NAER.EQ.5)THEN
            DO JK=1,NLEV
               DO JL=KSTART,KSTOP
                  ZSAER(JL,JK,1)=AMAX1(ZSAER(JL,JK,1),CAEROS)
                  ZSAER(JL,JK,2)=AMAX1(ZSAER(JL,JK,2),CAEROS)
                  ZSAER(JL,JK,3)=AMAX1(ZSAER(JL,JK,3),CAEROS)
                  ZSAER(JL,JK,4)=CAEROS
                  ZSAER(JL,JK,5)=AMAX1(ZSAER(JL,JK,5),CAEROS)
               ENDDO
            ENDDO
C
C       ACTIVATE COMMENT LINES PRECEDED BY CRVD
C       WHEN THE VALUE OF NAER IS INCREASED TO 11 OR 12
C
CRVD      ELSE IF(NAER.EQ.11.OR.NAER.EQ.12)THEN
C           IF (NAER.LT.JPAER) STOP ' INCREASE JPAER IN RADINT '
CRVD        DO 5381 JAER=1,NAER
CRVD        DO 5371 JK=1,NLEV
CRVD          JKL=NLEV+1-JK
CRVD          DO 5361 JL=KSTART,KSTOP
CRVD          ZSAER(JL,JK,JAER)=ZFAER(JL,JK,JAER)*CRAOGR*TM1(JL,JKL)*
CRVD     +                     ZDP(JL,JKL)/ZAP(JL,JKL)
CRVD 5361   CONTINUE
CRVD 5371   CONTINUE
CRVD 5381   CONTINUE
         ELSE IF (NAER.EQ.0) THEN
            DO JAER=1,5
               DO JK=1,NLEV
                  DO JL=KSTART,KSTOP
                     ZSAER(JL,JK,JAER)=CAEROS
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
C
C        SWITCH FOR SOLAR CLEAR SKY DIAGNOSTIC
C
         IF(LSOLC) THEN
            KMODE=1
         ELSE
            KMODE=0
         ENDIF
C
C        CFC SWITCH
C
         IF(LCFC) THEN
            KCFC=1
         ELSE
            KCFC=0
         ENDIF
C
C*   KCFC IS THE MULTIPLICATION FACTOR FOR CFC'S
C    KCFC=0 ; NO CFC'S INVOLVED IN RADIATION COMPUTATION
C
C
C
C
C*           5.3      CALL TO *RADLSW*
C
         CALL RADLSW(KLON , NLEV  , KSTART , KSTOP, KMODE, KAER , KCFC,
     &               ZSCT , LOLAND, LOGLAC , ZSAER, ZALSO, ZCLC , ZMU0,
     &               ZQOF , APHM1(1,NLEVP1), ZDP  , ZCLWA, ZQ   , ZQS ,
     &               APHM1, TM1   , ZHTI   , ZFLS , ZFLT , ZFLSC, ZFLTC,
     &               ZCLIC)
C
C***
C
C*         5.4     STORAGE OF THE OUTPUT.
C
         DO JK=1,NLEVP1
            DO JL=KSTART,KSTOP
               TRSOL(JL,JK)=ZFLS(JL,JK)/(ZSCT*ZMU0(JL))
               EMTER(JL,JK)=ZFLT(JL,JK)/(STBO*ZHTI(JL,JK)**4)
            ENDDO
         ENDDO
C
C
C     CLEAR SKY FLUXES
C
         DO JL=KSTART,KSTOP
            EMTEF(JL,1)=ZFLTC(JL,1)/(STBO*ZHTI(JL,1)**4)
            EMTEF(JL,2)=ZFLTC(JL,2)/(STBO*ZHTI(JL,NLEVP1)**4)
            TRSOF(JL,1)=ZFLSC(JL,1)/(ZSCT*ZMU0(JL))
            TRSOF(JL,2)=ZFLSC(JL,2)/(ZSCT*ZMU0(JL))
         ENDDO
C
C*         7.     NECESSARY COMPUTATIONS IF SUBROUTINE IS BY-PASSED.
C
C***
      ELSE
C***
         DO JK=1,NLEVP1
            DO JL=KSTART,KSTOP
               TRSOL(JL,JK)=0.
               EMTER(JL,JK)=0.
            ENDDO
         ENDDO
         DO JK=1,2
            DO JL=KSTART,KSTOP
               EMTEF(JL,JK)=0.
               TRSOF(JL,JK)=0.
            ENDDO
         ENDDO
C***
      ENDIF
C***
C
C
      RETURN
      END SUBROUTINE RADINT
