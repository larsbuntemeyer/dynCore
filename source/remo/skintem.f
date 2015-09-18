C
C     SUBROUTINE SKINTEM
C
CTS 250100
C
C**** *SKINTEM* CALCULATE SEA-ICE SKIN-TEMPERATURE
C
C     F.LUNKEIT                UNIHH          11.04.91
C     D. JACOB CHANGES FOR HIRHAM  1.12 93
C
C     PURPOSE.
C     --------
C     CALCULATE ICE SKIN-TEMPERATURE AS A PROGNOSTIC VARIABLE
C
C**   INTERFACE.
C     ----------
C     *CALL* *SKINTEM* CALLED FROM PHECHAM
C
C     METHOD.
C     -------
C     SURFACE FLUXES ARE LINEARISED AROUND SURFACE-TEMPERATURE
C
C     Q(TSKIN)=Q(TSM1M)+DQDT*(TSKIN-TSM1M)
C
C     THE HEAT-BALANCE EQUATION: Q(TSKIN)*DT=CP*(TSKIN-TSKINM)
C     LEEDS TO:
C
C     QC(TSM1M,TSKIN)=TSKIN*QD(TSM1M)  ==> TSKIN=QC/QD
C
C     QC=CONSTANT FLUX-TERMS ONLY DEPENDET ON TSM1M AND TSKIN
C     QD=DIRIVATIONS OF Q
C
C     EXTERNALS.
C     ----------
C     NONE.
C
C     REFERENCES.
C     -----------
C     NONE.
C
C ---------------------------------------------------------------------
      SUBROUTINE SKINTEM
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
     & (KIDIA  , KFDIA , KLON   ,KLEVP1, TWODT ,
     &  PEMTER , PTSM1M, PSICED , PTEFF,
     &  PALBEDO, PALSOI, PSRFL  , PTHFLI, PDHFTI,
     &  PTSI   ,PTSIM1M, PAHFICE, PQRES , PTSLIN, INFRL  )
C
      IMPLICIT NONE
C
      INCLUDE "COMCON"
      INCLUDE "COMPH2"
      INCLUDE "faktinf.h"
C-----------------------------------------------------------------------
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: KIDIA, KFDIA, KLON, KLEVP1
      REAL,    INTENT(IN)    :: TWODT
C
      REAL,    INTENT(IN)    :: PEMTER(KLON,KLEVP1),
     &                          PTSM1M(KLON)       ,
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
     &                          PSICED(KLON)       ,
     &                          PALBEDO(KLON)      ,
     &                          PALSOI(KLON)       ,
     &                          PSRFL(KLON)        ,
     &                          PTHFLI(KLON)       ,
     &                          PDHFTI(KLON)       ,
     &                          PTSIM1M(KLON)      
c   
      REAL,    INTENT(INOUT) :: PTEFF(KLON)        ,
     &                          PTSI(KLON)         ,
     &                          PAHFICE(KLON)      ,
     &                          PQRES(KLON)        ,
     &                          PTSLIN(KLON)
C
      INTEGER, INTENT(IN)    :: INFRL(KLON)
C-----------------------------------------------------------------------
C     Local Declarations
C
CTS 250100
C
      REAL    :: ZEMTERM(KLON)
      INTEGER :: JL
      REAL :: ZALPHA, ZCPCON, ZCPDT, ZCPICE, ZDIAGT, ZDICE, ZDICEFL, 
     &        ZDQICE, ZDQIDT, ZDSFLX, ZDTHFL, ZDTRFL, ZEMISS, ZICEFL, 
     &        ZQRES, ZRHOICE, ZSFLX, ZSIG, ZSOFL, ZTHFLI
      REAL :: ZTMELT, ZTMST, ZTRFL, ZTSEA
C-----------------------------------------------------------------------
C
C*    1 LOCATE AND POSITION SPACE
C     - ------ --- -------- -----
C
C
C*    2 SET UP CONSTANTS
C     - --- -- ---------
C
      ZTMELT  = TMELT
      ZALPHA  = 2.0
      ZCPICE  = 2090.
      ZRHOICE = 1000.
      ZSIG    = STBO
      ZTMST   = TWODT
      ZDIAGT  = 0.5*TWODT
      ZDICE   = 0.10
      ZEMISS  = 0.996
C
      ZCPCON = ZRHOICE*ZCPICE*ZDICE
      ZCPDT  = ZCPCON/ZTMST
C
C
C*    3 COMPUTE NEW SKIN-TEMPERATURE
C     - ------- --- ----------------
C
      CALL COPYRE(PEMTER(KIDIA,KLEVP1),ZEMTERM(KIDIA),KFDIA-KIDIA+1)
C
      DO JL=KIDIA,KFDIA
C
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
         IF(INFRL(JL).LT.NINT(FAKINF).AND.PSICED(JL).GT.ZDICE) THEN
            ZTSEA=CTFREEZ
            ZTRFL=ZEMTERM(JL)*ZSIG*PTSIM1M(JL)**4
     &             +4.*ZEMISS*ZSIG*PTSIM1M(JL)**4
            ZSOFL=(1.-PALSOI(JL))*PSRFL(JL)/(1.-PALBEDO(JL))
            ZTHFLI=PTHFLI(JL)-PDHFTI(JL)*PTSIM1M(JL)
            ZICEFL=ZALPHA*ZTSEA/PSICED(JL)
            ZDQICE=ZCPDT*PTSIM1M(JL)
            ZSFLX=-ZTRFL-ZTHFLI-ZSOFL-ZICEFL-ZDQICE
C
C*    DERIVATIONS
C
            ZDTRFL=-4.*ZEMISS*ZSIG*PTSIM1M(JL)**3
            ZDTHFL=PDHFTI(JL)
            ZDICEFL=-ZALPHA/PSICED(JL)
            ZDQIDT=-ZCPDT
            ZDSFLX=ZDTRFL+ZDTHFL+ZDICEFL+ZDQIDT
C
            PTSI(JL)=ZSFLX/ZDSFLX
C
            IF(PTSI(JL).GT.ZTMELT) THEN
               ZQRES=(ZALPHA/PSICED(JL)+ZCPDT)*(ZTMELT-PTSI(JL))
               PTSI(JL)=ZTMELT
            ELSE
               ZQRES=0.0
            ENDIF
            PQRES(JL)=PQRES(JL)+ZQRES*ZDIAGT
            PAHFICE(JL)=PAHFICE(JL)+ZALPHA*(PTSI(JL)-ZTSEA)
     &           *ZDIAGT
            PTSLIN(JL)=PTSLIN(JL)+ZALPHA*(PTSI(JL)-ZTSEA)/PSICED(JL)
     &           *ZDIAGT
         ENDIF
CTS 250100
C
C*     ACCUMULATE SKIN-TEMPERATURE
C
C
CTS CHANGE FOR LAND-SEA-SEAICE DIFFERENTIATION
         PTEFF(JL)=PTEFF(JL)+ZDIAGT*PTSM1M(JL)
CTS 250100
C
C
      ENDDO
C
      RETURN
      END SUBROUTINE SKINTEM
