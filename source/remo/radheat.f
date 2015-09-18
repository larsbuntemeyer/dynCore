C
C   SUBROUTINE RADHEAT
C----------------------------------------------------------------------
C - INPUT 2D
C         PAPHM1,   PAPM1,   PQM1,   PTM1
C----------------------------------------------------------------------
C - IN/OUTPUT 1D
C         PSCLF0,  PSCLFS
C         PSRAD0,  PSRAD0U, PSRADS,  PSRADSU, PSRAF0
C         PSRAFS,  PTCLF0,  PTCLFS,  PTRAD0,  PTRADS
C         PTRADSU, PTRAF0,  PTRAFS,  PTSM1M
C----------------------------------------------------------------------
C - OUTPUT 2D
C         PEMTEF,   PEMTER,  PTRSOF,   PTRSOL
C----------------------------------------------------------------------
C - OUTPUT 1D
C         PALBEDO, PALB
C         PSRFL
C---------------------------------------------------------------------
C - INPUT/OUTPUT 2D
C         PTE
C---------------------------------------------------------------------
C
C**** *RADHEAT* - COMPUTES TEMPERATURE CHANGES DUE TO RADIATION.
C
C     J.F.GELEYN     E.C.M.W.F.     03/06/82.
C
C     E. VAN MEIJGAARD  KNMI-DE BILT DEC 93  UPDATED FOR HIRHAM4/RACMO
C
C     PURPOSE.
C
C
C          THIS ROUTINE COMPUTES THE TENDENCIES OF THE ATMOSPHERE'S
C     TEMPERATURE DUE TO THE EFFECTS OF LONG WAVE AND SHORT WAVE
C     RADIATION. THE COMPUTATION IS DONE ON THE T-1 TIME LEVEL USING
C     VALUES OF ATMOSPHERIC TRANSMISIVITIES AND EMISSIVITIES THAT HAVE
C     BEEN STORED AT THE LAST FULL RADIATION TIME STEP. THE SURFACE
C     SOLAR FLUX LATER TO BE USED IN THE SOIL PROCESS CALCULATIONS IS
C     ALSO STORED.
C
C**   INTERFACE.
C
C
C          *RADHEAT* IS CALLED FROM *PHYSC*.
C          THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE: TS,
C     T AND Q AT T-1 AND P ON LEVEL BOUNDARIES AND IN THE MIDDLE OF THE
C     LAYERS AT THE SAME TIME. ALSO USED ARE THE FOUR PARAMERTERS
C     DESCRIBING THE SUN'S POSITION. THE ROUTINE RETURNS ITS OUTPUT TO
C     THE LONG TERM STORAGE: TENDENCIES OF T AND SURFACE SOLAR FLUX.
C
C     METHOD.
C
C
C          A CALL TO SUBROUTINE *SOLANG* GIVES FIELDS OF SOLAR ZENITH
C     ANGLES AND RELATIVE DAY LENGTH FROM WHICH AN EFFECTIVE SOLAR
C     INFLUX IS COMPUTED. THE RESULTS ARE OF COURSE DIFFERENT DEPENDING
C     ON THE SWITCH ON OR OFF OF THE DIURNAL CYCLE. PRODUCT OF SOLAR
C     INFLUX BY TRANSMISSIVITIES LEADS TO SOLAR FLUXES. THEN THE
C     TEMPERATURES ARE INTERPOLATED/EXTRAPOLATED TO THE LAYER BOUNDARIES
C     (AT THE BOTTOM ONE TAKES THE SURFACE TEMPERATURE) AND A PRODUCT BY
C     EMISSIVITIES OF SIGMA*T**4 GIVES THERMAL FLUXES. THE TWO FLUXES
C     ARE ADDED AND DIVERGENCES COMPUTED TO GIVE HEATING RATES.
C
C     EXTERNALS.
C
C
C          *SOLANG*.
C
C     REFERENCE.
C
C
C          SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION FOR DETAILS
C     ABOUT THE MATHEMATICS OF THIS ROUTINE.
C
C---------------------------------------------------------------------
      SUBROUTINE RADHEAT
     &    (KLON  , KLEV  , KLEVP1, KSTART,
     &     KSTOP , CONACC, TWODT , LRAD ,
C          -- PROGNOSTIC FIELDS AND TENDENCIES --
     &     PTM1, PTE, PQM1,
C          -- PRESSURE FIELDS --
     &     PAPM1, PAPHM1,
C          -- GRID COORDINATES --
     &     SINLAT, COSLAT, SINLON, COSLON ,
C          -- SURFACE FIELDS  --
     &     PTSM1M, PSRFL , PALBEDO,
C          -- EMISSIVITIES AND TRANSMISSIVITIES --
     &     PEMTER, PTRSOL, PEMTEF, PTRSOF,
C          -- REMAINING ELEMENTS IN *RADHEAT* --
     &     PSRADS, PSRADSU, PSRAD0, PSRAD0U,
     &     PTRADS, PTRADSU, PTRAD0,
     &     PSCLF0, PSCLFS , PSRAF0, PSRAFS,
     &     PTCLF0, PTCLFS , PTRAF0, PTRAFS)
C
      IMPLICIT NONE
C
      INCLUDE "COMCON"
      INCLUDE "COMPH1"
C
C     -----------------------------------------------
C     Dummy Arguments
C
      INTEGER, INTENT(IN) :: KLON, KLEV, KLEVP1, KSTART, KSTOP
      REAL,    INTENT(IN) :: CONACC, TWODT
      LOGICAL, INTENT(IN) :: LRAD
C           -- PROGNOSTIC FIELDS AND TENDENCIES --
      REAL,    INTENT(INOUT) :: 
     &     PTM1(KLON,KLEV), PTE(KLON,KLEV), PQM1(KLON,KLEV),
C          -- PRESSURE FIELDS --
     &     PAPM1(KLON,KLEV),PAPHM1(KLON,KLEVP1),
C          -- GRID COORDINATES --
     &     SINLAT(KLON), COSLAT(KLON), SINLON(KLON), COSLON(KLON),
C          -- SURFACE FIELDS --
     &     PTSM1M(KLON), PSRFL(KLON), PALBEDO(KLON),
C          -- EMISSIVITIES --
     &     PEMTER(KLON,KLEVP1), PTRSOL(KLON,KLEVP1),
     &     PEMTEF(KLON,2)     , PTRSOF(KLON,2),
C          -- REMAINING ELEMENTS IN *RADHEAT* --
     &     PSRADS(KLON), PSRADSU(KLON), PSRAD0(KLON), PSRAD0U(KLON),
     &     PTRADS(KLON), PTRADSU(KLON), PTRAD0(KLON),
     &     PSCLF0(KLON), PSCLFS(KLON) , PSRAF0(KLON), PSRAFS(KLON) ,
     &     PTCLF0(KLON), PTCLFS(KLON) , PTRAF0(KLON), PTRAFS(KLON)
C
C     --------------------------------------------------------------
C     Local Variables
C     WORK ARRAYS
C
      REAL :: ZTI(KLON,KLEVP1), ZAMU0(KLON), ZDTDT(KLON), ZEMTERM(KLON),
     &        ZFLB(KLON)      , ZFLT(KLON) , ZFSO(KLON) , ZFSOB(KLON)  ,
     &        ZFTE(KLON)      , ZI0(KLON)  , ZRDAYL(KLON) ,
     &        ZTRSOLM(KLON)
      INTEGER :: JK, JL
      REAL    :: ZALSO, ZCONS2, ZCONS3, ZDIAGT, ZEMISS, ZQEMISS, ZSR0U, 
     &           ZTIM1, ZTIM2, ZTIM3, ZTRSU, ZZFAC, ZSRSU
C
C*    COMPUTATIONAL CONSTANTS.
C
C
      ZDIAGT=CONACC*TWODT
C
      ZCONS2=CDISSE*SOLC
      ZCONS3=G/CPD

      DATA ZEMISS/0.996/
      ZQEMISS=(1.-ZEMISS)/ZEMISS
C
C
C
      IF (LRAD) THEN
C***
C
C
C
C*         2.     SOLAR ANGLE COMPUTATIONS.
C
C
C*         2.1     INTRODUCE THE LATITUDE DEPENDENCY.
C
         ZTIM1=CZEN1
         ZTIM2=-CZEN2
         ZTIM3=CZEN3
C
C*         2.2     CALL TO *SOLANG*.
C
C***
         CALL SOLANG (KLON  , KSTART, KSTOP , ZTIM1 , ZTIM2, ZTIM3,
     &                COSLAT, SINLAT, COSLON, SINLON, ZAMU0, ZRDAYL)
C***
C
C*         2.3     EARTH'S CURVATURE CORRECTION TO CREATE THE INFLUX.
C
         DO JL=KSTART,KSTOP
            ZI0(JL)=ZCONS2*ZAMU0(JL)*ZRDAYL(JL)
         ENDDO
C
C*         3.     TEMPERATURES AT LAYERS' BOUDARIES.
C
C
C*         3.2     INTERPOLATION PROPER.
C
         DO JK=2,KLEV
            DO JL=KSTART,KSTOP
               ZTI(JL,JK)=(PTM1(JL,JK-1)*PAPM1(JL,JK-1)
     &              *(PAPM1(JL,JK)-PAPHM1(JL,JK))
     &              +PTM1(JL,JK)*PAPM1(JL,JK)
     &              *(PAPHM1(JL,JK)-PAPM1(JL,JK-1)))
     &              /(PAPHM1(JL,JK)*(PAPM1(JL,JK)-PAPM1(JL,JK-1)))
            ENDDO
         ENDDO
C
C*        3.3     SURFACE AND TOP OF ATMOSPHERE TEMPERATURE.
C
         DO JL=KSTART,KSTOP
            ZTI(JL,KLEVP1)=PTSM1M(JL)
            ZTI(JL,1)=PTM1(JL,1)-PAPM1(JL,1)*(PTM1(JL,1)-ZTI(JL,2))
     &           /(PAPM1(JL,1)-PAPHM1(JL,2))
         ENDDO
C
C*         4.     FLUXES AND THEIR DIVERGENCES.
C
C*         4.1     TOP FLUX.
C
         CALL COPYRE(PTRSOL(KSTART,1),ZTRSOLM(KSTART),KLON)
         CALL COPYRE(PEMTER(KSTART,1),ZEMTERM(KSTART),KLON)
         DO JL=KSTART,KSTOP
            ZFSO(JL)=ZI0(JL)*ZTRSOLM(JL)
            ZFTE(JL)=STBO*ZTI(JL,1)**4*ZEMTERM(JL)
            ZFLT(JL)=ZFSO(JL)+ZFTE(JL)
            PSRAD0(JL)=PSRAD0(JL)+ZDIAGT*ZFSO(JL)
            PTRAD0(JL)=PTRAD0(JL)+ZDIAGT*ZFTE(JL)
            ZSR0U=-ZI0(JL)*(1.-ZTRSOLM(JL))
            PSRAD0U(JL)=PSRAD0U(JL)+ZDIAGT*ZSR0U
C
            PSRAF0(JL)=PSRAF0(JL)+ZDIAGT*ZI0(JL)*PTRSOF(JL,1)
            PTRAF0(JL)=PTRAF0(JL)+ZDIAGT*STBO*ZTI(JL,1)**4*PEMTEF(JL,1)
            PSCLF0(JL)=PSRAD0(JL)-PSRAF0(JL)
            PTCLF0(JL)=PTRAD0(JL)-PTRAF0(JL)
C
         ENDDO
C
C*         4.2     VERTICAL LOOP, BOTTOM FLUX AND HEATING RATE.
C
C***
         DO JK=1,KLEV
C***
            CALL COPYRE(PTRSOL(KSTART,JK+1),ZTRSOLM(KSTART),KLON)
            CALL COPYRE(PEMTER(KSTART,JK+1),ZEMTERM(KSTART),KLON)
            DO JL=KSTART,KSTOP
               ZFSO(JL)=ZI0(JL)*ZTRSOLM(JL)
               ZFTE(JL)=STBO*ZTI(JL,JK+1)**4*ZEMTERM(JL)
               ZFLB(JL)=ZFSO(JL)+ZFTE(JL)
               ZDTDT(JL)=-ZCONS3*(ZFLB(JL)-ZFLT(JL))/((PAPHM1(JL,JK+1)
     &              -PAPHM1(JL,JK))*(1.+VTMPC2*PQM1(JL,JK)))
               PTE(JL,JK)=PTE(JL,JK)+ZDTDT(JL)
               ZFSOB(JL)=ZFSO(JL)
            ENDDO
C
C*         4.4     FLUX SWAP, END OF THE LOOP AND SURFACE SOLAR FLUX.
C
            DO JL=KSTART,KSTOP
               ZFLT(JL)=ZFLB(JL)
            ENDDO
C***
         ENDDO
C***
         DO JL=KSTART,KSTOP
            PSRFL(JL)=ZFSO(JL)
            PSRADS(JL)=PSRADS(JL)+ZDIAGT*ZFSO(JL)
            PTRADS(JL)=PTRADS(JL)+ZDIAGT*ZFTE(JL)
            ZALSO=PALBEDO(JL)
            ZSRSU=-ZFSOB(JL)*(1./(1.-ZALSO)-1.)
            ZZFAC=STBO*ZTI(JL,KLEVP1)**4
            ZTRSU=-(1.+ZEMTERM(JL)*ZQEMISS)*ZZFAC
            PSRADSU(JL)=PSRADSU(JL)+ZDIAGT*ZSRSU
            PTRADSU(JL)=PTRADSU(JL)+ZDIAGT*ZTRSU
C
            PSRAFS(JL)=PSRAFS(JL)+ZDIAGT*ZI0(JL)*PTRSOF(JL,2)
            PTRAFS(JL)=PTRAFS(JL)+
     &           ZDIAGT*STBO*ZTI(JL,KLEVP1)**4*PEMTEF(JL,2)
            PSCLFS(JL)=PSRADS(JL)-PSRAFS(JL)
            PTCLFS(JL)=PTRADS(JL)-PTRAFS(JL)
C
         ENDDO
C
C
C*         5.     NECESSARY COMPUTATIONS IF SUBROUTINE IS BY-PASSED.
C
C***
      ELSE
C***
         DO JL=KSTART,KSTOP
            PSRFL(JL)=0.
         ENDDO
C
      ENDIF
C***
C
C
C
C*         6.     RETURN WORKSPACE.
C
C
      RETURN
      END SUBROUTINE RADHEAT