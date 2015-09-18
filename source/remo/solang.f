C
C
C**** *SOLANG* - FOR SOLAR ZENITH ANGLE AND RELATIVE DAYLENGTH.
C
C     J.F.GELEYN     E.C.M.W.F.     03/06/82.
C
C     MODIFICATION TO FIT INTO HIRLAM:
C     J.H. CHRISTENSEN    DMI       10/12/91.
C
C     MODIFICATION BY DANIELA JACOB FOR REMO:
C     17.08.94
C     REMOVE OF LDIUR SWITCH FOR REGIONAL MODEL, ROUTINE WORKS ONLY
C     IN THE WAY OF LDIUR = TRUE.
C
C     PURPOSE.
C     --------
C
C          THIS ROUTINE GIVES DIFFERENT RESULTS DEPENDING ON A LOGICAL
C     SWITCH. IF LDIUR IS TRUE ONE OBTAINS ACTUAL SOLAR ZENITH ANGLES
C     AND VALUES OF ONE OR ZERO DEPENDING ON THE SIGN OF THE FORMER. IF
C     LDIUR IS FALSE ONE GETS THE SAME ANSWERS AT ALL POINTS, I.E. MEAN
C     VALUE OF THE DAYTIME SOLAR ZENITH ANGLE AND RELATIVE LENGTH OF
C     THE DAY.
C
C**   INTERFACE.
C     ----------
C
C          *SOLANG* IS CALLED FROM *RADMOD* AND FROM *RADHEAT*.
C          THERE ARE THREE DUMMY ARGUMENTS: *PTIM1*, *PTIM2* AND *PTIM3*
C     ARE  PARAMETERS ABOUT THE SUN'S POSITION.
C          THE ROUTINE RETURNS SOLAR ZENITH ANGLES AND RELATIVE DAY
C     LENGTHS TO THE LONG TERM STORAGE.
C
C     METHOD.
C     -------
C
C          STAIGHTFORWARD IN THE CASE "ON". FOR THE CASE "OFF" THE
C     TYPE OF "ON" COMPUTATION IS REPEATED  WITH 128 POINTS AND THE
C     RELEVANT MEAN VALUES ARE COMPUTED AND STORED.
C
C     EXTERNALS.
C     ----------
C
C          NONE.
C
C     REFERENCE.
C     ----------
C
C          SEE RELEVANT PART OF THE DOCUMENTATION.
C
C
      SUBROUTINE SOLANG (KLON,KSTART,KSTOP,
     &     PTIM1,PTIM2,PTIM3,
     &     COSLAT,SINLAT,COSLON,SINLON,
     &     AMU0,RDAYL)
C
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: KLON,KSTART,KSTOP
      REAL,    INTENT(IN)    :: PTIM1,PTIM2,PTIM3
C     DIMENSION ZMU0(128),ZRDAYL(128)
      REAL,    INTENT(INOUT) :: AMU0(KLON), RDAYL(KLON)
      REAL,    INTENT(INOUT) :: COSLAT(KLON),SINLAT(KLON),
     &                          COSLON(KLON),SINLON(KLON)
C-----------------------------------------------------------------------
C     Local Declarations
C   
      LOGICAL :: LO
      INTEGER :: JL
      REAL    :: ZTIM1, ZTIM2, ZTIM3
C
C     ------------------------------------------------------------------
C
C*         1.     LOCATE AND POSITION SPACE.
C                 ------ --- -------- ------
C
      ZTIM1=PTIM1
      ZTIM2=PTIM2
      ZTIM3=PTIM3
C
C     ------------------------------------------------------------------
C
C*         2.     COMPUTATIONS IF DIURNAL CYCLE "ON".
C                 ------------ -- ------- ----- -----
C
C***
      DO JL=KSTART,KSTOP
         AMU0(JL)=ZTIM1*SINLAT(JL)+ZTIM2*COSLAT(JL)*COSLON(JL)+
     &        ZTIM3*COSLAT(JL)*SINLON(JL)
         LO=AMU0(JL).GE.0.
CRP   CVMGT ERSETZT
C     AMU0(JL)=CVMGT(AMU0(JL),0.,LO)
C     RDAYL(JL)=CVMGT(1.,0.,LO)
         IF (LO) THEN
            RDAYL(JL)=1.
         ELSE
            AMU0(JL)=0.
            RDAYL(JL)=0.
         ENDIF
      ENDDO
C
C     ------------------------------------------------------------------
C
C*         3.     COMPUTATIONS IF DIURNAL CYCLE "OFF".
C                 ------------ -- ------- ----- ------
C  REMOVED FOR REMO -- 17.08.94
C
C     ------------------------------------------------------------------
C
C*         4.     RETURN.
C                 -------
C
      RETURN
      END SUBROUTINE SOLANG
