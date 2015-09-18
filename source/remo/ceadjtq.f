      SUBROUTINE CEADJTQ
     &    (KIDIA, KFDIA, KLON,     KLEV,     KK,
     &     PP,       PT,       PQ,       LDFLAG,   KCALL)
C
      IMPLICIT NONE
C
C          M.TIEDTKE         E.C.M.W.F.     12/89
C          D.SALMOND         CRAY(UK))      12/8/91
C
C          PURPOSE.
C          --------
C          TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT
C
C          INTERFACE
C          ---------
C          THIS ROUTINE IS CALLED FROM SUBROUTINES:
C              *CEBASE*   (T AND Q AT CONDENSTION LEVEL)
C              *CEASC*    (T AND Q AT CLOUD LEVELS)
C              *CEINI*    (ENVIRONMENTAL T AND QS VALUES AT HALF LEVELS)
C          INPUT ARE UNADJUSTED T AND Q VALUES,
C          IT RETURNS ADJUSTED VALUES OF T AND Q
C          NOTE: INPUT PARAMETER KCALL DEFINES CALCULATION AS
C               KCALL=0    ENV. T AND QS IN*CEINI*
C               KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CEBASE, CEASC)
C               KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CEDLFS,CEDDRAF)
C
C          EXTERNALS
C          ---------
C          3 LOOKUP TABLES ( TLUCUA, TLUCUB, TLUCUC )
C          FOR CONDENSATION CALCULATIONS.
C          THE TABLES ARE INITIALISED IN *SETPHYS*.
C
      INCLUDE "COMCON"
      INCLUDE "YOTLUC"
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: KIDIA, KFDIA, KLON, KLEV, KK, KCALL
      REAL,    INTENT(IN)    :: PP(KLON)
      LOGICAL, INTENT(IN)    :: LDFLAG(KLON)
      REAL,    INTENT(INOUT) :: PT(KLON,KLEV), PQ(KLON,KLEV)
C
C     Local Variables
C
      INTEGER                :: ISUM,IT,JL
      REAL                   :: ZCOND1,ZCOR,ZQSAT
      REAL, DIMENSION(KLON)  :: ZCOND, ZQP
C
C----------------------------------------------------------------------
C
C     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
C                  -----------------------------------------------------
C
      IF (KCALL.EQ.1 ) THEN
C
         ISUM=0
!DIR$    IVDEP
         DO JL=KIDIA,KFDIA
            IF(LDFLAG(JL)) THEN
               ZQP(JL)=1./PP(JL)
               IT=INT(PT(JL,KK)*1000.)
C
               ZQSAT=TLUCUA(IT)*ZQP(JL)
C
               ZQSAT=MIN(0.5,ZQSAT)
               ZCOR=1./(1.-VTMPC1*ZQSAT)
               ZQSAT=ZQSAT*ZCOR
               ZCOND(JL)=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(IT))
               ZCOND(JL)=MAX(ZCOND(JL),0.)
               PT(JL,KK)=PT(JL,KK)+TLUCUC(IT)*ZCOND(JL)
               PQ(JL,KK)=PQ(JL,KK)-ZCOND(JL)
               IF(ZCOND(JL).NE.0.0) ISUM=ISUM+1
            ENDIF
         ENDDO
C
         IF(ISUM.NE.0) THEN
C
!DIR$    IVDEP
            DO JL=KIDIA,KFDIA
               IF(LDFLAG(JL)) THEN
                  IF(ZCOND(JL).NE.0.) THEN
                     IT=INT(PT(JL,KK)*1000.)
                     ZQSAT=TLUCUA(IT)*ZQP(JL)
                     ZQSAT=MIN(0.5,ZQSAT)
                     ZCOR=1./(1.-VTMPC1*ZQSAT)
                     ZQSAT=ZQSAT*ZCOR
                     ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(IT))
                     PT(JL,KK)=PT(JL,KK)+TLUCUC(IT)*ZCOND1
                     PQ(JL,KK)=PQ(JL,KK)-ZCOND1
                  ENDIF
               ENDIF
            ENDDO
C
         ENDIF
C
      ENDIF

      IF(KCALL.EQ.2) THEN
C
         ISUM=0
!DIR$    IVDEP
         DO JL=KIDIA,KFDIA
            IF(LDFLAG(JL)) THEN
               IT=INT(PT(JL,KK)*1000.)
               ZQP(JL)=1./PP(JL)
               ZQSAT=TLUCUA(IT)*ZQP(JL)
               ZQSAT=MIN(0.5,ZQSAT)
               ZCOR=1./(1.-VTMPC1*ZQSAT)
               ZQSAT=ZQSAT*ZCOR
               ZCOND(JL)=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(IT))
               ZCOND(JL)=MIN(ZCOND(JL),0.)
               PT(JL,KK)=PT(JL,KK)+TLUCUC(IT)*ZCOND(JL)
               PQ(JL,KK)=PQ(JL,KK)-ZCOND(JL)
               IF(ZCOND(JL).NE.0.0) ISUM=ISUM+1
            ENDIF
         ENDDO
C
         IF(ISUM.NE.0) THEN
C
!DIR$    IVDEP
            DO JL=KIDIA,KFDIA
               IF(LDFLAG(JL)) THEN
                  IF(ZCOND(JL).NE.0.) THEN
                     IT=INT(PT(JL,KK)*1000.)
                     ZQSAT=TLUCUA(IT)*ZQP(JL)
                     ZQSAT=MIN(0.5,ZQSAT)
                     ZCOR=1./(1.-VTMPC1*ZQSAT)
                     ZQSAT=ZQSAT*ZCOR
                     ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(IT))
                     PT(JL,KK)=PT(JL,KK)+TLUCUC(IT)*ZCOND1
                     PQ(JL,KK)=PQ(JL,KK)-ZCOND1
                  ENDIF
               ENDIF
            ENDDO
C
         ENDIF
C
      ENDIF

      IF(KCALL.EQ.0) THEN
C
         ISUM=0
!DIR$    IVDEP
         DO JL=KIDIA,KFDIA
            IT=INT(PT(JL,KK)*1000.)
            ZQP(JL)=1./PP(JL)
            ZQSAT=TLUCUA(IT)*ZQP(JL)
            ZQSAT=MIN(0.5,ZQSAT)
            ZCOR=1./(1.-VTMPC1*ZQSAT)
            ZQSAT=ZQSAT*ZCOR
            ZCOND(JL)=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(IT))
            PT(JL,KK)=PT(JL,KK)+TLUCUC(IT)*ZCOND(JL)
            PQ(JL,KK)=PQ(JL,KK)-ZCOND(JL)
            IF(ZCOND(JL).NE.0.0) ISUM=ISUM+1
         ENDDO
C
         IF(ISUM.NE.0) THEN
C
!DIR$    IVDEP
            DO JL=KIDIA,KFDIA
               IT=INT(PT(JL,KK)*1000.)
               ZQSAT=TLUCUA(IT)*ZQP(JL)
               ZQSAT=MIN(0.5,ZQSAT)
               ZCOR=1./(1.-VTMPC1*ZQSAT)
               ZQSAT=ZQSAT*ZCOR
               ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(IT))
               PT(JL,KK)=PT(JL,KK)+TLUCUC(IT)*ZCOND1
               PQ(JL,KK)=PQ(JL,KK)-ZCOND1
            ENDDO
C
         ENDIF
C
      ENDIF

      IF(KCALL.EQ.4) THEN

!DIR$    IVDEP
         DO JL=KIDIA,KFDIA
            IT=INT(PT(JL,KK)*1000.)
            ZQP(JL)=1./PP(JL)
C
            ZQSAT=TLUCUA(IT)*ZQP(JL)
C
            ZQSAT=MIN(0.5,ZQSAT)
            ZCOR=1./(1.-VTMPC1*ZQSAT)
            ZQSAT=ZQSAT*ZCOR
            ZCOND(JL)=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(IT))
            PT(JL,KK)=PT(JL,KK)+TLUCUC(IT)*ZCOND(JL)
            PQ(JL,KK)=PQ(JL,KK)-ZCOND(JL)
         ENDDO

!DIR$    IVDEP
         DO JL=KIDIA,KFDIA
            IT=INT(PT(JL,KK)*1000.)
C
            ZQSAT=TLUCUA(IT)*ZQP(JL)
C
            ZQSAT=MIN(0.5,ZQSAT)
            ZCOR=1./(1.-VTMPC1*ZQSAT)
            ZQSAT=ZQSAT*ZCOR
            ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*TLUCUB(IT))
            PT(JL,KK)=PT(JL,KK)+TLUCUC(IT)*ZCOND1
            PQ(JL,KK)=PQ(JL,KK)-ZCOND1
         ENDDO

      ENDIF
C
      RETURN
      END SUBROUTINE CEADJTQ
