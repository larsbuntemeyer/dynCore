C
C**** *LEGTRI* - *LEGENDRE FUNCTIONS FOR A TRIANGULAR TRUNCATION.
C
C     J.F.GELEYN     E.C.M.W.F.     03/06/82.
C
C     MODIFICATIONS TO FIT INTO HIRLAM
C     J.H. CHRISTENSEN     DMI      10/12/91.
C         - 00 -                    15/5/93. REWRITTEN
C
C     PURPOSE.
C     --------
C
C          THIS ROUTINE COMPUTES THE VALUES *PALP* FOR THE ARGUMENT
C     *PSIN* OF THE NORMALISED *LEGENDRE ASSOCIATED FUNCTIONS IN TH
C     ORDER ((JN1=JM1,KCP),JM1=1,KCP) FOR JN=JN1-1 AND JM=JM1-1 .
C
C**   INTERFACE.
C     ----------
C
C          *LEGTRI* IS CALLED FROM *RADMOD*.
C          THERE ARE THREE DUMMY ARGUMENTS: *PSIN* IS THE SINE OF
C     LATITUDE.
C                                           *KCP* IS ONE PLUS THE LIMIT
C     WAVE NUMBER.
C                                           *PALP* IS THE ARRAY OF THE
C     RESULTS.
C
C     METHOD.
C     -------
C
C          SIMPLE RECURENCE FORMULA.
C
C     EXTERNALS.
C     ----------
C
C          NONE.
C
C     REFERENCE.
C     ----------
C
C          NONE.
C
      SUBROUTINE LEGTRI(KLON,KSTART,KSTOP,KCP,
     &       SINLAT,COSLAT,PALP,ZF3M,JJ)
C
      IMPLICIT NONE
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: KLON,KSTART,KSTOP,KCP
      REAL,    INTENT(IN)    :: COSLAT(KLON),SINLAT(KLON)
      REAL,    INTENT(INOUT) :: PALP(KLON,(KCP+1)*KCP/2)
C
C*    DECLARATION OF WORK ARRAYS
C
      INTEGER, INTENT(INOUT) :: JJ(KLON)
      REAL,    INTENT(INOUT) :: ZF3M(KLON)
C
C     Local Variabales
C
      INTEGER :: IC, ICP, JL, JM, JM1, JM2, JN
      REAL    :: Z2M, ZE1, ZE2, ZF1M, ZF2M, ZM, ZN, ZN2, ZRE1
C
C
C     ------------------------------------------------------------------
C
C*         1.     PRELIMINARY SETTING.
C                 ----------- --------
C
      ICP=KCP
C
C     ------------------------------------------------------------------
C
C*         2.     COMPUTATIONS.
C                 -------------
C
      IC=ICP-1
      ZF1M=SQRT(3.)
      DO JL=KSTART,KSTOP
         JJ(JL)=2
         PALP(JL,1)=1.
         PALP(JL,2)=ZF1M*SINLAT(JL)
         ZF3M(JL)=ZF1M
      ENDDO

      DO JM1=1,ICP
         JM=JM1-1
         ZM=JM
         Z2M=ZM+ZM
         ZRE1=SQRT(Z2M+3.)
         ZE1=1./ZRE1
         IF(JM.NE.0) THEN
            IF(JM.LT.IC) THEN
               DO JL=KSTART,KSTOP
                  ZF2M=ZF3M(JL)*COSLAT(JL)/SQRT(Z2M)
                  ZF3M(JL)=ZF2M*ZRE1
                  JJ(JL)=JJ(JL)+1
                  PALP(JL,JJ(JL))=ZF2M
                  JJ(JL)=JJ(JL)+1
                  PALP(JL,JJ(JL))=ZF3M(JL)*SINLAT(JL)
               ENDDO
            ELSE
               DO JL=KSTART,KSTOP
                  ZF2M=ZF3M(JL)*COSLAT(JL)/SQRT(Z2M)
                  ZF3M(JL)=ZF2M*ZRE1
                  JJ(JL)=JJ(JL)+1
                  PALP(JL,JJ(JL))=ZF2M
               ENDDO
            ENDIF
            IF(JM.EQ.IC) CYCLE
            IF(JM1.EQ.IC) CYCLE
         ENDIF
         JM2=JM+2
         DO JN=JM2,IC
            ZN=JN
            ZN2=ZN**2
            ZE2=SQRT((4.*ZN2-1.)/(ZN2-ZM**2))
            DO JL=KSTART,KSTOP
               JJ(JL)=JJ(JL)+1
               PALP(JL,JJ(JL))=ZE2*(SINLAT(JL)*PALP(JL,JJ(JL)-1)
     &                                    -ZE1*PALP(JL,JJ(JL)-2))
            ENDDO
            ZE1=1./ZE2
         ENDDO
      ENDDO
C
C     ------------------------------------------------------------------
C
C*         3.     RETURN.
C                 -------
C
      RETURN
      END SUBROUTINE LEGTRI
