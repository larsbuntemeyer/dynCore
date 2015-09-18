C
C**** *LWC*   - LONGWAVE RADIATION, CLOUD EFFECTS
C
C     PURPOSE.
C     --------
C           INTRODUCES CLOUD EFFECTS ON LONGWAVE FLUXES OR
C           RADIANCES
C
C**   INTERFACE.
C     ----------
C      SUBROUTINE LWC ( KLON,KLEV
C    S  ,  PBINT,PBSUIN,PCLDLW,PCNTRB,PEMIS,PFDN,PFUP,PHFG
C    S  ,  PFLUX                                                )
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C     ==== INPUTS ===
C PBINT  : (KLON,0:KLEV)     ; HALF LEVEL PLANCK FUNCTION
C PBSUIN : (KLON)             ; SURFACE PLANCK FUNCTION
C PCLDLW : (KLON,KLEV)       ; CLOUD FRACTIONAL COVER
C PCNTRB : (KLON,0:KLEV,0:KLEV); CLEAR-SKY ENERGY EXCHANGE
C     ==== OUTPUTS ===
C PFLUX(KLON,2,KLEV)         ; RADIATIVE FLUXES :
C PFLUC(KLON,2,KLEV+1)         ; CLEAR SKY FLUX
C                     1  ==>  UPWARD   FLUX TOTAL
C                     2  ==>  DOWNWARD FLUX TOTAL
C
C        IMPLICIT ARGUMENTS :   NONE
C        --------------------
C
C     METHOD.
C     -------
C
C          1. INITIALIZES ALL FLUXES TO CLEAR-SKY VALUES
C          2. EFFECT OF ONE OVERCAST UNITY EMISSIVITY CLOUD LAYER
C          3. EFFECT OF SEMI-TRANSPARENT, PARTIAL OR MULTI-LAYERED
C     CLOUDS
C
C     EXTERNALS.
C     ----------
C
C          NONE
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "IN CORE MODEL"
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C-----------------------------------------------------------------------
      SUBROUTINE LWC ( KLON,KLEV
     &  , PBINT,PBSUIN,PCLDLW,PCNTRB,PFLUC
     &  , PFLUX                                                )
C
      IMPLICIT NONE
C
      INCLUDE "YOMRDU"
      INCLUDE "YOMRDI"
C-----------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
      INTEGER, INTENT(IN)    :: KLON, KLEV
C
      REAL,    INTENT(IN)    :: 
     &     PBINT(KLON,KLEV+1),PBSUIN(KLON),PCLDLW(KLON,KLEV)
     &  ,  PCNTRB(KLON,KLEV+1,KLEV+1)
     &  ,  PFLUC(KLON,2,KLEV+1)
C
      REAL,    INTENT(INOUT) :: PFLUX(KLON,2,KLEV+1)
C
C-----------------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS
C              ------------
C
      INTEGER :: IMX(KLON), IMXP(KLON)
C
      REAL    :: ZCLEAR(KLON),ZCLOUD(KLON),ZDNF(KLON,KLEV+1,KLEV+1)
     &         , ZFD(KLON), ZFU(KLON)
     &         , ZUPF(KLON,KLEV+1,KLEV+1)
      REAL    :: ZCLM(KLON,KLEV+1,KLEV+1)
      INTEGER :: IMAXC, IMX1, IMX2, JCLOUD, JK, JK1, JK2, JKC, JKCP1,
     &           JKJ, JKM1, JKP1, JL
      LOGICAL :: LO1
      REAL    :: ZCFRAC
C
!DIR$ NOBOUNDS PFLUX,PFLUC
C-----------------------------------------------------------------------
C
C     ------------------------------------------------------------------
C
C*         1.     INITIALIZATION
C                 --------------
C
      IMAXC = 0
C
      DO JL = 1 , KLON
         IMX(JL)=0
         IMXP(JL)=0
         ZCLOUD(JL) = 0.
      ENDDO
C
C*         1.1    SEARCH THE LAYER INDEX OF THE HIGHEST CLOUD
C                 -------------------------------------------
C
      DO JK = 1 , KLEV
         DO JL = 1 , KLON
            LO1=PCLDLW(JL,JK).GT.ZEPSC
            IMX1=IMX(JL)
            IMX2=JK
            IF (LO1) THEN
               IMXP(JL)=IMX2
            ELSE
               IMXP(JL)=IMX1
            ENDIF
            IMAXC=MAX0(IMXP(JL),IMAXC)
            IMX(JL)=IMXP(JL)
         ENDDO
      ENDDO
C
      PFLUX(:,:,:) = PFLUC(:,:,:)
C
C     ------------------------------------------------------------------
C
C*         2.      EFFECT OF CLOUDINESS ON LONGWAVE FLUXES
C                  ---------------------------------------
C
CRP   UM, UNABHAENGIG VON KLON, IMMER DIESELBEN ERGEBNISSE ZU ERHALTEN,
CRP   WIRD IMAXC=KLEV GESETZT.
      IMAXC=KLEV
      IF (IMAXC.GT.0) THEN
C
C
C*         2.0     INITIALIZE TO CLEAR-SKY FLUXES
C                  ------------------------------
C
         DO JK1=1,KLEV+1
            DO JK2=1,KLEV+1
               DO JL = 1 , KLON
                  ZUPF(JL,JK2,JK1)=PFLUC(JL,1,JK1)
                  ZDNF(JL,JK2,JK1)=PFLUC(JL,2,JK1)
               ENDDO
            ENDDO
         ENDDO
C
C*         2.1     FLUXES FOR ONE OVERCAST UNITY EMISSIVITY CLOUD
C                  ----------------------------------------------
C
         DO JKC = 1 , IMAXC
            JCLOUD=JKC
            JKCP1=JCLOUD+1
C
C*         2.1.1   ABOVE THE CLOUD
C                  ---------------
C
            DO JK=JKCP1,KLEV+1

               JKM1=JK-1
               DO JL = 1 , KLON
                  ZFU(JL)=0.
               ENDDO
               IF (JK .GT. JKCP1) THEN
                  DO JKJ=JKCP1,JKM1
                     DO JL = 1 , KLON
                        ZFU(JL) = ZFU(JL) + PCNTRB(JL,JK,JKJ)
                     ENDDO
                  ENDDO
               ENDIF
C
               DO JL = 1 , KLON
                  ZUPF(JL,JKCP1,JK)=PBINT(JL,JK)-ZFU(JL)
               ENDDO
            ENDDO
C
C*         2.1.2   BELOW THE CLOUD
C                  ---------------
C
            DO JK=1,JCLOUD
               JKP1=JK+1
               DO JL = 1 , KLON
                  ZFD(JL)=0.
               ENDDO
C
               IF (JK .LT. JCLOUD) THEN
                  DO JKJ=JKP1,JCLOUD
                     DO JL = 1 , KLON
                        ZFD(JL) = ZFD(JL) + PCNTRB(JL,JK,JKJ)
                     ENDDO
                  ENDDO
               ENDIF
               DO JL = 1 , KLON
                  ZDNF(JL,JKCP1,JK)=-PBINT(JL,JK)-ZFD(JL)
               ENDDO
            ENDDO
C
         ENDDO
C
C
C*         2.2     CLOUD COVER MATRIX
C                  ------------------
C
C*    ZCLM(JK1,JK2) IS THE OBSCURATION FACTOR BY CLOUD LAYERS BETWEEN
C     HALF-LEVELS JK1 AND JK2 AS SEEN FROM JK1
C
         DO JK1 = 1 , KLEV+1
            DO JK2 = 1 , KLEV+1
               DO JL = 1 , KLON
                  ZCLM(JL,JK1,JK2) = 0.
               ENDDO
            ENDDO
         ENDDO
C
C
C*         2.4     CLOUD COVER BELOW THE LEVEL OF CALCULATION
C                  ------------------------------------------
C
         DO JK1 = 2 , KLEV+1
            DO JL = 1 , KLON
               ZCLEAR(JL)=1.
               ZCLOUD(JL)=0.
            ENDDO
            DO JK = JK1 - 1 , 1 , -1
               DO JL = 1 , KLON
                  ZCLEAR(JL)=ZCLEAR(JL)*
     &                 (1.0-AMAX1(PCLDLW(JL,JK),ZCLOUD(JL)))
     &                 /(1.0-AMIN1(ZCLOUD(JL),1.-ZEPSEC))
                  ZCLM(JL,JK1,JK) = 1.0 - ZCLEAR(JL)
                  ZCLOUD(JL) = PCLDLW(JL,JK)
               ENDDO
            ENDDO
         ENDDO
C
C
C*         2.5     CLOUD COVER ABOVE THE LEVEL OF CALCULATION
C                  ------------------------------------------
C
         DO JK1 = 1 , KLEV
            DO JL = 1 , KLON
               ZCLEAR(JL)=1.
               ZCLOUD(JL)=0.
            ENDDO
            DO JK = JK1 , KLEV
               DO JL = 1 , KLON
                  ZCLEAR(JL)=ZCLEAR(JL)*
     &                 (1.0-AMAX1(PCLDLW(JL,JK),ZCLOUD(JL)))
     &                 /(1.0-AMIN1(ZCLOUD(JL),1.-ZEPSEC))
                  ZCLM(JL,JK1,JK) = 1.0 - ZCLEAR(JL)
                  ZCLOUD(JL) = PCLDLW(JL,JK)
               ENDDO
            ENDDO
         ENDDO
C
C
C*         3.      FLUXES FOR PARTIAL/MULTIPLE LAYERED CLOUDINESS
C                  ----------------------------------------------
C
C*         3.1     DOWNWARD FLUXES
C                  ---------------
C
         DO JL = 1 , KLON
            PFLUX(JL,2,KLEV+1) = 0.
         ENDDO
C
         DO JK1 = KLEV , 1 , -1
C
C*                 CONTRIBUTION FROM CLEAR-SKY FRACTION
C
            DO JL = 1 , KLON
               ZFD (JL) = (1. - ZCLM(JL,JK1,KLEV)) * ZDNF(JL,1,JK1)
            ENDDO
C
C*                 CONTRIBUTION FROM ADJACENT CLOUD
C
            DO JL = 1 , KLON
               ZFD(JL) = ZFD(JL) + ZCLM(JL,JK1,JK1) * ZDNF(JL,JK1+1,JK1)
            ENDDO
C
C*                 CONTRIBUTION FROM OTHER CLOUDY FRACTIONS
C
            DO JK = KLEV-1 , JK1 , -1
               DO JL = 1 , KLON
                  ZCFRAC = ZCLM(JL,JK1,JK+1) - ZCLM(JL,JK1,JK)
                  ZFD(JL) =  ZFD(JL) + ZCFRAC * ZDNF(JL,JK+2,JK1)
               ENDDO
            ENDDO
C
            DO JL = 1 , KLON
               PFLUX(JL,2,JK1) = ZFD (JL)
            ENDDO
C
         ENDDO
C
C
C
C
C*         3.2     UPWARD FLUX AT THE SURFACE
C                  --------------------------
C
         DO JL = 1 , KLON
            PFLUX(JL,1,1) = ZEMISS*PBSUIN(JL)-(1.-ZEMISS)*PFLUX(JL,2,1)
         ENDDO
C
C
C
C*         3.3     UPWARD FLUXES
C                  -------------
C
         DO JK1 = 2 , KLEV+1
C
C*                 CONTRIBUTION FROM CLEAR-SKY FRACTION
C
            DO JL = 1 , KLON
               ZFU (JL) = (1. - ZCLM(JL,JK1,1)) * ZUPF(JL,1,JK1)
            ENDDO
C
C*                 CONTRIBUTION FROM ADJACENT CLOUD
C
            DO JL = 1 , KLON
               ZFU(JL) = ZFU(JL) + ZCLM(JL,JK1,JK1-1) * ZUPF(JL,JK1,JK1)
            ENDDO
C
C*                 CONTRIBUTION FROM OTHER CLOUDY FRACTIONS
C
            DO JK = 2 , JK1-1
               DO JL = 1 , KLON
                  ZCFRAC = ZCLM(JL,JK1,JK-1) - ZCLM(JL,JK1,JK)
                  ZFU(JL) = ZFU(JL) + ZCFRAC * ZUPF(JL,JK  ,JK1)
               ENDDO
            ENDDO
C
            DO JL = 1 , KLON
               PFLUX(JL,1,JK1) = ZFU (JL)
            ENDDO
C
         ENDDO
C
C
      ENDIF
C
C
C*         2.3     END OF CLOUD EFFECT COMPUTATIONS
C
C
      RETURN
C
      END SUBROUTINE LWC
