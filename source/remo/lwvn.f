C
C
C**** *LWVN*   - L.W., VERTICAL INTEGRATION, NEARBY LAYERS
C
C     PURPOSE.
C     --------
C           CARRIES OUT THE VERTICAL INTEGRATION ON NEARBY LAYERS
C           TO GIVE LONGWAVE FLUXES OR RADIANCES
C
C**   INTERFACE.
C     ----------
C     *CALL*     LWVN ( KLON,KLEV,KUAER,KTRAER
C    S  , KAER,KCFC
C    R  , PABCU,PDBSL,PGA,PGB,PGC,PGD
C    S  , PADJD,PADJU,PCNTRB,PDBDT                              )
C    S  , PADJD,PADJU,PCNTRB,PDBDT,PDWFSU )
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C     ==== INPUTS ===
C PABCU : (KLON,NUA,3*KLEV+1); ABSORBER AMOUNTS
C PDBSL  : (KLON,KLEV*2)     ; SUB-LAYER PLANCK FUNCTION GRADIENT
C     ==== OUTPUTS ===
C PADJ.. : (KLON,KLEV+1)     ; CONTRIBUTION OF ADJACENT LAYERS
C PCNTRB : (KLON,KLEV+1,KLEV+1); CLEAR-SKY ENERGY EXCHANGE MATRIX
C PDBDT  : (KLON,NUA,KLEV)   ; LAYER PLANCK FUNCTION GRADIENT
C PDWFSU : (KLON,NINT)        ; DOWNWARD BAND FLUX AT SURFACE
C
C        IMPLICIT ARGUMENTS :   NONE
C        --------------------
C
C     METHOD.
C     -------
C
C          1. PERFORMS THE VERTICAL INTEGRATION CORRESPONDING TO THE
C     CONTRIBUTIONS OF THE ADJACENT LAYERS USING A GAUSSIAN QUADRATURE
C
C     EXTERNALS.
C     ----------
C
C          *LWTT*
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
C      MODIFIED M. GIORGETTA  U. SCHLESE  JUL-93
C-----------------------------------------------------------------------
      SUBROUTINE LWVN ( KLON,KLEV,KUAER
     &  , KAER,KCFC
     &  , PABCU,PDBSL,PGA,PGB,PGC,PGD
     &  , PADJD,PADJU,PCNTRB,PDBDT,PDWFSU )
C
      IMPLICIT NONE
C
      INCLUDE "YOMLW"
C-----------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
      INTEGER, INTENT(IN)    :: KLON,KLEV,KUAER
     &                        , KAER,KCFC
      REAL,    INTENT(IN)    :: PABCU(KLON,NUA,3*KLEV+1)
     &                       ,  PDBSL(KLON,NINT,KLEV*2)
     &                       ,  PGA(KLON,8,2,KLEV), PGB(KLON,8,2,KLEV)
     &                       ,  PGC(KLON,5,2,KLEV), PGD(KLON,5,2,KLEV)
C
      REAL,    INTENT(INOUT) :: PADJD(KLON,KLEV+1), PADJU(KLON,KLEV+1)
     &                       ,  PCNTRB(KLON,KLEV+1,KLEV+1)
     &                       ,  PDBDT(KLON,NINT,KLEV)
     &                       ,  PDWFSU(KLON,NINT)
C
C-----------------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS
C              ------------
C
C
      REAL :: ZGLAYD(KLON),   ZGLAYU(KLON)
     &     ,  ZTT(KLON,NTRA), ZUU(KLON,NUA)
C
      INTEGER :: IG, ITMP2, JK, JK1, JK2, JL, KBS, KDD, KM12, KMU, KND, 
     &           KNU, KXD, KXU
      REAL    :: ZWTR
C
C-----------------------------------------------------------------------
!DIR$ NOBOUNDS
C
C*         1.    INITIALIZATION
C                --------------
C
C*         1.1     INITIALIZE LAYER CONTRIBUTIONS
C                  ------------------------------
C
      PADJD(:,:) = 0.
      PADJU(:,:) = 0.
C
C     ------------------------------------------------------------------
C
C*         2.      VERTICAL INTEGRATION
C                  --------------------
C
C
C*         2.1     CONTRIBUTION FROM ADJACENT LAYERS
C                  ---------------------------------
C
      DO JK = 1 , KLEV
C
C*         2.1.1   DOWNWARD LAYERS
C                  ---------------
C
         KM12 = 2 * (JK - 1)
         KND = (JK - 1) * NG1P1 + 1
         KXD = KND
         KNU = JK * NG1P1 + 1
         KXU = KND
C
         DO JL = 1 , KLON
            ZGLAYD(JL) = 0.
            ZGLAYU(JL) = 0.
         ENDDO
C
         DO IG = 1 , NG1
            KBS = KM12 + IG
            KDD = KXD + IG
            ZUU(:,1:KUAER) = PABCU(:,1:KUAER,KND) - PABCU(:,1:KUAER,KDD)
C
C
            CALL LWTT (KLON,KAER,KCFC,PGA(1,1,1,JK),PGB(1,1,1,JK)
     &           , PGC(1,1,1,JK),PGD(1,1,1,JK), ZUU,ZTT)
C
            IF (JK.EQ.1) THEN
               DO JL=1,KLON
                  PDWFSU(JL,1)=PDWFSU(JL,1)+
     &            WG1(IG)*PDBSL(JL,1,KBS)*ZTT(JL,1)          *ZTT(JL,10)
                  PDWFSU(JL,2)=PDWFSU(JL,2)+
     &            WG1(IG)*PDBSL(JL,2,KBS)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
                  PDWFSU(JL,3)=PDWFSU(JL,3)+
     &            WG1(IG)*PDBSL(JL,3,KBS)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
                  PDWFSU(JL,4)=PDWFSU(JL,4)+
     &            WG1(IG)*PDBSL(JL,4,KBS)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
                  PDWFSU(JL,5)=PDWFSU(JL,5)+
     &            WG1(IG)*PDBSL(JL,5,KBS)*ZTT(JL,3)          *ZTT(JL,14)
                  PDWFSU(JL,6)=PDWFSU(JL,6)+
     &            WG1(IG)*PDBSL(JL,6,KBS)*ZTT(JL,6)          *ZTT(JL,15)
               ENDDO
            ENDIF
C
            DO JL = 1 , KLON
               ZWTR=PDBSL(JL,1,KBS)*ZTT(JL,1)          *ZTT(JL,10)
     &             +PDBSL(JL,2,KBS)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     &             +PDBSL(JL,3,KBS)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     &             +PDBSL(JL,4,KBS)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     &             +PDBSL(JL,5,KBS)*ZTT(JL,3)          *ZTT(JL,14)
     &             +PDBSL(JL,6,KBS)*ZTT(JL,6)          *ZTT(JL,15)
               ZGLAYD(JL)=ZGLAYD(JL)+ZWTR*WG1(IG)
            ENDDO
C
C*         2.1.2   UPWARD LAYERS
C                  ---------------
C
            KMU = KXU + IG
            DO itmp2 = 1, KUAER
               ZUU(:,itmp2) = PABCU(:,itmp2,KMU) - PABCU(:,itmp2,KNU)
            ENDDO
C
C
            CALL LWTT (KLON,KAER,KCFC,PGA(1,1,1,JK),PGB(1,1,1,JK)
     &           , PGC(1,1,1,JK),PGD(1,1,1,JK), ZUU,ZTT)
C
            DO JL = 1 , KLON
               ZWTR=PDBSL(JL,1,KBS)*ZTT(JL,1)          *ZTT(JL,10)
     &             +PDBSL(JL,2,KBS)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     &             +PDBSL(JL,3,KBS)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     &             +PDBSL(JL,4,KBS)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     &             +PDBSL(JL,5,KBS)*ZTT(JL,3)          *ZTT(JL,14)
     &             +PDBSL(JL,6,KBS)*ZTT(JL,6)          *ZTT(JL,15)
               ZGLAYU(JL)=ZGLAYU(JL)+ZWTR*WG1(IG)
            ENDDO
C
         ENDDO
C
         DO JL = 1 , KLON
            PADJD(JL,JK) = ZGLAYD(JL)
            PCNTRB(JL,JK,JK+1) = ZGLAYD(JL)
            PADJU(JL,JK+1) = ZGLAYU(JL)
            PCNTRB(JL,JK+1,JK) = ZGLAYU(JL)
            PCNTRB(JL,JK  ,JK) = 0.0
         ENDDO
C
      ENDDO
C
      DO JK = 1 , KLEV
         JK2 = 2 * JK
         JK1 = JK2 - 1
         PDBDT(:,:,JK) = PDBSL(:,:,JK1) + PDBSL(:,:,JK2)
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END SUBROUTINE LWVN
