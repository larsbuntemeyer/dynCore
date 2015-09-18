C
C
C**** *LWVD*   - L.W., VERTICAL INTEGRATION, DISTANT LAYERS
C
C     PURPOSE.
C     --------
C           CARRIES OUT THE VERTICAL INTEGRATION ON THE DISTANT LAYERS
C
C**   INTERFACE.
C     ----------
C     *CALL*     LWVD ( KLON,KLEV,KUAER,KTRAER
C    S  , KAER,KCFC
C    S  , KXDIA,KXT
C    S  , PABCU,PDBDT
C    R  , PGA,PGB,PGC,PGD
C    S  , PCNTRB,PDISD,PDISU                              )
C    S  , PCNTRB,PDISD,PDISU,PDWFSU )
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C     ==== INPUTS ===
C PABCU : (KLON,NUA,3*KLEV+1); ABSORBER AMOUNTS
C PDBDT  : (KLON,KLEV)       ; LAYER PLANCK FUNCTION GRADIENT
C     ==== OUTPUTS ===
C PDIS.. : (KLON,KLEV+1)     ; CONTRIBUTION BY DISTANT LAYERS
C PCNTRB : (KLON,KLEV+1,KLEV+1); ENERGY EXCHANGE MATRIX
C PDWFSU : (KLON,NINT)        ; DOWNWARD BAND FLUX AT SURFACE
C
C        IMPLICIT ARGUMENTS :   NONE
C        --------------------
C
C     METHOD.
C     -------
C
C          1. PERFORMS THE VERTICAL INTEGRATION CORRESPONDING TO THE
C     CONTRIBUTIONS OF THE DISTANT LAYERS USING TRAPEZOIDAL RULE
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

C       MODIFIED  U. SCHLESE, M. GIORGETTA   JUL-93
C-----------------------------------------------------------------------
      SUBROUTINE LWVD ( KLON,KLEV,KUAER,KTRAER
     &  , KAER,KCFC
     &  , PABCU,PDBDT
     &  , PGA,PGB,PGC,PGD
     &  , PCNTRB,PDISD,PDISU,PDWFSU )
C
      IMPLICIT NONE
C
      INCLUDE "YOMLW"
!DIR$ NOBOUNDS
C-----------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
      INTEGER, INTENT(IN)    :: KLON,KLEV,KUAER,KTRAER
     &                     ,    KAER,KCFC
C
      REAL,    INTENT(IN)    :: PABCU(KLON,NUA,3*KLEV+1)
     &                    ,     PDBDT(KLON,NINT,KLEV)
     &                    ,     PGA(KLON,8,2,KLEV) , PGB(KLON,8,2,KLEV)
     &                    ,     PGC(KLON,5,2,KLEV),PGD(KLON,5,2,KLEV)
C
      REAL,    INTENT(INOUT) :: PCNTRB(KLON,KLEV+1,KLEV+1)
     &                    ,     PDISD(KLON,KLEV+1), PDISU(KLON,KLEV+1)
     &                    ,     PDWFSU(KLON,NINT)
C
C-----------------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS
C              ------------
C
      REAL    ::
     &     ZTT(KLON,NTRA), ZTT1(KLON,NTRA), ZTT2(KLON,NTRA)
     &  ,  ZUU(KLON,NUA)
C
      INTEGER :: 
     &     JK, JKJ, JKL, JKP1, JL, JLK, KD1, KD2, KJ, KJP1, KM1, KN, 
     &     KU1, KU2
      REAL    ::
     &     ZT1, ZT10, ZT11, ZT12, ZT13, ZT14, ZT15, ZT2, ZT3, ZT4, ZT5, 
     &     ZT6, ZT7, ZT8, ZT9, ZWW
C
C-----------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------
C
C*         1.1     INITIALIZE LAYER CONTRIBUTIONS
C                  ------------------------------
C
      PDISD(:,:) = 0.
      PDISU(:,:) = 0.
C
C     ------------------------------------------------------------------
C
C*         2.      VERTICAL INTEGRATION
C                  --------------------
C
C
C*         2.2     CONTRIBUTION FROM DISTANT LAYERS
C                  ---------------------------------
C
C
C*         2.2.1   DISTANT AND ABOVE LAYERS
C                  ------------------------
C
C
C*         2.2.2   FIRST UPPER LEVEL
C                  -----------------
C
      DO JK = 1 , KLEV-1
         JKP1=JK+1
         KN=(JK-1)*NG1P1+1
         KD1= JK  *NG1P1+1
C
C
         ZUU(:,1:KUAER) = PABCU(:,1:KUAER,KN) - PABCU(:,1:KUAER,KD1)
C
C
         CALL LWTT (KLON,KAER,KCFC,PGA(1,1,1,JK),PGB(1,1,1,JK)
     &        , PGC(1,1,1,JK),PGD(1,1,1,JK), ZUU,ZTT1)
C
C
C*         2.2.3   HIGHER UP
C                  ---------
C
         DO JKJ=JKP1,KLEV
            KJP1=JKJ+1
            KD2= JKJ  *NG1P1+1
C
C
            ZUU(:,1:KUAER) = PABCU(:,1:KUAER,KN) - PABCU(:,1:KUAER,KD2)
C
C
            CALL LWTT (KLON,KAER,KCFC,PGA(1,1,1,JKJ),PGB(1,1,1,JKJ)
     &           , PGC(1,1,1,JKJ),PGD(1,1,1,JKJ), ZUU,ZTT2)
C
            ZTT(:,1:KTRAER) = (ZTT1(:,1:KTRAER)+ZTT2(:,1:KTRAER))*0.5
C
            IF (JK.EQ.1) THEN
               DO JL=1,KLON
                  PDWFSU(JL,1)=PDWFSU(JL,1)+
     &                 PDBDT(JL,1,JKJ)*ZTT(JL,1)          *ZTT(JL,10)
                  PDWFSU(JL,2)=PDWFSU(JL,2)+
     &                 PDBDT(JL,2,JKJ)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
                  PDWFSU(JL,3)=PDWFSU(JL,3)+
     &                 PDBDT(JL,3,JKJ)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
                  PDWFSU(JL,4)=PDWFSU(JL,4)+
     &                 PDBDT(JL,4,JKJ)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
                  PDWFSU(JL,5)=PDWFSU(JL,5)+
     &                 PDBDT(JL,5,JKJ)*ZTT(JL,3)          *ZTT(JL,14)
                  PDWFSU(JL,6)=PDWFSU(JL,6)+
     &                 PDBDT(JL,6,JKJ)*ZTT(JL,6)          *ZTT(JL,15)
               ENDDO
            ENDIF
C
C
            DO JL = 1 , KLON
               ZWW=PDBDT(JL,1,JKJ)*ZTT(JL,1)          *ZTT(JL,10)
     &            +PDBDT(JL,2,JKJ)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     &            +PDBDT(JL,3,JKJ)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     &            +PDBDT(JL,4,JKJ)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     &            +PDBDT(JL,5,JKJ)*ZTT(JL,3)          *ZTT(JL,14)
     &            +PDBDT(JL,6,JKJ)*ZTT(JL,6)          *ZTT(JL,15)
C
               PDISD(JL,JK)=PDISD(JL,JK)+ZWW
               PCNTRB(JL,JK,KJP1)=ZWW
            ENDDO
C
            ZTT1(:,1:KTRAER) = ZTT2(:,1:KTRAER)
C
         ENDDO
      ENDDO
C
C
C*         2.2.4   DISTANT AND BELOW LAYERS
C                  ------------------------
C
C
C*         2.2.5   FIRST LOWER LEVEL
C                  -----------------
C
      DO JK=3,KLEV+1
         KN=(JK-1)*NG1P1+1
         KM1=JK-1
         KJ=JK-2
         KU1= KJ  *NG1P1+1
C
C
         ZUU(:,1:KUAER) = PABCU(:,1:KUAER,KU1) - PABCU(:,1:KUAER,KN)
C
C
         CALL LWTT (KLON,KAER,KCFC,PGA(1,1,1,KM1),PGB(1,1,1,KM1)
     &        ,PGC(1,1,1,KM1),PGD(1,1,1,KM1),ZUU,ZTT1)
C
C*         2.2.6   DOWN BELOW
C                  ----------
C
         DO JLK=1,KJ
            JKL=KM1-JLK
            KU2=(JKL-1)*NG1P1+1
C
C
            ZUU(:,1:KUAER) = PABCU(:,1:KUAER,KU2) - PABCU(:,1:KUAER,KN)
C
C
            CALL LWTT (KLON,KAER,KCFC,PGA(1,1,1,JKL),PGB(1,1,1,JKL)
     &           , PGC(1,1,1,JKL),PGD(1,1,1,JKL), ZUU,ZTT2)
C
C
            DO JL = 1 , KLON
               ZT1 =(ZTT1(JL,1) +ZTT2(JL,1))*0.5
               ZT2 =(ZTT1(JL,2) +ZTT2(JL,2))*0.5
               ZT3 =(ZTT1(JL,3) +ZTT2(JL,3))*0.5
               ZT4 =(ZTT1(JL,4) +ZTT2(JL,4))*0.5
               ZT5 =(ZTT1(JL,5) +ZTT2(JL,5))*0.5
               ZT6 =(ZTT1(JL,6) +ZTT2(JL,6))*0.5
               ZT7 =(ZTT1(JL,7) +ZTT2(JL,7))*0.5
               ZT8 =(ZTT1(JL,8) +ZTT2(JL,8))*0.5
               ZT9 =(ZTT1(JL,9) +ZTT2(JL,9))*0.5
               ZT10=(ZTT1(JL,10)+ZTT2(JL,10))*0.5
               ZT11=(ZTT1(JL,11)+ZTT2(JL,11))*0.5
               ZT12=(ZTT1(JL,12)+ZTT2(JL,12))*0.5
               ZT13=(ZTT1(JL,13)+ZTT2(JL,13))*0.5
               ZT14=(ZTT1(JL,14)+ZTT2(JL,14))*0.5
               ZT15=(ZTT1(JL,15)+ZTT2(JL,15))*0.5

               ZWW=PDBDT(JL,1,JKL)*ZT1         *ZT10
     &            +PDBDT(JL,2,JKL)*ZT2*ZT7*ZT11
     &            +PDBDT(JL,3,JKL)*ZT4*ZT8*ZT12
     &            +PDBDT(JL,4,JKL)*ZT5*ZT9*ZT13
     &            +PDBDT(JL,5,JKL)*ZT3         *ZT14
     &            +PDBDT(JL,6,JKL)*ZT6         *ZT15
               PDISU(JL,JK)=PDISU(JL,JK)+ZWW
               PCNTRB(JL,JK,JKL)=ZWW
            ENDDO
C
            ZTT1(:,1:KTRAER) = ZTT2(:,1:KTRAER)
C
         ENDDO
      ENDDO
C
      RETURN
C
      END SUBROUTINE LWVD
