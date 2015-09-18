C
C
C**** *LWVB*   - L.W., VERTICAL INTEGRATION, EXCHANGE WITH BOUNDARIES
C
C     PURPOSE.
C     --------
C           INTRODUCES THE EFFECTS OF THE BOUNDARIES IN THE VERTICAL
C           INTEGRATION
C
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C     ==== INPUTS ===
C KX...    (KLON,...          ; TEMPERATURE INDICES
C PABCU : (KLON,NUA,3*KLEV+1); ABSORBER AMOUNTS
C PADJ.. : (KLON,KLEV+1)     ; CONTRIBUTION BY ADJACENT LAYERS
C PBINT  : (KLON,KLEV+1)     ; HALF-LEVEL PLANCK FUNCTIONS
C PBTOP  : (KLON,NINT)        ; T.O.A. SPECTRAL PLANCK FUNCTION
C PDIS.. : (KLON,KLEV+1)     ; CONTRIBUTION BY DISTANT LAYERS
C     ==== OUTPUTS ===
C  IF KMODE = 0, 1, 2
C PFLUC(KLON,2,KLEV)         ; RADIATIVE FLUXES CLEAR-SKY:
C                     1  ==>  UPWARD   FLUX TOTAL
C PDWFSU : (KLON,NINT)        ; DOWNWARD BAND FLUX AT SURFACE
C
C        IMPLICIT ARGUMENTS :   NONE
C        --------------------
C
C     METHOD.
C     -------
C
C          1. COMPUTES THE ENERGY EXCHANGE WITH TOP AND SURFACE OF THE
C     ATMOSPHERE
C          2. COMPUTES THE COOLING-TO-SPACE AND HEATING-FROM-GROUND
C     TERMS FOR THE APPROXIMATE COOLING RATE ABOVE 10 HPA
C          3. ADDS UP ALL CONTRIBUTIONS TO GET THE CLEAR-SKY FLUXES
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
C    MODIFIED   U. SCHLESE  M. GIORGETTA
C-----------------------------------------------------------------------
      SUBROUTINE LWVB ( KLON,KLEV,KUAER,KAER,KCFC
     &  , PABCU,PADJD,PADJU,PBINT,PBTOP
     &  , PDISD,PDISU
     &  , PGASUR,PGBSUR,PGATOP,PGBTOP
     &  , PGCSUR,PGDSUR,PGCTOP,PGDTOP
     &  , PFLUC,PDWFSU )
C
      IMPLICIT NONE
C
      INCLUDE "YOMLW"
      INCLUDE "YOMRDI"
C-----------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
C
      INTEGER, INTENT(IN)    :: KLON,KLEV,KUAER,KAER,KCFC
C
      REAL,    INTENT(IN)    :: PABCU(KLON,NUA,3*KLEV+1)
     &                       ,  PADJD(KLON,KLEV+1), PADJU(KLON,KLEV+1)
     &                       ,  PBINT(KLON,KLEV+1), PBTOP(KLON,NINT)
     &                       ,  PDISD(KLON,KLEV+1), PDISU(KLON,KLEV+1)
     &                       ,  PGASUR(KLON,8,2), PGBSUR(KLON,8,2)
     &                       ,  PGATOP(KLON,8,2), PGBTOP(KLON,8,2)
     &                       ,  PGCSUR(KLON,5,2),PGDSUR(KLON,5,2)
     &                       ,  PGCTOP(KLON,5,2),PGDTOP(KLON,5,2)

C
      REAL,    INTENT(INOUT) :: PFLUC(KLON,2,KLEV+1),  PDWFSU(KLON,NINT)
C
C-----------------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS
C              ------------
C
      REAL    :: ZTT(KLON,NTRA), ZUU(KLON,NUA)
C
      INTEGER :: JK, JL, KN
      REAL    :: ZCNSOL, ZCNTOP
C
C
!DIR$ NOBOUNDS
C     ------------------------------------------------------------------
C
C*         2.      VERTICAL INTEGRATION
C                  --------------------
C
C
C*         2.3     EXCHANGE WITH TOP OF THE ATMOSPHERE
C                  -----------------------------------
C
      DO JK = 1 , KLEV
         KN=(JK-1)*NG1P1+1
C
         ZUU(:,1:KUAER)=PABCU(:,1:KUAER,KN)
C
C
         CALL LWTT (KLON,KAER,KCFC,PGATOP(1,1,1),PGBTOP(1,1,1)
     &        ,        PGCTOP(1,1,1),PGDTOP(1,1,1),ZUU,ZTT)
C
         IF (JK.EQ.1) THEN
            DO JL=1,KLON
               PDWFSU(JL,1)=PDWFSU(JL,1)-
     &              PBTOP(JL,1)*ZTT(JL,1)          *ZTT(JL,10)
               PDWFSU(JL,2)=PDWFSU(JL,2)-
     &              PBTOP(JL,2)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
               PDWFSU(JL,3)=PDWFSU(JL,3)-
     &              PBTOP(JL,3)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
               PDWFSU(JL,4)=PDWFSU(JL,4)-
     &              PBTOP(JL,4)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
               PDWFSU(JL,5)=PDWFSU(JL,5)-
     &              PBTOP(JL,5)*ZTT(JL,3)          *ZTT(JL,14)
               PDWFSU(JL,6)=PDWFSU(JL,6)-
     &              PBTOP(JL,6)*ZTT(JL,6)          *ZTT(JL,15)
            ENDDO
         ENDIF
C
         DO JL = 1 , KLON
            ZCNTOP=PBTOP(JL,1)*ZTT(JL,1)          *ZTT(JL,10)
     &            +PBTOP(JL,2)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     &            +PBTOP(JL,3)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     &            +PBTOP(JL,4)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     &            +PBTOP(JL,5)*ZTT(JL,3)          *ZTT(JL,14)
     &            +PBTOP(JL,6)*ZTT(JL,6)          *ZTT(JL,15)
C
            PFLUC(JL,2,JK)=ZCNTOP-PBINT(JL,JK)-PDISD(JL,JK)-PADJD(JL,JK)
         ENDDO
C
      ENDDO
C
      JK = KLEV+1
      KN=(JK-1)*NG1P1+1
C
      DO JL = 1 , KLON
         ZCNTOP= PBTOP(JL,1)
     &         + PBTOP(JL,2)
     &         + PBTOP(JL,3)
     &         + PBTOP(JL,4)
     &         + PBTOP(JL,5)
     &         + PBTOP(JL,6)
         PFLUC(JL,2,JK)=ZCNTOP-PBINT(JL,JK)-PDISD(JL,JK)-PADJD(JL,JK)
      ENDDO
C
C
C*         2.5     EXCHANGE WITH LOWER LIMIT
C                  -------------------------
C
      JK = 1
      KN=(JK-1)*NG1P1+1
C
      DO JL = 1 , KLON
C
         ZCNSOL=PDWFSU(JL,1)+PDWFSU(JL,2)+PDWFSU(JL,3)
     &         +PDWFSU(JL,4)+PDWFSU(JL,5)+PDWFSU(JL,6)
         ZCNSOL=(1.-ZEMISS)*ZCNSOL
C
         PFLUC(JL,1,JK)=ZCNSOL+PBINT(JL,JK)-PDISU(JL,JK)-PADJU(JL,JK)
      ENDDO
C
      DO JK = 2 , KLEV+1
         KN=(JK-1)*NG1P1+1
C
C
         ZUU(:,1:KUAER)=PABCU(:,1:KUAER,1)-PABCU(:,1:KUAER,KN)
C
C
         CALL LWTT (KLON,KAER,KCFC,PGASUR(1,1,1),PGBSUR(1,1,1)
     &        ,        PGCSUR(1,1,1),PGDSUR(1,1,1),ZUU,ZTT)
C
         DO JL = 1 , KLON
C
            ZCNSOL=PDWFSU(JL,1)*ZTT(JL,1)          *ZTT(JL,10)
     &            +PDWFSU(JL,2)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     &            +PDWFSU(JL,3)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     &            +PDWFSU(JL,4)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     &            +PDWFSU(JL,5)*ZTT(JL,3)          *ZTT(JL,14)
     &            +PDWFSU(JL,6)*ZTT(JL,6)          *ZTT(JL,15)
            ZCNSOL=(1.-ZEMISS)*ZCNSOL
C
            PFLUC(JL,1,JK)=ZCNSOL+PBINT(JL,JK)-PDISU(JL,JK)-PADJU(JL,JK)
         ENDDO
C
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END SUBROUTINE LWVB
