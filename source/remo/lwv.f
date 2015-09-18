C
C
C**** *LWV*   - LONGWAVE RADIATION, VERTICAL INTEGRATION
C
C     PURPOSE.
C     --------
C           CARRIES OUT THE VERTICAL INTEGRATION TO GIVE LONGWAVE
C           FLUXES OR RADIANCES

C
C**   INTERFACE.
C     ----------
C     *CALL*     LWV ( KLON,KLEV,KMODE,KFLUX,KRAD,KUAER,KTRAER
C    I  , KAER,KCFC
C    S  , KXDIA,KXT,KXTSU,KXTTP
C    S  , PABCU,PBINT,PBSUIN,PBSUR,PBTOP,PDBSL,PEMIS
C    R  , PGA,PGB,PGASUR,PGBSUR,PGATOP,PGBTOP
C    R  , PGC,PGD,PGCSUR,PGDSUR,PGCTOP,PGDTOP
C    S  , PCNTRB,PCTS,PHFG,PEMD,PEMU,PCOLC,PFLUC                 )
C

C        EXPLICIT ARGUMENTS :
C        --------------------
C     ==== INPUTS ===
C KX...    (KLON,...          ; TEMPERATURE INDICES
C PABCU : (KLON,NUA,3*KLEV+1); ABSORBER AMOUNTS
C PBINT  : (KLON,0:KLEV)     ; HALF-LEVEL PLANCK FUNCTION
C PBTOP  : (KLON,NINT)        ; T.O.A. SPECTRAL PLANCK FUNCTION
C PDBSL  : (KLON,KLEV*2)     ; SUB-LAYER PLANCK FUNCTION GRADIENT
C     ==== OUTPUTS ===
C PCNTRB : (KLON,KLEV+1,KLEV+1); CLEAR-SKY ENERGY EXCHANGE MATRIX
C  IF KMODE = 0, 1, 2
C PFLUC(KLON,2,KLEV)         ; RADIATIVE FLUXES CLEAR-SKY:
C                     1  ==>  UPWARD   FLUX TOTAL
C
C        IMPLICIT ARGUMENTS :   NONE
C        --------------------
C
C     METHOD.
C     -------
C
C          1. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING BETWEEN
C     CONTRIBUTIONS BY -  THE NEARBY LAYERS
C                      -  THE DISTANT LAYERS
C                      -  THE BOUNDARY TERMS
C          2. COMPUTES THE CLEAR-SKY DOWNWARD AND UPWARD EMISSIVITIES.
C
C     EXTERNALS.
C     ----------
C
C          *LWVN*, *LWVD*, *LWVB*
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
      SUBROUTINE LWV ( KLON,KLEV,KUAER,KTRAER
     &  , KAER,KCFC
     &  , PABCU,PBINT,PBTOP,PDBSL
     &  , PGA,PGB,PGASUR,PGBSUR,PGATOP,PGBTOP
     &  , PGC,PGD,PGCSUR,PGDSUR,PGCTOP,PGDTOP
     &  , PCNTRB,PFLUC                 )
C
      IMPLICIT NONE
C
      INCLUDE "YOMLW"
C-----------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
      INTEGER, INTENT(IN)    :: KLON,KLEV,KUAER,KTRAER
     &                        , KAER,KCFC
      REAL,    INTENT(IN)    :: PABCU(KLON,NUA,3*KLEV+1)
     &                       ,  PBINT(KLON,KLEV+1)
     &                       ,  PBTOP(KLON,NINT)
     &                       ,  PDBSL(KLON,NINT,KLEV*2)
     &                       ,  PGA(KLON,8,2,KLEV), PGB(KLON,8,2,KLEV)
     &                       ,  PGASUR(KLON,8,2), PGBSUR(KLON,8,2)
     &                       ,  PGATOP(KLON,8,2), PGBTOP(KLON,8,2)
     &                       ,  PGC(KLON,5,2,KLEV),PGD(KLON,5,2,KLEV)
     &                       ,  PGCSUR(KLON,5,2),PGDSUR(KLON,5,2)
     &                       ,  PGCTOP(KLON,5,2),PGDTOP(KLON,5,2)
C
      REAL,    INTENT(INOUT) :: PCNTRB(KLON,KLEV+1,KLEV+1)
     &                       ,  PFLUC(KLON,2,KLEV+1)
C
C-----------------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS
C              ------------
C
C
      REAL ZADJD(KLON,KLEV+1), ZADJU(KLON,KLEV+1)
     &  ,  ZDBDT(KLON,NINT,KLEV)
     &  ,  ZDISD(KLON,KLEV+1), ZDISU(KLON,KLEV+1)
     &  ,  ZDWFSU(KLON,NINT)
C
!DIR$ NOBOUNDS
C-----------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------
C
C*         1.1     INITIALIZE LAYER CONTRIBUTIONS
C                  ------------------------------
C
      ZADJD(:,:)=0.
      ZADJU(:,:)=0.
      ZDISD(:,:)=0.
      ZDISU(:,:)=0.
C
C
      ZDWFSU(:,:)=0.
C
C     ------------------------------------------------------------------
C
C*         2.      VERTICAL INTEGRATION
C                  --------------------
C
C     ------------------------------------------------------------------
C
C*         2.1     CONTRIBUTION FROM ADJACENT LAYERS
C                  ---------------------------------
C
      CALL LWVN ( KLON,KLEV,KUAER
     &  , KAER,KCFC
     &  , PABCU,PDBSL,PGA,PGB,PGC,PGD
     &  , ZADJD,ZADJU, PCNTRB,ZDBDT,ZDWFSU )
C
C     ------------------------------------------------------------------
C
C*         2.2     CONTRIBUTION FROM DISTANT LAYERS
C                  ---------------------------------
C
      CALL LWVD ( KLON,KLEV,KUAER,KTRAER
     &  , KAER,KCFC
     &  , PABCU,ZDBDT,PGA,PGB,PGC,PGD
     &  , PCNTRB,ZDISD,ZDISU,ZDWFSU )
C
C     ------------------------------------------------------------------
C
C*         2.3     EXCHANGE WITH THE BOUNDARIES
C                  ----------------------------
C
      CALL LWVB ( KLON,KLEV,KUAER,KAER,KCFC
     &  , PABCU,ZADJD,ZADJU,PBINT,PBTOP
     &  , ZDISD,ZDISU
     &  , PGASUR,PGBSUR,PGATOP,PGBTOP
     &  , PGCSUR,PGDSUR,PGCTOP,PGDTOP
     &  , PFLUC,ZDWFSU )
C
      RETURN
      END SUBROUTINE LWV
