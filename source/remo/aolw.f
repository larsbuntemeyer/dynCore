      SUBROUTINE AOLW ( KLON, KLEV, KEWAER, KAERH, KAER, KCFC
     &  , PCCO2,PCLDLW,PDP,PPMB,PQOF,PTL
     &  , PAER,PTAVE,PWV
     &  , PFLUX,PFLUC,PCFCABS )
      !
      IMPLICIT NONE
      !
C
C**** *LW*   - ORGANIZES THE LONGWAVE CALCULATIONS
C
C     PURPOSE.
C     --------
C           DEPENDING ON KMODE, COMPUTES LONGWAVE FLUXES AND/OR
C           RADIANCES
C
C**   INTERFACE.
C     ----------
C     *CALL*     LW ( KLON, KLEV, KAER, KCFC, KMODE, KFLUX, KRAD
C    &  , PCCO2,PCLDLW,PDP,PDT0,PEMIS,PPMB,PPSOL,PQOF,PTL
C    &  , PAER,PTAVE,PVIEW,PWV
C    &  , PEMD,PEMU,PCOLR,PCOLC,PFLUX,PFLUC             )
C

C        *LW* IS CALLED FROM *RADLSW*
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C PCCO2  :                     ; CONCENTRATION IN CO2 (PA/PA)
C PQOF   : (KLON,KLEV)       ; CONCENTRATION IN OZONE (PA/PA)
C PTAVE  : (KLON,KLEV)       ; TEMPERATURE
C PTL    : (KLON,KLEV+1)     ; HALF LEVEL TEMPERATURE
C PPMB   : (KLON,KLEV+1)     ; HALF LEVEL PRESSURE
C PWV    : (KLON,KLEV)       ; SPECIFIC HUMIDITY PA/PA
C PCLDLW : (KLON,KLEV)       ; CLOUD FRACTIONAL COVER
C PAER   : (KDLON,KFLEV,5+KEWAER); AEROSOL OPTICAL THICKNESS (1,..,5)
C                                  TANRE ET AL., 1984
C                                  AEROSOL MASS MIXING RATIO (KG/KG)
C                                  (6,...,5+KEWAER) COMPUTED IN ECHAM4
C     ==== OUTPUTS ===
C  IF KMODE = 0, 1, 2
C PFLUX(KLON,2,KLEV)         ; RADIATIVE FLUXES :
C                     1  ==>  UPWARD   FLUX TOTAL
C                     2  ==>  DOWNWARD FLUX TOTAL
C PFLUC(KLON,2,KLEV)         ; RADIATIVE FLUXES CLEAR SKY:
C                     1  ==>  UPWARD   FLUX TOTAL
C                     2  ==>  DOWNWARD FLUX TOTAL
C
C        IMPLICIT ARGUMENTS :   NONE
C        --------------------
C
C     METHOD.
C     -------
C
C          1. COMPUTES THE PRESSURE AND TEMPERATURE WEIGHTED AMOUNTS OF
C     ABSORBERS.
C          2. COMPUTES THE PLANCK FUNCTIONS ON THE INTERFACES AND THE
C     GRADIENT OF PLANCK FUNCTIONS IN THE LAYERS.
C          3. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING THE CON-
C     TRIBUTIONS OF THE ADJACENT AND DISTANT LAYERS AND THOSE FROM THE
C     BOUNDARIES.
C          4. COMPUTES THE CLEAR-SKY DOWNWARD AND UPWARD EMISSIVITIES.
C          5. INTRODUCES THE EFFECTS OF THE CLOUDS ON THE FLUXES.
C
C     EXTERNALS.
C     ----------
C
C          *LWU*, *LWB*, *LWV*, *LWC*
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
C       MODIFIED: ROB VAN DORLAND, KNMI, 95-05-10
C-----------------------------------------------------------------------
      !IMPLICIT LOGICAL (L)
C
      INCLUDE "YOMLW"
      INCLUDE "YOMRDU"
C-----------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
      !
      ! Dummy Arguments
      !
      INTEGER, INTENT(IN)  :: KLON,KLEV
      INTEGER, INTENT(IN)  :: KEWAER, KAER, KCFC
      REAL,    INTENT(IN)  :: PCCO2
      INTEGER, INTENT(IN)  :: KAERH(KLON,KLEV)
      !
      REAL,    INTENT(IN)  :: PCLDLW(KLON,KLEV),PDP(KLON,KLEV)
     &                     ,  PPMB(KLON,KLEV+1)
     &                     ,  PQOF(KLON,KLEV), PTL(KLON,KLEV+1)
     &                     ,  PAER(KLON,KLEV,5+KEWAER), PTAVE(KLON,KLEV)
     &                     ,  PWV(KLON,KLEV)
      REAL,    INTENT(IN)  :: PCFCABS(4)
      !
      REAL,    INTENT(OUT) :: PFLUX(KLON,2,KLEV+1)
     &                     ,  PFLUC(KLON,2,KLEV+1)
C
C
      

C-----------------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS
C              ------------
C
      REAL ZABCU(KLON,NUA,3*KLEV+1),ZB(KLON,NINT,KLEV+1)
     &  ,  ZBINT(KLON,KLEV+1)
     &  ,  ZBSUI(KLON), ZBTOP(KLON,NINT)
     &  ,  ZCNTRB(KLON,KLEV+1,KLEV+1)
     &  ,  ZDBSL(KLON,NINT,KLEV*2)
     &  ,  ZGA(KLON,8,2,KLEV), ZGB(KLON,8,2,KLEV)
     &  ,  ZGASUR(KLON,8,2)   , ZGBSUR(KLON,8,2)
     &  ,  ZGATOP(KLON,8,2)   , ZGBTOP(KLON,8,2)
     &  ,   ZGC(KLON,5,2,KLEV), ZGD(KLON,5,2,KLEV)
     &  ,   ZGCSUR(KLON,5,2),ZGDSUR(KLON,5,2)
     &  ,   ZGCTOP(KLON,5,2),ZGDTOP(KLON,5,2)
C
C     ------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------
C
C     ------------------------------------------------------------------
C
C*         1.1   COMPUTES ABSORBER AMOUNTS
C                -------------------------
C
      CALL AOLWU ( KLON,KLEV, KEWAER, KAERH
     &  ,  PAER,PCCO2,PDP,PPMB,PQOF,PTAVE,PWV
     &  ,  ZABCU,PCFCABS   )
C
C     ------------------------------------------------------------------
C
C*         2.    COMPUTES PLANCK FUNCTIONS
C                -------------------------
C
       CALL LWB ( KLON,KLEV
     &          , PTAVE,PTL
     &          , ZB,ZBINT,ZBSUI,ZBTOP,ZDBSL
     &          , ZGA,ZGB,ZGASUR,ZGBSUR,ZGATOP,ZGBTOP
     &          , ZGC,ZGD,ZGCSUR,ZGDSUR,ZGCTOP,ZGDTOP  )
C
C     ------------------------------------------------------------------
C
C*         3.    PERFORMS THE VERTICAL INTEGRATION
C                ---------------------------------
C
      CALL LWV ( KLON,KLEV,NUAER,NTRAER
     &  , KAER,KCFC
     &  , ZABCU,ZBINT,ZBTOP,ZDBSL
     &  , ZGA,ZGB,ZGASUR,ZGBSUR,ZGATOP,ZGBTOP
     &  , ZGC,ZGD,ZGCSUR,ZGDSUR,ZGCTOP,ZGDTOP
     &  , ZCNTRB,PFLUC                     )
C
C
C     ------------------------------------------------------------------
C
C*         4.    INTRODUCES THE EFFECTS OF CLOUDS
C                --------------------------------
C
       CALL LWC ( KLON,KLEV
     &  , ZBINT,ZBSUI,PCLDLW,ZCNTRB,PFLUC
     &  ,  PFLUX                                                )
C
      RETURN
      END SUBROUTINE AOLW
