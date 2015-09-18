C
C     SUBROUTINE INBOXS
C
C**** INBOXS   -   UP: HORIZONTAL TRANSPORTS INTO GRID BOXES
C
C**   AUFRUF   :   CALL INBOXS IN HP *INBOXS*
C**   ENTRIES  :      ---
C**   ZWECK    :
C     COMPUTATION OF THE HORIZONTAL TRANSPORTS OF QD, QW, KINETIC ENERGY,
C     SENSIBLE HEAT, AND POTENTIAL ENERGY INTO EACH GRID BOX DURING
C     THE TIME BETWEEN TO CALLS OF DIAGNOS.
C     OUTPUT VARIABLES CONTAIN CHANGES OF WATER CONTENTS AND ENERGY, ALREADY
C     INTEGRATED OVER EACH GRID BOX.
C     ALL WATER VARIABLES Q*BOXS ARE ACCUMULATIONS. THEY ARE TIME INTEGRATED
C     USING DT.
C     ALL ENERGY VARIABLES ARE ACCUMULATIONS, BUT DIVIDED BY 3600. THEY ARE
C     TIME INTEGRATED USING DTDEH. TO OBTAIN THE TIME AVERAGE, ONLY DIVISION
C     BY THE NUMBER OF HOURS IS REQUIRED, AS IN OTHER TIME AVERAGE
C     VARIABLES OF THE MODEL.

C     ATTENTION: BOXES AT THE BOUNDARIES OF THE COMPUTATIONAL DOMAIN
C                ARE NOT CONSIDERED!
C

C**   VERSIONS-
C**   DATUM    :   20.02.97
C**
C     INPUT PARAMETERS
C     ----------------
C
C     U,V      - WIND COMPONENTS   [M/S]
C     T        - TEMPERATURE       [K]
C     QD       - SPECIFIC HUMIDITY              [KG/KG]
C     QW       - SPECIFIC CLOUD WATER CONTENT   [KG/KG]
C     QI       - SPECIFIC CLOUD ICE CONTENT     [KG/KG]
C     FI       - GEOPOTENTIAL                   [
C     FIB      - GEOPOTENTIAL AT THE SURFACE    [
C     PS       - SURFACE PRESSSURE  [PA]
C     DAK,DBK  - COEFFICIENTS FOR ETA VERTICAL COORDINATE
C     CPHI     - COS(PHI)
C     IE,JE,KE  - NUMBERS OF GRIDPOINTS (ZONAL, MERIDIONAL, VERTICAL)
C     DLAM   - MESH WIDTH IN LAMBDA DIRECTION  [GRAD]
C     DPHI   - MESH WIDTH IN PHI DIRECTION     [GRAD]
C     DTDEH    - TIME STEP DIVIDED BY 3600       [S]
C     G,RERD,WCP - GRAV. ACC., MEAN EARTH RADIUS, HEAT CAPACITY OF DRY AIR
C
C
C     OUTPUT PARAMETERS
C     -----------------
C
C     QDBOXS - TIME INTEGRAL OF HORIZONTAL TRANSPORT OF WATER VAPOUR  [KG/M**2]
C     QWBOXS - TIME INTEGRAL OF HORIZONTAL TRANSPORT OF CLOUD WATER   [KG/M**2]
C     QIBOXS - TIME INTEGRAL OF HORIZONTAL TRANSPORT OF CLOUD ICE   [KG/M**2]
C     EKBOXS - TIME INTEGRAL OF HORIZONTAL TRANSPORT OF KINETIC ENERGY[J/M**2]
C              DIVIDED BY 3600.
C     FHBOXS - TIME INTEGRAL OF HORIZONTAL TRANSPORT OF SENSIBLE HEAT [J/M**2]
C              DIVIDED BY 3600.
C     FIBOXS - TIME INTEGRAL OF HORIZONTAL TRANSPORT OF POT.   ENERGY [J/M**2]
C              DIVIDED BY 3600.
C
C**   COMMON-
C**   BLOECKE  :   PHYKON, HIGKON, ORG, GRID
C**
C**   METHODE  :
C**   FEHLERBE-
C**   HANDLUNG :      ---
C**   VERFASSER:   RALPH NOLTE-HOLBE
C**                WORKSTATION VERSION: BURKHARDT ROCKEL
C
      SUBROUTINE INBOXS(
     &                  U     , V     , T     , QD    , QW  , PS    ,
     &                  FI    , FIB   , DAK   , DBK   , CPHI, QDBOXS,
     &                  QWBOXS, EKBOXS, FHBOXS, FIBOXS, QI  , QIBOXS)
C
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
      INCLUDE "phykon.h"
      INCLUDE "higkon.h"
      INCLUDE "grid.h"
      INCLUDE "org.h"
C
C     Dummy Arguments
C
      REAL, INTENT(IN) ::
     &           T      (IE,JE,KE,3)       ,
     &           QD     (IE,JE,KE,3)       ,
     &           QW     (IE,JE,KE,3)       ,
     &           U      (IE,JE,KE,3)       ,
     &           V      (IE,JE,KE,3)       ,
     &           FI     (IE,JE,KE,2)       ,
     &           PS     (IE,JE,3)          ,
     &           FIB    (IE,JE)            ,
     &           DAK    (KE)               ,
     &           DBK    (KE)               ,
     &           CPHI   (JE,2)             ,
     &           QI     (IE,JE,KE,3)       
      REAL, INTENT(INOUT) ::
     &           QDBOXS (IE,JE)            ,
     &           QWBOXS (IE,JE)            ,
     &           EKBOXS (IE,JE)            ,
     &           FHBOXS (IE,JE)            ,
     &           FIBOXS (IE,JE)            ,
     &           QIBOXS (IE,JE)
C
C     Local Variables
C
      REAL :: ZFHW,ZFHS,ZFHN,ZFHE,ZFH0,ZDT,ZDQW,ZDQI,ZDQD,ZDPIPJ,ZDPIMJ,
     &        ZDPIJP,ZDPIJM,ZDPIJ,ZDPHI,ZDLAM,ZDFI,ZDFH,ZDEK,ZABOX,
     &        ZFI0, ZFI0P, ZFIE, ZFIN, ZFIS, ZFIW, ZPI, ZQD0, ZQDE,  
     &        ZQDS, ZQDW, ZQI0, ZQIE, ZQIN, ZQIS, ZQIW, ZQW0, ZQWE, 
     &        ZQWS, ZQWW, ZU2IJ, ZU2IJM, ZU2IJP, ZU2IMJ, ZU2IMJM,  
     &        ZUV2E, ZUV2N, ZUV2S, ZUV2W, ZV2IJ, ZV2IJM, ZV2IMJ,  
     &        ZV2IPJ, ZV2IPJM, ZQDN, ZU2IMJP, ZV2IMJM, ZQWN
C
      INTEGER :: NX2,NX,K,J,I
C
      NX = NJ
      NX2= NJ2
      ZPI=4.*ATAN(1.)
      ZDPHI = DPHI/180.*ZPI
      ZDLAM = DLAM/180.*ZPI
C *** THIS K-LOOP IS SPLITTED BECAUSE OF    '' FI(KE+1) = FIB '':
      DO K = 1,KE-1
         DO J = 2,JE-1
C **        AREA OF ONE GRID BOX:
            ZABOX = (RERD**2)*CPHI(J,1)*ZDPHI*ZDLAM
            DO I = 2,IE-1
C**            DPS IN BOX (I,J,K) AND ITS NEIGHBOURS:
               ZDPIJ  = DAK(K) + DBK(K)*PS(I,  J  ,NX)
               ZDPIPJ = DAK(K) + DBK(K)*PS(I+1,J  ,NX)
               ZDPIMJ = DAK(K) + DBK(K)*PS(I-1,J  ,NX)
               ZDPIJP = DAK(K) + DBK(K)*PS(I  ,J+1,NX)
               ZDPIJM = DAK(K) + DBK(K)*PS(I  ,J-1,NX)
C**            INTERPOLATION TO BOX BOUNDARY VALUES: EAST, WEST, NORTH, SOUTH
               ZQD0 =      QD(I  ,J,K,NX)
               ZQDE = 0.5*(QD(I+1,J,K,NX)+ZQD0)
               ZQDW = 0.5*(QD(I-1,J,K,NX)+ZQD0)
               ZQDN = 0.5*(QD(I,J+1,K,NX)+ZQD0)
               ZQDS = 0.5*(QD(I,J-1,K,NX)+ZQD0)

               ZQW0 =      QW(I  ,J,K,NX)
               ZQWE = 0.5*(QW(I+1,J,K,NX)+ZQW0)
               ZQWW = 0.5*(QW(I-1,J,K,NX)+ZQW0)
               ZQWN = 0.5*(QW(I,J+1,K,NX)+ZQW0)
               ZQWS = 0.5*(QW(I,J-1,K,NX)+ZQW0)

               ZFH0 = T(I,J,K,NX)
               ZFHE = 0.5*( T(I+1,J,K,NX)+ZFH0)
               ZFHW = 0.5*( T(I-1,J,K,NX)+ZFH0)
               ZFHN = 0.5*( T(I,J+1,K,NX)+ZFH0)
               ZFHS = 0.5*( T(I,J-1,K,NX)+ZFH0)

               ZQI0 =      QI(I  ,J,K,NX)
               ZQIE = 0.5*(QI(I+1,J,K,NX)+ZQI0)
               ZQIW = 0.5*(QI(I-1,J,K,NX)+ZQI0)
               ZQIN = 0.5*(QI(I,J+1,K,NX)+ZQI0)
               ZQIS = 0.5*(QI(I,J-1,K,NX)+ZQI0)

C ***          GEOPOTENTIAL AT BOX BOUNDARIES:
               ZFI0 = FI(I,J ,K,  NX2)
               ZFI0P= FI(I,J ,K+1,NX2)
               ZFIE = 0.25*( ZFI0 +FI(I+1,J  ,K  ,NX2  )
     &                      +ZFI0P+FI(I+1,J  ,K+1,NX2) )
               ZFIN = 0.25*( ZFI0 +FI(I  ,J+1,K  ,NX2  )
     &                      +ZFI0P+FI(I  ,J+1,K+1,NX2) )
               ZFIW = 0.25*( ZFI0 +FI(I-1,J  ,K  ,NX2  )
     &                      +ZFI0P+FI(I-1,J  ,K+1,NX2) )
               ZFIS = 0.25*( ZFI0 +FI(I  ,J-1,K  ,NX2  )
     &                      +ZFI0P+FI(I  ,J-1,K+1,NX2) )


C**       -------------  HORIZONTAL FLUXES QD, QW, WCP*RHO*T  ---------
            ZDQD=(-U(I  ,J  ,K,NX)*ZQDE*            ZDPHI*(ZDPIJ+ZDPIPJ)
     &            -V(I  ,J  ,K,NX)*ZQDN*CPHI(J  ,2)*ZDLAM*(ZDPIJ+ZDPIJP)
     &            +U(I-1,J  ,K,NX)*ZQDW*            ZDPHI*(ZDPIJ+ZDPIMJ)
     &            +V(I  ,J-1,K,NX)*ZQDS*CPHI(J-1,2)*ZDLAM*(ZDPIJ+ZDPIJM)
     &           )*0.5*RERD/G

            ZDQW=(-U(I  ,J  ,K,NX)*ZQWE*            ZDPHI*(ZDPIJ+ZDPIPJ)
     &            -V(I  ,J  ,K,NX)*ZQWN*CPHI(J  ,2)*ZDLAM*(ZDPIJ+ZDPIJP)
     &            +U(I-1,J  ,K,NX)*ZQWW*            ZDPHI*(ZDPIJ+ZDPIMJ)
     &            +V(I  ,J-1,K,NX)*ZQWS*CPHI(J-1,2)*ZDLAM*(ZDPIJ+ZDPIJM)
     &           )*0.5*RERD/G

            ZDFH=(-U(I  ,J  ,K,NX)*ZFHE*            ZDPHI*(ZDPIJ+ZDPIPJ)
     &            -V(I  ,J  ,K,NX)*ZFHN*CPHI(J  ,2)*ZDLAM*(ZDPIJ+ZDPIJP)
     &            +U(I-1,J  ,K,NX)*ZFHW*            ZDPHI*(ZDPIJ+ZDPIMJ)
     &            +V(I  ,J-1,K,NX)*ZFHS*CPHI(J-1,2)*ZDLAM*(ZDPIJ+ZDPIJM)
     &           )*0.5*RERD/G*WCP
            ZDFI=(-U(I  ,J  ,K,NX)*ZFIE*            ZDPHI*(ZDPIJ+ZDPIPJ)
     &            -V(I  ,J  ,K,NX)*ZFIN*CPHI(J  ,2)*ZDLAM*(ZDPIJ+ZDPIJP)
     &            +U(I-1,J  ,K,NX)*ZFIW*            ZDPHI*(ZDPIJ+ZDPIMJ)
     &            +V(I  ,J-1,K,NX)*ZFIS*CPHI(J-1,2)*ZDLAM*(ZDPIJ+ZDPIJM)
     &           )*0.5*RERD/G
            ZDQI=(-U(I  ,J  ,K,NX)*ZQIE*            ZDPHI*(ZDPIJ+ZDPIPJ)
     &            -V(I  ,J  ,K,NX)*ZQIN*CPHI(J  ,2)*ZDLAM*(ZDPIJ+ZDPIJP)
     &            +U(I-1,J  ,K,NX)*ZQIW*            ZDPHI*(ZDPIJ+ZDPIMJ)
     &            +V(I  ,J-1,K,NX)*ZQIS*CPHI(J-1,2)*ZDLAM*(ZDPIJ+ZDPIJM)
     &           )*0.5*RERD/G
C**       ----------  KINETIC ENERGY FLUXES  -------------------------
C**            VELOCITY SQUARED ON BOX BOUNDARIES:
               ZU2IJ   = U(I  ,J  ,K,NX)**2
               ZU2IMJ  = U(I-1,J  ,K,NX)**2
               ZU2IJP  = U(I  ,J+1,K,NX)**2
               ZU2IJM  = U(I  ,J-1,K,NX)**2
               ZU2IMJP = U(I-1,J+1,K,NX)**2
               ZU2IMJM = U(I-1,J-1,K,NX)**2
               ZV2IJ   = V(I  ,J  ,K,NX)**2
               ZV2IPJ  = V(I+1,J  ,K,NX)**2
               ZV2IMJ  = V(I-1,J  ,K,NX)**2
               ZV2IJM  = V(I  ,J-1,K,NX)**2
               ZV2IPJM = V(I+1,J-1,K,NX)**2
               ZV2IMJM = V(I-1,J-1,K,NX)**2
C**            EAST BOUNDARY:
               ZUV2E = ZU2IJ + (CPHI(J,  2)*(ZV2IJ+ZV2IPJ) +
     &                      CPHI(J-1,2)*(ZV2IJM+ZV2IPJM))/(4*CPHI(J,1))
C**            NORTH BOUNDARY:
               ZUV2N = ZV2IJ + (ZU2IMJ + ZU2IJ + ZU2IMJP + ZU2IJP)*0.25
C**            WEST BOUNDARY:
               ZUV2W = ZU2IMJ+ (CPHI(J,2)*(ZV2IMJ+ZV2IJ) +
     &                      CPHI(J-1,2)*(ZV2IMJM+ZV2IJM))/(4*CPHI(J,1))
C**            SOUTH BOUNDARY:
               ZUV2S = ZV2IJM+ (ZU2IMJ + ZU2IJ + ZU2IMJM + ZU2IJM)*0.25
C**            ADD 1/2  RHO ZU2 * (NORMAL VELOCITY ON THE BOUNDARY)
C**            FOR EACH BOUNDARY:
          ZDEK=(-ZUV2E*U(I  ,J  ,K,NX)*ZDPHI            *(ZDPIJ+ZDPIPJ)
     &          -ZUV2N*V(I  ,J  ,K,NX)*CPHI(J  ,2)*ZDLAM*(ZDPIJ+ZDPIJP)
     &          +ZUV2W*U(I-1,J  ,K,NX)*ZDPHI            *(ZDPIJ+ZDPIMJ)
     &          +ZUV2S*V(I  ,J-1,K,NX)*CPHI(J-1,2)*ZDLAM*(ZDPIJ+ZDPIJM)
     &         )*0.25*RERD/G
C***        ---------- END OF KINETIC ENERGY COMPUTATION  ---------------
C***           ADD THE CHANGES TO THE OUTPUT VARIABLES AND NORMALIZE TO
C***           VALUES PER AREA:
               ZDT=3600.*DTDEH
               QDBOXS(I,J) = QDBOXS(I,J) + ZDT*ZDQD/ZABOX
               QWBOXS(I,J) = QWBOXS(I,J) + ZDT*ZDQW/ZABOX
               EKBOXS(I,J) = EKBOXS(I,J) + DTDEH*ZDEK/ZABOX
               FHBOXS(I,J) = FHBOXS(I,J) + DTDEH*ZDFH/ZABOX
               FIBOXS(I,J) = FIBOXS(I,J) + DTDEH*ZDFI/ZABOX
               QIBOXS(I,J) = QIBOXS(I,J) + ZDT*ZDQI/ZABOX
            ENDDO
         ENDDO
      ENDDO ! DO K = 1,KE-1



C***  AND NOW THE ABOVE FOR K=KE:

      K = KE
      DO J = 2,JE-1
C **     AREA OF ONE GRID BOX:
         ZABOX = (RERD**2)*CPHI(J,1)*ZDPHI*ZDLAM
         DO I = 2,IE-1
C**         DPS IN BOX (I,J,K) AND ITS NEIGHBOURS:
            ZDPIJ  = DAK(K) + DBK(K)*PS(I,  J  ,NX)
            ZDPIPJ = DAK(K) + DBK(K)*PS(I+1,J  ,NX)
            ZDPIMJ = DAK(K) + DBK(K)*PS(I-1,J  ,NX)
            ZDPIJP = DAK(K) + DBK(K)*PS(I  ,J+1,NX)
            ZDPIJM = DAK(K) + DBK(K)*PS(I  ,J-1,NX)
C**           INTERPOLATION TO BOX BOUNDARY VALUES: EAST, WEST, NORTH, SOUTH
            ZQD0 =      QD(I  ,J,K,NX)
            ZQDE = 0.5*(QD(I+1,J,K,NX)+ZQD0)
            ZQDW = 0.5*(QD(I-1,J,K,NX)+ZQD0)
            ZQDN = 0.5*(QD(I,J+1,K,NX)+ZQD0)
            ZQDS = 0.5*(QD(I,J-1,K,NX)+ZQD0)

            ZQW0 =      QW(I  ,J,K,NX)
            ZQWE = 0.5*(QW(I+1,J,K,NX)+ZQW0)
            ZQWW = 0.5*(QW(I-1,J,K,NX)+ZQW0)
            ZQWN = 0.5*(QW(I,J+1,K,NX)+ZQW0)
            ZQWS = 0.5*(QW(I,J-1,K,NX)+ZQW0)

            ZFH0 = T(I,J,K,NX)
            ZFHE = 0.5*( T(I+1,J,K,NX)+ZFH0)
            ZFHW = 0.5*( T(I-1,J,K,NX)+ZFH0)
            ZFHN = 0.5*( T(I,J+1,K,NX)+ZFH0)
            ZFHS = 0.5*( T(I,J-1,K,NX)+ZFH0)

            ZQI0 =      QI(I  ,J,K,NX)
            ZQIE = 0.5*(QI(I+1,J,K,NX)+ZQI0)
            ZQIW = 0.5*(QI(I-1,J,K,NX)+ZQI0)
            ZQIN = 0.5*(QI(I,J+1,K,NX)+ZQI0)
            ZQIS = 0.5*(QI(I,J-1,K,NX)+ZQI0)

C ***       GEOPOTENTIAL AT BOX BOUNDARIES:
            ZFI0 = FI (I,J,K,NX2)
            ZFI0P= FIB(I,J      )
            ZFIE=0.25*( ZFI0  + FI (I+1,J  ,K,NX2)
     &                 +ZFI0P + FIB(I+1,J        ) )
            ZFIN=0.25*( ZFI0  + FI (I  ,J+1,K,NX2)
     &                 +ZFI0P + FIB(I  ,J+1      ) )
            ZFIW=0.25*( ZFI0  + FI (I-1,J  ,K,NX2)
     &                 +ZFI0P + FIB(I-1,J        ) )
            ZFIS=0.25*( ZFI0  + FI (I  ,J-1,K,NX2)
     &                 +ZFI0P + FIB(I  ,J-1      ) )

C**       -------------  HORIZONTAL FLUXES QD, QW, WCP*RHO*T  --------
            ZDQD=(-U(I  ,J  ,K,NX)*ZQDE*            ZDPHI*(ZDPIJ+ZDPIPJ)
     &            -V(I  ,J  ,K,NX)*ZQDN*CPHI(J  ,2)*ZDLAM*(ZDPIJ+ZDPIJP)
     &            +U(I-1,J  ,K,NX)*ZQDW*            ZDPHI*(ZDPIJ+ZDPIMJ)
     &            +V(I  ,J-1,K,NX)*ZQDS*CPHI(J-1,2)*ZDLAM*(ZDPIJ+ZDPIJM)
     &           )*0.5*RERD/G

            ZDQW=(-U(I  ,J  ,K,NX)*ZQWE*            ZDPHI*(ZDPIJ+ZDPIPJ)
     &            -V(I  ,J  ,K,NX)*ZQWN*CPHI(J  ,2)*ZDLAM*(ZDPIJ+ZDPIJP)
     &            +U(I-1,J  ,K,NX)*ZQWW*            ZDPHI*(ZDPIJ+ZDPIMJ)
     &            +V(I  ,J-1,K,NX)*ZQWS*CPHI(J-1,2)*ZDLAM*(ZDPIJ+ZDPIJM)
     &           )*0.5*RERD/G

            ZDFH=(-U(I  ,J  ,K,NX)*ZFHE*            ZDPHI*(ZDPIJ+ZDPIPJ)
     &            -V(I  ,J  ,K,NX)*ZFHN*CPHI(J  ,2)*ZDLAM*(ZDPIJ+ZDPIJP)
     &            +U(I-1,J  ,K,NX)*ZFHW*            ZDPHI*(ZDPIJ+ZDPIMJ)
     &            +V(I  ,J-1,K,NX)*ZFHS*CPHI(J-1,2)*ZDLAM*(ZDPIJ+ZDPIJM)
     &           )*0.5*RERD/G*WCP
            ZDFI=(-U(I  ,J  ,K,NX)*ZFIE*            ZDPHI*(ZDPIJ+ZDPIPJ)
     &            -V(I  ,J  ,K,NX)*ZFIN*CPHI(J  ,2)*ZDLAM*(ZDPIJ+ZDPIJP)
     &            +U(I-1,J  ,K,NX)*ZFIW*            ZDPHI*(ZDPIJ+ZDPIMJ)
     &            +V(I  ,J-1,K,NX)*ZFIS*CPHI(J-1,2)*ZDLAM*(ZDPIJ+ZDPIJM)
     &           )*0.5*RERD/G
            ZDQI=(-U(I  ,J  ,K,NX)*ZQIE*            ZDPHI*(ZDPIJ+ZDPIPJ)
     &            -V(I  ,J  ,K,NX)*ZQIN*CPHI(J  ,2)*ZDLAM*(ZDPIJ+ZDPIJP)
     &            +U(I-1,J  ,K,NX)*ZQIW*            ZDPHI*(ZDPIJ+ZDPIMJ)
     &            +V(I  ,J-1,K,NX)*ZQIS*CPHI(J-1,2)*ZDLAM*(ZDPIJ+ZDPIJM)
     &           )*0.5*RERD/G
C**       ----------  KINETIC ENERGY FLUXES  -------------------------
C**              VELOCITY SQUARED ON BOX BOUNDARIES:
            ZU2IJ   = U(I  ,J  ,K,NX)**2
            ZU2IMJ  = U(I-1,J  ,K,NX)**2
            ZU2IJP  = U(I  ,J+1,K,NX)**2
            ZU2IJM  = U(I  ,J-1,K,NX)**2
            ZU2IMJP = U(I-1,J+1,K,NX)**2
            ZU2IMJM = U(I-1,J-1,K,NX)**2
            ZV2IJ   = V(I  ,J  ,K,NX)**2
            ZV2IPJ  = V(I+1,J  ,K,NX)**2
            ZV2IMJ  = V(I-1,J  ,K,NX)**2
            ZV2IJM  = V(I  ,J-1,K,NX)**2
            ZV2IPJM = V(I+1,J-1,K,NX)**2
            ZV2IMJM = V(I-1,J-1,K,NX)**2
C**              EAST BOUNDARY:
            ZUV2E = ZU2IJ + (CPHI(J,  2)*(ZV2IJ+ZV2IPJ) +
     &                      CPHI(J-1,2)*(ZV2IJM+ZV2IPJM))/(4*CPHI(J,1))
C**              NORTH BOUNDARY:
            ZUV2N = ZV2IJ + (ZU2IMJ + ZU2IJ + ZU2IMJP + ZU2IJP)*0.25
C**              WEST BOUNDARY:
            ZUV2W = ZU2IMJ+ (CPHI(J,2)*(ZV2IMJ+ZV2IJ) +
     &                      CPHI(J-1,2)*(ZV2IMJM+ZV2IJM))/(4*CPHI(J,1))
C**              SOUTH BOUNDARY:
            ZUV2S = ZV2IJM+ (ZU2IMJ + ZU2IJ + ZU2IMJM + ZU2IJM)*0.25
C**              ADD 1/2  RHO ZU2 * (NORMAL VELOCITY ON THE BOUNDARY)
C**              FOR EACH BOUNDARY:
          ZDEK=(-ZUV2E*U(I  ,J  ,K,NX)*ZDPHI            *(ZDPIJ+ZDPIPJ)
     &          -ZUV2N*V(I  ,J  ,K,NX)*CPHI(J  ,2)*ZDLAM*(ZDPIJ+ZDPIJP)
     &          +ZUV2W*U(I-1,J  ,K,NX)*ZDPHI            *(ZDPIJ+ZDPIMJ)
     &          +ZUV2S*V(I  ,J-1,K,NX)*CPHI(J-1,2)*ZDLAM*(ZDPIJ+ZDPIJM)
     &         )*0.25*RERD/G
C**         ---------- END OF KINETIC ENERGY COMPUTATION  ---------------
C**         ADD THE CHANGES TO THE OUTPUT VARIABLES AND NORMALIZE TO
C**         VALUES PER AREA:
            ZDT=3600.*DTDEH
            QDBOXS(I,J) = QDBOXS(I,J) + ZDT*ZDQD/ZABOX
            QWBOXS(I,J) = QWBOXS(I,J) + ZDT*ZDQW/ZABOX
            EKBOXS(I,J) = EKBOXS(I,J) + DTDEH*ZDEK/ZABOX
            FHBOXS(I,J) = FHBOXS(I,J) + DTDEH*ZDFH/ZABOX
            FIBOXS(I,J) = FIBOXS(I,J) + DTDEH*ZDFI/ZABOX
            QIBOXS(I,J) = QIBOXS(I,J) + ZDT*ZDQI/ZABOX
         ENDDO
      ENDDO
C***   END OF K=KE
      RETURN
      END SUBROUTINE INBOXS
