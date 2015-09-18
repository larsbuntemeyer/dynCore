C
C**** *LWU* - LONGWAVE EFFECTIVE ABSORBER AMOUNTS
C
C     PURPOSE.
C     --------
C           COMPUTES ABSORBER AMOUNTS INCLUDING PRESSURE AND
C           TEMPERATURE EFFECTS
C
C**   INTERFACE.
C     ----------
C     *CALL*     LWU ( KLON,KLEV,KAER,KMODE,KFLUX,KRAD
C    1  ,  PAER,PCCO2,PDP,PPMB,PPSOL,PQOF,PTAVE,PVIEW,PWV
C    2  ,  KXDIA,KXT,KXTSU,KXTTP, PABCU                      )
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C     ==== INPUTS ===
C PAER   : (KLON,KLEV,5+KEWAER); AEROSOL OPTICAL THICKNESS (1,..,5)
C                                  TANRE ET AL., 1984
C                                  AEROSOL MASS MIXING RATIO (KG/KG)
C                                  (6,...,5+KEWAER) COMPUTED IN ECHAM4
C PCCO2  :                     ; CONCENTRATION IN CO2 (PA/PA)
C PDP    : (KLON,KLEV)       ; LAYER PRESSURE THICKNESS (PA)
C PPMB   : (KLON,0:KLEV)     ; HALF LEVEL PRESSURE
C PQOF   : (KLON,KLEV)       ; CONCENTRATION IN OZONE (PA/PA)
C PTAVE  : (KLON,KLEV)       ; TEMPERATURE
C PWV    : (KDLSUR,KLEV)      ; SPECIFIC HUMIDITY PA/PA
C     ==== OUTPUTS ===
C KX...  : (KLON,...          ; TEMPERATURE INDICES
C PABCU  :(KLON,NUA,3*KLEV+1); EFFECTIVE ABSORBER AMOUNTS
C
C        IMPLICIT ARGUMENTS :   NONE
C        --------------------
C
C     METHOD.
C     -------
C
C          1. COMPUTES THE PRESSURE AND TEMPERATURE WEIGHTED AMOUNTS OF
C     ABSORBERS.
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
C         MODIFIED : 92-05-25   ROB VAN DORLAND  *KNMI*
C       MODIFIED: ROB VAN DORLAND, KNMI, 95-05-10
C-----------------------------------------------------------------------
      SUBROUTINE AOLWU ( KLON, KLEV, KEWAER, KAERH
     &  ,  PAER,PCCO2,PDP,PPMB,PQOF,PTAVE,PWV
     &  ,  PABCU,PCFCABS     )
C
      IMPLICIT NONE 
C
      INCLUDE "YOMAER"
      INCLUDE "YOMLW"
      INCLUDE "YOMRDU"
      INCLUDE "YOMCFC"
      INCLUDE "YOMRDI"
      INCLUDE "COMCON"
C-----------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
      INTEGER, INTENT(IN) :: KLON, KLEV, KEWAER
      INTEGER, INTENT(IN) :: KAERH(KLON,KLEV)
      REAL,    INTENT(IN) :: PCCO2
      REAL,    INTENT(IN) :: 
     &     PAER(KLON,KLEV,5+KEWAER), PDP(KLON,KLEV)
     &  ,  PPMB(KLON,KLEV+1), PQOF(KLON,KLEV)
     &  ,  PTAVE(KLON,KLEV),  PWV(KLON,KLEV)
C
      REAL,    INTENT(INOUT) :: PABCU(KLON,NUA,3*KLEV+1)
      REAL,    INTENT(IN)    :: PCFCABS(4)
C
C-----------------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS
C              ------------
      REAL ZABLY(KLON,NUA,3*KLEV+1)
     &  ,  ZDUC(KLON, 3*KLEV+1)
      REAL ZPHIO(KLON),ZPSC2(KLON),ZPSC3(KLON),ZPSH1(KLON)
     &  ,  ZPSH2(KLON),ZPSH3(KLON),ZPSH4(KLON),ZPSH5(KLON)
     &  ,  ZPSH6(KLON),ZPSIO(KLON),ZTCON(KLON)
     &  ,  ZPHIO2(KLON),ZPSIO2(KLON)
      REAL ZPSN2(KLON),ZPSN3(KLON),ZPSN6(KLON)
     &  ,  ZPSM3(KLON),ZPSM6(KLON)
      REAL ZSSIG(KLON,3*KLEV+1)
     &  ,  ZUAER(KLON,NINT), ZXOZ(KLON), ZXWV(KLON)
C
      INTEGER :: IAE, ICAE, IG1, IH, JAE1, JAE2, JAE3, JAER, JC, JCP1,
     &           JJ, JJPN, JK, JKI, JKIP1, JKJ, JKJP, JKJPN, JKJR
      INTEGER :: JKK, JKL, JL, JI 
      REAL :: ZALUP, ZCAC8, ZCAH1, ZCAH2, ZCAH3, ZCAH4, ZCAH5, ZCAH6, 
     &        ZCAM4, ZCAM5, ZCAN1, ZCAN2, ZCAN3, ZCBC8, ZCBH1, ZCBH2, 
     &        ZCBH3, ZCBH4, ZCBH5, ZCBH6
      REAL :: ZCBM4, ZCBM5, ZCBN1, ZCBN2, ZCBN3, ZCH4UP, ZDIFF, ZDPM, 
     &        ZFPPW, ZN2OUP, ZTAVI, ZTX, ZTX2, ZU6, ZUP, ZUPM, ZUPMCO2, 
     &        ZUPMH2O, ZUPMO3, ZZABLY
      REAL :: ZZABME, ZZABNI
C
C
!DIR$ NOBOUNDS PABCU
C-----------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------
C
      ZDIFF=DIFF
C
C-----------------------------------------------------------------------
C
C
C*         2.    PRESSURE OVER GAUSS SUB-LEVELS
C                ------------------------------
C
      DO JL = 1 , KLON
         ZSSIG(JL, 1 ) = PPMB(JL,1) * 100.
      ENDDO
C
      DO JK = 1 , KLEV
         JKJ=(JK-1)*NG1P1+1
         JKJR = JKJ
         JKJP = JKJ + NG1P1
         DO JL = 1 , KLON
            ZSSIG(JL,JKJP)=PPMB(JL,JK+1)* 100.
         ENDDO
         DO IG1=1,NG1
            JKJ=JKJ+1
            DO JL = 1 , KLON
               ZSSIG(JL,JKJ)= (ZSSIG(JL,JKJR) + ZSSIG(JL,JKJP)) * 0.5
     &           + RT1(IG1) * (ZSSIG(JL,JKJP) - ZSSIG(JL,JKJR)) * 0.5
            ENDDO
         ENDDO
      ENDDO
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C
C*         4.    PRESSURE THICKNESS AND MEAN PRESSURE OF SUB-LAYERS
C                --------------------------------------------------
C
      DO JKI=1,3*KLEV
         JKIP1=JKI+1
         DO JL = 1 , KLON
            ZABLY(JL,5,JKI)=(ZSSIG(JL,JKI)+ZSSIG(JL,JKIP1))*0.5
            ZABLY(JL,3,JKI)=(ZSSIG(JL,JKI)-ZSSIG(JL,JKIP1))
     &           /(10.*G)
         ENDDO
      ENDDO
C
      DO JK = 1 , KLEV
         JKL = KLEV+1 - JK
         DO JL = 1 , KLON
            ZXWV(JL) = AMAX1 (PWV(JL,JKL) , ZEPSCQ )
            ZXOZ(JL) = AMAX1 (PQOF(JL,JKL) / PDP(JL,JKL) , ZEPSCO )
         ENDDO
         JKJ=(JK-1)*NG1P1+1
         JKJPN=JKJ+NG1
         DO JKK=JKJ,JKJPN
            DO JL = 1 , KLON
               ZDPM = ZABLY(JL,3,JKK)
               ZUPM = ZABLY(JL,5,JKK)             * ZDPM / 101325.
               ZUPMH2O = ( ZABLY(JL,5,JKK) + PVGH2O ) * ZDPM / 101325.
               ZUPMCO2 = ( ZABLY(JL,5,JKK) + PVGCO2 ) * ZDPM / 101325.
               ZUPMO3  = ( ZABLY(JL,5,JKK) + PVGO3  ) * ZDPM / 101325.
               ZDUC(JL,JKK)=ZDPM
               ZABLY(JL,12,JKK)=ZXOZ(JL)*ZDPM
               ZABLY(JL,13,JKK)=ZXOZ(JL)*ZUPMO3
               ZABLY(JL,28,JKK)=ZXOZ(JL) * ZDPM
               ZABLY(JL,29,JKK)=ZXOZ(JL) * ZUPMO3
               ZU6=ZXWV(JL) * ZUPM
               ZFPPW= 1.6078 *ZXWV(JL)/(1.+0.608*ZXWV(JL))
               ZABLY(JL, 6,JKK)=ZXWV(JL)*ZUPMH2O
               ZABLY(JL,11,JKK)=ZU6*ZFPPW
               ZABLY(JL,10,JKK)=ZU6*(1.-ZFPPW)
               ZABLY(JL, 9,JKK)=PCCO2   *ZUPMCO2
               ZABLY(JL,8,JKK)= PCCO2 * ZDPM
               ZABLY(JL,16,JKK)= ZNITOX * ZUPM
               ZABLY(JL,18,JKK)= ZMETHA * ZUPM
               ZABLY(JL,22,JKK)= ZDPM
            ENDDO
         ENDDO
      ENDDO
C
C-----------------------------------------------------------------------
C
C
C*         5.    CUMULATIVE ABSORBER AMOUNTS FROM TOP OF ATMOSPHERE
C                --------------------------------------------------
C
      DO JI = 1, KLON*NUA
         PABCU(JI,1,3*KLEV+1) = 0.
      ENDDO
C
      DO JK = 1 , KLEV
         JJ=(JK-1)*NG1P1+1
         JJPN=JJ+NG1
         JKL=KLEV+1-JK
C
C
C*         5.1  CUMULATIVE AEROSOL AMOUNTS FROM TOP OF ATMOSPHERE
C               --------------------------------------------------
C
         JAE1=3*KLEV+1-JJ
         JAE2=3*KLEV+1-(JJ+1)
         JAE3=3*KLEV+1-JJPN
         DO IAE=1,5
            DO JL = 1 , KLON
               ZUAER(JL,IAE) = (CAER(IAE,1)*PAER(JL,JKL,1)
     &      +CAER(IAE,2)*PAER(JL,JKL,2)+CAER(IAE,3)*PAER(JL,JKL,3)
     &      +CAER(IAE,4)*PAER(JL,JKL,4)+CAER(IAE,5)*PAER(JL,JKL,5))
     &      /(ZDUC(JL,JAE1)+ZDUC(JL,JAE2)+ZDUC(JL,JAE3))
            ENDDO
C           CONTRIBUTION GADS AEROSOLS
            DO JAER=1,KEWAER
               DO JL = 1 , KLON
                  ICAE=NDFAER(JAER)
                  IH=KAERH(JL,JK)
                  ZUAER(JL,IAE)=ZUAER(JL,IAE)+CAERN(IH,IAE,ICAE)
     &                 *FCVAER(ICAE)*PAER(JL,JKL,5+JAER)
               ENDDO
            ENDDO
         ENDDO
C
C
C
C*         5.2  INTRODUCES TEMPERATURE EFFECTS ON ABSORBER AMOUNTS
C               --------------------------------------------------
C
         DO JL = 1 , KLON
            ZTAVI=PTAVE(JL,JKL)
            ZTCON(JL)=EXP(6.08*(296./ZTAVI-1.))
            ZTX=ZTAVI-TREF
            ZTX2=ZTX*ZTX
            ZZABLY = ZABLY(JL,6,JAE1)+ZABLY(JL,6,JAE2)+ZABLY(JL,6,JAE3)
            ZUP=AMIN1( AMAX1( 0.5*C10E*ALOG( ZZABLY ) + 5., 0.), 6.0)
            ZCAH1=AT(1,1)+ZUP*(AT(1,2)+ZUP*AT(1,3))
            ZCBH1=BT(1,1)+ZUP*(BT(1,2)+ZUP*BT(1,3))
            ZPSH1(JL)=EXP(ZCAH1*ZTX+ZCBH1*ZTX2)
            ZCAH2=AT(2,1)+ZUP*(AT(2,2)+ZUP*AT(2,3))
            ZCBH2=BT(2,1)+ZUP*(BT(2,2)+ZUP*BT(2,3))
            ZPSH2(JL)=EXP(ZCAH2*ZTX+ZCBH2*ZTX2)
            ZCAH3=AT(3,1)+ZUP*(AT(3,2)+ZUP*AT(3,3))
            ZCBH3=BT(3,1)+ZUP*(BT(3,2)+ZUP*BT(3,3))
            ZPSH3(JL)=EXP(ZCAH3*ZTX+ZCBH3*ZTX2)
            ZCAH4=AT(4,1)+ZUP*(AT(4,2)+ZUP*AT(4,3))
            ZCBH4=BT(4,1)+ZUP*(BT(4,2)+ZUP*BT(4,3))
            ZPSH4(JL)=EXP(ZCAH4*ZTX+ZCBH4*ZTX2)
            ZCAH5=AT(5,1)+ZUP*(AT(5,2)+ZUP*AT(5,3))
            ZCBH5=BT(5,1)+ZUP*(BT(5,2)+ZUP*BT(5,3))
            ZPSH5(JL)=EXP(ZCAH5*ZTX+ZCBH5*ZTX2)
            ZCAH6=AT(6,1)+ZUP*(AT(6,2)+ZUP*AT(6,3))
            ZCBH6=BT(6,1)+ZUP*(BT(6,2)+ZUP*BT(6,3))
            ZPSH6(JL)=EXP(ZCAH6*ZTX+ZCBH6*ZTX2)
C
            ZZABLY = ZABLY(JL,9,JAE1)+ZABLY(JL,9,JAE2)+ZABLY(JL,9,JAE3)
            ZALUP = C10E * ALOG ( ZZABLY )
            ZUP   = AMAX1( 0.0 , 5.0 + 0.5 * ZALUP )
            ZPSC2(JL) = (ZTAVI/TREF) ** ZUP
            ZCAC8=AT(8,1)+ZUP*(AT(8,2)+ZUP*(AT(8,3)))
            ZCBC8=BT(8,1)+ZUP*(BT(8,2)+ZUP*(BT(8,3)))
            ZPSC3(JL)=EXP(ZCAC8*ZTX+ZCBC8*ZTX2)
            ZPHIO(JL)=EXP(OCT(1)*ZTX+OCT(2)*ZTX2)
            ZPSIO(JL)=EXP(2.*(OCT(3)*ZTX+OCT(4)*ZTX2))
            ZPHIO2(JL)=EXP(ODT(1)*ZTX+ODT(2)*ZTX2)
            ZPSIO2(JL)=EXP(2.*(ODT(3)*ZTX+ODT(4)*ZTX2))
            ZZABNI=ZABLY(JL,16,JAE1)+ZABLY(JL,16,JAE2)+ZABLY(JL,16,JAE3)
            ZN2OUP = AMAX1(-1.5,5.0 + 0.5*C10E*ALOG(ZZABNI))
            ZCAN1 = CT(1,1)+ZN2OUP*(CT(1,2)+ZN2OUP*CT(1,3))
            ZCBN1 = DT(1,1)+ZN2OUP*(DT(1,2)+ZN2OUP*DT(1,3))
            ZPSN2(JL) = EXP(ZCAN1*ZTX+ZCBN1*ZTX2)
            ZCAN2 = CT(2,1)+ZN2OUP*(CT(2,2)+ZN2OUP*CT(2,3))
            ZCBN2 = DT(2,1)+ZN2OUP*(DT(2,2)+ZN2OUP*DT(2,3))
            ZPSN3(JL) = EXP(ZCAN2*ZTX+ZCBN2*ZTX2)
            ZCAN3 = CT(3,1)+ZN2OUP*(CT(3,2)+ZN2OUP*CT(3,3))
            ZCBN3 = DT(3,1)+ZN2OUP*(DT(3,2)+ZN2OUP*DT(3,3))
            ZPSN6(JL) = EXP(ZCAN3*ZTX+ZCBN3*ZTX2)
            ZZABME=ZABLY(JL,18,JAE1)+ZABLY(JL,18,JAE2)+ZABLY(JL,18,JAE3)
            ZCH4UP = AMAX1(-1.5,5.0 + 0.5*C10E*ALOG(ZZABME))
            ZCAM4 = CT(4,1)+ZCH4UP*(CT(4,2)+ZCH4UP*CT(4,3))
            ZCBM4 = DT(4,1)+ZCH4UP*(DT(4,2)+ZCH4UP*DT(4,3))
            ZPSM3(JL) = EXP(ZCAM4*ZTX+ZCBM4*ZTX2)
            ZCAM5 = CT(5,1)+ZCH4UP*(CT(5,2)+ZCH4UP*CT(5,3))
            ZCBM5 = DT(5,1)+ZCH4UP*(DT(5,2)+ZCH4UP*DT(5,3))
            ZPSM6(JL) = EXP(ZCAM5*ZTX+ZCBM5*ZTX2)
         ENDDO
C
         DO JKK=JJ,JJPN
            JC=3*KLEV+1-JKK
            JCP1=JC+1
            DO JL = 1 , KLON
               PABCU(JL,10,JC)=PABCU(JL,10,JCP1)
     &              +ZABLY(JL,10,JC)           *ZDIFF
               PABCU(JL,11,JC)=PABCU(JL,11,JCP1)
     &              +ZABLY(JL,11,JC)*ZTCON(JL)*ZDIFF
C
               PABCU(JL,12,JC)=PABCU(JL,12,JCP1)
     &              +ZABLY(JL,12,JC)*ZPHIO(JL)*ZDIFF
               PABCU(JL,13,JC)=PABCU(JL,13,JCP1)
     &              +ZABLY(JL,13,JC)*ZPSIO(JL)*ZDIFF
C
               PABCU(JL,28,JC)=PABCU(JL,28,JCP1)
     &              +ZABLY(JL,28,JC)*ZPHIO2(JL)*ZDIFF
               PABCU(JL,29,JC)=PABCU(JL,29,JCP1)
     &              +ZABLY(JL,29,JC)*ZPSIO2(JL)*ZDIFF
C
               PABCU(JL,7,JC)=PABCU(JL,7,JCP1)
     &              +ZABLY(JL,9,JC)*ZPSC2(JL)*ZDIFF
               PABCU(JL,8,JC)=PABCU(JL,8,JCP1)
     &              +ZABLY(JL,9,JC)*ZPSC3(JL)*ZDIFF
               PABCU(JL,9,JC)=PABCU(JL,9,JCP1)
     &              +ZABLY(JL,9,JC)*ZPSC3(JL)*ZDIFF
C
               PABCU(JL,1,JC)=PABCU(JL,1,JCP1)
     &              +ZABLY(JL,6,JC)*ZPSH1(JL)*ZDIFF
               PABCU(JL,2,JC)=PABCU(JL,2,JCP1)
     &              +ZABLY(JL,6,JC)*ZPSH2(JL)*ZDIFF
               PABCU(JL,3,JC)=PABCU(JL,3,JCP1)
     &              +ZABLY(JL,6,JC)*ZPSH5(JL)*ZDIFF
               PABCU(JL,4,JC)=PABCU(JL,4,JCP1)
     &              +ZABLY(JL,6,JC)*ZPSH3(JL)*ZDIFF
               PABCU(JL,5,JC)=PABCU(JL,5,JCP1)
     &              +ZABLY(JL,6,JC)*ZPSH4(JL)*ZDIFF
               PABCU(JL,6,JC)=PABCU(JL,6,JCP1)
     &              +ZABLY(JL,6,JC)*ZPSH6(JL)*ZDIFF
C
               PABCU(JL,14,JC)=PABCU(JL,14,JCP1)
     &              +ZABLY(JL,16,JC)*ZPSN2(JL)*ZDIFF
               PABCU(JL,15,JC)=PABCU(JL,15,JCP1)
     &              +ZABLY(JL,16,JC)*ZPSN3(JL)*ZDIFF
               PABCU(JL,16,JC)=PABCU(JL,16,JCP1)
     &              +ZABLY(JL,16,JC)*ZPSN6(JL)*ZDIFF
C
               PABCU(JL,17,JC)=PABCU(JL,17,JCP1)
     &              +ZABLY(JL,18,JC)*ZPSM3(JL)*ZDIFF
               PABCU(JL,18,JC)=PABCU(JL,18,JCP1)
     &              +ZABLY(JL,18,JC)*ZPSM6(JL)*ZDIFF
C
               PABCU(JL,19,JC)=PABCU(JL,19,JCP1)
     &              +ZABLY(JL,22,JC)*PCFCABS(1)*ZDICFC
               PABCU(JL,20,JC)=PABCU(JL,20,JCP1)
     &              +ZABLY(JL,22,JC)*PCFCABS(2)*ZDICFC
               PABCU(JL,21,JC)=PABCU(JL,21,JCP1)
     &              +ZABLY(JL,22,JC)*PCFCABS(3)*ZDICFC
               PABCU(JL,22,JC)=PABCU(JL,22,JCP1)
     &              +ZABLY(JL,22,JC)*PCFCABS(4)*ZDICFC
C
               PABCU(JL,23,JC)=PABCU(JL,23,JCP1)
     &              +ZUAER(JL,1)    *ZDUC(JL,JC)*ZDIFF
               PABCU(JL,24,JC)=PABCU(JL,24,JCP1)
     &              +ZUAER(JL,2)    *ZDUC(JL,JC)*ZDIFF
               PABCU(JL,25,JC)=PABCU(JL,25,JCP1)
     &              +ZUAER(JL,3)    *ZDUC(JL,JC)*ZDIFF
               PABCU(JL,26,JC)=PABCU(JL,26,JCP1)
     &              +ZUAER(JL,4)    *ZDUC(JL,JC)*ZDIFF
               PABCU(JL,27,JC)=PABCU(JL,27,JCP1)
     &              +ZUAER(JL,5)    *ZDUC(JL,JC)*ZDIFF
            ENDDO
         ENDDO
C
      ENDDO
C
C
      RETURN
C
      END
