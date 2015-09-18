C
C
C**** *OZONE* - FOR THE YEARLY CYCLE OF THE OZONE DISTIBUTION.
C
C     J.F.GELEYN     E.C.M.W.F.     03/06/82.
C
C     PURPOSE.
C     --------
C
C          THIS ROUTINE COMPUTES INSTANTANEOUS VALUES OF A T5 SPETRAL
C     DISTRIBUTION FOR TWO OZONE PARAMETERS (TOTAL QUANTITY AND PRESSURE
C     AT THE MAXIMUM OF CONCENTRATION,BOTH IN *PASCAL) FROM THE TIME
C     OF THE YEAR (SEE *ORBIT*).
C
C**   INTERFACE.
C     ----------
C
C          *OZONE* IS CALLED FROM *PHYSC* AT THE FIRST LATITUDE ROW AT
C     THE TIME OF A FULL RADIATION COMPUTATION.
C          THERE ARE FIVE DUMMY ARGUMENTS: *PYTIME* IS THE TIME OF THE
C     YEAR (IN RADIANS).
C                                          *POZQC*, *POZQS*, *POZHC* AND
C     *POZHS* ARE ARRAYS FOR THE T5 DISTRIBUTIONS (*Q FOR QUANTITY, *H
C     FOR HEIGHT, *C FOR COSINE AND *S FOR SINE).
C
C     METHOD.
C     -------
C
C          STAIGHTFORWARD, A SECOND ORDER *FOURIER DEVELOPMENT FOR THE
C     TIME OF THE YEAR.
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
      SUBROUTINE OZONE (PYTIME,POZQC,POZQS,POZHC,POZHS)
C
      IMPLICIT NONE
C
C     Dummy Arguments
C
      REAL, INTENT(IN)    :: PYTIME
      REAL, INTENT(INOUT) :: POZQC(21),POZQS(15),POZHC(21),POZHS(15)
C
C     Local Variables
C
      INTEGER :: JMN
      REAL    :: ZC1YT, ZC2YT, ZS1YT, ZS2YT, ZYTIME
C
C*    DATA STATEMENTS.
C     ---- -----------
C
C          *ZOZ Q/H C/S N* (N=0,4) CORRESPONDS TO THE *POZ Q/H C/S*
C     (SEE ABOVE) AND TO THE FIVE TERMS OF THE *FOURIER DEVELOPMENT.
C
      REAL ::   ZOZQC0(21),ZOZQC1(21),ZOZQC2(21),ZOZQC3(21),ZOZQC4(21),
     &          ZOZQS0(15),ZOZQS1(15),ZOZQS2(15),ZOZQS3(15),ZOZQS4(15),
     &          ZOZHC0(21),ZOZHC1(21),ZOZHC2(21),ZOZHC3(21),ZOZHC4(21),
     &          ZOZHS0(15),ZOZHS1(15),ZOZHS2(15),ZOZHS3(15),ZOZHS4(15)
C
      DATA ZOZQC0/
     &+.6012E-01,+.1887E-02,+.7410E-02,+.9950E-03,-.1426E-02,-.2072E-03,
     &           -.4954E-03,+.7955E-05,-.3701E-03,+.4116E-04,-.4163E-04,
     &                      -.2933E-03,+.2154E-04,-.2849E-03,-.1604E-03,
     &                                 -.1054E-03,+.4974E-03,+.1047E-03,
     &                                            +.8323E-04,+.2874E-03,
     &                                                       +.1333E-03/
      DATA ZOZQS0/
     &           +.4210E-03,-.9591E-03,+.2811E-03,-.2257E-03,-.1713E-03,
     &                      -.3538E-03,+.1095E-03,-.4390E-03,-.5605E-05,
     &                                 +.1478E-03,+.2849E-03,+.3430E-03,
     &                                            +.8248E-04,+.1442E-03,
     &                                                       -.1375E-04/
      DATA ZOZHC0/
     &+.3166E+04,+.8663E+02,+.9401E+03,+.1999E+02,-.3530E+03,-.3311E+02,
     &           -.4903E+02,-.4015E+00,-.1333E+02,+.5675E+01,+.7221E+01,
     &                      -.3001E+02,+.7570E+01,-.1142E+02,-.1365E+02,
     &                                 -.1502E+02,+.4911E+02,+.1425E+02,
     &                                            +.8983E+01,+.3064E+02,
     &                                                       +.1693E+02/
      DATA ZOZHS0/
     &           +.4231E+02,-.7391E+02,+.1273E+02,+.2086E+02,-.1597E+02,
     &                      -.3591E+02,+.1059E+02,-.2779E+02,-.6923E+01,
     &                                 +.1397E+02,+.2387E+02,+.2883E+02,
     &                                            +.8626E+01,+.1607E+02,
     &                                                       -.2676E+01/
      DATA ZOZQC1/
     &+.7090E-04,+.4930E-05,+.6829E-03,+.1897E-03,+.7226E-04,-.2807E-03,
     &           +.4970E-04,-.1753E-03,-.7843E-04,-.1649E-03,-.1037E-03,
     &                      -.4830E-04,-.6304E-04,-.1100E-03,-.7952E-04,
     &                                 +.1326E-04,+.2599E-04,+.9926E-05,
     &                                            -.9247E-05,-.3521E-05,
     &                                                       -.1780E-04/
      DATA ZOZQS1/
     &           +.6333E-04,+.1145E-03,+.1192E-03,+.4934E-04,+.2699E-04,
     &                      +.3684E-04,-.2395E-05,+.2045E-04,-.8684E-04,
     &                                 +.5301E-04,-.4176E-05,+.4103E-04,
     &                                            +.2783E-04,+.1754E-04,
     &                                                       +.1116E-04/
      DATA ZOZHC1/
     &-.3450E+02,+.2148E+03,+.3376E+02,+.6535E+02,-.1564E+02,-.4273E+02,
     &           +.9553E+01,-.4647E+01,-.6129E+01,-.6727E+01,-.6761E+01,
     &                      -.2467E+01,-.2181E+01,-.5361E+01,-.2395E+01,
     &                                 +.5952E+00,+.2106E+01,-.1367E+01,
     &                                            -.2349E+01,+.3532E+00,
     &                                                       -.3169E+01/
      DATA ZOZHS1/
     &           +.3977E+01,+.5032E+01,+.6226E+01,-.3625E+00,-.1373E+01,
     &                      +.4600E+01,+.4312E+01,+.2882E+01,-.6351E+01,
     &                                 +.5731E+01,-.2574E+01,+.3235E+00,
     &                                            +.2806E+01,+.8133E+00,
     &                                                       +.2032E+01/
      DATA ZOZQC2/
     &+.8571E-03,+.3086E-02,+.9287E-03,+.2787E-03,+.1826E-03,-.1006E-03,
     &           +.1092E-03,-.1266E-03,+.5372E-04,-.1188E-03,-.3285E-04,
     &                      -.1783E-04,-.3018E-05,-.8709E-04,-.8707E-04,
     &                                 +.8633E-04,+.3530E-04,+.4863E-04,
     &                                            +.3917E-05,-.3252E-04,
     &                                                       -.1936E-06/
      DATA ZOZQS2/
     &           -.8822E-04,+.1341E-03,+.3095E-04,+.8230E-04,+.2735E-04,
     &                      +.1714E-04,-.9406E-04,+.1912E-04,-.5402E-04,
     &                                 +.3571E-04,+.3897E-04,+.4487E-04,
     &                                            +.3079E-04,+.3196E-04,
     &                                                       -.2391E-05/
      DATA ZOZHC2/
     &+.5216E+02,+.1613E+03,+.3284E+02,-.7670E+02,-.9548E+01,+.1608E+02,
     &           +.1023E+02,-.1090E+02,+.2748E+01,-.3846E+01,-.4135E+01,
     &                      +.1255E+01,-.3301E-01,-.5273E+01,-.7247E+01,
     &                                 +.1387E+02,+.4184E+01,+.6495E+01,
     &                                            +.2944E+01,-.1947E+01,
     &                                                       +.1132E+01/
      DATA ZOZHS2/
     &           -.1968E+02,+.1192E+02,-.1194E+01,+.1084E+01,+.2946E+01,
     &                      +.2630E+01,-.1256E+02,+.1395E+01,-.2222E+01,
     &                                 +.4864E+01,+.6450E+01,+.5568E+01,
     &                                            +.5292E+01,+.4876E+01,
     &                                                       -.7579E+00/
      DATA ZOZQC3/
     &-.2759E-03,-.2781E-03,-.1087E-03,-.1633E-03,-.3627E-04,-.4242E-04,
     &           +.6045E-05,-.1703E-04,+.4562E-04,-.1009E-04,+.2663E-04,
     &                      -.1786E-04,+.1550E-04,-.9135E-06,+.2372E-04,
     &                                 +.1100E-05,+.2299E-04,+.4659E-05,
     &                                            +.2423E-05,+.7321E-05,
     &                                                       +.8852E-05/
      DATA ZOZQS3/
     &           -.3678E-04,-.2219E-04,-.3911E-04,-.4398E-04,-.1142E-04,
     &                      -.9121E-05,-.2011E-04,+.4711E-06,-.3775E-05,
     &                                 +.3866E-05,+.2400E-04,+.2043E-04,
     &                                            -.1824E-05,-.5550E-05,
     &                                                       +.2506E-05/
      DATA ZOZHC3/
     &-.1534E+03,-.2095E+02,-.1006E+03,-.7385E+01,+.5203E+01,+.9434E+00,
     &           -.3814E+00,-.3175E+01,+.3366E+01,+.3378E+00,+.2740E+00,
     &                      -.2669E+01,+.8452E+00,+.3498E+00,+.2192E+01,
     &                                 -.4024E+00,+.1544E+01,-.4588E+00,
     &                                            +.6998E+00,+.6263E+00,
     &                                                       +.1228E+01/
      DATA ZOZHS3/
     &           -.3588E+01,+.2076E+00,-.2088E+01,-.4159E+01,+.2244E+00,
     &                      -.7751E+00,-.2749E+01,+.7234E+00,+.4390E+00,
     &                                 -.1646E+00,+.1700E+01,+.1046E+01,
     &                                            -.7856E+00,-.1644E+01,
     &                                                       +.2648E+00/
      DATA ZOZQC4/
     &-.1460E-03,+.3422E-03,-.3529E-04,+.1791E-03,-.1917E-03,-.2558E-04,
     &           +.6547E-04,+.6401E-04,+.4823E-04,+.7084E-05,+.2895E-04,
     &                      -.1561E-04,+.8179E-06,+.1028E-04,-.7667E-05,
     &                                 -.4347E-05,+.7293E-05,-.5735E-05,
     &                                            +.7838E-05,-.2933E-05,
     &                                                       +.3686E-05/
      DATA ZOZQS4/
     &           -.4560E-05,-.5292E-04,-.1252E-04,+.1850E-04,-.2273E-04,
     &                      +.6552E-05,+.1422E-04,-.6545E-05,+.7998E-06,
     &                                 +.2845E-04,+.2497E-04,+.2844E-04,
     &                                            +.3855E-06,-.1487E-04,
     &                                                       +.1954E-05/
      DATA ZOZHC4/
     &+.9260E+01,-.9055E+01,+.5460E+01,-.7603E+01,-.3329E+02,-.1048E+02,
     &           +.9328E+01,+.4597E+01,+.3827E+01,-.3201E+01,+.1708E+01,
     &                      -.1548E+01,-.5323E+00,+.3039E+01,+.5740E+00,
     &                                 +.1353E+00,-.2354E+01,+.2818E+00,
     &                                            +.1113E+01,-.1891E+01,
     &                                                       -.3074E+00/
      DATA ZOZHS4/
     &           -.2446E+01,+.4199E+01,-.2571E+01,+.8194E+01,+.4206E+00,
     &                      +.3856E+01,+.1159E+01,+.2547E+01,-.1314E+01,
     &                                 +.2331E+01,+.1144E+01,-.4408E+00,
     &                                            -.6797E+00,-.2598E+01,
     &                                                       +.8953E+00/
C
C     ------------------------------------------------------------------
C
C*         1.     PRELIMINARY SETTING.
C                 ----------- --------
C
      ZYTIME=PYTIME
C
C     ------------------------------------------------------------------
C
C*         2.     COMPUTATIONS.
C                 -------------
C
      ZC1YT=COS(ZYTIME)
      ZS1YT=SIN(ZYTIME)
      ZC2YT=ZC1YT**2-ZS1YT**2
      ZS2YT=2.*ZS1YT*ZC1YT
      DO JMN=1,21
         POZQC(JMN)=ZOZQC0(JMN)+2.*(ZOZQC1(JMN)*ZC1YT+ZOZQC2(JMN)*ZS1YT
     &        +ZOZQC3(JMN)*ZC2YT+ZOZQC4(JMN)*ZS2YT)
         POZHC(JMN)=ZOZHC0(JMN)+2.*(ZOZHC1(JMN)*ZC1YT+ZOZHC2(JMN)*ZS1YT
     &        +ZOZHC3(JMN)*ZC2YT+ZOZHC4(JMN)*ZS2YT)
      ENDDO
      DO JMN=1,15
         POZQS(JMN)=ZOZQS0(JMN)+2.*(ZOZQS1(JMN)*ZC1YT+ZOZQS2(JMN)*ZS1YT
     &        +ZOZQS3(JMN)*ZC2YT+ZOZQS4(JMN)*ZS2YT)
         POZHS(JMN)=ZOZHS0(JMN)+2.*(ZOZHS1(JMN)*ZC1YT+ZOZHS2(JMN)*ZS1YT
     &        +ZOZHS3(JMN)*ZC2YT+ZOZHS4(JMN)*ZS2YT)
      ENDDO
C
C     ------------------------------------------------------------------
C
C*         3.     RETURN.
C                 -------
C
      RETURN
      END SUBROUTINE OZONE
