C                                                                       
C*** *INILUC*  PRESET CONSTANTS IN *COMCON*.                            
C                                                                       
C             R. PODZUN    DKRZ     19/01/1995.                         
C                                                                       
C     PURPOSE.                                                          
C     --------                                                          
C                                                                       
C             PRESET CONSTANTS IN *COMCON*, *YOTLUC*.                  
C                                                                       
C**   INTERFACE.                                                        
C     ----------                                                        
C                                                                       
C             *INILUC* IS CALLED FROM *INIT*.                           
C                                                                       
C     EXTERNALS.                                                        
C     ----------                                                        
C                                                                       
C             NONE.                                                     
C    
      SUBROUTINE INILUC
C
      IMPLICIT NONE
C
      INCLUDE "COMCON"
      INCLUDE "YOTLUC"
C
      INTEGER :: I, IT
      REAL    :: TT, Z5ALSCP, Z5ALVCP, ZALSDCP, ZALVDCP, ZAVI1, ZAVI2,  
     &           ZAVI4, ZAVI5, ZAVL1, ZAVL2, ZCVM4, ZCVM5, ZLDCP, ZT,
     &           ZAVL3, ZAVL4, ZAVL5, ZAVM1, ZAVM2, ZAVM3, ZAVM4, ZAVM5,
     &           ZAVI3
C
C     INITIALISE LOOKUP TABLES FOR CUADJTQ                              
C                                                                       
      Z5ALVCP=C5LES*ALV/CPD                                             
      Z5ALSCP=C5IES*ALS/CPD                                             
      ZALVDCP=ALV/CPD
      ZALSDCP=ALS/CPD
      ZAVL1=-6096.9385
      ZAVL2=21.2409642
      ZAVL3=-2.711193
      ZAVL4=1.673952
      ZAVL5=2.433502
      ZAVI1=-6024.5282
      ZAVI2=29.32707
      ZAVI3=1.0613868
      ZAVI4=-1.3198825
      ZAVI5=-0.49382577
      TT=100.
      DO I=JPTLUCU1,JPTLUCU2 
         ZT=TT*1000.
         IT=NINT(ZT)
         IF(TT-TMELT.GT.0.) THEN
            ZCVM4=C4LES
            ZCVM5=Z5ALVCP
            ZLDCP=ZALVDCP
            ZAVM1=ZAVL1
            ZAVM2=ZAVL2
            ZAVM3=ZAVL3
            ZAVM4=ZAVL4
            ZAVM5=ZAVL5
         ELSE
            ZCVM4=C4IES
            ZCVM5=Z5ALSCP
            ZLDCP=ZALSDCP
            ZAVM1=ZAVI1
            ZAVM2=ZAVI2
            ZAVM3=ZAVI3
            ZAVM4=ZAVI4
            ZAVM5=ZAVI5
         ENDIF
         TLUCUC(IT)=ZLDCP
         TLUCUA(IT)=EXP((ZAVM1/TT+ZAVM2+ZAVM3*TT/100.
     &        +ZAVM4*TT*TT/1.E5+ZAVM5*LOG(TT)))*RD/RV
         TLUCUB(IT)=ZCVM5*(1./(TT-ZCVM4))**2  
         TLUCUAW(IT)=EXP((ZAVL1/TT+ZAVL2+ZAVL3*TT/100.
     &        +ZAVL4*TT*TT/1.E5+ZAVL5*LOG(TT)))*RD/RV
         TT=TT+0.001
      ENDDO
C
      RETURN
      END SUBROUTINE INILUC
