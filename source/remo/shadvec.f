C
C     SUBROUTINE SHADVEC
C
C**   AUFRUF   :   CALL SHADVEC IN UP *CHEM*
C**   ZWECK    :   MODIFIKATION VON SPURENSTOFFKONZENTRATIONEN 
C**                DURCH HORIZONTALADVECTION
C**   DATUM    :   OKT.96
C**                
C**
C**   EXTERNALS:   -
C**
C**   EINGABE-
C**   PARAMETER:   K:   INDEX FUER VERTIKALE SCHICHT 
C**                ZDT: ZEITSCHRITT FUER VORWAERTSEXTRAPOLATION
C**                IXY: IXY=1: ZUERST X-RICHTUNG, DANN Y-RICHTUNG
C**                     IXY=2: ZUERST Y-RICHTUNG, DANN X-RICHTUNG
C**   AUSGABE-
C**   PARAMETER:   KEINE
C**
C**   COMMON-
C**   BLOECKE  :   PARCHEM, ORG, HIGKON, PHYKON, CHPS3D, PRC
C**
C**   METHODE  :   RAEUMLICH: SMOLARKIEWITZ, 2 KORREKTUREN
C**                ZEITLICH:  EULER VORWAERTS
C**   FEHLERBE-
C**   HANDLUNG :   ---
C**   VERFASSER:   B.LANGMANN
C
C
      SUBROUTINE SHADVEC (K, ZDT, IXY, L, 
     &                    U, V, TRAC, ACPHIR, CPHI,NSPRED,NX,
     &                    IOFF,IEOFF,JOFF,JEOFF,ION,IEON,JON,JEON)
C
      IMPLICIT NONE
C
      INCLUDE "parorg.h"
      INCLUDE "higkon.h"
C-----------------------------------------------------------------------
C     Dummy Arguments
C
      INTEGER, INTENT(IN) :: K,NX,IOFF,IEOFF,JOFF,JEOFF,L,
     &                       ION,IEON,JON,JEON,NSPRED,IXY
      REAL,    INTENT(IN) :: ZDT
C
C     ATMOSPHAERENFELDER                                                 
      REAL, INTENT(IN) :: U(IE,JE,KE,3), V(IE,JE,KE,3)  

C     KONZENTRATIONSFELDER                                                 
      REAL, INTENT(INOUT) :: TRAC(IE,JE,KE,NSPRED)

C     EXTERNE PARAMETER
      REAL, INTENT(IN) :: ACPHIR(JE,2), CPHI(JE,2)
C-----------------------------------------------------------------------
C     Local Declarations
C
C     LOKALE FELDER
      INTEGER, PARAMETER :: NCORR=1
C
      REAL ::   WT(IEOFF,JEOFF,NCORR+1),DTODX(JEOFF),DTODY(JEOFF),
     &          UCORR(IEOFF,JEOFF,NCORR+1)
C
      INTEGER :: I, II, J, KCB
C
C.....FUNKTIONEN 
C     TODO: Find out, why this macro statement gives slightly different
C           results than the function statement below.
      REAL F,UD,CI,CIP1,UIH,DT0DX
C
      F(CI,CIP1,UIH)=(((ABS(UIH)+UIH)/2.)*CI  ) +
     &     (((UIH-ABS(UIH))/2.)*CIP1)
      UD(CI,CIP1,UIH,DT0DX)=(ABS(UIH)-DT0DX*UIH*UIH)
     &     *(CIP1-CI)/(CI+CIP1+1.E-35)
C-----------------------------------------------------------------------

      DO J=JOFF,JEOFF
         DTODX(J)=ZDT*ACPHIR(J,1)*EDDLAM
         DTODY(J)=ZDT*ACPHIR(J,1)*CPHI(J,2)*EDDPHI
      ENDDO

C.....INITIALISIERUNG VON WT
C
      DO KCB=1,NCORR+1
         DO J=JOFF,JEOFF
            DO I=IOFF,IEOFF
               WT(I,J,KCB)=TRAC(I,J,K,L)
            ENDDO
         ENDDO
      ENDDO

C***********************************************************************
      IF(IXY.EQ.1) THEN
C***********************************************************************
C
C........ZUNAECHST ADVEKTION IN X-RICHTUNG
C
C        UPSTREAM
         DO J=JOFF,JEOFF
            DO I=IOFF+1,IEOFF-1
               WT(I,J,1)=TRAC(I,J,K,L)-DTODX(J)* (
     &              F(TRAC(I,J,K,L),TRAC(I+1,J,K,L),U(I,J,K,NX))-
     &              F(TRAC(I-1,J,K,L),TRAC(I,J,K,L),U(I-1,J,K,NX)) )
            ENDDO
         ENDDO
C
C........ANTIDIFFUSIONSKORREKTUR IN X-RICHTUNG
C       
C        INITIALISIERUNG DER KORREKTURGESCHWINDIGKEIT 
         DO J=JOFF,JEOFF
            DO I=IOFF,IEOFF
               UCORR(I,J,1)= U(I,J,K,NX)
            ENDDO
         ENDDO
C
         DO II=1,NCORR

C           UCORR: INNERHALB DES I-GEBIETES

            DO J=JOFF,JEOFF
               DO I=IOFF+1,IEOFF-2
                  UCORR(I,J,II+1)= UD( WT(I,J,II), WT(I+1,J,II),
     &                 UCORR(I,J,II), DTODX(J) )
               ENDDO
            ENDDO

C           WESTLICHE RANDBEDINGUNG FUER I=IOFF

            IF(ION .EQ. 3) THEN

               DO J=JOFF,JEOFF
                  IF(U(IOFF,J,K,NX).LT.0.) THEN
                     WT(IOFF,J,II)=WT(IOFF+1,J,II)-U(IOFF+1,J,K,NX)
     &                    /U(IOFF,J,K,NX)* 
     &                    (WT(IOFF+2,J,II)-WT(IOFF+1,J,II))
                  ENDIF

                  UCORR(IOFF,J,II+1)=UD(WT(IOFF,J,II),WT(IOFF+1,J,II),
     &                 UCORR(IOFF,J,II),DTODX(J))

                  IF (U(IOFF,J,K,NX).LT.0.) THEN
                     WT(IOFF+1,J,II+1)=WT(IOFF+1,J,II)
                     WT(IOFF+2,J,II+1)=WT(IOFF+2,J,II)- DTODX(J)*(
     &                    F(WT(IOFF+2,J,II),WT(IOFF+3,J,II),
     &                    UCORR(IOFF+2,J,II+1))-
     &                    UCORR(IOFF+1,J,II+1)*WT(IOFF+2,J,II) )
                  ELSE
                     DO I=IOFF+1,IOFF+2
                        WT(I,J,II+1)= WT(I,J,II) - DTODX(J)*(
     &                    F(WT(I,J,II),WT(I+1,J,II),UCORR(I,J,II+1))-
     &                    F(WT(I-1,J,II),WT(I,J,II),UCORR(I-1,J,II+1)) )
                     ENDDO
                  ENDIF
               ENDDO

            ELSE

               DO J=JOFF,JEOFF
                  UCORR(IOFF,J,II+1)=UD(WT(IOFF,J,II),WT(IOFF+1,J,II),
     &                 UCORR(IOFF,J,II),DTODX(J))
               ENDDO

            ENDIF

C           OESTLICHE RANDBEDINGUNG FUER  I=IE-1

            IF(IEON .EQ. 3) THEN

               DO J=JOFF,JEOFF
                  IF (U(IEOFF-1,J,K,NX).GT.0.) THEN
                     WT(IEOFF,J,II)=WT(IEOFF-1,J,II)-U(IEOFF-2,J,K,NX)/
     &                    U(IEOFF-1,J,K,NX)*
     &                    (WT(IEOFF-2,J,II)-WT(IEOFF-1,J,II))
                  ENDIF

                  UCORR(IEOFF-1,J,II+1)=UD(WT(IEOFF-1,J,II),
     &                 WT(IEOFF,J,II),UCORR(IEOFF-1,J,II),DTODX(J))

                  IF (U(IEOFF-1,J,K,NX).GT.0.) THEN
                     WT(IEOFF-1,J,II+1)=WT(IEOFF-1,J,II)
                     WT(IEOFF-2,J,II+1)=WT(IEOFF-2,J,II)-DTODX(J)*(
     &                    UCORR(IEOFF-2,J,II+1)*WT(IEOFF-2,J,II)-
     &                    F(WT(IEOFF-3,J,II),WT(IEOFF-2,J,II),
     &                    UCORR(IEOFF-3,J,II+1)))
                  ELSE
                     DO I=IEOFF-2,IEOFF-1
                        WT(I,J,II+1)= WT(I,J,II) - DTODX(J)*(
     &                    F(WT(I,J,II),WT(I+1,J,II),UCORR(I,J,II+1))-
     &                    F(WT(I-1,J,II),WT(I,J,II),UCORR(I-1,J,II+1)) )
                     ENDDO
                  ENDIF
               ENDDO

            ELSE

               DO J=JOFF,JEOFF
                  UCORR(IEOFF-1,J,II+1)=UD(WT(IEOFF-1,J,II),
     &                 WT(IEOFF,J,II),UCORR(IEOFF-1,J,II),DTODX(J))
               ENDDO

            ENDIF

C           KORREKTURN INNERHALB DES I-GEBIETES

            DO J=JOFF,JEOFF
               DO I=IOFF+ION,IEOFF-IEON
                  WT(I,J,II+1)= WT(I,J,II) - DTODX(J)*(
     &                 F(WT(I,J,II),WT(I+1,J,II),UCORR(I,J,II+1))-
     &                 F(WT(I-1,J,II),WT(I,J,II),UCORR(I-1,J,II+1)) )

               ENDDO
            ENDDO

         ENDDO

C
C........DANN ADVEKTION IN Y-RICHTUNG
C
C        ZURUECK ZU TRAC 
C        INITIALISIERUNG DER KORREKTURGESCHWINDIGKEIT
         DO J=JOFF,JEOFF
            DO I=IOFF,IEOFF
               TRAC(I,J,K,L)=WT(I,J,ncorr+1)
               UCORR(I,J,1)=V(I,J,K,NX)
               WT(I,J,1)=TRAC(I,J,K,L)
            ENDDO
         ENDDO

C        UPSTREAM 
         DO J=JOFF+1,JEOFF-1
            DO I=IOFF,IEOFF
               WT(I,J,1)=TRAC(I,J,K,L)-DTODY(J)* (
     &              F(TRAC(I,J,K,L),TRAC(I,J+1,K,L),V(I,J,K,NX))-
     &              F(TRAC(I,J-1,K,L),TRAC(I,J,K,L),V(I,J-1,K,NX)) )
            ENDDO
         ENDDO
C
C........ANTIDIFFUSIONSKORREKTUR IN Y-RICHTUNG
C       
         DO II=1,NCORR

C        UCORR: INNERHALB DES J-GEBIETES

            DO I=IOFF,IEOFF
               DO J=JOFF+1,JEOFF-2
                  UCORR(I,J,II+1)=UD( WT(I,J,II), WT(I,J+1,II),
     &                 UCORR(I,J,II), DTODY(J))
               ENDDO
            ENDDO

C           SUEDLICHE RANDBEDINGUNG FUER J=JOFF

            IF(JON .EQ. 3) THEN

               DO I=IOFF,IEOFF
                  IF(V(I,JOFF,K,NX).LT.0.) THEN
                     WT(I,JOFF,II)=WT(I,JOFF+1,II)-V(I,JOFF+1,K,NX)/
     &                    V(I,JOFF,K,NX)* 
     &                    (WT(I,JOFF+2,II)-WT(I,JOFF+1,II))
                  ENDIF

                  UCORR(I,JOFF,II+1)=UD(WT(I,JOFF,II),WT(I,JOFF+1,II),
     &                 UCORR(I,JOFF,II),DTODY(JOFF))

                  IF (V(I,JOFF,K,NX).LT.0.) THEN
                     WT(I,JOFF+1,II+1)=WT(I,JOFF+1,II)
                     WT(I,JOFF+2,II+1)=WT(I,JOFF+2,II)-DTODY(JOFF+2)*(
     &                    F(WT(I,JOFF+2,II),WT(I,JOFF+3,II),
     &                    UCORR(I,JOFF+2,II+1)) -
     &                    UCORR(I,JOFF+1,II+1)*WT(I,JOFF+2,II) )
                  ELSE
                     DO J=JOFF+1,JOFF+2
                        WT(I,J,II+1)= WT(I,J,II) - DTODY(J)*(
     &                    F(WT(I,J,II),WT(I,J+1,II),UCORR(I,J,II+1))-
     &                    F(WT(I,J-1,II),WT(I,J,II),UCORR(I,J-1,II+1)) )
                     ENDDO
                  ENDIF
               ENDDO

            ELSE

               DO I=IOFF,IEOFF
                  UCORR(I,JOFF,II+1)=UD(WT(I,JOFF,II),WT(I,JOFF+1,II),
     &                 UCORR(I,JOFF,II),DTODY(JOFF))
               ENDDO
            ENDIF

C           NOERDLICHE RANDBEDINGUNG FUER J=JE-1

            IF(JEON .EQ. 3) THEN

               DO I=IOFF,IEOFF
                  IF (V(I,JEOFF-1,K,NX).GT.0.) THEN
                     WT(I,JEOFF,II) = WT(I,JEOFF-1,II)-V(I,JEOFF-2,K,NX)
     &                    /V(I,JEOFF-1,K,NX)*(WT(I,JEOFF-2,II)
     &                    -WT(I,JEOFF-1,II))
                  ENDIF

                  UCORR(I,JEOFF-1,II+1)=UD(WT(I,JEOFF-1,II),
     &                WT(I,JEOFF,II),UCORR(I,JEOFF-1,II),DTODY(JEOFF-1))

                  IF (V(I,JEOFF-1,K,NX).GT.0.) THEN
                     WT(I,JEOFF-1,II+1)=WT(I,JEOFF-1,II)
                     WT(I,JEOFF-2,II+1)=WT(I,JEOFF-2,II)- DTODY(JEOFF-2)
     &                    *(UCORR(I,JEOFF-2,II+1)*WT(I,JEOFF-2,II)-
     &                    F(WT(I,JEOFF-3,II),WT(I,JEOFF-2,II),
     &                    UCORR(I,JEOFF-3,II+1)))
                  ELSE
                     DO J=JEOFF-2,JEOFF-1
                        WT(I,J,II+1)= WT(I,J,II) - DTODY(J)*(
     &                    F(WT(I,J,II),WT(I,J+1,II),UCORR(I,J,II+1))-
     &                    F(WT(I,J-1,II),WT(I,J,II),UCORR(I,J-1,II+1)) )
                     ENDDO
                  ENDIF
               ENDDO

            ELSE

               DO I=IOFF,IEOFF
                  UCORR(I,JEOFF-1,II+1)=UD(WT(I,JEOFF-1,II),
     &                WT(I,JEOFF,II),UCORR(I,JEOFF-1,II),DTODY(JEOFF-1))
               ENDDO
            ENDIF

C           KORREKTUREN INNERHALB DES J-GEBIETES

            DO I=IOFF,IEOFF
               DO J=JOFF+JON,JEOFF-JEON
                  WT(I,J,II+1)= WT(I,J,II) - DTODY(J)*(
     &                 F(WT(I,J,II),WT(I,J+1,II),UCORR(I,J,II+1))-
     &                 F(WT(I,J-1,II),WT(I,J,II),UCORR(I,J-1,II+1)))
               ENDDO
            ENDDO

         ENDDO
C
C***********************************************************************
      ELSE
C***********************************************************************
C
C........ZUNAECHST ADVEKTION IN Y-RICHTUNG
C
C        UPSTREAM 
         DO J=JOFF+1,JEOFF-1
            DO I=IOFF,IEOFF
               WT(I,J,1)=TRAC(I,J,K,L)-DTODY(J)* (
     &              F(TRAC(I,J,K,L),TRAC(I,J+1,K,L),V(I,J,K,NX))-
     &              F(TRAC(I,J-1,K,L),TRAC(I,J,K,L),V(I,J-1,K,NX)) )
            ENDDO
         ENDDO
C
C........ANTIDIFFUSIONSKORREKTUR IN Y-RICHTUNG
C           
C        INITIALISIERUNG DER KORREKTURGESCHWINDIGKEIT 
         DO J=JOFF,JEOFF
            DO I=IOFF,IEOFF
               UCORR(I,J,1)=V(I,J,K,NX)
            ENDDO
         ENDDO

         DO II=1,NCORR

C           UCORR: INNERHALB DES J-GEBIETES

            DO I=IOFF,IEOFF
               DO J=JOFF+1,JEOFF-2
                  UCORR(I,J,II+1)=UD( WT(I,J,II), WT(I,J+1,II),
     &                 UCORR(I,J,II), DTODY(J))
               ENDDO
            ENDDO

C           SUEDLICHE RANDBEDINGUNG FUER J=JOFF

            IF(JON .EQ. 3) THEN

               DO I=IOFF,IEOFF
                  IF(V(I,JOFF,K,NX).LT.0.) THEN
                     WT(I,JOFF,II)=WT(I,JOFF+1,II)-V(I,JOFF+1,K,NX)/
     &                  V(I,JOFF,K,NX)*(WT(I,JOFF+2,II)-WT(I,JOFF+1,II))
                  ENDIF

                  UCORR(I,JOFF,II+1)=UD(WT(I,JOFF,II),WT(I,JOFF+1,II),
     &                 UCORR(I,JOFF,II),DTODY(JOFF))

                  IF (V(I,JOFF,K,NX).LT.0.) THEN
                     WT(I,JOFF+1,II+1)=WT(I,JOFF+1,II)
                     WT(I,JOFF+2,II+1)=WT(I,JOFF+2,II)-DTODY(JOFF+2)*(
     &                    F(WT(I,JOFF+2,II),WT(I,JOFF+3,II),
     &                    UCORR(I,JOFF+2,II+1)) -
     &                    UCORR(I,JOFF+1,II+1)*WT(I,JOFF+2,II) )
                  ELSE
                     DO J=JOFF+1,JOFF+2
                        WT(I,J,II+1)= WT(I,J,II) - DTODY(J)*(
     &                    F(WT(I,J,II),WT(I,J+1,II),UCORR(I,J,II+1))-
     &                    F(WT(I,J-1,II),WT(I,J,II),UCORR(I,J-1,II+1)) )
                     ENDDO
                  ENDIF
               ENDDO

            ELSE

               DO I=IOFF,IEOFF
                  UCORR(I,JOFF,II+1)=UD(WT(I,JOFF,II),WT(I,JOFF+1,II),
     &                 UCORR(I,JOFF,II),DTODY(JOFF))
               ENDDO

            ENDIF

C           NOERDLICHE RANDBEDINGUNG FUER J=JE-1

            IF(JEON .EQ. 3) THEN

               DO I=IOFF,IEOFF
                  IF (V(I,JEOFF-1,K,NX).GT.0.) THEN
                     WT(I,JEOFF,II) = WT(I,JEOFF-1,II)-V(I,JEOFF-2,K,NX)
     &                    /V(I,JEOFF-1,K,NX)*(WT(I,JEOFF-2,II)
     &                    -WT(I,JEOFF-1,II))
                  ENDIF

                  UCORR(I,JEOFF-1,II+1)=UD(WT(I,JEOFF-1,II),
     &                WT(I,JEOFF,II),UCORR(I,JEOFF-1,II),DTODY(JEOFF-1))

                  IF (V(I,JEOFF-1,K,NX).GT.0.) THEN
                     WT(I,JEOFF-1,II+1)=WT(I,JEOFF-1,II)
                     WT(I,JEOFF-2,II+1)=WT(I,JEOFF-2,II)- DTODY(JEOFF-2)
     &                    *(UCORR(I,JEOFF-2,II+1)*WT(I,JEOFF-2,II)-
     &                    F(WT(I,JEOFF-3,II),WT(I,JEOFF-2,II),
     &                    UCORR(I,JEOFF-3,II+1)))
                  ELSE
                     DO J=JEOFF-2,JEOFF-1
                        WT(I,J,II+1)= WT(I,J,II) - DTODY(J)*(
     &                    F(WT(I,J,II),WT(I,J+1,II),UCORR(I,J,II+1))-
     &                    F(WT(I,J-1,II),WT(I,J,II),UCORR(I,J-1,II+1)) )
                     ENDDO
                  ENDIF
               ENDDO

            ELSE

               DO I=IOFF,IEOFF
                  UCORR(I,JEOFF-1,II+1)=UD(WT(I,JEOFF-1,II),
     &                WT(I,JEOFF,II),UCORR(I,JEOFF-1,II),DTODY(JEOFF-1))
               ENDDO
            ENDIF

C           KORREKTUREN INNERHALB DES J-GEBIETES

            DO I=IOFF,IEOFF
               DO J=JOFF+JON,JEOFF-JEON
                  WT(I,J,II+1)= WT(I,J,II) - DTODY(J)*(
     &                 F(WT(I,J,II),WT(I,J+1,II),UCORR(I,J,II+1))-
     &                 F(WT(I,J-1,II),WT(I,J,II),UCORR(I,J-1,II+1)))
               ENDDO
            ENDDO

         ENDDO
C
C........DANN ADVEKTION IN X-RICHTUNG
C
C        ZURUECK ZU TRAC
C        INITIALISIERUNG DER KORREKTURGESCHWINDIGKEIT 
         DO J=JOFF,JEOFF
            DO I=IOFF,IEOFF
               TRAC(I,J,K,L)=WT(I,J,ncorr+1)
               UCORR(I,J,1)= U(I,J,K,NX)
               WT(I,J,1)=TRAC(I,J,K,L)
            ENDDO
         ENDDO

         DO J=JOFF,JEOFF
            DO I=IOFF+1,IEOFF-1
               WT(I,J,1)=TRAC(I,J,K,L)-DTODX(J)* (
     &              F(TRAC(I,J,K,L),TRAC(I+1,J,K,L),U(I,J,K,NX))-
     &              F(TRAC(I-1,J,K,L),TRAC(I,J,K,L),U(I-1,J,K,NX)) )
            ENDDO
         ENDDO
C
C........ANTIDIFFUSIONSKORREKTUR IN X-RICHTUNG
C            
         DO II=1,NCORR

C           UCORR: INNERHALB DES I-GEBIETES

            DO J=JOFF,JEOFF
               DO I=IOFF+1,IEOFF-2
                  UCORR(I,J,II+1)= UD( WT(I,J,II), WT(I+1,J,II),
     &                 UCORR(I,J,II), DTODX(J) )
               ENDDO
            ENDDO

C           WESTLICHE RANDBEDINGUNG FUER I=IOFF

            IF(ION .EQ. 3) THEN

               DO J=JOFF,JEOFF
                  IF(U(IOFF,J,K,NX).LT.0.) THEN
                     WT(IOFF,J,II)=WT(IOFF+1,J,II)-U(IOFF+1,J,K,NX)
     &                    /U(IOFF,J,K,NX)*
     &                    (WT(IOFF+2,J,II)-WT(IOFF+1,J,II))
                  ENDIF

                  UCORR(IOFF,J,II+1)=UD(WT(IOFF,J,II),WT(IOFF+1,J,II),
     &                 UCORR(IOFF,J,II),DTODX(J))

                  IF (U(IOFF,J,K,NX).LT.0.) THEN
                     WT(IOFF+1,J,II+1)=WT(IOFF+1,J,II)
                     WT(IOFF+2,J,II+1)=WT(IOFF+2,J,II)- DTODX(J)*(
     &                    F(WT(IOFF+2,J,II),WT(IOFF+3,J,II),
     &                    UCORR(IOFF+2,J,II+1))-
     &                    UCORR(IOFF+1,J,II+1)*WT(IOFF+2,J,II) )
                  ELSE
                     DO I=IOFF+1,IOFF+2
                        WT(I,J,II+1)= WT(I,J,II) - DTODX(J)*(
     &                    F(WT(I,J,II),WT(I+1,J,II),UCORR(I,J,II+1))-
     &                    F(WT(I-1,J,II),WT(I,J,II),UCORR(I-1,J,II+1)) )
                     ENDDO
                  ENDIF

               ENDDO

            ELSE

               DO J=JOFF,JEOFF
                  UCORR(IOFF,J,II+1)=UD(WT(IOFF,J,II),WT(IOFF+1,J,II),
     1                 UCORR(IOFF,J,II),DTODX(J))
               ENDDO

            ENDIF

C           OESTLICHE RANDBEDINGUNG FUER  I=IE-1

            IF(IEON .EQ. 3) THEN

               DO J=JOFF,JEOFF
                  IF (U(IEOFF-1,J,K,NX).GT.0.) THEN
                     WT(IEOFF,J,II)=WT(IEOFF-1,J,II)-U(IEOFF-2,J,K,NX)/
     &                    U(IEOFF-1,J,K,NX)* (WT(IEOFF-2,J,II)
     &                    -WT(IEOFF-1,J,II))
                  ENDIF

                  UCORR(IEOFF-1,J,II+1)=UD(WT(IEOFF-1,J,II),
     &                 WT(IEOFF,J,II),UCORR(IEOFF-1,J,II),DTODX(J))

                  IF (U(IEOFF-1,J,K,NX).GT.0.) THEN
                     WT(IEOFF-1,J,II+1)=WT(IEOFF-1,J,II)
                     WT(IEOFF-2,J,II+1)=WT(IEOFF-2,J,II)-DTODX(J)*(
     &                    UCORR(IEOFF-2,J,II+1)*WT(IEOFF-2,J,II)-
     &                    F(WT(IEOFF-3,J,II),WT(IEOFF-2,J,II),
     &                    UCORR(IEOFF-3,J,II+1)))
                  ELSE
                     DO I=IEOFF-2,IEOFF-1
                        WT(I,J,II+1)= WT(I,J,II) - DTODX(J)*(
     &                    F(WT(I,J,II),WT(I+1,J,II),UCORR(I,J,II+1))-
     &                    F(WT(I-1,J,II),WT(I,J,II),UCORR(I-1,J,II+1)) )
                     ENDDO
                  ENDIF
               ENDDO

            ELSE

               DO J=JOFF,JEOFF
                  UCORR(IEOFF-1,J,II+1)=UD(WT(IEOFF-1,J,II),
     &                 WT(IEOFF,J,II),UCORR(IEOFF-1,J,II),DTODX(J))
               ENDDO

            ENDIF

C           KORREKTURN INNERHALB DES I-GEBIETES

            DO J=JOFF,JEOFF
               DO I=IOFF+ION,IEOFF-IEON
                  WT(I,J,II+1)= WT(I,J,II) - DTODX(J)*(
     &                 F(WT(I,J,II),WT(I+1,J,II),UCORR(I,J,II+1))-
     &                 F(WT(I-1,J,II),WT(I,J,II),UCORR(I-1,J,II+1)) )

               ENDDO
            ENDDO

         ENDDO
C***********************************************************************
      ENDIF
C***********************************************************************

C
C.....ZURUECK ZU TRAC
C
      DO J=JOFF,JEOFF
         DO I=IOFF,IEOFF
            TRAC(I,J,K,L)=WT(I,J,NCORR+1)
         ENDDO
      ENDDO
C
      RETURN
C
C      CONTAINS
CC
C        REAL FUNCTION F(CI,CIP1,UIH)
C          IMPLICIT NONE
C          REAL :: CI,CIP1,UIH
C          F=(((ABS(UIH)+UIH)/2.)*CI  ) +
C     &      (((UIH-ABS(UIH))/2.)*CIP1)
C        END FUNCTION F
CCC
C        REAL FUNCTION UD(CI,CIP1,UIH,DT0DX)
C          IMPLICIT NONE
C          REAL :: CI,CIP1,UIH,DT0DX
C          UD=(ABS(UIH)-DT0DX*UIH*UIH)
C     &      *(CIP1-CI)/(CI+CIP1+1.E-35)
C        END FUNCTION UD
C     
      END SUBROUTINE SHADVEC
