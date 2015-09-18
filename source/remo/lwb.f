C
C**** *LWB*   - COMPUTES BLACK-BODY FUNCTIONS FOR LONGWAVE CALCULATIONS
C
C     PURPOSE.
C     --------
C           COMPUTES PLANCK FUNCTIONS
C
C**   INTERFACE.
C     ----------
C      SUBROUTINE LWB ( KLON,KLEV,KAER,KMODE,KFLUX,KRAD
C    S  , PDT0,PTAVE,PTL
C    S  , PBINT,PBSUIN,PBSUR,PBTOP,PDBSL                   )
C    S  , PB,PBINT,PBSUIN,PBSUR,PBTOP,PDBSL
C    S  , PGA,PGB,PGASUR,PGBSUR,PGATOP,PGBTOP
C    S  , PGC,PGD,PGCSUR,PGDSUR,PGCTOP,PGDTOP    )
C        EXPLICIT ARGUMENTS :
C        --------------------
C     ==== INPUTS ===
C PTAVE  : (KLON,KLEV)       ; TEMPERATURE
C PTL    : (KLON,0:KLEV)     ; HALF LEVEL TEMPERATURE
C     ==== OUTPUTS ===
C PBINT  : (KLON,0:KLEV)     ; HALF LEVEL PLANCK FUNCTION
C PBSUIN : (KLON)             ; SURFACE PLANCK FUNCTION
C PBTOP  : (KLON,NINT)        ; TOP SPECTRAL PLANCK FUNCTION
C PDBSL  : (KLON,NINT,KLEV*2); SUB-LAYER PLANCK FUNCTION GRADIENT
C PGA    : (KLON,8,2,KLEV); DB/DT-WEIGHTED LAYER PADE APPROXIMANTS
C PGB    : (KLON,8,2,KLEV); DB/DT-WEIGHTED LAYER PADE APPROXIMANTS
C PGC   :(KLON,5,2,KLEV); DB/DT-WEIGHTED LAYER PADE APPROXIMANTS
C PGD   :(KLON,5,2,KLEV); DB/DT-WEIGHTED LAYER PADE APPROXIMANTS
C
C        IMPLICIT ARGUMENTS :   NONE
C        --------------------
C
C     METHOD.
C     -------
C
C          1. COMPUTES THE PLANCK FUNCTION ON ALL LEVELS AND HALF LEVELS
C     FROM A POLYNOMIAL DEVELOPMENT OF PLANCK FUNCTION
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
C-----------------------------------------------------------------------
      SUBROUTINE LWB ( KLON,KLEV
     &  , PTAVE,PTL
     &  , PB,PBINT,PBSUIN,PBTOP,PDBSL
     &  , PGA,PGB,PGASUR,PGBSUR,PGATOP,PGBTOP
     &  , PGC,PGD,PGCSUR,PGDSUR,PGCTOP,PGDTOP    )
C
      IMPLICIT NONE
C
      INCLUDE "YOMLW"
C-----------------------------------------------------------------------
C        Dummy Arguments
C*       0.1   ARGUMENTS
C              ---------
C
      INTEGER, INTENT(IN)    :: KLON, KLEV
      REAL,    INTENT(IN)    :: PTAVE(KLON,KLEV),PTL(KLON,KLEV+1)
C
      REAL,    INTENT(INOUT) :: 
     &     PB(KLON,NINT,KLEV+1), PBINT(KLON,KLEV+1)
     &  ,  PBSUIN(KLON)
     &  ,  PBTOP(KLON,NINT), PDBSL(KLON,NINT,KLEV*2)
     &  ,  PGA(KLON,8,2,KLEV), PGB(KLON,8,2,KLEV)
     &  ,  PGASUR(KLON,8,2), PGBSUR(KLON,8,2)
     &  ,  PGATOP(KLON,8,2), PGBTOP(KLON,8,2)
     &  ,  PGC(KLON,5,2,KLEV), PGD(KLON,5,2,KLEV)
     &  ,  PGCSUR(KLON,5,2),PGDSUR(KLON,5,2)
     &  ,  PGCTOP(KLON,5,2),PGDTOP(KLON,5,2)
C
C-----------------------------------------------------------------------
C        Local Variables
C*       0.2   LOCAL ARRAYS
C              ------------
      INTEGER :: INDB(KLON),INDS(KLON)
      REAL    :: ZBLAY(KLON,KLEV),ZBLEV(KLON,KLEV+1)
C
      INTEGER :: INDSU, INDT, INDTO, INDTP, IXTOX, IXTX, JF,  
     &           JG, JK, JK1, JK2, JL, JNU
      REAL    :: ZDST1, ZDSTO1, ZDSTOX, ZDSTX, ZRES, ZRES2, ZTI, ZTI2
C
C     ------------------------------------------------------------------
C
C
C*         1.0     PLANCK FUNCTIONS AND GRADIENTS
C                  ------------------------------
C
      DO JK = 1 , KLEV+1
         DO JL = 1 , KLON
            PBINT(JL,JK) = 0.
         ENDDO
      ENDDO
      DO JL = 1 , KLON
         PBSUIN(JL) = 0.
      ENDDO
C
      DO JNU=1,NINT
C
C
C*         1.1   LEVELS FROM SURFACE TO KLEV

C                ----------------------------
C
         DO JK = 1 , KLEV
            DO JL = 1 , KLON
               ZTI=(PTL(JL,JK)-TSTAND)/TSTAND
               ZRES = XP(1,JNU)+ZTI*(XP(2,JNU)+ZTI*(XP(3,JNU)
     &          +ZTI*(XP(4,JNU)+ZTI*(XP(5,JNU)+ZTI*(XP(6,JNU) )))))
               PBINT(JL,JK)=PBINT(JL,JK)+ZRES
               PB(JL,JNU,JK)= ZRES
               ZBLEV(JL,JK) = ZRES
               ZTI2=(PTAVE(JL,JK)-TSTAND)/TSTAND
               ZRES2=XP(1,JNU)+ZTI2*(XP(2,JNU)+ZTI2*(XP(3,JNU)
     &        +ZTI2*(XP(4,JNU)+ZTI2*(XP(5,JNU)+ZTI2*(XP(6,JNU) )))))
               ZBLAY(JL,JK) = ZRES2
            ENDDO
         ENDDO
C
C
C*         1.2   TOP OF THE ATMOSPHERE AND SURFACE
C                ---------------------------------
C
         DO JL = 1 , KLON
            ZTI=(PTL(JL,KLEV+1)-TSTAND)/TSTAND
            ZTI2 = (PTL(JL,1) - TSTAND) / TSTAND
            ZRES = XP(1,JNU)+ZTI*(XP(2,JNU)+ZTI*(XP(3,JNU)
     &       +ZTI*(XP(4,JNU)+ZTI*(XP(5,JNU)+ZTI*(XP(6,JNU) )))))
            ZRES2 = XP(1,JNU)+ZTI2*(XP(2,JNU)+ZTI2*(XP(3,JNU)
     &       +ZTI2*(XP(4,JNU)+ZTI2*(XP(5,JNU)+ZTI2*(XP(6,JNU) )))))
            PBINT(JL,KLEV+1) = PBINT(JL,KLEV+1)+ZRES
            PB(JL,JNU,KLEV+1)= ZRES
            ZBLEV(JL,KLEV+1) = ZRES
            PBTOP(JL,JNU) = ZRES
            PBSUIN(JL) = PBSUIN(JL) + ZRES2
         ENDDO
C
C
C*         1.3   GRADIENTS IN SUB-LAYERS
C                -----------------------
C
         DO JK = 1 , KLEV
            JK2 = 2 * JK
            JK1 = JK2 - 1
            DO JL = 1 , KLON
               PDBSL(JL,JNU,JK1) = ZBLAY(JL,JK  ) - ZBLEV(JL,JK)
               PDBSL(JL,JNU,JK2) = ZBLEV(JL,JK+1) - ZBLAY(JL,JK)
            ENDDO
         ENDDO
C
      ENDDO
C
C*         2.0   CHOOSE THE RELEVANT SETS OF PADE APPROXIMANTS
C                ---------------------------------------------
C
      DO JL=1,KLON
         ZDSTO1 = (PTL(JL,KLEV+1)-TINTP(1)) / TSTP
         IXTOX = MAX0( 1, MIN0( MXIXT, INT( ZDSTO1 + 1. ) ) )
         ZDSTOX = (PTL(JL,KLEV+1)-TINTP(IXTOX))/TSTP
CRP   CVMGT ERSETZT
C     INDTO = CVMGT(IXTOX,IXTOX+1, ZDSTOX.LT.0.5 )
         IF (ZDSTOX.LT.0.5) THEN
            INDTO = IXTOX
         ELSE
            INDTO = IXTOX+1
         ENDIF
         INDB(JL)=INDTO
         ZDST1 = (PTL(JL,1)-TINTP(1)) / TSTP
         IXTX = MAX0( 1, MIN0( MXIXT, INT( ZDST1 + 1. ) ) )
         ZDSTX = (PTL(JL,1)-TINTP(IXTX))/TSTP
CRP   CVMGT ERSETZT
C     INDT = CVMGT(IXTX,IXTX+1, ZDSTX.LT.0.5 )
         IF (ZDSTX.LT.0.5) THEN
            INDT = IXTX
         ELSE
            INDT = IXTX+1
         ENDIF
         INDS(JL)=INDT
      ENDDO
C
      DO JF=1,2
         DO JG=1, 8
            DO JL=1,KLON
               INDSU=INDS(JL)
               PGASUR(JL,JG,JF)=GA(INDSU,2*JG-1,JF)
               PGBSUR(JL,JG,JF)=GB(INDSU,2*JG-1,JF)
               INDTP=INDB(JL)
               PGATOP(JL,JG,JF)=GA(INDTP,2*JG-1,JF)
               PGBTOP(JL,JG,JF)=GB(INDTP,2*JG-1,JF)
            ENDDO
         ENDDO
      ENDDO
C
      DO JF=1,2
         DO JG=1,5
            DO JL=1,KLON
               INDSU=INDS(JL)
               PGCSUR(JL,JG,JF)=GC(INDSU,2*JG-1,JF)
               PGDSUR(JL,JG,JF)=GD(INDSU,2*JG-1,JF)
               INDTP=INDB(JL)
               PGCTOP(JL,JG,JF)=GC(INDTP,2*JG-1,JF)
               PGDTOP(JL,JG,JF)=GD(INDTP,2*JG-1,JF)
            ENDDO
         ENDDO
      ENDDO
C
      DO JK=1,KLEV
         DO JL=1,KLON
            ZDST1 = (PTAVE(JL,JK)-TINTP(1)) / TSTP
            IXTX = MAX0( 1, MIN0( MXIXT, INT( ZDST1 + 1. ) ) )
            ZDSTX = (PTAVE(JL,JK)-TINTP(IXTX))/TSTP
CRP   CVMGT ERSETZT
C     INDT = CVMGT(IXTX,IXTX+1, ZDSTX.LT.0.5 )
            IF (ZDSTX.LT.0.5) THEN
               INDT = IXTX
            ELSE
               INDT = IXTX+1
            ENDIF
            INDB(JL)=INDT
         ENDDO
C
         DO JF=1,2
            DO JG=1, 8
               DO JL=1,KLON
                  INDT=INDB(JL)
                  PGA(JL,JG,JF,JK)=GA(INDT,2*JG,JF)
                  PGB(JL,JG,JF,JK)=GB(INDT,2*JG,JF)
               ENDDO
            ENDDO
         ENDDO
         DO JF=1,2
            DO JG=1,5
               DO JL=1,KLON
                  INDT=INDB(JL)
                  PGC(JL,JG,JF,JK)=GC(INDT,2*JG,JF)
                  PGD(JL,JG,JF,JK)=GD(INDT,2*JG,JF)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C
C     ------------------------------------------------------------------
C
      RETURN
      END SUBROUTINE LWB
