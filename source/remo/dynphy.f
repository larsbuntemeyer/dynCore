      SUBROUTINE DYNPHY(FI,NLON,NLEV,NSTART,NSTOP,NTASK)
      IMPLICIT NONE
C
      INTEGER, INTENT(IN)     ::   NLON,NLEV,NTASK
      INTEGER, INTENT(IN)     ::   NSTART(NTASK),NSTOP(NTASK)
      REAL,    INTENT(INOUT)  ::   FI(NLON,NLEV)   
C
      REAL    :: ZH(NLON*NLEV)
      INTEGER :: I,J,K,IC
C
      IC=0
      DO K=1,NTASK
         DO J=1,NLEV
            DO I=NSTART(K),NSTOP(K)
               IC=IC+1
               ZH(IC)=FI(I,J)
            ENDDO
         ENDDO
      ENDDO
C
      IC=0
      DO J=1,NLEV
         DO I=1,NLON
            IC=IC+1
            FI(I,J)=ZH(IC)
         ENDDO
      ENDDO
C
      RETURN
      END SUBROUTINE DYNPHY
