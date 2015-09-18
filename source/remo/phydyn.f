      SUBROUTINE PHYDYN(FI,NLON,NLEV,NSTART,NSTOP,NTASK)
      !
      ! This subroutine is not used.
      !
      IMPLICIT NONE
      !
      ! Dummy arguments
      !
      INTEGER, INTENT(IN)     :: NLON,NLEV,NTASK
      INTEGER, INTENT(IN)     :: NSTART(NTASK),NSTOP(NTASK)
      REAL, INTENT(INOUT)     :: FI(NLON*NLEV)
      !
      ! Local variables
      !
      INTEGER     ::     I,J,K,IC
      REAL        ::     ZH(NLON,NLEV)
      !
      IC=0
      DO K=1,NTASK
         DO J=1,NLEV
            DO I=NSTART(K),NSTOP(K)
               IC=IC+1
               ZH(I,J)=FI(IC)
            ENDDO
         ENDDO
      ENDDO
      !
      IC=0
      DO J=1,NLEV
         DO I=1,NLON
            IC=IC+1
            FI(IC)=ZH(I,J)
         ENDDO
      ENDDO
      !
      RETURN
      END SUBROUTINE PHYDYN
