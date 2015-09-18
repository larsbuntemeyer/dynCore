      SUBROUTINE MINVLAP(A,N)
      IMPLICIT NONE
C     Dummy Arguments
      INTEGER, INTENT(IN)     :: N
      REAL,    INTENT(INOUT) :: A(N,N)
C     Local Variables
      REAL    :: WORK(2*N), WORK1(N,N)
      INTEGER :: IPIV(N)
      INTEGER :: I,J,INFO
C
      DO I=1,N
         DO J=1,N
            WORK1(I,J)=A(I,J)
         ENDDO
      ENDDO
C
      CALL DGETRF(N,N,WORK1,N,IPIV,INFO)
      CALL DGETRI(N,WORK1,N,IPIV,WORK,2*N,INFO)
C
      DO I=1,N
         DO J=1,N
            A(I,J)=WORK1(I,J)
         ENDDO
      ENDDO
C
      RETURN
      END SUBROUTINE MINVLAP
