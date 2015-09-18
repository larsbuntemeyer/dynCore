!
!     COMPUTE THE ALPHA TERM FOR THE LATERAL BOUNDARY TREATMENT FROM RMY
!     AND FAKRMY (SEE EM/DM MANUAL CHAPTER II.4)
!
!     WRITTEN BY KEVIN SIECK
!
      SUBROUTINE GETALPHA(RMY, FAKRMY, ALPHABOUND)
!
      IMPLICIT NONE
!
      INCLUDE "parorg.h"
!
!     DUMMY ARGUMENTS
!
      REAL, INTENT(IN) :: RMY(IE,JE,3)
      REAL, INTENT(IN) :: FAKRMY
      REAL, INTENT(OUT):: ALPHABOUND(IE,JE,3)
!
!     LOCAL VARIABLES
!
      INTEGER :: I, J
!
!     COMPUTE ALPHA VALUE
!     INDEX 1: MASS POINT VARIABLES
!     INDEX 2: U VELOCITY COMPONENT
!     INDEX 3: V VELOCITY COMPONENT
!
      DO J = 1, JE
        DO I = 1, IE
          ALPHABOUND(I,J,1) = FAKRMY*RMY(I,J,1)/
     &         (1. + FAKRMY*RMY(I,J,1))
          ALPHABOUND(I,J,2) = FAKRMY*RMY(I,J,2)/
     &         (1. + FAKRMY*RMY(I,J,2))
          ALPHABOUND(I,J,3) = FAKRMY*RMY(I,J,3)/
     &         (1. + FAKRMY*RMY(I,J,3))
        END DO
      END DO
!
!     SET THE BOUNDARIES TO 1. (ONLY NECESSARY FOR MASS POINT VARIABLES)
!
      IF (NEIGHBOR(1) .EQ. -1) THEN
        ALPHABOUND(1:2,:,1) = 1.0
      END IF

      IF (NEIGHBOR(2) .EQ. -1) THEN
        ALPHABOUND(:,JE-1:JE,1) = 1.0
      END IF

      IF (NEIGHBOR(3) .EQ. -1) THEN
        ALPHABOUND(IE-1:IE,:,1) = 1.0
      END IF

      IF (NEIGHBOR(4) .EQ. -1) THEN
        ALPHABOUND(:,1:2,1) = 1.0
      END IF

      RETURN
      END SUBROUTINE GETALPHA