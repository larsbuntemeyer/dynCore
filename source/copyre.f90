!

! Code converted using TO_F90 by Alan Miller
! Date: 2015-09-07  Time: 14:27:29

!     ----------------------------------------------------------------

SUBROUTINE copyre(pa,pb,klen)

IMPLICIT NONE

!**   *SUBROUTINE* *COPYRE*  COPY *KLEN* REAL VALUES FROM
!                            *PA* TO *PB*.


REAL, INTENT(IN)           :: pa(klen)
REAL, INTENT(OUT)          :: pb(klen)
INTEGER, INTENT(IN)        :: klen

!     Dummy Arguments





!     Local Variables

INTEGER :: jk
!     ----------------------------------------------------------------

!*        1.       COPY PA INTO PB.
!                  ---- -- ---- ---
DO jk=1,klen
  pb(jk)=pa(jk)
END DO

!     ----------------------------------------------------------------

RETURN
END SUBROUTINE copyre
