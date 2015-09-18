SUBROUTINE remark(ytext)

! Code converted using TO_F90 by Alan Miller
! Date: 2015-09-07  Time: 14:27:18

IMPLICIT NONE
CHARACTER (LEN=*), intent(in) :: ytext
WRITE(*,'(A)') ytext
RETURN

END SUBROUTINE remark
