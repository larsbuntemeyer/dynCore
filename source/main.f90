!
!
!
program main 
   !
   use Driver
   !
   implicit none
   !
   call write_banner
   !
   call Driver_init
   call Driver_evolve
   call Driver_finish
   !
   write(*,*) ''   
   write(*,*) 'finished DynCore'   
   write(*,*) ''   
   !
contains 
   !
   subroutine write_banner
      !
      implicit none 
      !
      write(*,*) ''
      write(*,*) ''
      write(*,*) '====================================================='
      write(*,*) '                      DynCore                        '
      write(*,*) '====================================================='
      write(*,*) ''
      write(*,*) ''
      !
   end subroutine write_banner
   !
end program main 
!
