!
module Io
   !
   use Grid 
   use Database
   !
   implicit none
   !
   public :: Io_init, Io_write_to_file
   !
   character(len=80),save   :: outdata_base
   character(len=80),save   :: indata_base
   integer          ,save   :: outdata_unit
   !
   !
   contains 
   !
subroutine Io_init
   !
   implicit none
   !
   write(*,*) '----- Io_init ---------------------'
   !
   outdata_unit = 2 
   outdata_base = 'output.dat'
   open(unit=outdata_unit,file=outdata_base)
   !
   write(*,*) 'output data file:',outdata_base
   write(*,*) '----- Io_init done ----------------'
   ! 
end subroutine Io_init
   !
   subroutine Io_write_to_file
      !
      implicit none
      integer :: i,j,k
      !
      do i=ib,ie
         do j=jb,je
            do k=kb,ke
               !
               !write(outdata_unit,'(3I5,3F13.8)') i,j,k,dens(i,j,k),u(i,j,k),pres(i,j,k)
               write(outdata_unit,'(4F13.8)') xcCoord(i),dens(i,j,k),u(i,j,k),pres(i,j,k)
               !
            enddo  
         enddo  
      enddo  
      !       
   end subroutine Io_write_to_file
   !
end module Io
!
