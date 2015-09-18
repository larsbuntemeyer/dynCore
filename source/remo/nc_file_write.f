
      subroutine nc_file_write( varname, tndx, twod_in, threed_in )
c-----------------------------------------------------------------------
c     ... write netcdf output file
c-----------------------------------------------------------------------

      include 'org.h'
      include 'parorg.h'

c-----------------------------------------------------------------------
c     ... dummy arguments
c-----------------------------------------------------------------------
      integer, intent(in)          :: tndx
      character(len=*), intent(in) :: varname
      real, optional, intent(in)   :: twod_in(:,:)
      real, optional, intent(in)   :: threed_in(:,:,:)
      
c-----------------------------------------------------------------------
c     ... local variables
c-----------------------------------------------------------------------
      integer :: k, zcnt
      integer :: ncid, nc_stat, varid
      integer :: dims(4), strt_ndx(4), length(4)
      real :: twod_out(moie,moje)
      real :: threed_out(moie,moje,moke)

      character(len=132) :: mess

      include 'netcdf.inc'

c---------------------------------------------------------------------
c	... check for variable
c---------------------------------------------------------------------
      if( .not. present( twod_in) .and.
     1    .not. present( threed_in ) ) then
        return
      endif
c---------------------------------------------------------------------
c	... gather variable to root process
c---------------------------------------------------------------------
      if( present( twod_in ) ) then
        call collectdata( twod_out, twod_in, 200, 0 )
      elseif( present( threed_in ) ) then
        select case( trim(varname) )
          case( 'DWDT', 'T' )
            zcnt = ke
          case( 'PINT', 'W' )
            zcnt = ke1
        end select
        do k = 1,zcnt
          call collectdata( threed_out(1,1,k), threed_in(1,1,k),
     1                      300, k )
        end do
      endif
          
c---------------------------------------------------------------------
c	... only the root process writes
c---------------------------------------------------------------------
      if( myid == 0 ) then
        if( trim(varname) == 'APRL' ) then
        write(*,*) 'nc_file_write: diagnostics'
        write(*,*) 'nc_file_write: twod_in dims = ',
     1     size(twod_in,dim=1),size(twod_in,dim=2)
        write(*,*) 'nc_file_write: iah,ieh,jah,jeh = ',
     1       iah,ieh,jah,jeh
        write(*,*) 'nc_file_write: moie,moje = ',moie,moje
          write(*,*) 'nc_file_write: twod_in'
          do k = 2,ieh,5
            write(*,'(1p5g15.7)') twod_in(k:min(k+4,ieh),3)
          end do
          write(*,*) 'nc_file_write: twod_out dims = ',
     1       size(twod_out,dim=1),size(twod_out,dim=2)
          write(*,*) 'nc_file_write: twod_out'
          do k = 1,ieh-1,5
            write(*,'(1p5g15.7)') twod_out(k:min(k+4,ieh-1),2)
          end do
          write(*,*) 
     1  'min,max APRL = ',minval(twod_out),maxval(twod_out)
        endif
c---------------------------------------------------------------------
c	... open file
c---------------------------------------------------------------------
        nc_stat = nf_open( 'remo_out.nc', nf_write, ncid )
        if( nc_stat /= nf_noerr ) then
           mess = 'nc_file_init: failed to open remo_out.nc'
           call handle_ncerr
        endif

        nc_stat = nf_inq_varid( ncid, 'times', varid )
        nc_stat = nf_put_vara_int( ncid, varid, (/tndx/), (/1/), tndx )
c---------------------------------------------------------------------
c	... write the variable
c---------------------------------------------------------------------
        select case( trim(varname) )
          case( 'PS', 'APRL', 'TSRF', 'PINT_SRF' )
            nc_stat = nf_inq_varid( ncid, trim(varname), varid )
            if( nc_stat /= nf_noerr ) then
              mess = 'nc_file_init: failed to get variable '
     1               // trim(varname) // ' id'
              call handle_ncerr
            endif
            strt_ndx(:3) = (/ 1, 1, tndx /)
            length(:3) = (/ moie, moje, 1 /)
            nc_stat = nf_put_vara_double( ncid, varid, strt_ndx(:3),
     1                                  length(:3), twod_out )
            if( nc_stat /= nf_noerr ) then
              mess = 'nc_file_init: failed to write variable '
     1                // trim(varname)
              call handle_ncerr
            endif
          case( 'DWDT', 'PINT', 'T', 'W' )
            nc_stat = nf_inq_varid( ncid, trim(varname), varid )
            if( nc_stat /= nf_noerr ) then
              mess = 'nc_file_init: failed to get variable '
     1               // trim(varname) // ' id'
              call handle_ncerr
            endif
            strt_ndx(:) = (/ 1, 1, 1, tndx /)
            length(:) = (/ moie, moje, zcnt, 1 /)
            nc_stat = nf_put_vara_double( ncid, varid, strt_ndx,
     1                                  length, threed_out )
            if( nc_stat /= nf_noerr ) then
              mess = 'nc_file_init: failed to write '
     1               // trim(varname)
              call handle_ncerr
            endif
        end select

c---------------------------------------------------------------------
c     	... close the file
c---------------------------------------------------------------------
        nc_stat = nf_close( ncid )
        if( nc_stat /= nf_noerr ) then
           mess = 'nc_file_init: failed to close remo_out.nc'
           call handle_ncerr
        endif
      endif

      CONTAINS

      subroutine handle_ncerr
c---------------------------------------------------------------------
c	... netcdf error handling routine
c---------------------------------------------------------------------

      write(*,*) 'handle_ncerr: ' // trim(mess)
      write(*,*) nf_strerror( ret )
      call mpi_finalize( ret )
      stop 'netcdf error'

      end subroutine handle_ncerr

      end subroutine nc_file_write
